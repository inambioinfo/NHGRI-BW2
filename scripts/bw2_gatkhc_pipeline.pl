#!/usr/bin/perl -w

# $Id: bw2_gatkhc_pipeline.pl 7656 2016-03-17 00:36:18Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use NHGRI::BW2::Session;

our %Opt;

$ENV{PATH} = "/home/nhansen/projects/gtb/perl/scripts:$ENV{PATH}";

=head1 NAME

bw2_gatkhc_pipeline.pl - Run a set of fastq sequence files throught the GATK germline variant calling pipeline.

=head1 SYNOPSIS

  bw2_gatkhc_pipeline.pl --fastq_fof file_of_fastq_file_pairs --merged_bam name_of_analysis_ready_bam --ucsc_ref ref_to_align_to --analysis_name name_of_analysis [options]

For complete documentation, run C<bw2_gatkhc_pipeline.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of pairs of fastq files (or single fastq files for single end reads), and copies them to the biowulf2 cluster, aligning them with bwa mem, adding readtags, marking duplicates, realigning around known indels and recalibrating base qualities, and then finally running the GATK HaplotypeCaller to call diploid germline variants and genotypes.  Resulting VCF files, and the analysis-ready BAM file, are copied back to the local server.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $analysis_name = $Opt{analysis_name};
if (!$analysis_name) { # create a unique one
    my $pid = $$;
    my $time = time();
    $analysis_name = "gatkhcpipeline.$pid.$time";
}

my $nisc_dir = $Opt{nisc_dir};
if ($nisc_dir) { # process input files from a NISC MPG whole exome directory
    create_fastq_files_from_exome_directory($nisc_dir, $analysis_name);
}

my $fastq_fof = $Opt{fastq_fof};
my $tags_file = $Opt{tags_file};
my $ucsc_ref = $Opt{ucsc_ref};
my $merged_bam = $Opt{merged_bam};

## STEP 1 - align with BWA ##
my $rmdup_bamfile_fof = run_bwa($fastq_fof, $tags_file, $ucsc_ref, $analysis_name);

## STEP 2 - mark PCR duplicates ##
my $indel_target_bamfile_fof = mark_duplicates($rmdup_bamfile_fof, $ucsc_ref, $analysis_name);

## STEP 3 - realign around known indels ##
my $recal_bamfile_fof = realign_indels($indel_target_bamfile_fof, $ucsc_ref, $analysis_name);

## STEP 4 - recalibrate base qualities ##
my $merge_bamfile_fof = recalibrate_qualities($recal_bamfile_fof, $ucsc_ref, $analysis_name);

## STEP 5 - merge to one BAM file ##
merge_bamfiles($merge_bamfile_fof, $merged_bam, $analysis_name);

## STEP 6 - call GATK HaplotypeCaller
call_gatkhc($merged_bam, $ucsc_ref, $analysis_name);

## STEP 7 - check output files ##
check_output();

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( ucsc_ref => 'hg19', merged_bam => 'sample.bam'
           );
    GetOptions(\%Opt, qw( fastq_fof=s tags_file=s nisc_dir=s ucsc_ref=s merged_bam=s 
           analysis_name=s skip_bam2fastq skip_bwa skip_markdup skip_realign 
           skip_recalibrate skip_merge variantonly manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_gatkhc_pipeline.pl, ", q$Revision: 7656 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --fastq_fof!") if (!$Opt{nisc_dir} && !$Opt{fastq_fof});
    pod2usage("Missing required field --tags_file!") if (!$Opt{nisc_dir} && !$Opt{tags_file});
}

sub create_fastq_files_from_exome_directory {
    my $nisc_dir = shift;
    my $analysis_name = shift;

    # find logfile:
    opendir NISC, $nisc_dir
        or die "Couldn\'t read directory $nisc_dir!\n";
    my @potential_logfiles = grep /\.log$/, readdir NISC;
    closedir NISC;

    my $logfile;
    if (@potential_logfiles == 1) {
        $logfile = "$nisc_dir/$potential_logfiles[0]";
    }
    else {
        die "Unable to find logfile in directory $nisc_dir!\n";
    }
    my $sample = $logfile;
    $sample =~ s:.*/::;
    $sample =~ s:\.log$::;

    my $solexa_read_string = `grep 'File is' $logfile | grep 'solexa' | awk '{print \$6}' | sed 's:,::'`;
    my @solexa_read_bams = split /\n/, $solexa_read_string;
    
    if (!$Opt{'skip_bam2fastq'}) {

        my $bamfile_fof = "$analysis_name.lanebams.fof";
        my $bamfile_fh = Open($bamfile_fof, "w");
        foreach my $solexa_bam (@solexa_read_bams) {
            print $bamfile_fh "$solexa_bam\n";
        }
    
        close $bamfile_fh;
    
        # run bam2fastq on the NISC lane bams:
        system("bw2_bam2fastq.pl --bamfile_fof $bamfile_fof --analysis_name $analysis_name --skip_retrieve") == 0
            or die "Couldn\'t run bam2fastq on NISC lane bam files for analysis $analysis_name";
    }

    my $fastq_fof = "$analysis_name.bwfastqs.fof";
    $Opt{fastq_fof} = $fastq_fof;
    my $fastq_fh = Open($fastq_fof, "w");
    my $tag_fof = "$analysis_name.tags.fof";
    $Opt{tags_file} = $tag_fof;
    my $tag_fh = Open($tag_fof, "w");
    my $session_obj = NHGRI::BW2::Session->new();
    foreach my $solexa_bam (@solexa_read_bams) {
        my $bam_base = $solexa_bam;
        $bam_base =~ s:.*/::;
        $bam_base =~ s:\.bam$::;
        my $fastq_string = $session_obj->remote_command(-command => "ls /data/$ENV{USER}/$analysis_name/output_files/$bam_base.*.fq.gz");

        my @fastq_files = split /\n/, $fastq_string;

        my ($fastq0_file) = grep /$bam_base\.0\.fq.gz$/, @fastq_files;
        my ($fastq1_file) = grep /$bam_base\.1\.fq.gz$/, @fastq_files;
        my ($fastq2_file) = grep /$bam_base\.2\.fq.gz$/, @fastq_files;

        #print $fastq_fh "$fastq1_file\t$fastq2_file\n";
        my ($id, $lib, $fc);
        if ($solexa_bam =~ m:/([^/]+)\.(\d+)\.(\d+)\.bam:) {
            $id = "$1.$2";
            $lib = $3;
            $fc = $1;
        }
        if ($fastq0_file) {
            print $fastq_fh "$fastq0_file\n";
            print $tag_fh "ID:$id\tSM:$sample\tLB:$lib\tFC:$fc\n";
        }
        if ($fastq1_file && $fastq2_file) {
            print $fastq_fh "$fastq1_file\t$fastq2_file\n";
            print $tag_fh "ID:$id\tSM:$sample\tLB:$lib\tFC:$fc\n";
        }
    }
    
    close $fastq_fh;
    close $tag_fh;

    return 0;
}

sub run_bwa {
    my $fastq_fof = shift;
    my $tags_file = shift;
    my $ref = shift;
    my $analysis_name = shift;

    my $markdup_file = "$analysis_name.markdupbams.fof";

    if (!$Opt{skip_bwa}) {
        my $bwfastq_string = ($Opt{'nisc_dir'}) ? ' --bwfastqs' : '';
        system("bw2_bwa_align.pl --fastq_fof $fastq_fof --tags_file $tags_file --ref $ref --analysis_name $analysis_name$bwfastq_string --skip_retrieve") == 0
            or die "Couldn\'t perform BWA alignment step for analysis $analysis_name";

        # write a file of BAM file locations for the mark duplicate step:
    
        my $out_fh = Open($markdup_file, "w");
        my $in_fh = Open($fastq_fof);
    
        while (<$in_fh>) {
            chomp;
            my @fields = split /\t/, $_;
            my $rgbam = $fields[0];
            $rgbam =~ s:.*/:/data/$ENV{USER}/$analysis_name/output_files/:;
            $rgbam =~ s:\.[01]\.f(ast){0,1}q(\.gz){0,1}:.bam:;
            print $out_fh "$rgbam\n";
        }

        close $in_fh;
        close $out_fh;
    }

    return $markdup_file;
}

sub mark_duplicates {
    my $bamfile_fof = shift;
    my $ref = shift;
    my $analysis_name = shift;

    my $realign_file = "$analysis_name.indelrealign.fof";

    if (!$Opt{skip_markdup}) {
        system("bw2_mark_duplicates.pl --bam_paths $bamfile_fof --bwbams --analysis_name $analysis_name --skip_retrieve --index") == 0
            or die "Couldn\'t mark duplicates for analysis $analysis_name";
    
        # write a file of BAM file locations for the indel realignment step:
    
        my $out_fh = Open($realign_file, "w");
        my $in_fh = Open($bamfile_fof);
    
        while (<$in_fh>) {
            chomp;
            my @fields = split /\t/, $_;
            my $realignbam = $fields[0];
            $realignbam =~ s:\.bam$:.markdup.bam:;
            print $out_fh "$realignbam\n";
        }
    
        close $in_fh;
        close $out_fh;
    }

    return $realign_file;
}

sub realign_indels {
    my $bamfile_fof = shift;
    my $ref = shift;
    my $analysis_name = shift;

    my $realign_file = "$analysis_name.targetinfo.fof";
    my $recalbams_fof = "$analysis_name.recalbams.fof";

    if (!$Opt{skip_realign}) {
        system("bw2_create_indeltargets.pl --bam_paths $bamfile_fof --ref $ref --known_resource gold --bwbams --analysis_name $analysis_name --skip_retrieve") == 0
            or die "Couldn\'t create targets for indel realignment for analysis $analysis_name";

        # write a file of BAM file locations for the indel realignment step:
    
        my $out_fh = Open($realign_file, "w");
        my $in_fh = Open($bamfile_fof);
    
        while (<$in_fh>) {
            chomp;
            my @fields = split /\t/, $_;
            my $realignbam = $fields[0];
            my $intervalfile = $realignbam;
            $intervalfile =~ s/\.bam/.intervals/;
            print $out_fh "$realignbam\t$intervalfile\n";
        }
    
        close $in_fh;
        close $out_fh;
    
        system("bw2_realign_indels.pl --bam_info $realign_file --ref $ref --known_resource gold --bwbams --analysis_name $analysis_name --skip_retrieve") == 0
            or die "Couldn\'t realign around known indels for analysis $analysis_name";
    
        # write a file of BAM file locations for the recalibration step:
    
        $out_fh = Open($recalbams_fof, "w");
        $in_fh = Open($bamfile_fof);
    
        while (<$in_fh>) {
            chomp;
            my @fields = split /\t/, $_;
            my $realignedbam = $fields[0];
            $realignedbam =~ s/\.bam/.realigned.bam/;
            print $out_fh "$realignedbam\n";
        }
    
        close $in_fh;
        close $out_fh;
    }
    
    return $recalbams_fof;
}

sub recalibrate_qualities {
    my $recal_fof = shift;
    my $ref = shift;
    my $analysis_name = shift;

    my $mergebams_fof = "$analysis_name.mergebams.fof";

    if (!$Opt{skip_recalibrate}) {
        system("bw2_recalibrate_bases.pl --bam_paths $recal_fof --ref $ref --bwbams --analysis_name $analysis_name --skip_retrieve") == 0
        or die "Couldn\'t recalibrate qualities for analysis $analysis_name";
    
        # write a file of BAM file locations for the merge step:
    
        my $out_fh = Open($mergebams_fof, "w");
        my $in_fh = Open($recal_fof);
    
        while (<$in_fh>) {
            chomp;
            my @fields = split /\t/, $_;
            my $recalibratedbam = $fields[0];
            $recalibratedbam =~ s/\.bam/.recal.bam/;
            print $out_fh "$recalibratedbam\n";
        }
    
        close $in_fh;
        close $out_fh;

    }

    return $mergebams_fof;
}

sub merge_bamfiles {
    my $mergebams_fof = shift;
    my $merged_bam = shift;
    my $analysis_name = shift;

    if (!$Opt{skip_merge}) {
        system("bw2_merge_bamfiles.pl --bam_paths $mergebams_fof --merged_bam $merged_bam --analysis_name $analysis_name --skip_retrieve --index --bwbams") == 0
            or die "Couldn\'t merge bams for analysis $analysis_name";
    }

    return 0;

}

sub call_gatkhc {
    my $bamfile = shift;
    my $ref = shift;
    my $analysis_name = shift;

    # BAM file resides in the output directory on the biowulf server:
    $bamfile = "/data/$ENV{USER}/$analysis_name/output_files/$bamfile";

    my $bw_gatk_string = "bw2_call_gatkhc.pl --bam $bamfile --ref $ref --analysis_name $analysis_name --bwbams";
    if ($Opt{variantonly}) {
        $bw_gatk_string .= " --variantonly";
    }
    system( $bw_gatk_string ) == 0
        or die "Couldn\'t call GATK HaplotypeCaller for analysis $analysis_name";

    return 0;

}

sub check_output {

    opendir JOBACCT, "job_accounting"
        or die "Couldn\'t open directory job_accounting to read sacct files!\n";
    foreach my $sacct_file (grep /sacct/, readdir JOBACCT) {
        my $sacct_fh = Open("job_accounting/$sacct_file");
        while (<$sacct_fh>) {
            next if (/^\s*AllocCPUS/); # header line
            next if (/\-{10}/); # header line
            chomp;
            my @fields = split /\s+/, $_; 
        }
        close $sacct_fh;
    }
    closedir JOBACCT;
}

__END__

=head1 OPTIONS

=over 4

=item B<--fastq_fof>

A file containing paths to local fastq files to be aligned with bwa mem.  For paired end reads, two tab-delimited fastq files should be included per line, and for single-end reads, one per line.

=item B<--tags_file>

A file containing read tags (in format "ID:140612_OPTIMUS_C46BJANXX.7\tLB:6935221\tPL:illumina\tSM:AP24 PU:140612_OPTIMUS_C46BJANXX") for each fastq file or file pair to be aligned with bwa mem.

=item B<--ucsc_ref>

A reference name to be found in the UCSC directory on biowulf2.

=item B<--analysis_name>

Optional name for the analysis subdirectory (of /data/$user) on biowulf2.

=item B<--variantonly>

Option that will be passed on to the script bw2_call_gatkhc.pl, and will request GATK to produce only a smaller VCF file containing variant sites only, rather than a full .g.vcf file.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software is a United States Government Work. Anyone may use the software on a worldwide and royalty-free
 basis for any purpose and anyone may reproduce and prepare derivative works without restriction. Although all
 reasonable efforts have been taken to ensure the accuracy and reliability of the software, the National Human Genome
 Research Institute (NHGRI), National Institutes of Health (NIH) and the U.S. Government do not and cannot warrant the
 performance or any results that may be obtained by using this software. NHGRI, NIH and the U.S. Government disclaim
 all warranties as to performance, merchantability or fitness for any particular purpose.  No indemnification is
 intended or provided by the US government.

In any work produced using or in any product derived from this material, the authors would appreciate being
 acknowledged per scientific custom as the source of the software.

=cut

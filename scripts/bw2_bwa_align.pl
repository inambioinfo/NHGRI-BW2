#!/usr/bin/perl -w

# $Id: bw2_bwa_align.pl 7656 2016-03-17 00:36:18Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use NHGRI::BW2::Session;
use NHGRI::BW2::Analysis;
use NHGRI::BW2::SwarmJob;
use GTB::File qw(Open);

our %Opt;

our $BWAVERSION = '0.7.12';
our $SAMTOOLSVERSION = '0.1.19';

=head1 NAME

bw2_bwa_align.pl - Submit a BWA mem alignment to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_bwa_align.pl --fastq_fof file_of_fastq_file_pairs --ref ref_to_align_to [options]

For complete documentation, run C<bw2_bwa_align.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of pairs of fastq files (or single fastq files for single end reads), and copies them to the biowulf2 cluster, aligning them with bwa mem.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $fastq_fof = $Opt{fastq_fof};
my $analysis_name = $Opt{analysis_name};

my ($ra_job_fastqs, $ra_input_files) = read_fastq_files($fastq_fof);

my $ra_tags; # read from file specified by "--tags_file" option
my $tags_file = $Opt{'tags_file'};
if ($tags_file) {
    $ra_tags = read_tags_file($tags_file);
}

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $analysis_name,
                      -analysis_basename => 'bwamem',
                      -input_files => $ra_input_files);

my $analysis_dir = $analysis_obj->create_analysis_directory();
unless($Opt{bwfastqs}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_bwa_swarm_string($ra_job_fastqs, $ra_tags, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/bwamemswarm.cmd");

$swarm_job->modules(["bwa/$BWAVERSION", "samtools/$SAMTOOLSVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
$swarm_job->threads_per_process(4);
$swarm_job->gb_per_process(24);
$swarm_job->walltime_per_process(12);
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

unless($Opt{skip_retrieve}) {
    $analysis_obj->retrieve_output_files(-local_dir => 'bwamem_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'bwamem_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/bwamem.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( ref => 'hg19'
           );
    GetOptions(\%Opt, qw( fastq_fof=s ref=s skip_retrieve tags_file=s bwfastqs
                analysis_name=s manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_bwa_align.pl, ", q$Revision: 7656 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --fastq_fof!") if (!$Opt{fastq_fof});
}

sub read_fastq_files {
    my $fileoffiles = shift;

    my @fastq_files = ();
    my @job_fastqs = ();
    my $fh = Open($fileoffiles);
    while (<$fh>) {
        if (/^(\S+)\t(\S+)$/) {
            my ($read1file, $read2file) = ($1, $2);
            push @fastq_files, $read1file;
            push @fastq_files, $read2file;
            push @job_fastqs, [$read1file, $read2file];
        }
        elsif (/^(\S+)$/) {
            my $singlereadfile = $1;
            push @fastq_files, $singlereadfile;
            push @job_fastqs, [$singlereadfile];

        }
    }
    close $fh;

    return ([@job_fastqs], [@fastq_files]);
}

sub read_tags_file {
    my $tag_file = shift;

    my @tag_strings = ();
    my $fh = Open($tag_file);
    while (<$fh>) {
        chomp;
        s/\s/\\t/g;
        push @tag_strings, $_;
    }
    close $fh;

    return [@tag_strings];
}

sub create_bwa_swarm_string {
    my $ra_job_fastqs = shift;
    my $ra_tags = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my $output_directory = $analysis_obj->output_dir();

    my $swarm_string = '';
    my $ref_string = find_biowulf_reference($analysis_obj->session()); 
    my @tag_strings = ($ra_tags) ? @{$ra_tags} : ();
    foreach my $ra_fastqs (@{$ra_job_fastqs}) {
        my @remote_fastq_files = @{$ra_fastqs};

        if (!$Opt{bwfastqs}) {
            @remote_fastq_files = map {$_=~ s:.*/:$input_directory/:; $_} @remote_fastq_files;
        }
        my $fastq_string = (@remote_fastq_files==2) ? join ' ', @remote_fastq_files : $remote_fastq_files[0];
        my $bam_name = $remote_fastq_files[0];
        $bam_name =~ s:.*/::;
        $bam_name =~ s/\.[01]\.f(ast){0,1}q(\.gz){0,1}//;

        my $tag_string = '';
        if (@tag_strings) {
            my $next_string = shift @tag_strings;
            if ($next_string !~ /PL:(\S+)/) { # add platform
                $next_string .= "\\tPL:illumina";
            }
            $tag_string = " -R \'\@RG\\t$next_string\'";
        }

        $swarm_string .= "bwa mem -t \$SLURM_CPUS_PER_TASK -M$tag_string $ref_string $fastq_string | samtools view -uS - | samtools sort -m 2000000000 - $output_directory/$bam_name\n";
        
    }

    return $swarm_string;
}

sub find_biowulf_reference {
    my $session_obj = shift;
    
    my $bw_ref_string = $session_obj->remote_command(-command => "ls /fdb/igenomes/Homo_sapiens/UCSC/$Opt{ref}/Sequence/BWAIndex/version0.6.0/genome.fa");
    if ($bw_ref_string =~ m!(/fdb/igenomes/\S+)!) {
        return $1;
    }
    else {
        die "Couldn\'t find reference files for $Opt{ref}!\n";
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--fastq_fof>

A file containing paths to local fastq files to be aligned with bwa mem.

=item B<--ref>

A reference name to be found in the UCSC directory on biowulf2.

=item B<--analysis_name>

Optional name for the analysis subdirectory (of /data/$user) on biowulf2.

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

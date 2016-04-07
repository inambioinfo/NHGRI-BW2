#!/usr/bin/perl -w

# $Id: bw2_call_gatkhc.pl 7648 2016-03-12 00:29:04Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use NHGRI::BW2::Session;
use NHGRI::BW2::Analysis;
use NHGRI::BW2::SwarmJob;
use GTB::File qw(Open);

our %Opt;

our $GATKVERSION = '3.5-0';

=head1 NAME

bw2_call_gatkhc.pl - Submit GATK "HaplotypeCaller" jobs to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_call_gatkhc.pl --bam bam_file_path --ref ucsc_reference_name

For complete documentation, run C<bw2_call_gatkhc.pl -man>

=head1 DESCRIPTION

This script submits a set of swarm jobs to run GATK's haplotype caller on each chromosome of a specified bam file (which resides either locally, or if the --bwbams option is specified, on the biowulf2 server). 

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $bamfile = $Opt{bam};
my $ref = $Opt{ref}; # UCSC reference name, for now
my $analysis_name = $Opt{analysis_name};

my $session_obj = NHGRI::BW2::Session->new();

# check for bam index, if file is local:

my $bamindex;
if (!$Opt{bwbams}) {
    if (-e "$bamfile.bai") {
        $bamindex = "$bamfile.bai";
    }
    else {
        $bamindex = $bamfile;
        $bamindex =~ s:.bam$:.bai:;
        if (!(-e $bamindex)) {
            die "Unable to find bam index for file $bamfile--please index.\n";
        }
    }
}

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $analysis_name,
                      -analysis_basename => 'gatkhc',
                      -input_files => [$bamfile, $bamindex]);

my $analysis_dir = $analysis_obj->create_analysis_directory();

unless($Opt{bwbams}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_haplotypecaller_swarm_string($bamfile, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/gatkhcswarm.cmd");

$swarm_job->modules(["GATK/$GATKVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
$swarm_job->gb_per_process(10);
$swarm_job->gres("lscratch:100");
my $walltime = $Opt{'walltime'};
$swarm_job->walltime_per_process($walltime);
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

if (!$Opt{skip_retrieve}) {
   $analysis_obj->retrieve_output_files(-local_dir => 'gatkhc_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'gatkhc_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/gatkhc.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( ref => 'hg19', walltime => 8
           );
    GetOptions(\%Opt, qw( bam=s bwbams skip_retrieve ref=s chroms=s variantonly
                          analysis_name=s walltime=s manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_call_gatkhc.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --bam!") if (!$Opt{bam});
}

sub create_haplotypecaller_swarm_string {
    my $bam_file = shift;
    my $analysis_obj = shift;

    if (!$Opt{bwbams}) {
        my $input_directory = $analysis_obj->input_dir();
        $bam_file =~ s:.*/::;
        $bam_file = "$input_directory/$bam_file";
    }

    my $output_directory = $analysis_obj->output_dir();

    my $reffasta = find_biowulf_reference($analysis_obj->session()); 
    my @all_chroms = ($Opt{chroms}) ? split /,/, $Opt{chroms} :
              get_biowulf_chromosomes($analysis_obj->session(), $reffasta);

    my $swarm_string = '';

    foreach my $chrom (@all_chroms) {
        my $output_file = $bam_file;
        my $file_ext = ($Opt{'variantonly'}) ? '.vcf.gz' : '.g.vcf.gz';
        $output_file =~ s/\.bam$/.$chrom$file_ext/;
        $output_file =~ s:.*/::;
        $output_file = "$output_directory/$output_file";

        # Note: if you change java memory or cpu allocation below, change
        # values requested for job too (via the swarm_job object calls to 
        # gb_per_process and threads_per_process in this script)
        
        $swarm_string .= "GATK -m 8g HaplotypeCaller -R $reffasta -L $chrom -I $bam_file -o $output_file";

        if (!$Opt{'variantonly'}) {
            $swarm_string .= ' --emitRefConfidence GVCF';
        }
        $swarm_string .= "\n";
    }

    return $swarm_string;
}

sub find_biowulf_reference {
    my $session_obj = shift;
    
    my $bw_ref_string = $session_obj->remote_command(-command => "ls /fdb/igenomes/*/UCSC/$Opt{ref}/Sequence/WholeGenomeFasta/genome.fa");
    if ($bw_ref_string =~ m!(/fdb/igenomes/\S+)!) {
        return $1;
    }
    else {
        die "Couldn\'t find reference files for $Opt{ref}!\n";
    }
}

sub get_biowulf_chromosomes {
    my $session_obj = shift;
    my $bw_ref_fasta = shift;
    
    my $bw_chrom_string = $session_obj->remote_command(-command => "awk \'{print \$1}\' $bw_ref_fasta.fai");
    my @chroms = split /\n/, $bw_chrom_string;

    if (!@chroms) {
        die "Couldn\'t find reference files for $Opt{ref}! Command returned string:\n$bw_chrom_string";
    }

    return @chroms;
}

__END__

=head1 OPTIONS

=over 4

=item B<--bam_paths>

A file containing paths to local or remote bam files for which indel targets will be calculated by GATK.

=item B<--ref>

A reference name to be found in the UCSC directory on biowulf2.

=item B<--walltime>

Value to pass to slurm requesting time for each job.  As an integer, it will be interpreted as number of hours.  As a string, it is formatted as HH:MI:SS (default 8 hours).

=item B<--chroms>

Comma-delimited list of chromosomes to process (default is to analyze all chromosomes in the provided reference).

=item B<--variantonly>

Write a VCF with only variant locations (by default, HaplotypeCaller is run with the --emitRefConfidence GVCF option, so will write GVCF format).

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

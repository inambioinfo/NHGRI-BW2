#!/usr/bin/perl -w

# $Id: bw2_mark_duplicates.pl 7648 2016-03-12 00:29:04Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use NHGRI::BW2::Session;
use NHGRI::BW2::Analysis;
use NHGRI::BW2::SwarmJob;
use GTB::File qw(Open);

our %Opt;

our $PICARDVERSION = '1.139';

=head1 NAME

bw2_mark_duplicates.pl - Submit a Picard MarkDuplicates job to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_mark_duplicates.pl --bam_paths file_of_bam_locations

For complete documentation, run C<bw2_mark_duplicates.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of bam file locations (either locally, or if the --bwbams option is specified, on the biowulf2 server), and runs Picard's MarkDuplicates command on them. 

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $bamfile_fof = $Opt{bam_paths};
my $analysis_name = $Opt{analysis_name};

my $ra_bamfiles = read_bam_paths($bamfile_fof);

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $analysis_name,
                      -analysis_basename => 'markdups',
                      -input_files => $ra_bamfiles);

my $analysis_dir = $analysis_obj->create_analysis_directory();

unless($Opt{bwbams}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_markdups_swarm_string($ra_bamfiles, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/markdupswarm.cmd");

$swarm_job->modules(["picard/$PICARDVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
$swarm_job->threads_per_process(5);
$swarm_job->gb_per_process(12);
$swarm_job->gres("lscratch:32");
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

if (!$Opt{skip_retrieve}) {
   $analysis_obj->retrieve_output_files(-local_dir => 'markdup_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'markdup_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/markdups.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( 
           );
    GetOptions(\%Opt, qw( bam_paths=s bwbams skip_retrieve index analysis_name=s
                          manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_mark_duplicates.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --bam_paths!") if (!$Opt{bam_paths});
}

sub read_bam_paths {
    my $fileoffiles = shift;

    my @bam_files = ();
    my $fh = Open($fileoffiles);
    while (<$fh>) {
        # bampath is only column
        if (/^(\S+)$/) {
            my $bamfile = $1;
            push @bam_files, $bamfile;
        }
        else {
            die "File of BAM files should contain one column per row with a BAM file location!\n";
        }
    }
    close $fh;

    return [@bam_files];
}

sub create_markdups_swarm_string {
    my $ra_bam_files = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my $output_directory = $analysis_obj->output_dir();

    my $swarm_string = '';
    foreach my $bamfile (@{$ra_bam_files}) {
        my $input_file = $bamfile;
        if (!$Opt{bwbams}) {
            $input_file =~ s:.*/::;
            $input_file = "$input_directory/$input_file";
        }
        my $output_file = $bamfile;
        $output_file =~ s/\.bam$/.markdup.bam/;
        $output_file =~ s:.*/::;
        $output_file = "$output_directory/$output_file";
        my $metrics_file = $output_file;
        $metrics_file =~ s/\.bam$/.metrics/;

        # Note: if you change java memory or cpu allocation below, change
        # values requested for job too (via the swarm_job object calls to 
        # gb_per_process and threads_per_process in this script)

        my $index_string = ($Opt{index}) ? ' CREATE_INDEX=true' : '';
        
        $swarm_string .= "java -Xmx10g -XX:ParallelGCThreads=5 -jar \$PICARDJARPATH/picard.jar MarkDuplicates INPUT=$input_file OUTPUT=$output_file METRICS_FILE=$metrics_file ASSUME_SORTED=true TMP_DIR=/lscratch/\$SLURM_JOBID$index_string\n";
        
    }

    return $swarm_string;
}

__END__

=head1 OPTIONS

=over 4

=item B<--bam_paths>

A file containing paths to local or remote bam files within which duplicates will be marked, one path per line.

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

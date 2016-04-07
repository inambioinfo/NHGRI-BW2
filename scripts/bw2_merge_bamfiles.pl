#!/usr/bin/perl -w

# $Id: bw2_merge_bamfiles.pl 7648 2016-03-12 00:29:04Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use NHGRI::BW2::Session;
use NHGRI::BW2::Analysis;
use NHGRI::BW2::SbatchJob;
use GTB::File qw(Open);

our %Opt;

our $PICARDVERSION = '1.139';

=head1 NAME

bw2_merge_bamfiles.pl - Submit a Picard MergeSamFiles job to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_merge_bamfiles.pl --bam_paths file_of_bam_locations --merged_bam name_of_output_file

For complete documentation, run C<bw2_merge_bamfiles.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of bam file locations (either locally, or if the --bwbams option is specified, on the biowulf2 server), and runs Picard's MergeSamFiles command on them. 

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $bamfile_fof = $Opt{bam_paths};
my $merged_bam = $Opt{merged_bam};
my $analysis_name = $Opt{analysis_name};

my $ra_bamfiles = read_bam_paths($bamfile_fof);

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $analysis_name,
                      -analysis_basename => 'mergebams',
                      -input_files => $ra_bamfiles);

my $analysis_dir = $analysis_obj->create_analysis_directory();

unless($Opt{bwbams}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $sbatch_job = NHGRI::BW2::SbatchJob->new();
$sbatch_job->session($session_obj);
my $sbatch_string = create_mergebams_sbatch_string($ra_bamfiles, $merged_bam, $analysis_obj);
$sbatch_job->sbatch_command_string($sbatch_string);

my $remote_sbatch_path = $analysis_obj->scripts_dir();
$sbatch_job->modules(["picard/$PICARDVERSION"]);
$sbatch_job->logfile_dir($analysis_obj->logfile_dir());
$sbatch_job->walltime_per_process(30);
$sbatch_job->gb_per_process(6);

$sbatch_job->write_and_transfer_sbatch_file(
                      -remotepath => "$remote_sbatch_path/mergebamssbatch.cmd");

$sbatch_job->submit_sbatch_command();
$sbatch_job->wait_for_sbatch();

if (!$Opt{skip_retrieve}) {
   $analysis_obj->retrieve_output_files(-local_dir => 'mergebams_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'mergebams_logfiles');

my $jobid = $sbatch_job->jobid();
my $sacct_info = $sbatch_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/mergebams.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( 
           );
    GetOptions(\%Opt, qw( bam_paths=s merged_bam=s bwbams skip_retrieve index analysis_name=s
                          manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_merge_bamfiles.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --bam_paths!") if (!$Opt{bam_paths});
    pod2usage("Missing required field --merged_bam (name for merged file)!") if (!$Opt{merged_bam});
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

sub create_mergebams_sbatch_string {
    my $ra_bam_files = shift;
    my $output_file = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my $output_directory = $analysis_obj->output_dir();

    $output_file = "$output_directory/$output_file";

    my $sbatch_string = '';
    my @input_bams = ();
    if (!$Opt{bwbams}) {
        @input_bams = map { s:.*::; "$input_directory/$_" } @{$ra_bam_files};
    }
    else {
        @input_bams = @{$ra_bam_files};
    }
    @input_bams = map { "INPUT=$_" } @input_bams;

    # Note: if you change java memory or cpu allocation below, change
    # values requested for job too (via the swarm_job object calls to 
    # gb_per_process and threads_per_process in this script)

    my $index_string = ($Opt{index}) ? ' CREATE_INDEX=true' : '';
    my $sort = 'coordinate'; # until we implement options
   
    my $input_string = join ' ', @input_bams; 
    $sbatch_string .= "java -Xmx5g -jar \$PICARDJARPATH/picard.jar MergeSamFiles $input_string OUTPUT=$output_file SORT_ORDER=coordinate ASSUME_SORTED=true$index_string\n";
        

    return $sbatch_string;
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

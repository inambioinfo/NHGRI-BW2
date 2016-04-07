#!/usr/bin/perl -w

# $Id: bw2_bam2fastq.pl 7648 2016-03-12 00:29:04Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use NHGRI::BW2::Session;
use NHGRI::BW2::Analysis;
use NHGRI::BW2::SwarmJob;
use GTB::File qw(Open);

our %Opt;

our $SAMTOOLSVERSION = '0.1.19';

=head1 NAME

bw2_bam2fastq.pl - Submit a swarm of bam2fastq jobs to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_bam2fastq.pl --bamfile_fof file_of_bam_file_paths [options]

For complete documentation, run C<bw2_bam2fastq.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of local bam file paths (or paths on the biowulf2 server if the --bwbams option is specified), and copies them to the biowulf2 cluster, converting them to fastq files with bam2fastq.pl.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $bamfile_fof = $Opt{bamfile_fof};
my $analysis_name = $Opt{analysis_name};

my $ra_input_files = read_bamfiles($bamfile_fof);

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $analysis_name,
                      -analysis_basename => 'bam2fastq',
                      -input_files => $ra_input_files);

my $analysis_dir = $analysis_obj->create_analysis_directory();
unless($Opt{bwbams}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_bam2fastq_swarm_string($ra_input_files, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/bam2fastqswarm.cmd");

$swarm_job->modules(["samtools/$SAMTOOLSVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
$swarm_job->threads_per_process(2);
$swarm_job->gb_per_process(4);
$swarm_job->gres("lscratch:20");
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

unless($Opt{skip_retrieve}) {
    $analysis_obj->retrieve_output_files(-local_dir => 'bam2fastq_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'bam2fastq_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( 
           );
    GetOptions(\%Opt, qw( bamfile_fof=s skip_retrieve bwbams
                analysis_name=s manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_bam2fastq.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --bamfile_fof!") if (!$Opt{bamfile_fof});
}

sub read_bamfiles {
    my $fileoffiles = shift;

    my @bamfiles = ();
    my $fh = Open($fileoffiles);
    while (<$fh>) {
        if (/^(\S+)$/) {
            my $bamfile = $1;
            push @bamfiles, $bamfile;
        }
    }
    close $fh;

    return [@bamfiles];
}

sub create_bam2fastq_swarm_string {
    my $ra_bamfiles = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my $output_directory = $analysis_obj->output_dir();

    my $swarm_string = '';
    foreach my $bamfile (@{$ra_bamfiles}) {
        my $remote_bamfile = $bamfile;

        if (!$Opt{bwbams}) {
            $remote_bamfile =~ s:.*/:$input_directory/:;
        }
        my $bam_base = $bamfile;
        $bam_base =~ s:.*/::;
        $bam_base =~ s:\.bam::;

        $swarm_string .= "export TMPDIR=/lscratch/\$SLURM_JOBID; export PERL5LIB=/home/nhansen/INSTALL/lib/site_perl/5.12.1:\$PERL5LIB; /home/nhansen/INSTALL/bin/bam2fastq.pl $remote_bamfile --prefix $output_directory/$bam_base\n";
        
    }

    return $swarm_string;
}

__END__

=head1 OPTIONS

=over 4

=item B<--bamfile_fof>

A file containing paths to local bam files to be converted to fastq files (or, if --bwbams option is specified, paths to bam files on the biowulf2 server).

=item B<--bwbams>

Flag option to specify that bamfiles to be converted already reside on the biowulf2 server.

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

#!/usr/bin/perl -w

# $Id: bw2_star_align.pl 7648 2016-03-12 00:29:04Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use NHGRI::BW2::Session;
use NHGRI::BW2::Analysis;
use NHGRI::BW2::SwarmJob;
use GTB::File qw(Open);

our %Opt;

our $STARVERSION = '2.5.1b';
our $SAMTOOLSVERSION = '1.2';
our $BW_REFERENCE_PATH = "/data/$ENV{USER}/references/STAR";

=head1 NAME

bw2_star_align.pl - Submit a RNA STAR alignment swarm to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_star_align.pl --fastq_fof file_of_fastq_file_pairs --genome_dir preindexed_star_genome_directory [options]

For complete documentation, run C<bw2_star_align.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of pairs of fastq files (or single fastq files for single end reads), and copies them to the biowulf2 cluster, aligning them with star mem.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $fastq_fof = $Opt{fastq_fof};
my $genome_dir = $Opt{genome_dir};

my ($ra_job_fastqs, $ra_input_files) = read_fastq_files($fastq_fof);

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $Opt{analysis_name},
                      -analysis_basename => 'staralign',
                      -input_files => $ra_input_files);

my $analysis_dir = $analysis_obj->create_analysis_directory();
unless ($Opt{bwfastqs}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_star_swarm_string($ra_job_fastqs, $genome_dir, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/staralignswarm.cmd");

$swarm_job->modules(["STAR/$STARVERSION", "samtools/$SAMTOOLSVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
$swarm_job->threads_per_process($Opt{threads_per_process});
$swarm_job->walltime_per_process($Opt{walltime});
$swarm_job->gb_per_process($Opt{gb_per_process});
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

unless($Opt{skip_retrieve}) {
    $analysis_obj->retrieve_output_files(-local_dir => 'staralign_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'staralign_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/staralign.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( threads_per_process => 12, gb_per_process => 36, walltime => 4, readlength => 101
           );
    GetOptions(\%Opt, qw( fastq_fof=s genome_dir=s skip_retrieve threads_per_process=i 
                gb_per_process=i walltime=s readlength=i bwfastqs twopass 
                analysis_name=s manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_star_align.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --fastq_fof!") if (!$Opt{fastq_fof});
    pod2usage("Missing required field --genome_dir!") if (!$Opt{genome_dir});
}

sub read_fastq_files {
    my $fileoffiles = shift;

    my @fastq_files = ();
    my @job_fastqs = ();
    my $fh = Open($fileoffiles);
    while (<$fh>) {
        next if (/^#/); # comment lines allowed
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

sub create_star_swarm_string {
    my $ra_job_fastqs = shift;
    my $genome_dir = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my $output_directory = $analysis_obj->output_dir();

    my $swarm_string = '';
    foreach my $ra_fastqs (@{$ra_job_fastqs}) {
        my @remote_fastq_files = @{$ra_fastqs};

        unless ($Opt{bwfastqs}) {
            @remote_fastq_files = map {$_=~ s:.*/:$input_directory/:; $_} @remote_fastq_files;
        }
        my $fastq_string = (@remote_fastq_files==2) ? join ' ', @remote_fastq_files : $remote_fastq_files[0];
        my $bam_name = $remote_fastq_files[0];
        $bam_name =~ s:.*/::;
        $bam_name =~ s/\.f(ast){0,1}q(\.gz){0,1}//;

        my $cat_prog = ($remote_fastq_files[0] =~ /\.gz$/) ? 'zcat' : 'cat';

        my $overhang = $Opt{readlength} - 1;

        my $twopassopt = ($Opt{twopass}) ? ' --twopassMode Basic' : '';
           
        $swarm_string .= "STAR --runThreadN \$SLURM_CPUS_PER_TASK --genomeDir $BW_REFERENCE_PATH/$genome_dir --readFilesIn $fastq_string --readFilesCommand $cat_prog --sjdbOverhang $overhang$twopassopt --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $output_directory/$bam_name.\n";
    }

    return $swarm_string;
}

__END__

=head1 OPTIONS

=over 4

=item B<--fastq_fof>

A file containing paths to local fastq files to be aligned with star mem, or, if --bwfastqs option has been specified, paths to fastq files on biowulf2.

=item B<--genome_dir>

The path to an indexed STAR genome directory on biowulf2, formatted for the correct readlength (use bw2_index_star_genome.pl).

=item B<--readlength>

The readlength of reads to be aligned (in a better world, this would be determined by examining the files).

=item B<--twopass>

This flag causes STAR to be run with the option "--twopassMode Basic", which will run STAR twice, extracting discovered junctions after the first run and inserting them into the genome index before re-mapping in a second alignment step.

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

#!/usr/bin/perl -w

# $Id: bw2_index_star_genome.pl 7648 2016-03-12 00:29:04Z nhansen $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use NHGRI::BW2::Session;
use NHGRI::BW2::Analysis;
use NHGRI::BW2::SbatchJob;
use GTB::File qw(Open);

our %Opt;

our $STARVERSION = '2.5.1b';
our $BW_REFERENCE_PATH = "/data/$ENV{USER}/references/STAR";

=head1 NAME

bw2_index_star_genomes.pl - Submit a STAR job to the biowulf2 cluster to call with the --genomeGenerate option to index a genome with optional annotation.

=head1 SYNOPSIS

  bw2_index_star_genomes.pl --ref_fasta fasta file of reference contigs --gtf gtf_file_of_annotations --readlength length_of_reads_to_align --genome_dir subdirectory_name_on_biowulf

For complete documentation, run C<bw2_index_star_genomes.pl -man>

=head1 DESCRIPTION

This script takes a reference fasta file and an optional GTF-formatted file of gene annotations, and calls STAR on the biowulf2 cluster to format an index using its "genomeGenerate" command.  Users can specify the read length of reads to be aligned with --readlength, and the name of the subdirectory in which to store the index with --genome_dir.  If the reference fasta already resides on biowulf2, it can be specified with the --bw_ucsc_reference or --bw_reference_path options, and if the gtf file already resides on biowulf, its path can be specified with the --bw_gtf option.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $ref_fasta = $Opt{ref_fasta};
my $bw_ucsc_reference = $Opt{bw_ucsc_reference};
my $bw_reference_path = $Opt{bw_reference_path};
my $gtf = $Opt{gtf};
my $bw_gtf = $Opt{bw_gtf};
my $readlength = $Opt{readlength};
my $genome_dir = $Opt{genome_dir};

my $session_obj = NHGRI::BW2::Session->new();

my $ra_input_files = [];

if ($ref_fasta) {
    push @{$ra_input_files}, $ref_fasta;
}

if ($gtf) {
    push @{$ra_input_files}, $gtf;
}

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_basename => 'stargenerate',
                      -input_files => $ra_input_files);

my $analysis_dir = $analysis_obj->create_analysis_directory();

if (@{$ra_input_files}) {
    $analysis_obj->transfer_input_files();
}

# Determine biowulf2 reference fasta path from options:
my ($ref_path, $gtf_path);
if ($ref_fasta) { # local fasta transferred to this analysis's biowulf2 "input_directory"
    my $input_directory = $analysis_obj->input_dir();
    $ref_path = $ref_fasta;
    $ref_path =~ s:.*/::;
    $ref_path = "$input_directory/$ref_path";
}
elsif ($bw_reference_path) {
    $ref_path = $bw_reference_path;
}
elsif ($bw_ucsc_reference) { # find UCSC reference on biowulf2
    $ref_path = find_biowulf_ucsc_reference_fasta($bw_ucsc_reference, $session_obj);
}

# Now determine biowulf2 gtf path (if there is one)

if ($gtf) { # local gtf transferred to this analysis's biowulf2 "input_directory"
    my $input_directory = $analysis_obj->input_dir();
    $gtf_path = $gtf;
    $gtf_path =~ s:.*/::;
    $gtf_path = "$input_directory/$gtf_path";
}
elsif ($bw_gtf) {
    $gtf_path = $bw_gtf;
}

# job submission and checking here
my $sbatch_job = NHGRI::BW2::SbatchJob->new();
$sbatch_job->session($session_obj);
my $sbatch_string = create_stargenerate_sbatch_string($ref_path, $gtf_path, $genome_dir, $analysis_obj);
$sbatch_job->sbatch_command_string($sbatch_string);

my $remote_sbatch_path = $analysis_obj->scripts_dir();
$sbatch_job->modules(["STAR/$STARVERSION"]);
$sbatch_job->logfile_dir($analysis_obj->logfile_dir());
$sbatch_job->walltime_per_process($Opt{walltime});
$sbatch_job->gb_per_process($Opt{gb_per_process});
$sbatch_job->threads_per_process($Opt{threads_per_process});

$sbatch_job->write_and_transfer_sbatch_file(
                      -remotepath => "$remote_sbatch_path/stargeneratesbatch.cmd");

$session_obj->make_remote_directory(-dirname => "$BW_REFERENCE_PATH/$genome_dir");

$sbatch_job->submit_sbatch_command();
$sbatch_job->wait_for_sbatch();

if (!$Opt{skip_retrieve}) {
   $analysis_obj->retrieve_output_files(-local_dir => 'stargenerate_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'stargenerate_logfiles');

my $jobid = $sbatch_job->jobid();
my $sacct_info = $sbatch_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/stargenerate.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( gb_per_process => 35, walltime => 4, threads_per_process => 4
           );
    GetOptions(\%Opt, qw( ref_fasta=s bw_reference_path=s bw_ucsc_reference=s gtf=s bw_gtf=s genome_dir=s
                          readlength=i gb_per_process=s walltime=s threads_per_process=i manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_index_star_genomes.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Must specify a reference path with --ref_fasta, --bw_reference_path, or --bw_ucsc_reference!") if (!$Opt{ref_fasta} && !$Opt{bw_reference_path} && !$Opt{bw_ucsc_reference});
    pod2usage("Missing required parameter --genome_dir (name for STAR reference subdirectory to store indexes in)!") if (!$Opt{genome_dir});
}

sub create_stargenerate_sbatch_string {
    my $reference_path = shift;
    my $gtf_path = shift;
    my $genome_dir = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my $readlengthminusone = $Opt{readlength} - 1;

    my $sbatch_string .= "STAR --runThreadN \$SLURM_CPUS_PER_TASK --genomeDir $BW_REFERENCE_PATH/$genome_dir --runMode genomeGenerate --genomeFastaFiles $reference_path --sjdbGTFfile $gtf_path --sjdbOverhang $readlengthminusone";
        

    return $sbatch_string;
}

sub find_biowulf_ucsc_reference_fasta {
    my $ucsc_ref_name = shift;
    my $session_obj = shift;
    
    my $bw_ref_string = $session_obj->remote_command(-command => "ls /fdb/igenomes/Homo_sapiens/UCSC/$ucsc_ref_name/Sequence/WholeGenomeFasta/genome.fa");
    if ($bw_ref_string =~ m!(/fdb/igenomes/\S+)!) {
        return $1;
    }
    else {
        die "Couldn\'t find reference files for $ucsc_ref_name!\n";
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--bam_paths>

A file containing paths to local or remote bam files within which duplicates will be marked, one path per line.

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

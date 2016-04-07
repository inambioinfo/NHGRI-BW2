#!/usr/bin/perl -w

# $Id: bw2_gatk_asereadcounter.pl 7648 2016-03-12 00:29:04Z nhansen $

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

bw2_gatk_asereadcounter.pl - Submit a set of GATK "ASEReadCounter" jobs to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_gatk_asereadcounter.pl --bam_fof file_of_bam_file_paths --sites_fof file_of_vcf_files --ucsc_ref ucsc_reference_name

For complete documentation, run C<bw2_gatk_asereadcounter.pl -man>

=head1 DESCRIPTION

This script submits a set of swarm jobs to run GATK's ASEReadCounter command on each of a specified list of bam files (which reside either locally, or if the --bwbams option is specified, on the biowulf2 server), interrogating sites in each of a set of specified VCF files.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $bamfof = $Opt{bam_fof};
my $sitesfof = $Opt{sites_fof};
my $ucsc_ref = $Opt{ucsc_ref}; # UCSC ucsc_reference name, for now

my @input_files = ();
my $ra_bamfiles = read_file_of_files($bamfof);
my $ra_sitesfiles = read_file_of_files($sitesfof);

if (!$Opt{bwbams}) {
    push @input_files, @{$ra_bamfiles};
}

if (!$Opt{bwsites}) {
    push @input_files, @{$ra_sitesfiles};
    my @index_files = map { "$_.tbi" } @{$ra_sitesfiles};
    push @input_files, @index_files;
}

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_basename => 'asereadcounter');

my $analysis_dir = $analysis_obj->create_analysis_directory();

if (@input_files) { # transfer input files if necessary
    $analysis_obj->input_files([@input_files]);
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_asereadcounter_swarm_string($ra_bamfiles, $ra_sitesfiles, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/asereadcounterswarm.cmd");

$swarm_job->modules(["GATK/$GATKVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
my $gb_per_process = $Opt{'gb_per_process'};
$swarm_job->gb_per_process($gb_per_process);
my $walltime = $Opt{'walltime'};
$swarm_job->walltime_per_process($walltime);
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

if (!$Opt{skip_retrieve}) {
   $analysis_obj->retrieve_output_files(-local_dir => 'asereadcounter_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'asereadcounter_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/asereadcounter.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( ref => 'hg19', walltime => 8, gb_per_process => 6
           );
    GetOptions(\%Opt, qw( bam_fof=s bwbams sites_fof=s bwsites skip_retrieve ucsc_ref=s
                          walltime=s gb_per_process=i manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_gatk_asereadcounter.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --bam_fof!") if (!$Opt{bam_fof});
    pod2usage("Missing required field --sites_fof!") if (!$Opt{sites_fof});
}

sub create_asereadcounter_swarm_string {
    my $ra_bams = shift;
    my $ra_sites = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my @bwbampaths = ($Opt{bwbams}) ? @{$ra_bams} :  map { s:.*/::; "$input_directory/$_" } @{$ra_bams};
    my @bwsitespaths = ($Opt{bwsites}) ? @{$ra_sites} :  map { s:.*/::; "$input_directory/$_" } @{$ra_sites};

    my $output_directory = $analysis_obj->output_dir();

    my $reffasta = find_biowulf_reference($analysis_obj->session()); 

    my $swarm_string = '';

    foreach my $bamfile (@bwbampaths) {
        foreach my $sites_file (@bwsitespaths) {
            my $output_file1 = $bamfile;
            $output_file1 =~ s/\.bam$//;
            $output_file1 =~ s:.*/::;
            $output_file1 = "$output_directory/$output_file1";

            my $output_file2 = $sites_file;
            $output_file2 =~ s/\.vcf(\.gz){0,1}$//;
            $output_file2 =~ s:.*/::;

            my $output_file = "$output_file1.$output_file2.csv";

            # Note: if you change java memory or cpu allocation below, change
            # values requested for job too (via the swarm_job object calls to 
            # gb_per_process and threads_per_process in this script)
            
            $swarm_string .= "GATK -m $Opt{gb_per_process}g ASEReadCounter -R $reffasta -sites $sites_file -I $bamfile -o $output_file -U ALLOW_N_CIGAR_READS\n";
        }
    }

    return $swarm_string;
}

sub read_file_of_files {
    my $file = shift;

    my @return_list = ();
    my $fh = Open($file);
    while (<$fh>) {
        chomp;
        next if (/^#/);
        if (/^(\S+)$/) {
            push @return_list, $_;
        }
        else {
            print STDERR "Skipping line $_--should be only one field per line!\n";
        }
    }

    return [@return_list];
}

sub find_biowulf_reference {
    my $session_obj = shift;
    
    my $bw_ref_string = $session_obj->remote_command(-command => "ls /fdb/igenomes/*/UCSC/$Opt{ucsc_ref}/Sequence/WholeGenomeFasta/genome.fa");
    if ($bw_ref_string =~ m!(/fdb/igenomes/\S+)!) {
        return $1;
    }
    else {
        die "Couldn\'t find reference files for $Opt{ucsc_ref}!\n";
    }
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

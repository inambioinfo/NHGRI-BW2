#!/usr/bin/perl -w

# $Id: bw2_create_indeltargets.pl 7648 2016-03-12 00:29:04Z nhansen $

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

bw2_create_indeltargets.pl - Submit a GATK "RealignerTargetCreator" job to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_create_indeltargets.pl --bam_paths file_of_bam_locations

For complete documentation, run C<bw2_create_indeltargets.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of bam file locations (either locally, or if the --bwbams option is specified, on the biowulf2 server), and runs GATK's RealignerTargetCreator command on them. 

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $bamfile_fof = $Opt{bam_paths};
my $ref = $Opt{ref};
my $analysis_name = $Opt{analysis_name};

my $ra_bamfiles = read_bam_paths($bamfile_fof);

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $analysis_name,
                      -analysis_basename => 'targetcreator',
                      -input_files => $ra_bamfiles);

my $analysis_dir = $analysis_obj->create_analysis_directory();

unless($Opt{bwbams}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_targetcreator_swarm_string($ra_bamfiles, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/targetcreatorswarm.cmd");

$swarm_job->modules(["GATK/$GATKVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
$swarm_job->gb_per_process(8);
$swarm_job->gres("lscratch:5");
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

if (!$Opt{skip_retrieve}) {
   $analysis_obj->retrieve_output_files(-local_dir => 'targetcreator_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'targetcreator_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/targetcreator.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( ref => 'hg19', known_resource => 'gold'
           );
    GetOptions(\%Opt, qw( bam_paths=s bwbams skip_retrieve ref=s known_resource=s analysis_name=s
                          manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_create_indeltargets.pl, ", q$Revision: 7648 $, "\n"; }
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

sub create_targetcreator_swarm_string {
    my $ra_bam_files = shift;
    my $analysis_obj = shift;

    my $input_directory = $analysis_obj->input_dir();
    my $output_directory = $analysis_obj->output_dir();

    my $reffasta = find_biowulf_reference($analysis_obj->session()); 
    my $known_indel_vcf = find_known_indel_file($analysis_obj->session()); 

    my $swarm_string = '';
    foreach my $bamfile (@{$ra_bam_files}) {
        my $input_file = $bamfile;
        if (!$Opt{bwbams}) {
            $input_file =~ s:.*/::;
            $input_file = "$input_directory/$input_file";
        }
        my $output_file = $bamfile;
        $output_file =~ s/\.bam$/.intervals/;
        $output_file =~ s:.*/::;
        $output_file = "$output_directory/$output_file";

        # Note: if you change java memory or cpu allocation below, change
        # values requested for job too (via the swarm_job object calls to 
        # gb_per_process and threads_per_process in this script)
        
        $swarm_string .= "GATK -m 7g RealignerTargetCreator -R $reffasta -I $input_file -o $output_file -known $known_indel_vcf\n";
        
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

sub find_known_indel_file {
    my $session_obj = shift;
    
    my $bw_ref_string = $session_obj->remote_command(-command => "ls /fdb/GATK_resource_bundle/$Opt{ref}/*.indels*.vcf*");
    if (($Opt{known_resource} =~ /gold/i) && ($bw_ref_string =~ m!(/fdb/GATK_resource_bundle/$Opt{ref}/\S+\_gold\_standard\.indels.*vcf(\.gz){0,1})!)) {
        print STDERR "Using known indel file $1\n";
        return $1;
    }
    elsif (($Opt{known_resource} =~ /1(000|k)g/i) && ($bw_ref_string =~ m!(/fdb/GATK_resource_bundle/$Opt{ref}/1000G.*\.indels.*vcf(\.gz){0,1})!)) {
        print STDERR "Using known indel file $1\n";
        return $1;
    }
    else {
        die "Couldn\'t find known indel files for $Opt{ref}!\n";
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--bam_paths>

A file containing paths to local or remote bam files for which indel targets will be calculated by GATK.

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

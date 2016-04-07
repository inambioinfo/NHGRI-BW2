#!/usr/bin/perl -w

# $Id: bw2_add_readtags.pl 7648 2016-03-12 00:29:04Z nhansen $

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

bw2_add_readtags.pl - Submit a Picard AddOrReplaceReadGroups job to the biowulf2 cluster.

=head1 SYNOPSIS

  bw2_add_readtags.pl --bam_info file_of_bam_tags --pl platform [options]

For complete documentation, run C<bw2_add_readtags.pl -man>

=head1 DESCRIPTION

This script takes a file containing a list of bam files with five tab-delimited columns per line: path to BAM file (locally, or if --bwbams option was specified, on biowulf2), sample name, library name, flowcell, and id.  It then adds the tags "RGSM", "RGLB", "RGPU", and "RGID" to each BAM file using the Picard function "AddOrReplaceReadGroups".  In addition, it adds the "RGPL" tag, with default value "illumina", or optionally a different value specified by the --pl option.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $bamfile_fof = $Opt{bam_info};

my $ra_bam_hashes = read_bam_info($bamfile_fof);
my @bam_paths = ($Opt{'bwbams'}) ? () : map {$_->{bam}} @{$ra_bam_hashes};

my $session_obj = NHGRI::BW2::Session->new();

my $analysis_obj = NHGRI::BW2::Analysis->new(
                      -session => $session_obj,
                      -analysis_name => $Opt{analysis_name},
                      -analysis_basename => 'addreadtags',
                      -input_files => \@bam_paths);

my $analysis_dir = $analysis_obj->create_analysis_directory();

unless($Opt{bwbams}) {
    $analysis_obj->transfer_input_files();
}

# job submission and checking here
my $swarm_job = NHGRI::BW2::SwarmJob->new();
$swarm_job->session($session_obj);
my $swarm_string = create_readtag_swarm_string($ra_bam_hashes, $analysis_obj);
$swarm_job->swarm_command_string($swarm_string);
my $remote_swarm_path = $analysis_obj->scripts_dir();
$swarm_job->write_and_transfer_swarm_file(
                      -remotepath => "$remote_swarm_path/addtagsswarm.cmd");

$swarm_job->modules(["picard/$PICARDVERSION"]);
$swarm_job->logfile_dir($analysis_obj->logfile_dir());
$swarm_job->threads_per_process(5);
$swarm_job->gb_per_process(6);
$swarm_job->submit_swarm_command();
$swarm_job->wait_for_swarm();

if (!$Opt{skip_retrieve}) {
   $analysis_obj->retrieve_output_files(-local_dir => 'addreadtag_output');
}

$analysis_obj->retrieve_log_files(-local_dir => 'addreadtag_logfiles');

my $jobid = $swarm_job->jobid();
my $sacct_info = $swarm_job->all_sacct_info();
mkdir "job_accounting";
open SACCT, ">job_accounting/addreadtags.sacct.$jobid";
print SACCT "$sacct_info";
close SACCT;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( platform => 'illumina'
           );
    GetOptions(\%Opt, qw( bam_info=s bwbams skip_retrieve platform=s index
                          analysis_name=s manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "bw2_add_readtags.pl, ", q$Revision: 7648 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    # Add criteria to test option arguments below
    pod2usage("Missing required field --bam_info!") if (!$Opt{bam_info});
}

sub read_bam_info {
    my $fileoffiles = shift;

    my @bam_files = ();
    my $fh = Open($fileoffiles);
    while (<$fh>) {
        # bampath sample library flowcell id
        if (/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)$/) {
            my ($bamfile, $sample, $library, $flowcell, $id) = ($1, $2, $3, $4, $5);
            push @bam_files, {'bam' => $bamfile,
                              'sample' => $sample,
                              'library' => $library,
                              'flowcell' => $flowcell,
                              'id' => $id};
        }
        else {
            die "File of BAM files should contain five tab-separated columns: BAM file, sample name, library id, flowcell, and id!\n";
        }
    }
    close $fh;

    return [@bam_files];
}

sub create_readtag_swarm_string {
    my $ra_bam_info = shift;
    my $readtag_analysis_obj = shift;

    my $input_directory = $readtag_analysis_obj->input_dir();
    my $output_directory = $readtag_analysis_obj->output_dir();

    my $index_string = ($Opt{index}) ? ' CREATE_INDEX=true' : '';

    my $swarm_string = '';
    foreach my $rh_baminfo (@{$ra_bam_info}) {
        my $bamfile = $rh_baminfo->{bam};
        my $input_file = $bamfile;
        if (!$Opt{bwbams}) {
            $input_file =~ s:.*/::;
            $input_file = "$input_directory/$input_file";
        }
        my $output_file = $bamfile;
        $output_file =~ s/\.bam$/.rg.bam/;
        $output_file =~ s:.*/::;
        $output_file = "$output_directory/$output_file";

        my $sample = $rh_baminfo->{sample};
        my $library = $rh_baminfo->{library};
        my $flowcell = $rh_baminfo->{flowcell};
        my $id = $rh_baminfo->{id};

        # Note: if you change java memory or cpu allocation below, change
        # values requested for job too (via the swarm_job object calls to 
        # gb_per_process and threads_per_process in this script)
        
        $swarm_string .= "java -Xmx4g -XX:ParallelGCThreads=5 -jar \$PICARDJARPATH/picard.jar AddOrReplaceReadGroups INPUT=$input_file OUTPUT=$output_file RGPL=$Opt{platform} RGSM=$sample RGLB=$library RGPU=$flowcell RGID=$id$index_string\n";
        
    }

    return $swarm_string;
}

__END__

=head1 OPTIONS

=over 4

=item B<--bam_info>

A file containing paths to local or remote bam files to be tagged, with extra colums specifying sample name, library id, flowcell, and id.

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

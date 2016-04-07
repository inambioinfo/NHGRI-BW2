package NHGRI::BW2::SbatchJob;
#######################################################################

=head1 NAME

SbatchJob.pm - A Perl module representing a slurm sbatch job on the new 
biowulf2 cluster.

=head1 DESCRIPTION

  This module is meant to be a parent class for a variety of different types
  of slurm sbatch jobs on the biowulf2 cluster.  It handles job submission,
  logging, and monitoring.

=head1 DATE

 October 20, 2015

=head1 AUTHORS

 Nancy Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=over

=cut

#######################################################################
require 5.005_62;
use strict;
use warnings;

use File::Temp qw/ :mktemp /;

use NHGRI::BW2::Session;
use GTB::File qw(Open);

our $VERSION  = '0.01';
our $REVISION = '$Id:$';

use vars qw( $DEFAULT_THREADS_PER_PROCESS $DEFAULT_GB_PER_PROCESS $DEFAULT_WALLTIME $DEFAULT_WAIT_SECS );
$DEFAULT_THREADS_PER_PROCESS = 2;
$DEFAULT_GB_PER_PROCESS = 1.5;
$DEFAULT_WAIT_SECS = 300; # number of seconds to wait between checking job status for submitted jobs
$DEFAULT_WALLTIME = '4:00:00';

#######################################################################

=item new()

  This method creates a new SbatchJob object.

  Input:  -session passes a NHGRI::BW2::Session object that connects to the
              biowulf2 cluster
          -job_name will give a name to the SbatchJob analysis
          -job_id is the job id assigned to this job by slurm
          -partition is the partition to which the job will be 
              submitted (norm, b1, etc.)
          -extra_sbatch_opts is a string specifying additional
              options to be specified when the job is submitted
              with sbatch.
       
  Output: New SbatchJob object

=cut

#######################################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $session = $params{-session};
    my $job_name = $params{-job_name};
    my $job_id = $params{-job_id};
    my $partition = $params{-partition};
    my $extra_sbatch_opts = $params{-extra_sbatch_opts};
 
    my $self = { job_name => $job_name,
                 job_id => $job_id,
                 session => $session,
                 partition => $partition,
                 extra_sbatch_opts => $extra_sbatch_opts };

    bless $self, $class;
    return $self;

} ## end new

#######################################################################

=item session()

  The method returns the current NHGRI::BW2::Session object for this run, 
  creating it if it's not yet defined.

  Input: None.
  Output: NHGRI::BW2::Session object.

=cut

#######################################################################
sub session {
    my $self  = shift;

    if (!defined ($self->{session}))
    {
        $self->{session} = NHGRI::BW2::Session->new();
    }

    return $self->{session};

} ## end session

#######################################################################

=item logfile_dir()

  The method gets or sets the name of the logs directory.  When a job
  is submitted, this path is specified with the "--logdir" option.

  Input: Optional argument sets value.
  Output: Path to the logfile directory  (scalar string).

=cut

#######################################################################
sub logfile_dir {
    my $self  = shift;

    if (defined (my $new_logfile_dir = shift))
    {
        $self->{logfile_dir} = $new_logfile_dir;
    }

    return $self->{logfile_dir};

} ## end logfile_dir

#######################################################################

=item modules()

  The method gets or sets the modules needed for this job.  When a job
  is submitted, this comma-delimited list of modules is specified with 
  the "--module" option.

  Input: Optional argument sets value.
  Output: Comma-delimited list of required modules (scalar string).

=cut

#######################################################################
sub modules {
    my $self  = shift;

    if (defined (my $new_modules = shift))
    {
        $self->{modules} = $new_modules;
    }

    return $self->{modules};

} ## end modules

#######################################################################

=item threads_per_process()

  The method gets or sets the threads_per_process needed for this job.  
  When a job is submitted, this value of threads_per_process is passed
  to sbatch with the "--cpus-per-task" option.

  Note: Within the swarm command file, it is best to specify threads per
  process to the programs you are running with the $SLURM_CPUS_PER_TASK
  environment variable.  This way, you assure that your program won't 
  try to use more threads than have been requested.

  Input: Optional argument sets value.
  Output: Number of threads per process to run each job on (scalar number).

=cut

#######################################################################
sub threads_per_process {
    my $self  = shift;

    if (defined (my $new_threads_per_process = shift))
    {
        $self->{threads_per_process} = $new_threads_per_process;
    }

    return $self->{threads_per_process};

} ## end threads_per_process

#######################################################################

=item gb_per_process()

  The method gets or sets the gb_per_process needed for this job.  
  When a job is submitted, this value of gb_per_process is passed
  to swarm with the "-g" option, or to sbatch with the "--mem=" 
  option.

  Input: Optional argument sets value.
  Output: Number of Gb memory per process allocated to each subjob 
      (scalar number).

=cut

#######################################################################
sub gb_per_process {
    my $self  = shift;

    if (defined (my $new_gb_per_process = shift))
    {
        $self->{gb_per_process} = $new_gb_per_process;
    }

    return $self->{gb_per_process};

} ## end gb_per_process

#######################################################################

=item walltime_per_process()

  The method gets or sets the walltime_per_process needed for this job.  
  When a job is submitted, this value of walltime_per_process is passed
  to swarm or sbatch with the "-t" option.

  Input: Optional argument sets value.
  Output: Number of hours of walltime required for each job in the swarm
      (scalar number).

=cut

#######################################################################
sub walltime_per_process {
    my $self  = shift;

    if (defined (my $new_walltime_per_process = shift))
    {
        $self->{walltime_per_process} = $new_walltime_per_process;
    }

    return $self->{walltime_per_process};

} ## end walltime_per_process

#######################################################################

=item gres()

  The method gets or sets the gres (general resources) that should be 
  requested when this job is submitted.

  Input: Optional argument sets value.
  Output: Value of gres (scalar string).

=cut

#######################################################################
sub gres {
    my $self  = shift;

    if (defined (my $new_gres = shift))
    {
        $self->{gres} = $new_gres;
    }

    return $self->{gres};

} ## gres

#######################################################################

=item partition()

  The method gets or sets the partition that should be requested when 
  this job is submitted (with the --partition option of the swarm
  command).

  Input: Optional argument sets value.
  Output: Value of partition (scalar string).

=cut

#######################################################################
sub partition {
    my $self  = shift;

    if (defined (my $new_partition = shift))
    {
        $self->{partition} = $new_partition;
    }

    return $self->{partition};

} ## partition

#######################################################################

=item job_name()

  The method gets or sets the value of "job_name", which is the 
  job_name for the run.

  Input: Optional argument sets value.
  Output: Name of the job (scalar string).

=cut

#######################################################################
sub job_name {
    my $self  = shift;

    if (defined (my $new_job_name = shift))
    {
        $self->{job_name} = $new_job_name;
    }
    return $self->{job_name};

} ## end job_name

#######################################################################

=item jobid()

  The method gets or sets the value of "jobid", which is the job id for 
  the sbatch job on the biowulf2 cluster.

  Input: Optional argument sets value.
  Output: Id of the job (scalar integer).

=cut

#######################################################################
sub jobid {
    my $self  = shift;

    if (defined (my $new_jobid = shift))
    {
        $self->{jobid} = $new_jobid;
    }
    
    return $self->{jobid};

} ## end jobid

#######################################################################

=item extra_sbatch_opts()

  The method gets or sets the extra_sbatch_opts that should be used when this
  job is submitted.

  Input: Optional argument sets value.
  Output: Value of extra_sbatch_opts (scalar string).

=cut

#######################################################################
sub extra_sbatch_opts {
    my $self  = shift;

    if (defined (my $new_extra_sbatch_opts = shift))
    {
        $self->{extra_sbatch_opts} = $new_extra_sbatch_opts;
    }

    return $self->{extra_sbatch_opts};

} ## extra_sbatch_opts

#######################################################################

=item sbatch_command_file()

  The method gets or sets the value of "sbatch_command_file", which is the 
  path to the sbatch command file for the run.

  Input: Optional argument sets value.
  Output: Path to the sbatch command file (scalar string).

=cut

#######################################################################
sub sbatch_command_file {
    my $self  = shift;

    if (defined (my $new_sbatch_command_file = shift))
    {
        $self->{sbatch_command_file} = $new_sbatch_command_file;
    }
    return $self->{sbatch_command_file};

} ## end sbatch_command_file

#######################################################################

=item write_and_transfer_sbatch_file()

  The method will write a temporary local sbatch file, then transfer the
  file to the specified remote path on the biowulf2 server, setting the
  object's "sbatch_command_file" path to the new remote location.

  Input: -remotepath specifies the path to the new sbatch file on the 
         biowulf2 cluster.
  Output: Path to the remote sbatch command file (scalar string).

=cut

#######################################################################
sub write_and_transfer_sbatch_file {
    my $self  = shift;
    my %params = @_;

    my $session = $self->session();

    my $remotepath = $params{-remotepath} or
        die "Must specify a remote path with -remotepath in write_and_transfer_sbatch_file!\n";

    my $local_sbatch = mktemp( "sbatchXXXXXX" );
    $self->write_sbatch_file(-file => $local_sbatch);
    $session->transfer_file_to_biowulf(
                 -localpath => $local_sbatch,
                 -remotepath => $remotepath);

    $self->sbatch_command_file($remotepath);

    unlink $local_sbatch;

    return $remotepath;

} ## end write_and_transfer_sbatch_file

#######################################################################

=item write_sbatch_file()

  The method writes a local file at the path specified by the "-file"
  param and print the value of "sbatch_command_string" to that file,
  along with necessary slurm directives and commands to load modules.

  Input: -file specifies the local path to which the sbatch file should
         be written.
  Output: 1 if successful (dies otherwise)

=cut

#######################################################################
sub write_sbatch_file {
    my $self  = shift;
    my %params = @_;
    my $file = $params{'-file'};

    if (!defined ($file)) {
        die "Must pass a file path to write_sbatch_file with -file parameter!\n";
    }

    my $sbatch_string = $self->{sbatch_command_string};
    if (!$sbatch_string) {
        die "Must set the value of sbatch_command_string before calling write_sbatch_file method!\n";
    }

    my $sbatch_fh = Open($file, "w");
    print $sbatch_fh "#!/bin/bash\n";
    my $ra_modules = $self->modules() || [];
    if (@{$ra_modules}) {
        foreach my $module (@{$ra_modules}) {
            print $sbatch_fh "module load $module || exit 1\n";
        }
    }
    print $sbatch_fh "\n";
    print $sbatch_fh $self->{sbatch_command_string};
    close $sbatch_fh;

    return 1;

} ## end write_sbatch_file

#######################################################################

=item sbatch_command_string()

  The method gets or sets the value of "sbatch_command_string", which is 
  the string used to submit this sbatch job.

  Input: Optional argument sets value.
  Output: sbatch submit command (scalar string).

=cut

#######################################################################
sub sbatch_command_string {
    my $self  = shift;

    if (defined (my $new_sbatch_command_string = shift))
    {
        $self->{sbatch_command_string} = $new_sbatch_command_string;
    }
    return $self->{sbatch_command_string};

} ## end sbatch_command_string

#######################################################################

=item all_sacct_info()

  The method calls the sacct command on the biowulf2 cluster for the
  job and returns the string returned by the command.

  Input: None.
  Output: The output of the sacct command

=cut

#######################################################################
sub all_sacct_info {
    my $self = shift;
    my %params = @_;
    my $jobid = $self->jobid();
    if (!$jobid) {
         die "Can\'t call sacct on a job that hasn\'t been submitted yet!\n";
    }

    my $command = "sacct -j $jobid --format=ALL";

    my $session = $self->session();

    my $return_string = $session->remote_command(-command => $command);

    return $return_string;
    
} ## end all_sacct_info

#######################################################################

=item jobhist_info()

  The method calls the sacct command on the biowulf2 cluster for the
  job and returns the string returned by the command.

  Input: None.
  Output: The output of the sacct command

=cut

#######################################################################
sub jobhist_info {
    my $self = shift;
    my %params = @_;
    my $jobid = $self->jobid();
    if (!$jobid) {
         die "Can\'t call sacct on a job that hasn\'t been submitted yet!\n";
    }

    my $command = "jobhist $jobid";

    my $session = $self->session();
    my $return_string = $session->remote_command(-command => $command);

    return $return_string;
    
} ## end jobhist_info

#######################################################################

=item wait_for_sbatch()

  This method checks every $DEFAULT_WAIT_SECS seconds for a job to
  finish running.

  Input: -jobid passes the job id for which we are waiting.
         -wait specifies a number of seconds to wait between checks.
  Output: Number of job components with a "State" value other than
      "COMPLETED".

=cut

#######################################################################
sub wait_for_sbatch {
    my $self  = shift;
    my %params = @_;
    my $jobid = $params{-jobid} || $self->jobid()
        or die "Must pass a job id value to \"wait_for_sbatch\" method!";
    my $wait_seconds = $params{-wait} || $DEFAULT_WAIT_SECS;

    my $finished = 0;

    my $session = $self->session();
    while (!$finished)
    {
        my $sacct = $session->remote_command(-command => "sacct -j $jobid --format=state\%-30");
        print "SACCT output: $sacct\n";

        my %status_counts = ('COMPLETED' => 0,
                             'OTHER' => 0);

        foreach my $line (split /\n/, $sacct) {
            next if ($line =~ /State/);
            next if ($line =~ /^[\-\s]+$/);
 
            if ($line =~ /^PENDING\s*/) {
                $finished = 0;
                last;
            }
            elsif ($line =~ /^RUNNING\s*/) {
                $finished = 0;
                last;
            }
            elsif ($line =~ /^COMPLETED\s*/) {
                $status_counts{'COMPLETED'}++;
            }
            else {
                $status_counts{'OTHER'}++;
            }
            $finished = 1;
        }

        if (!$finished) {
            print "Waiting for job(s) to finish\n";
            sleep $wait_seconds;
        }
        else {
            return $status_counts{'OTHER'};
        }
    }

} ## wait_for_sbatch

#######################################################################

=item submit_sbatch_command()

  The method submits the job with all of its parameters, and sets the 
  "jobid" value of the object.

  Input: -extra_opts => string to be used as an option to sbatch
         Note: without requested resources, sbatch on
         biowulf2 currently requests 2 CPUs and 4 GB memory per 
         process.  This option overrides the object's value of 
         "extra_sbatch_opts".
         -dry_run => if true, will generate the submit command string,
         store it in "sbatch_command_string", but not submit the job.
  Output: Output from sbatch

=cut

#######################################################################
sub submit_sbatch_command {
    my $self  = shift;
    my %params = @_;

    my $optstring = $params{'-extra_opts'} || $self->extra_sbatch_opts();

    my $run_string = "sbatch";

    if (my $tpp = $self->threads_per_process()) {
        $run_string .= " --cpus-per-task=$tpp"; 
    }
    else {
        $self->threads_per_process($DEFAULT_THREADS_PER_PROCESS);
    }

    if (my $gpp = $self->gb_per_process()) {
        $run_string .= " --mem=$gpp"."g"; 
    }
    else {
        $self->gb_per_process($DEFAULT_GB_PER_PROCESS);
    }

    if (my $wt = $self->walltime_per_process()) {
        if ($wt =~ /^\d+$/) {
            $wt = "$wt:00:00";
        }
        $run_string .= " --time=$wt"; 
    }
    else {
        $self->walltime_per_process($DEFAULT_WALLTIME);
    }

    if (my $gres = $self->gres()) {
        $run_string .= " -gres=$gres"; 
    }

    if (my $logfile_dir = $self->logfile_dir()) {
        $run_string .= " -o $logfile_dir/sbatch_\%A_\%a.o -e $logfile_dir/sbatch_\%A_\%a.e"; 
    }

    if (my $partition = $self->partition()) {
        $run_string .= " --partition=$partition"; 
    }

    if ($optstring) {
        $run_string .= " $optstring";
    }

    my $sbatch_command_file = $self->sbatch_command_file();
    $run_string .= " $sbatch_command_file";

    $self->sbatch_command_string($run_string);

    unless ($params{'-quiet'}) {
        print "$run_string\n";
    }

    my $session = $self->session();
    my $response = ($params{-dry_run}) ? '' : $session->remote_command( -command => $run_string);

    if ($response =~ /(\d+)/) {
        $self->jobid($1);
    }
    else {
        die "Unable to parse job id from server response $response!\n";
    }

    return $self->jobid();

} ## end submit_sbatch_command

#######################################################################

=back

=cut

1;
__END__

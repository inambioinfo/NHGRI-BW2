package NHGRI::BW2::Biowulf2Job;
#######################################################################

=head1 NAME

Biowulf2Job.pm - A Perl module representing a slurm job on the new 
biowulf2 cluster.

=head1 DESCRIPTION

  This module is meant to be a parent class for swarm (SwarmJob) or
  sbatch (SbatchJob) objects represenging jobs on the biowulf2 cluster.
  It handles job submission, logging, and monitoring.

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

our $VERSION  = '0.01';
our $REVISION = '$Id:$';
our $PICARDVERSION = '1.139';

use vars qw( $DEFAULT_THREADS_PER_PROCESS $DEFAULT_GB_PER_PROCESS $DEFAULT_WALLTIME $DEFAULT_WAIT_SECS );
$DEFAULT_THREADS_PER_PROCESS = 2;
$DEFAULT_GB_PER_PROCESS = 1.5;
$DEFAULT_WAIT_SECS = 300; # number of seconds to wait between checking job status for submitted jobs
$DEFAULT_WALLTIME = '4:00:00';

#######################################################################

=item new()

  This method creates a new Biowulf2Job object.

  Input:  -session passes a NHGRI::BW2::Session object that connects to the
              biowulf2 cluster
          -job_name will give a name to the Biowulf2Job analysis
          -job_id is the job id assigned to this job by slurm
       
  Output: New Biowulf2Job object

=cut

#######################################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $session = $params{-session};
    my $job_name = $params{-job_name};
    my $job_id = $params{-job_id};
    my $partition = $params{-partition};
 
    my $self = { job_name => $job_name,
                 job_id => $job_id,
                 session => $session,
                 partition => $partition };

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
  to swarm with the "-t" option.

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
  the swarm job on the biowulf2 cluster.

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

=item retrieve_log_files()

  This method scp's all contents of a job's log file directory on biowulf2
  to a directory on the localhost.

  Input: -local_dir specifies local directory in which to deposit files.  
      Otherwise, output files are written to the current working directory.
  Output: 1 if successful.

=cut

#######################################################################
sub retrieve_log_files {
    my $self = shift;
    my %params = @_;
    my $local_dir = $params{-local_dir} || '.';
    my $session = $self->session();
    my $logfile_dir = $self->logfile_dir();

    $session->transfer_directory_from_biowulf(-remotepath => $logfile_dir,
                                              -localpath => $local_dir);

} ## retrieve_log_files

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

=back

=cut

1;
__END__

package NHGRI::BW2::SwarmJob;
#######################################################################

=head1 NAME

SwarmJob.pm - A Perl module representing a slurm swarm job on the new 
biowulf2 cluster.

=head1 DESCRIPTION

  This module is meant to be a parent class for a variety of different types
  of slurm swarm jobs on the biowulf2 cluster.  It handles job submission
  and monitoring.

=head1 DATE

 October 8, 2015

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

  This method creates a new SwarmJob object.

  Input:  -session passes a NHGRI::BW2::Session object that connects to the
              biowulf2 cluster
          -job_name will give a name to the SwarmJob analysis
          -job_id is the job id assigned to this job by slurm
          -partition is the partition to which the job will be
              submitted (norm, b1, etc.)
          -extra_swarm_opts is a string specifying additional
              options to be specified when the job is submitted
              with swarm.

       
  Output: New SwarmJob object

=cut

#######################################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $session = $params{-session};
    my $job_name = $params{-job_name};
    my $job_id = $params{-job_id};
    my $partition = $params{-partition};
    my $extra_swarm_opts = $params{-extra_swarm_opts};
 
    my $self = { job_name => $job_name,
                 job_id => $job_id,
                 session => $session,
                 partition => $partition,
                 extra_swarm_opts => $extra_swarm_opts };

    bless $self, $class;
    return $self;

} ## end new

#######################################################################

=item extra_swarm_opts()

  The method gets or sets the extra_swarm_opts that should be used when this
  job is submitted.

  Input: Optional argument sets value.
  Output: Value of extra_swarm_opts (scalar string).

=cut

#######################################################################
sub extra_swarm_opts {
    my $self  = shift;

    if (defined (my $new_extra_swarm_opts = shift))
    {
        $self->{extra_swarm_opts} = $new_extra_swarm_opts;
    }

    return $self->{extra_swarm_opts};

} ## extra_swarm_opts

#######################################################################

=item swarm_command_file()

  The method gets or sets the value of "swarm_command_file", which is the 
  path to the swarm command file for the run.

  Input: Optional argument sets value.
  Output: Path to the swarm command file (scalar string).

=cut

#######################################################################
sub swarm_command_file {
    my $self  = shift;

    if (defined (my $new_swarm_command_file = shift))
    {
        $self->{swarm_command_file} = $new_swarm_command_file;
    }
    return $self->{swarm_command_file};

} ## end swarm_command_file

#######################################################################

=item write_and_transfer_swarm_file()

  The method will write a temporary local swarm file, then transfer the
  file to the specified remote path on the biowulf2 server, setting the
  object's "swarm_command_file" path to the new remote location.

  Input: -remotepath specifies the path to the new swarm file on the 
         biowulf2 cluster.
  Output: Path to the remote swarm command file (scalar string).

=cut

#######################################################################
sub write_and_transfer_swarm_file {
    my $self  = shift;
    my %params = @_;

    my $session = $self->session();

    my $remotepath = $params{-remotepath} or
        die "Must specify a remote path with -remotepath in write_and_transfer_swarm_file!\n";

    my $local_swarm = mktemp( "swarmXXXXXX" );
    $self->write_swarm_file(-file => $local_swarm);
    $session->transfer_file_to_biowulf(
                 -localpath => $local_swarm,
                 -remotepath => $remotepath);

    $self->swarm_command_file($remotepath);

    unlink $local_swarm;

    return $remotepath;

} ## end write_and_transfer_swarm_file

#######################################################################

=item write_swarm_file()

  The method writes a local file at the path specified by the "-file"
  param and print the value of "swarm_command_string" to that file.

  Input: -file specifies the local path to which the swarm file should
         be written.
  Output: 1 if successful (dies otherwise)

=cut

#######################################################################
sub write_swarm_file {
    my $self  = shift;
    my %params = @_;
    my $file = $params{'-file'};

    if (!defined ($file)) {
        die "Must pass a file path to write_swarm_file with -file parameter!\n";
    }

    my $swarm_string = $self->{swarm_command_string};
    if (!$swarm_string) {
        die "Must set the value of swarm_command_string before calling write_swarm_file method!\n";
    }

    my $swarm_fh = Open($file, "w");
    print $swarm_fh $self->{swarm_command_string};
    close $swarm_fh;

    return 1;

} ## end write_swarm_file

#######################################################################

=item swarm_command_string()

  The method gets or sets the value of "swarm_command_string", which is 
  the string used to submit this swarm job.

  Input: Optional argument sets value.
  Output: Swarm submit command (scalar string).

=cut

#######################################################################
sub swarm_command_string {
    my $self  = shift;

    if (defined (my $new_swarm_command_string = shift))
    {
        $self->{swarm_command_string} = $new_swarm_command_string;
    }
    return $self->{swarm_command_string};

} ## end swarm_command_string

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

=item number_of_jobs()

  The method gets or sets the value of "number_of_jobs", which is the 
  number of commands in the swarm file for this SwarmJob job.

  Input: Optional argument sets value.
  Output: Number of jobs (scalar number).

=cut

#######################################################################
sub number_of_jobs {
    my $self  = shift;

    if (defined (my $new_number_of_jobs = shift))
    {
        $self->{number_of_jobs} = $new_number_of_jobs;
    }
    
    return $self->{number_of_jobs};

} ## end number_of_jobs

#######################################################################

=item wait_for_swarm()

  This method checks every $DEFAULT_WAIT_SECS seconds for all commands 
  in a swarm to finish running.

  Input: -jobid passes the job id for which we are waiting.
         -wait specifies a number of seconds to wait between checks.
  Output: Number of swarm components with a "State" value other than
      "COMPLETED".

=cut

#######################################################################
sub wait_for_swarm {
    my $self  = shift;
    my %params = @_;
    my $jobid = $params{-jobid} || $self->jobid()
        or die "Must pass a job id value to \"wait_for_swarm\" method!";
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

} ## wait_for_swarm

#######################################################################

=item submit_swarm_command()

  The method submits the job with all of its parameters, and sets the 
  "jobid" value of the object.

  Input: -extra_opts => string to be used as an option to swarm
         Note: without requested resources, sbatch on
         biowulf2 currently requests 2 CPUs and 1.5GB memory per 
         process.  This option overrides the object's value of 
         "extra_swarm_opts".
         -dry_run => if true, will generate the submit command string,
         store it in "swarm_command_string", but not submit the swarm
         command.
  Output: Output from swarm

=cut

#######################################################################
sub submit_swarm_command {
    my $self  = shift;
    my %params = @_;

    my $optstring = $params{'-extra_opts'} || $self->extra_swarm_opts();

    my $swarm_command_file = $self->swarm_command_file();

    my $run_string = "swarm -f $swarm_command_file";

    if (my $tpp = $self->threads_per_process()) {
        $run_string .= " -t $tpp"; 
    }
    else {
        $self->threads_per_process($DEFAULT_THREADS_PER_PROCESS);
    }

    if (my $gpp = $self->gb_per_process()) {
        $run_string .= " -g $gpp"; 
    }
    else {
        $self->gb_per_process($DEFAULT_GB_PER_PROCESS);
    }

    if (my $wt = $self->walltime_per_process()) {
        if ($wt =~ /^\d+$/) {
            $wt = "$wt:00:00";
        }
        $run_string .= " --time $wt"; 
    }
    else {
        $self->walltime_per_process($DEFAULT_WALLTIME);
    }

    if (my $gres = $self->gres()) {
        $run_string .= " -gres=$gres"; 
    }

    if (my $logfile_dir = $self->logfile_dir()) {
        $run_string .= " --logdir $logfile_dir"; 
    }

    if (my $ra_modules = $self->modules()) {
        my $module_string = join ',', @{$ra_modules};
        $run_string .= " --module $module_string"; 
    }

    if (my $partition = $self->partition()) {
        $run_string .= " --partition $partition"; 
    }

    if ($optstring) {
        $run_string .= " $optstring";
    }

    $self->swarm_command_string($run_string);

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

    my $submit_file = $self->jobid()."_swarm_submit_command.txt";
    my $submit_fh = Open($submit_file, "w");
    print $submit_fh "$run_string\n";
    close $submit_fh;

    my $submit_remotepath = $swarm_command_file;
    $submit_remotepath =~ s:[^/]+$::; # remove filename and keep directory
    $submit_remotepath .= $submit_file;

    $session->transfer_file_to_biowulf(
                 -localpath => $submit_file,
                 -remotepath => $submit_remotepath);

    #$self->submit_return_when_finished(-job_id => $self->jobid());

    return $self->jobid();

} ## end submit_swarm_command

#######################################################################

=item submit_return_when_finished()

  The method submits the job with all of its parameters, and sets the 
  "jobid" value of the object.

  Input: -job_id => required job id.  After completion of this job
         on the biowulf2 cluster (--dependency=afterany:JobID), a
         file will be created in the current working directory on the 
         local server using the "touch" command.
  Output: Output from swarm

=cut

#######################################################################
sub submit_return_when_finished {
    my $self  = shift;
    my %params = @_;
    my $jobid = $params{-job_id}
        or die "Must supply a job id with the -job_id parameter to submit_return_when_finished.\n";

    my $submit_string = "sbatch --dependency=afterany:$jobid\n";

} # end submit_return_when_finished

#######################################################################

=back

=cut

1;
__END__

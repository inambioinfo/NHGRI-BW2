package NHGRI::BW2::Analysis;
#######################################################################

=head1 NAME

Analysis.pm - A Perl module to track and contain information about an
  analysis consisting of jobs run on the biowulf2 cluster.

=head1 DESCRIPTION

  This module contains methods to set up an analysis directory on the
  biowulf2 server, transfer input files from the local host to the 
  remote server, run and monitor jobs that are part of the analysis,
  and transfer output files and job history back from the remote to
  the local server.

=head1 DATE

 October 9, 2015

=head1 AUTHORS

 Nancy Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=over

=cut

#######################################################################
require 5.005_62;
use strict;
use warnings;

use NHGRI::BW2::Session;

our $VERSION  = '0.01';
our $REVISION = '$Id:$';

#######################################################################

=item new()

  This method creates a new Analysis object.

  Input:  -session passes an NHGRI::BW2::Session object to be used for
              communication with the biowulf2 server.  If not passed in
              the constructor, a new Session object will be created if
              needed. (Optional)
          -analysis_name contains the name of this analysis.  If not 
              specified, the analysis name is constructed from the 
              analysis "base name", followed by the process id, 
              followed by the Linux time in "seconds since the Epoch". (Optional)
          -analysis_basename contains a base name that is used to 
              construct the "analysis_name" if it is not specified
              (see above) Default: "analysis". (Optional)
          -analysis_dir contains the path to a directory on the biowulf2
              server in which input, output, scripts, and logfiles
              are stored for the analysis.  By default, it is a
              directory with name "analysis_name" within the user's
              data directory. (Optional)
          -input_dir contains the path to a directory on the biowulf2
              server in which input files are stored for the analysis.
              By default, it is a directory with name "input_files" 
              within the analysis directory. (Optional)
          -input_files contains a reference to a list of paths of files
              on the local host which will be tranfserred to the 
              biowulf2 host's input_dir for use in the analysis. (Optional)
          -temporary_dir contains the path to a directory on the biowulf2
              server in which temporary files are stored for the analysis.
              By default, it is a directory with name "temporary_files" 
              within the analysis directory. (Optional)
          -output_dir contains the path to a directory on the biowulf2
              server in which output files are stored for the analysis.
              By default, it is a directory with name "output_files" 
              within the analysis directory. (Optional)
          -scripts_dir contains the path to a directory on the biowulf2
              server in which batch and swarm scripts are stored for the 
              analysis.  By default, it is a directory with name "scripts" 
              within the analysis directory. (Optional)
          -logfile_dir contains the path to a directory on the biowulf2
              server in which log files are stored for the analysis.
              By default, it is a directory with name "logfiles" within 
              the analysis directory. (Optional)
          -use_scratch if true, will create the analysis directory in 
              the user's scratch directory with name "analysis_name"
              (only relevant if -analysis_dir is not specified). (Optional)
       
  Output: New Analysis object

=cut

#######################################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $session = $params{-session};
    my $analysis_name = $params{-analysis_name};
    my $analysis_basename = $params{-analysis_basename} || 'analysis';
    my $analysis_dir = $params{-analysis_dir};
    my $input_dir = $params{-input_dir};
    my $input_files = $params{-input_files} || [];
    my $temporary_dir = $params{-temporary_dir};
    my $output_dir = $params{-output_dir};
    my $logfile_dir = $params{-logfile_dir};
    my $use_scratch = $params{-use_scratch} || 0;
 
    my $self = { session => $session,
                 analysis_name => $analysis_name,
                 analysis_basename => $analysis_basename,
                 analysis_dir => $analysis_dir,
                 input_dir => $input_dir,
                 input_files => $input_files,
                 temporary_dir => $temporary_dir,
                 output_dir => $output_dir,
                 logfile_dir => $logfile_dir,
                 use_scratch => $use_scratch };

    bless $self, $class;
    return $self;

} ## end new

#######################################################################

=item analysis_name()

  The method gets or sets the value of "analysis_name", which is a 
  string giving a unique name to this analysis.

  Input: Optional argument sets value.
  Output: Value of analysis_name (scalar string)

=cut

#######################################################################
sub analysis_name {
    my $self  = shift;

    if (defined (my $new_analysis_name = shift)) {
        $self->{analysis_name} = $new_analysis_name;
    }

    if (!(defined($self->{analysis_name}))) {
        my $pid = $$;
        my $timeslug = time();

        my $basename = $self->analysis_basename();
        $self->{analysis_name} = "$basename.$pid.$timeslug";
    }

    return $self->{analysis_name};

} ## end analysis_name

#######################################################################

=item session()

  The method gets or sets the value of "session", which is an 
  NHGRI::BW2::Session object used for communication with the biowulf2
  host.

  Input: Optional argument sets value.  If undefined, a new object
      will be created and returned.
  Output: Value of session (NHGRI::BW2::Session object)

=cut

#######################################################################
sub session {
    my $self  = shift;

    if (defined (my $new_session = shift)) {
        $self->{session} = $new_session;
    }

    if (!(defined($self->{session}))) {
        $self->{session} = NHGRI::BW2::Session->new();
    }

    return $self->{session};

} ## end session

#######################################################################

=item analysis_basename()

  The method gets or sets the value of "analysis_basename", which is a 
  string from which the analysis name is constructed if it was not
  specified by the user.

  Input: Optional argument sets value.
  Output: Value of analysis_basename (scalar string)

=cut

#######################################################################
sub analysis_basename {
    my $self  = shift;

    if (defined (my $new_analysis_basename = shift)) {
        $self->{analysis_basename} = $new_analysis_basename;
    }

    return $self->{analysis_basename};

} ## end analysis_basename

#######################################################################

=item analysis_dir()

  The method gets or sets the value of "analysis_dir", which is a 
  string giving a unique dir to this analysis.  If undefined, the
  value of analysis_dir will be set to be a subdirectory of the user's
  data or scratch directory (the latter if the -use_scratch option is
  passed) named with the object's value of "analysis_name".

  Input: Optional argument sets value.
  Output: Value of analysis_dir (scalar string)

=cut

#######################################################################
sub analysis_dir {
    my $self  = shift;

    if (defined (my $new_analysis_dir = shift)) {
        $self->{analysis_dir} = $new_analysis_dir;
    }

    if (!(defined($self->{analysis_dir}))) {
        my $session_obj = $self->session();
        my $username = $session_obj->username();
        my $analysis_name = $self->analysis_name();
        my $data_or_scratch_dir = ($self->use_scratch()) ? 
                 $session_obj->scratch_dir() : $session_obj->data_dir();
        $self->{analysis_dir} = "$data_or_scratch_dir/$username/$analysis_name";
    }

    return $self->{analysis_dir};

} ## end analysis_dir

#######################################################################

=item create_analysis_directory()

  The method creates the directory "analysis_dir" on the biowulf2 server,
  as well as the subdirectories "input_dir", "output_dir", "logfile_dir",
  and "temporary_dir".

  Input: Optional argument sets value.
  Output: Value of analysis_dir (scalar string)

=cut

#######################################################################
sub create_analysis_directory {
    my $self  = shift;

    my $session = $self->session();
    my $analysis_dir = $self->analysis_dir();
    my $input_dir = $self->input_dir();
    my $output_dir = $self->output_dir();
    my $scripts_dir = $self->scripts_dir();
    my $logfile_dir = $self->logfile_dir();
    my $temporary_dir = $self->temporary_dir();

    foreach my $dir ( $analysis_dir, $input_dir, $output_dir, $scripts_dir, $logfile_dir, $temporary_dir) {
        $session->make_remote_directory(-dirname => $dir);
    }

    return $analysis_dir;

} ## end create_analysis_directory

#######################################################################

=item input_dir()

  The method gets or sets the value of "input_dir", which is a 
  string giving the biowulf2 path to a directory to be used for input
  files that will be transferred from the local host for analysis.  
  If undefined, the value of input_dir will be set to be a subdirectory
  of the analysis directory named "input_files".

  Input: Optional argument sets value.
  Output: Value of input_dir (scalar string)

=cut

#######################################################################
sub input_dir {
    my $self  = shift;

    if (defined (my $new_input_dir = shift)) {
        $self->{input_dir} = $new_input_dir;
    }

    if (!(defined($self->{input_dir}))) {
        my $analysis_dir = $self->analysis_dir();
        $self->{input_dir} = "$analysis_dir/input_files";
    }

    return $self->{input_dir};

} ## end input_dir

#######################################################################

=item input_files()

  The method gets or sets the value of "input_files", which is a 
  reference to a list of local file paths of files which should
  be transferred to the biowulf2 "input_dir" before launching jobs.

  Input: Optional argument sets value.
  Output: Value of input_files (reference to a list of scalar strings)

=cut

#######################################################################
sub input_files {
    my $self  = shift;

    if (defined (my $new_input_files = shift)) {
        $self->{input_files} = $new_input_files;
    }

    return $self->{input_files};

} ## end input_files

#######################################################################

=item temporary_dir()

  The method gets or sets the value of "temporary_dir", which is a 
  string giving the biowulf2 path to a directory to be used for temporary
  files that will be transferred from the local host for analysis.  
  If undefined, the value of temporary_dir will be set to be a subdirectory
  of the analysis directory named "temporary_files".

  Input: Optional argument sets value.
  Output: Value of temporary_dir (scalar string)

=cut

#######################################################################
sub temporary_dir {
    my $self  = shift;

    if (defined (my $new_temporary_dir = shift)) {
        $self->{temporary_dir} = $new_temporary_dir;
    }

    if (!(defined($self->{temporary_dir}))) {
        my $analysis_dir = $self->analysis_dir();
        $self->{temporary_dir} = "$analysis_dir/temporary_files";
    }

    return $self->{temporary_dir};

} ## end temporary_dir

#######################################################################

=item output_dir()

  The method gets or sets the value of "output_dir", which is a 
  string giving the biowulf2 path to a directory to be used for output
  files that will be transferred from the local host for analysis.  
  If undefined, the value of output_dir will be set to be a subdirectory
  of the analysis directory named "output_files".

  Input: Optional argument sets value.
  Output: Value of output_dir (scalar string)

=cut

#######################################################################
sub output_dir {
    my $self  = shift;

    if (defined (my $new_output_dir = shift)) {
        $self->{output_dir} = $new_output_dir;
    }

    if (!(defined($self->{output_dir}))) {
        my $analysis_dir = $self->analysis_dir();
        $self->{output_dir} = "$analysis_dir/output_files";
    }

    return $self->{output_dir};

} ## end output_dir

#######################################################################

=item logfile_dir()

  The method gets or sets the value of "logfile_dir", which is a 
  string giving the biowulf2 path to a directory to be used for logfile
  files that will be transferred from the local host for analysis.  
  If undefined, the value of logfile_dir will be set to be a subdirectory
  of the analysis directory named "logfiles".

  Input: Optional argument sets value.
  Output: Value of logfile_dir (scalar string)

=cut

#######################################################################
sub logfile_dir {
    my $self  = shift;

    if (defined (my $new_logfile_dir = shift)) {
        $self->{logfile_dir} = $new_logfile_dir;
    }

    if (!(defined($self->{logfile_dir}))) {
        my $analysis_dir = $self->analysis_dir();
        $self->{logfile_dir} = "$analysis_dir/logfiles";
    }

    return $self->{logfile_dir};

} ## end logfile_dir

#######################################################################

=item use_scratch()

  The method gets or sets the value of "use_scratch", which is 1 or 0
  depending on whether this analysis should use the "scratch" directory
  on biowulf2 to store its files.

  Input: Optional argument sets value.
  Output: Value of use_scratch (scalar 0 or 1)

=cut

#######################################################################
sub use_scratch {
    my $self  = shift;

    if (defined (my $new_use_scratch = shift)) {
        $self->{use_scratch} = $new_use_scratch;
    }

    return $self->{use_scratch};

} ## end use_scratch

#######################################################################

=item scripts_dir()

  The method gets or sets the value of "scripts_dir", which is a 
  string giving the biowulf2 path to a directory to be used for scripts
  files that will be transferred from the local host for analysis.  
  If undefined, the value of scripts_dir will be set to be a subdirectory
  of the analysis directory named "scripts".

  Input: Optional argument sets value.
  Output: Value of scripts_dir (scalar string)

=cut

#######################################################################
sub scripts_dir {
    my $self  = shift;

    if (defined (my $new_scripts_dir = shift)) {
        $self->{scripts_dir} = $new_scripts_dir;
    }

    if (!(defined($self->{scripts_dir}))) {
        my $analysis_dir = $self->analysis_dir();
        $self->{scripts_dir} = "$analysis_dir/scripts";
    }

    return $self->{scripts_dir};

} ## end scripts_dir

#######################################################################

=item jobs()

  The method gets or sets the value of "jobs", which is a reference to 
  a list of Job (NHGRI::BW2::SwarmJob) objects for this analysis.

  Input: Optional argument sets value.
  Output: Reference to a list of Job objects.

=cut

#######################################################################
sub jobs {
    my $self  = shift;

    if (defined (my $new_jobs = shift))
    {
        $self->{jobs} = $new_jobs;
    }

    return $self->{jobs};

} ## end jobs

#######################################################################

=item transfer_input_files()

  This method uses the parameters specified to transfer this analysis's
  input files to the remote host's "input_dir" directory.

  Input: -md5sum_check if true, has the Session module do a check of 
         md5sums after the transfer to be sure files transferred 
         successfully (default is no check)
  Output: Full path of the input subdirectory (scalar string).

=cut

#######################################################################
sub transfer_input_files {
    my $self  = shift;
    my %params = @_;
    my $md5sum_check = $params{'-md5sum_check'} || 0;

    my $session = $self->session();

    my $input_dir = $self->input_dir();
    my $ra_input_files = $self->input_files();
    $session->transfer_files_to_biowulf_directory(
                       -localfiles => $ra_input_files,
                       -remotepath => $input_dir);

} ## end transfer_input_files

#######################################################################

=item retrieve_output_files()

  This method scp's all output files from the biowulf2 "output_dir" 
  directory to the localhost.

  Input: Analysis object, -local_dir specifies local directory in which to 
      deposit files.  Otherwise, output files are written to the current 
      working directory.
  Output: 1 if successful.

=cut

#######################################################################
sub retrieve_output_files {
    my $self = shift;
    my %params = @_;
    my $local_dir = $params{-local_dir} || '.';
    my $session = $self->session();
    my $output_dir = $self->output_dir();

    $session->transfer_directory_from_biowulf(
                       -localpath => $local_dir,
                       -remotepath => $output_dir);

} ## retrieve_output_files

#######################################################################

=item retrieve_log_files()

  This method scp's the analysis's "logfiles" directory on biowulf2 to the 
  localhost's logs directory.

  Input: Analysis object, -local_dir specifies local directory in which to 
      deposit files.  Otherwise, output files are written to the current 
      working directory.
  Output: 1 if successful.

=cut

#######################################################################
sub retrieve_log_files {
    my $self = shift;
    my %params = @_;
    my $local_dir = $params{-local_dir} || '.';
    my $session = $self->session();
    my $logfile_dir = $self->logfile_dir();

    $session->transfer_directory_from_biowulf(
                       -localpath => $local_dir,
                       -remotepath => $logfile_dir);

} ## retrieve_log_files

#######################################################################

=item remove_temporary_directory()

  This method removes the remote temporary directory for this analysis.

  Input: None.
  Output: 1 if successful.

=cut

#######################################################################
sub remove_work_directory {
    my $self  = shift;
    my $temporary_directory = $self->temporary_directory();
    my $session = $self->session();

    if (!$temporary_directory || $temporary_directory !~ m:^/:) {
        die "Refusing to issue a remove command without a valid directory!\n";
    }

    $session->remote_command("rm -r $temporary_directory");

} ## remove_temporary_directory

#######################################################################

=back

=cut

1;
__END__

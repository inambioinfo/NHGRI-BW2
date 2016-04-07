package NHGRI::BW2::Session;
#######################################################################

=head1 NAME

Session.pm - A Perl module to handle transfer of files and commands to 
          and from the CIT "biowulf2" cluster, as well as handle some 
          job control.

=head1 DESCRIPTION

  This module can create connections to CIT's biowulf2 server, transfer 
  files to and from the user's account directories, and handle some job 
  submission and monitoring.

  In order to use this module, the user must have a working account on 
  biowulf2, and have set up ssh authentication using an empty 
  passphrase.

=head1 DATE

 September 30, 2015

=head1 AUTHORS

 Nancy Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=over

=cut

#######################################################################
require 5.005_62;
use strict;
use warnings;

use IPC::Open2;
use Net::SSH qw(sshopen3);
use POSIX ":sys_wait_h";
use File::Temp qw/ tempdir /;
use File::Spec;

our $VERSION  = '0.01';
our $REVISION = '$Id:$';

use vars qw( $HOSTNAME $CPU_LIMIT $DATADIR $HOMEDIR $SCRATCHDIR );

$HOSTNAME = 'biowulf2.nih.gov';
$DATADIR = '/data'; # will be followed by a username
$SCRATCHDIR = '/scratch'; # will be followed by a username
$HOMEDIR = '/home'; # will be followed by a username

#######################################################################

=item new()

  This method creates a new Session object.

  Input:  -hostname can pass the name of a beowulf login node (cluster 
              must be running slurm), but if not specified, the default 
              value of $HOSTNAME is used.
          -username can pass the username to be used on the remote 
              host.  By default, the username of the user running the 
              script is used.
  Output: New Session object

=cut

#######################################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $hostname = $params{-hostname} || $HOSTNAME;
    my $username = $params{-username} || (getpwuid($<))[0] || '';
 
    my $self = { hostname => $hostname,
                 username => $username };

    bless $self, $class;
    return $self;

} ## end new

#######################################################################

=item hostname()

  The method gets or sets the value of "hostname", 
  which is a string specifying the name of the cluster
  login node.

  Input: Optional argument sets value.
  Output: Value of hostname (scalar string)

=cut

#######################################################################
sub hostname {
    my $self  = shift;

    if (defined (my $new_hostname = shift))
    {
        $self->{hostname} = $new_hostname;
    }

    return $self->{hostname};

} ## end hostname

#######################################################################

=item username()

  The method gets or sets the value of "username", 
  which is a string specifying a username on the cluster.

  Input: Optional argument sets value.
  Output: Value of username (scalar string)

=cut

#######################################################################
sub username {
    my $self  = shift;

    if (defined (my $new_username = shift))
    {
        $self->{username} = $new_username;
    }

    return $self->{username};

} ## end username

#######################################################################

=item cpu_limit()

  The method returns the scalar constant $CPU_LIMIT,
  which is the maximum number of cpus we can request
  on the biowulf cluster.

  Input: None.
  Output: Value of $CPU_LIMIT (scalar number)

=cut

#######################################################################
sub cpu_limit {

    return $CPU_LIMIT;

} ## end cpu_limit

#######################################################################

=item data_dir()

  The method returns the scalar constant $DATADIR, which is the path
  to the "data" directory on biowulf2.

  Input: None.
  Output: Value of $DATADIR (scalar string)

=cut

#######################################################################
sub data_dir {

    return $DATADIR;

} ## end data_dir

#######################################################################

=item scratch_dir()

  The method returns the scalar constant $SCRATCHDIR, which is the path
  to the "scratch" directory on biowulf2.

  Input: None.
  Output: Value of $SCRATCHDIR (scalar string)

=cut

#######################################################################
sub scratch_dir {

    return $SCRATCHDIR;

} ## end scratch_dir

#######################################################################

=item transfer_directory_to_biowulf()

  The method uses scp to copy a directory to the biowulf host.

  Input: -localpath specifies the path to the directory on the local 
          host.  With or without a trailing "/", this directory and its
          contents will be transferred to the biowulf host.
         -remotepath specifies the path to a directory on the biowulf
          host.  If the directory exists, the new directory will be 
          placed inside this directory with a name that is the same as 
          the local directory specified in "-localpath".  If the remote
          directory does not exist, the local directory will be copied
          to a new directory with the path specified by "-remotepath".
         -as_files, if true, will cause the files in the localpath 
          directory to be deposited within the remote directory specified
          by remotepath, even if remotepath already exists (if false, a
          subdirectory with the local directory's name will be created
          within that directory).
         -md5sum_check if true, will perform an md5sum check by
          calculating md5sum values on localhost and biowulf and 
          comparing.  If they are not the same, program will die.
          WARNING: md5sum_check is not yet implemented for directories!!!
  Output: 1 if successful, 0 otherwise

=cut

#######################################################################
sub transfer_directory_to_biowulf {
    my $self  = shift;
    my %params = @_;

    my $localpath = $params{-localpath};
    if ($params{-as_files}) {
        $localpath =~ s:/{0,1}$:/*:;
    }
    else {
        $localpath =~ s:/$::;
    }
    my $remotepath = $params{-remotepath};
    
    my $remotehost = $self->hostname();
    my $username = $self->username();

    my $command = "scp -r $localpath $username\@$remotehost:$remotepath 2> /dev/null";

    unless ($params{-quiet})
    {
        print "$command\n";
    }

    my $result = system("$command");

    if ($result)
    {
        die "Something went wrong transferring $localpath to $remotehost\n";
    }

    ####TODO####
    #if ($params{-md5sum_check}) {
        #$self->check_md5sums_from_transfer(-localpath => $localpath, 
                                           #-remotepath => $remotepath);
    #}

} ## end transfer_directory_to_biowulf

#######################################################################

=item transfer_directory_from_biowulf()

  The method uses scp to copy a directory from the biowulf host to the 
  localhost.

  Input: -remotepath specifies the path to the directory on the remote
          host.  With or without a trailing "/", this directory and its
          contents will be transferred to the local host.
         -localpath specifies the path to a directory on the local host.
          If the directory exists, the new directory will be placed inside 
          this directory with a name that is the same as the remote
          directory specified in "-remotepath".  If the local directory
          does not exist, the remote directory will be copied to a new 
          directory with the path specified by "-localpath".
         -as_files, if true, will cause the files in the remotepath 
          directory to be deposited within the local directory specified
          by localpath, even if localpath already exists (if false, a
          subdirectory with the remote directory's name will be created
          within that directory).
         -md5sum_check if true, will perform an md5sum check by calculating 
          md5sum values on localhost and biowulf and comparing.  If they are 
          not the same, program will die.
          WARNING: md5sum_check is not yet implemented for directories!!!
         
  Output: 1 if successful, 0 otherwise

=cut

#######################################################################
sub transfer_directory_from_biowulf {
    my $self  = shift;
    my %params = @_;

    my $localpath = $params{-localpath};
    if ($params{-as_files}) {
        $localpath =~ s:/{0,1}$:/*:;
    }
    else {
        $localpath =~ s:/$::;
    }
    my $remotepath = $params{-remotepath};
    
    my $remotehost = $self->hostname();
    my $username = $self->username();

    my $command = "scp -r $username\@$remotehost:$remotepath $localpath 2> /dev/null";

    unless ($params{-quiet})
    {
        print "$command\n";
    }

    my $result = system("$command");

    if ($result)
    {
        die "Something went wrong transferring $remotepath from $remotehost\n";
    }
    
    ####TODO####
    #if ($params{-md5sum_check}) {
        #$self->check_md5sums_from_transfer(-localpath => $localpath, 
                                           #-remotepath => $remotepath);
    #}

} ## end transfer_directory_from_biowulf

#######################################################################

=item transfer_file_to_biowulf()

  The method uses scp to copy a file to the biowulf host.

  Input: -localpath specifies the path to the file on the local host
         -remotepath specifies the path to the file on the biowulf host
         -md5sum_check if true, will perform an md5sum check by
          calculating md5sum values on localhost and biowulf and 
          comparing.  If they are not the same, program will die.
  Output: 1 if successful, 0 otherwise

=cut

#######################################################################
sub transfer_file_to_biowulf {
    my $self  = shift;
    my %params = @_;

    my $localpath = $params{-localpath};
    my $remotepath = $params{-remotepath};
    
    my $remotehost = $self->hostname();
    my $username = $self->username();

    my $command = "scp $localpath $username\@$remotehost:$remotepath 2> /dev/null";

    unless ($params{-quiet})
    {
        print "$command\n";
    }

    my $result = system("$command");

    if ($result)
    {
        die "Something went wrong transferring $localpath to $remotehost\n";
    }

    if ($params{-md5sum_check}) {
        $self->check_md5sums_from_transfer(-localpath => $localpath, 
                                           -remotepath => $remotepath);
    }

} ## end transfer_file_to_biowulf

#######################################################################

=item transfer_file_from_biowulf()

  The method uses scp to copy a file from the biowulf host to the 
  localhost.

  Input: -remotepath specifies the path to the file on the biowulf host
         -localpath specifies the (new) path to the file on the local host
         -md5sum_check if true, will perform an md5sum check by
          calculating md5sum values on localhost and biowulf and 
          comparing.  If they are not the same, program will die.
         
  Output: 1 if successful, 0 otherwise

=cut

#######################################################################
sub transfer_file_from_biowulf {
    my $self  = shift;
    my %params = @_;

    my $localpath = $params{-localpath};
    my $remotepath = $params{-remotepath};
    
    my $remotehost = $self->hostname();
    my $username = $self->username();

    my $command = "scp $username\@$remotehost:$remotepath $localpath 2> /dev/null";

    unless ($params{-quiet})
    {
        print "$command\n";
    }

    my $result = system("$command");

    if ($result)
    {
        die "Something went wrong transferring $remotepath from $remotehost\n";
    }
    
    if ($params{-md5sum_check}) {
        $self->check_md5sums_from_transfer(-localpath => $localpath, 
                                           -remotepath => $remotepath);
    }

} ## end transfer_file_from_biowulf

#######################################################################

=item transfer_files_to_biowulf_directory()

  The method uses scp to copy a set of files to a directory on the biowulf 
  host.  If the remote directory does not exist on biowulf, it will be
  created.

  Input: -localfiles is a reference to a list of file paths on the local 
          host
         -remotepath specifies the path to the directory on the biowulf
          host (will be created with specified mode if it does not exist)
         -mode specifies the mode with which a new directory will be 
          created on biowulf, if necessary (default 775)
         -md5sum_check if true, will perform an md5sum check by
          calculating md5sum values on localhost and biowulf and 
          comparing.  If they are not the same, program will die.
          ###WARNING - md5sum_check is not yet implemented for this
          method!!!!
  Output: 1 if successful, 0 otherwise

=cut

#######################################################################
sub transfer_files_to_biowulf_directory {
    my $self  = shift;
    my %params = @_;
    my $mode = $params{'-mode'} || 775;

    my $ra_localfiles = $params{-localfiles};
    my $remotepath = $params{-remotepath};
    $remotepath =~ s:/$::;

    my $run_string = "mkdir -p -m $mode $remotepath";
    $self->remote_command(-command => $run_string);
    
    # create a local directory with symlinks

    my $localdir = tempdir( CLEANUP => 1 );
    foreach my $localpath (@{$ra_localfiles}) {
        my $basename = $localpath;
        $basename =~ s:.*/::;
        $localpath = File::Spec->rel2abs( $localpath );
        symlink $localpath, "$localdir/$basename";
    }
    
    my $remotehost = $self->hostname();
    my $username = $self->username();

    my $command = "scp $localdir/* $username\@$remotehost:$remotepath/ 2> /dev/null";

    unless ($params{-quiet})
    {
        print "$command\n";
    }

    my $result = system("$command");

    if ($result)
    {
        die "Something went wrong transferring $localdir to $remotehost\n";
    }

    #####TODO#####
    #if ($params{-md5sum_check}) {
        #$self->check_md5sums_from_transfer(-localpath => $localpath, 
                                           #-remotepath => $remotepath);
    #}

} ## end transfer_files_to_biowulf_directory

#######################################################################

=item check_md5sums_from_transfer()

  The method calculates md5sums for the local and the remote copies
  of a file, compares them, and returns 0 if they are the same.  If
  the md5sums differ, the program dies.

  Input: -localpath passes the path to the file on the localhost
         -remotepath passes the path to the file on biowulf2
  Output: 0 if the two md5sums match, dies otherwise.

=cut

#######################################################################
sub check_md5sums_from_transfer {
    my $self  = shift;
    my %params = @_;
    my $localpath = $params{'-localpath'};
    my $remotepath = $params{'-remotepath'};

    my $local_md5sum = `md5sum $localpath`;
    $local_md5sum =~ s/\n.*//s; # only first line
    $local_md5sum =~ s/\s.*//s; # only first field
    my $session = $self->session();
    my $remote_md5sum = $session->remote_command(-command => "md5sum $remotepath");
    $remote_md5sum =~ s/\n.*//s;
    $remote_md5sum =~ s/\s.*//s;
    chomp $remote_md5sum;

    print "$localpath vs. $remotepath:\n$local_md5sum\n$remote_md5sum\n";

    if ($local_md5sum eq $remote_md5sum) {
        return 0;
    }
    else {
        die "md5sum values for local file $localpath and remote file $remotepath do not match!(local $local_md5sum remote $remote_md5sum)\n";
    }

} ## end check_md5sums_from_transfer

#######################################################################

=item make_user_data_directory()

  The method uses ssh to execute a command on the BW2 host to create a 
  subdirectory in the user's "data" directory.

  Input: -dirname specifies the name of the subdirectory 
         -mode gives optional mode value (default 775)
  Output: The full path (on the biowulf2 server) of the new directory.

=cut

#######################################################################
sub make_user_data_directory {
    my $self  = shift;
    my %params = @_;
    my $subdir_name = $params{'-dirname'};
    my $mode = $params{'-mode'} || '775';

    my $remotehost = $self->hostname();
    my $username = $self->username();
    my $newdir = "$DATADIR/$username/$subdir_name";

    print "Creating directory $newdir on $remotehost.\n";

    my $run_string = "mkdir -p -m $mode $newdir";

    unless ($params{-quiet})
    {
        print "$run_string\n";
    }

    my $response = $self->remote_command(-command => $run_string);

    return $newdir;

} ## end make_user_data_directory

#######################################################################

=item make_user_scratch_directory()

  The method uses ssh to execute a command on the BW2 host to create a 
  subdirectory in the user's "scratch" directory.

  Input: -dirname specifies the name of the subdirectory 
         -mode gives optional mode value (default 775)
  Output: The full path (on the biowulf2 server) of the new directory.

=cut

#######################################################################
sub make_user_scratch_directory {
    my $self  = shift;
    my %params = @_;
    my $subdir_name = $params{'-dirname'};
    my $mode = $params{'-mode'} || '775';

    my $remotehost = $self->hostname();
    my $username = $self->username();
    my $newdir = "$SCRATCHDIR/$username/$subdir_name";

    print "Creating directory $newdir on $remotehost.\n";

    my $run_string = "mkdir -p -m $mode $newdir";

    unless ($params{-quiet})
    {
        print "$run_string\n";
    }

    my $response = $self->remote_command( -command => $run_string);

    return $newdir;

} ## end make_user_scratch_directory

#######################################################################

=item make_remote_directory()

  The method uses ssh to execute a command on the BW2 host to create a 
  directory on the remote host.

  Input: -dirname specifies the path to the new directory
         -mode gives optional mode value (default 775)
  Output: The full path (on the biowulf2 server) of the new directory.

=cut

#######################################################################
sub make_remote_directory {
    my $self  = shift;
    my %params = @_;
    my $directory = $params{'-dirname'};
    my $mode = $params{'-mode'} || '775';

    my $remotehost = $self->hostname();
    my $username = $self->username();

    print "Creating directory $directory on $remotehost.\n";

    my $run_string = "mkdir -p -m $mode $directory";

    unless ($params{-quiet})
    {
        print "$run_string\n";
    }

    my $response = $self->remote_command(-command => $run_string);

    return $directory;

} ## end make_remote_directory

#######################################################################


=item remote_command()

  The method executes the command passed in the argument
  on the remote host, and returns the output.

  Input: -command passes the command to be executed
  Output: The output of the command when executed on the
      remote server (does not currently include STDERR).

=cut

#######################################################################
sub remote_command {
    my $self  = shift;
    my %params = @_;
    my $cmd = $params{'-command'};

    print "In remote_command with value $cmd!\n";

    my $remotehost = $self->hostname();
    my $username = $self->username();

    unless ($params{-quiet})
    {
        print "Executing $cmd\n";
    }

    my ($writer, $reader, $error);
    my $pid = sshopen3("$username\@$remotehost", *WRITER, *READER, *ERROR, $cmd);
    my $response = '';

    while(<READER>) {
        $response .= $_;
    }

    close READER;

    while(<ERROR>) {
        #$response .= $_;
    }

    close ERROR;
    close WRITER;

    waitpid($pid, 0); # will this kill it?

    return ($response) ? $response : '';

} ## end remote_command

#######################################################################

=item move_biowulf_file()

  The method moves a source to destination file on the biowulf2 server.

  Input: -source passes the path to the source file on biowulf2
         -destination passes the path to the destination file on biowulf2
  Output: The output of the command when executed on the
      remote server.

=cut

#######################################################################
sub move_biowulf_file {
    my $self  = shift;
    my %params = @_;
    my $source = $params{'-source'};
    my $destination = $params{'-destination'};

    my $response = $self->remote_command( -command => "mv $source $destination" );

    return ($response) ? $response : '';

} ## end move_biowulf_file

#######################################################################

=item remote_content()

  The method cats the remote file given in the argument
  and returns the return string.

  Input: -file specifies the remote path to the file
  Output: The output of the cat command.

=cut

#######################################################################
sub remote_content {
    my $self = shift;
    my %params = @_;
    my $remote_file = $params{'-file'};

    unless ($params{-quiet})
    {
        print "Executing cat $remote_file.\n";
    }

    my $return_string = $self->remote_command(-command => "cat $remote_file");

    return $return_string;
    
} ## end remote_content

#######################################################################

=item remove_data_directory()

  This method removes the remote data subdirectory 
  specified in the argument.

  Input: -dirname specifies the subdirectory name to remove
  Output: 1 if successful.

=cut

#######################################################################
sub remove_data_directory {
    my $self  = shift;
    my %params = @_;
    my $subdir = $params{'-dirname'};
    my $username = $self->username();
    my $remotehost = $self->hostname();
    $subdir = "$DATADIR/$username/$subdir";

    unless ($params{'-quiet'}) {
        print "Removing data subdirectory $subdir on $remotehost.\n";
    }

    my $run_string = "rm -rf $subdir";

    my $response = $self->remote_command( -command => $run_string);

    return ($response) ? $response : '';

} ## remove_data_directory

#######################################################################

=back

=cut

1;
__END__

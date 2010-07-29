# $Id: Process.pm,v 1.1.1.1 2004/10/07 19:44:13 cegrant Exp $
# Paul Pavlidis 2000
# Functions for managing processes. Specifically, for handling background processes, errors in same, and cancelations by users.

package Process;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
	     cancelJob
	     forkOK
	     isDone
	     readSemaphore
	     setUpSpawn
	     writeSemaphore
	     );

use strict;
use MetaGlobals;
use CGIutil;

#-----------------------------------------------------------------
# End a background process which was canceled by the user
#-----------------------------------------------------------------
sub cancelJob($$) {
  my ($pid, $uid) = @_;
  my $runningModPerl = $ENV{MOD_PERL};

  # Make sure the pid we are handed agrees with the one the user's job is running under.
  warn "Trying to kill $pid for $uid. I am process $$, $ENV{MOD_PERL}\n";
  
  my $sema = readSemaphore($uid);
  if ($sema) {
    $sema =~ m/^([0-9]+)\s.*?([0-9]+)\s/;
    my $userspid = $2;
    die "Hey! That isn't the PID I expected! ($userspid != $pid, $sema)\n" unless $userspid == $pid;

  # Be very paranoid about the pid. (this might be redundant)
    sanity($pid);
    $pid =~ m/^([0-9]+)$/;
    $pid = $1;
    
    # Is it already dead?
    if (kill 0 => $pid) {
      warn "$pid is still alive\n";
    } else {
      warn "$pid is already dead, probably of natural causes.\n";
      return 1;
    }
    
  # Try to kill it.
    warn "Looks okay to kill. Sending 'QUIT' to $pid\n";
    kill 'QUIT' => $pid; # we've registered a QUIT handler to 'die'.
    
    # Did it die? Wait briefly to make sure change registers. (This is debug code).
    sleep 1; 
    if ( kill 0 => $pid ) {
      warn "Seems I failed to kill $pid: $!\n";
      return 0;
    } else {
      warn "Seems I have successfully killed $pid (or it died of natural causes)\n";
      return 1;
    }
  } else {
    return 0;
  }
}  # sub cancelJob


#-----------------------------------------------------------------
# Figure out if it is okay to fork on this system.
#-----------------------------------------------------------------
sub forkOK() {
  # Fork is okay under unix but not windows unless perl is v5.6
  # (supposedly) We don't want to fork at all under modperl
  # (performance killer), but use another trick (post_connection)
  my $runningModPerl = $ENV{MOD_PERL};
  my $forkOK = ( $^O !~ /MSWin32/i && !$runningModPerl) ? 1 : 0; # don't try to fork under windows. It don't work (?)
  return $forkOK;
} # sub forkOK


#-----------------------------------------------------------------
# Figure out if a given job is complete
#------------------------------------------------------------------
sub isDone($$$) {
  my ($pid, $uid, $sema) = @_;

#  warn "isDone got PID:$pid UID:$uid $sema\n";

  # be optimistic: is it done? (note that the order of these is very
  # important!)
  if ($sema && $sema =~ /DONE/) {
    return 'finished';
  } elsif ($sema && $sema =~ /ABORT/) {
    return 'aborted';
  } elsif ($sema && $sema =~ /RUN/) {
    return 'running';
  } elsif ($sema && $sema =~ /WAITING/) {
    return 'waiting';
  } else {
    # not in a readily recognizable state.
  }

  # Be very paranoid about the pid. (this might be redundant)
  sanity($pid);
  $pid =~ m/^([0-9]+)$/;
  $pid = $1;

  # get a pid that was really the one that we should be using, as a check.
  $sema =~ m/^([0-9]+)\s([0-9]+)\s/;
  my $userspid = $2;
  die "Hey! That isn't the PID I expected! ($userspid != $pid)\n" unless $userspid == $pid;

  my $runningModPerl = $ENV{MOD_PERL};  

  # Hmm. Not done. Is the process still alive? (if not, we have a potential problem)
  if (! kill 0 => $pid) {
    warn "process $pid died\n";
    return 'died';
  } else {
    return 'running';
  }
} # sub isDone


#-------------------------------------------------------------------
# Read line from semaphore file.
#-------------------------------------------------------------------
sub readSemaphore($) {
  my ($uid) = @_;
  my $semaphoreFile =  "$UPLOAD_DIR/$uid/${uid}.$SEMAPHORE_SUFFIX";
#  print STDERR "Reading from $semaphoreFile...\n";
  my $sema;
  if (-e $semaphoreFile) {
    eval {  open (SEMA, "<$semaphoreFile");
	    flock(SEMA, LOCK_SH() );
	    $sema = <SEMA>;
	    close SEMA;
#	    print STDERR "Read '$sema' from $semaphoreFile\n";
	  };
    die "Error reading semaphore file: $@\n" if $@;
  } else {
    warn "$semaphoreFile does not exist\n";
  }
  return $sema;
} # sub readSemaphore


#-------------------------------------------------------------------
# Set up environment so we can keep track of background processes
#-------------------------------------------------------------------
sub setUpSpawn($) {
  my ($forkOK) = @_;
  # establish some details of how we deal with errant children, etc.
  if ($forkOK) { # what to do when certain signals are received.
    $SIG{CHLD} = 'IGNORE'; # Autoreap zombies. child process sends this signal to the parent when it quits. See pg 415 of camel 3rd ed.
    $SIG{QUIT} = sub { die "\nAction canceled by user.\n" };
    setpgrp(0, 0); # processid (0=current), process group (default) # not be platform-independent.
  }
} # sub setUpSpawn


#-------------------------------------------------------------------
# Write something to a semaphore file.
#-------------------------------------------------------------------
sub writeSemaphore($$$) {
  my ($uid, $text, $overwrite) = @_;
  $overwrite = $overwrite ? ">" : ">>";
  my $semaphoreFile =  "$UPLOAD_DIR/$uid/${uid}.$SEMAPHORE_SUFFIX";
#  print STDERR "writing to $semaphoreFile\n";
  eval {  open (SEMA, "$overwrite$semaphoreFile");
	  flock(SEMA, LOCK_EX() );
	  $|++;
	  print SEMA $text;
	  close SEMA;
	};
  die "Error reading semaphore file: $@\n" if $@;
  chmod (0666, $semaphoreFile);
#  print STDERR "Successfully wrote $text to $semaphoreFile\n";
} # sub writeSemaphore


1;

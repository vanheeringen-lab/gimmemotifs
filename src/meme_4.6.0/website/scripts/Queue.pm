# FILE: Queue.pm
# PROJECT: svm web server
# AUTHOR: Paul Pavlidis, based on code by Phan Lu, Darrin Lewis, Andrew Liu, Ilan Wapinski.
# CREATED: 12/02
# $Id: Queue.pm,v 1.1.1.1 2004/10/07 19:44:13 cegrant Exp $
package Queue;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(enqueue dequeue);

use MetaGlobals;
use Log;
use Process;
use strict;

############################################################################
##### QUEUEING THE JOB. Place the job in the queue. Only returns when
##### the queue is free or we give up.
############################################################################
sub enqueue {
  my ($uid, $log) = @_;

  # create a queue file
  my $queue_file = "$QUEUE_DIR/${uid}.queue";
  open (Q, ">$queue_file") || die("Cannot open queue file: $!");
  print Q time, "\n";
  $log->log("Enqueued $uid");
  $log->debug("Enqueued $uid");
  close Q;

  my $t = 0;
  my %ages;
  my $old_pos = 10000000000;
  $log->log("Waiting...");
  writeSemaphore($uid, "WAITING|", 0);

  cleanqueue();

  while ( $t < $PATIENCE ) {
    $^T = time; # what happens to this?
    undef %ages;
    opendir( QDIR, "$QUEUE_DIR" );
    while ( defined ($_ = readdir(QDIR))) {
      if (/queue$/) {
	my $age = -M "$QUEUE_DIR/$_";
	$ages{$_} = $age;
      }
    }
    closedir(QDIR);

    my $qfile;
    my $pos = 0;
    foreach $qfile (keys %ages) {
      if ( $ages{$qfile} > $ages{"${uid}.queue"} ) {
	$pos++; # how many are ahead of me.
      }
    }

    # sanity check
    if ($MAXJOBS < $SIMULTANEOUSJOBS) {
      return("ABORT|CONFIGERROR");
      $log->log("Notify the administrator that the MAXJOBS settings is too low given the SIMULTANEOUSJOBS setting");
    }

    # process job: ABORT or RUN or WAIT
    if ( $pos > $MAXJOBS ) {
      $log->log("Sorry, too many jobs waiting now. Try later."); # Should we do this? What does it hurt to have a long queue?
      return( "ABORT|MAXJOBS" );
    } else {
      if ( $pos < $SIMULTANEOUSJOBS ) {
	$log->log("It's your turn.");
	return( "RUN" );
      } else {
	if ($pos < $old_pos) {
	  $log->log("$pos jobs are ahead of yours."); # todo: show how many are actually running. Look for old jobs.
	}
	$old_pos = $pos;
	sleep $WAIT;
	$t += $WAIT;
      }
    }
  } # while waiting

  if ( $t > $PATIENCE ) {
    $log->log("Could not run the job. Ran out of patience.");
    return( "ABORT|PATIENCE" );
  }
}

############################################################################
# Clean up some 'zombies' (internal function)
############################################################################
sub cleanqueue {
  # At least one job in the queue should be running. If not, kill the oldest one.
  # To figure out which one is running is a bit complicated....maybe too tricky.
  return 1;
}

############################################################################
##### DEQUEUEING THE JOB
############################################################################
sub dequeue {
  my($uid, $log) = @_;
  my $success = 0;
  if (-e "$QUEUE_DIR/${uid}.queue") {
    $success = unlink "$QUEUE_DIR/${uid}.queue";
  } else {
    # Can happen if it already was cleaned up.
    #    $log->debug("Tried to dequeue $QUEUE_DIR/${uid}.queue but the file is not there");
    $success = 1;
  }

  if ($success < 1) {
    $log->log("Dequeue could not remove the queue file. Please contact the system administrator!!!");
  } else {
    $log->log("Dequeued $uid");
  }
}

1;

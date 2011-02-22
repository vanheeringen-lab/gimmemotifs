#!/bin/csh

set now=`date +%s`
if ($#argv == 1) then
  set start=$1
  @ elapsed = $now - $start
  @ ds = $elapsed % 60
  @ dm = ($elapsed / 60) % 60
  @ dh = $elapsed / 3600
  if ($dh > 0) then
    printf "%d hour(s), %d minute(s) and %d second(s)" $dh $dm $ds
  else if ($dm > 0) then
    printf "%d minute(s) and %d second(s)" $dm $ds
  else
    printf "%d second(s)" $ds
  endif
else
  echo $now
endif


#!INSERTCSHCOMMAND   
#   Bootstrap
#   Copyright (C) 2000 Christopher J. Basten 
# Usage governed by the terms of the 
# GNU  General Public License,  version 2 or higher
#  See the file COPYING in this directory
#
#   This file was meant as an example.  You  will need to edit it 
#   to work on your particular system with your data files.
#
#  Start by setting the variables needed.  
# 
set stem=qtlcart                        #  filename stem
set hypo=30                             #  hypothesis for SSupdate
set model=3                             #  analysis model
set reps=1000                           #  number of bootstraps
set email=basten\@statgen.ncsu.edu      #  email address for notice
set templog=temp.log                    #  temporary log file
set qbin=/home/basten/qtlcart/bin       #  where are the QTL Cart binaries
set bin=/usr/bin                        #  where are the system programs
#
$bin/rm -f $templog
echo "Bootstrap experiment started " > $templog
$bin/date >>  $templog
$bin/echo "Stem: " $stem >> $templog
$bin/echo "Model: " $model >> $templog
$bin/echo "Reps: " $reps >> $templog
$bin/echo "Email: " $email >> $templog
$bin/mv $stem.log $stem.logsave
$bin/rm -f $stem.z${model}a
$qbin/SSupdate.pl -I $hypo < $stem.z > $stem.z$model.boot
$bin/mv $stem.z $stem.zsave
set i=1
while ( $i <= $reps )
$qbin/Prune -A -V -i $stem.cro -b 1  >>&  $templog
$bin/nice $qbin/Zmapqtl -A -V -M $model -i $stem.crb >>&  $templog
$qbin/SSupdate.pl -I $hypo -f $stem.z$model.boot < $stem.z > $stem.z$model.new
$bin/mv $stem.z$model.new $stem.z$model.boot
$bin/rm $stem.z
@ i++
end
$qbin/SSupdate.pl -I $hypo -c -f $stem.z$model.boot > $stem.z$model.booted
$bin/mv $stem.logsave $stem.log
$bin/mv $stem.zsave $stem.z
$bin/echo "Bootstrap experiment ended " >> $templog
$bin/date >>  $templog
/usr/ucb/mail $email <  $templog

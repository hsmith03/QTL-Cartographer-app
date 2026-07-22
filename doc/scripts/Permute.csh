#!INSERTCSHCOMMAND 
#           Permute 
#   Copyright (C) 2000 Christopher J. Basten 
# Usage governed by the terms of the GNU  General Public License,  version 2 or higher
#  See the file COPYING in this directory
#
#   This file was meant as an example.  You  will need to edit it 
#   to work on your particular system with your data files.
#
#  Start by setting the variables needed.  
# 
set stem=qtlcart                        #  filename stem
set column=4                            #  LR column to process
set model=3                             #  analysis model
set reps=1000                           #  number of bootstraps
set email=basten\@statgen.ncsu.edu      #  email address for notice
set templog=temp.log                    #  temporary log file
set qbin=/home/basten/qtlcart/bin       #  where are the QTL Cart binaries
set bin=/usr/bin                        #  where are the system programs
#
#   Should only need to change what is above.
#
$bin/rm -f $templog
echo "Permutation test started " > $templog
$bin/date >>  $templog
$bin/echo "Stem: " $stem >> $templog
$bin/echo "Model: " $model >> $templog
$bin/echo "Reps: " $reps >> $templog
$bin/echo "Email: " $email >> $templog
$bin/mv $stem.log $stem.logsave
$bin/mv $stem.z $stem.zsave
$bin/rm -f $stem.z{$model}e
set i=1
$qbin/GetMaxLR.pl -i < $stem.zsave > $stem.z{$model}.ewt
$qbin/CWTupdate.pl -C $column < $stem.zsave > $stem.z{$model}.cwt
while ( $i <= $reps )
$qbin/Prune -A -V -i $stem.cro -b 2  >>&  $templog
$bin/nice $qbin/Zmapqtl -A -V -M $model -i $stem.crb  >>&  $templog
$qbin/GetMaxLR.pl -r $i -C $column < $stem.z  >> $stem.z{$model}.ewt
$qbin/CWTupdate.pl -f $stem.z{$model}.cwt -C $column < $stem.z  > $stem.z{$model}.newcwt
$bin/mv $stem.z{$model}.newcwt $stem.z{$model}.cwt 
$bin/rm -f $stem.z
@ i++
end
$bin/echo "Now your can run EWThreshold.pl on $stem.z$model.ewt" >> $templog
$bin/mv $stem.zsave $stem.z
$bin/mv $stem.logsave $stem.log
$bin/date >>  $templog
/usr/ucb/mail $email <  $templog

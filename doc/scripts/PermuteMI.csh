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
set column=5                            #  LR column to process
set ic=6                                #  Information criterion code
set model=SmprtSec                      #  analysis model
set reps=1000                           #  number of bootstraps
set email=basten\@statgen.ncsu.edu      #  email address for notice
set templog=temp.log                    #  temporary log file
set qbin=/home/basten/qtlcart/bin       #  where are the QTL Cart binaries
set bin=/usr/bin                        #  where are the system programs
#
#   Should only need to change what is above.
#
$bin/rm -f $templog
$bin/echo "Permutation test started " > $templog
$bin/date >>  $templog
$bin/echo "Stem: " $stem >> $templog
$bin/echo "Work: " $model >> $templog
$bin/echo "IC: " $ic >> $templog
$bin/echo "Reps: " $reps >> $templog
$bin/echo "Email: " $email >> $templog
$bin/mv $stem.log $stem.logsave
$bin/mv -f $stem.mim $stem.mim.save
set i=1
$bin/echo "Experimentwise test maxima " > $stem.mim.ewt
$bin/date >>  $stem.mim.ewt
$bin/echo "Stem: " $stem >> $stem.mim.ewt
$bin/echo "Work: " $model >> $stem.mim.ewt
$bin/echo "IC: " $ic >> $stem.mim.ewt
$bin/echo "Reps: " $reps >> $stem.mim.ewt
$bin/echo "Email: " $email >> $stem.mim.ewt
$bin/echo "-start" >> $stem.mim.ewt
while ( $i <= $reps )
$qbin/Prune -A -V -i $stem.cro -b 2  >>&  $templog
$bin/nice $qbin/MImapqtl -A -V -I $model -S $ic -i $stem.crb  >>&  $templog
$qbin/GetMaxLR.pl -r $i -C $column < $stem.mim  >> $stem.mim.ewt
$bin/rm -f $stem.mim
@ i++
end
$bin/echo "-end" >> $stem.mim.ewt
$bin/echo "Now your can run EWThreshold.pl on $stem.mim.ewt" >> $templog
$bin/mv $stem.mim.save $stem.mim
$bin/mv $stem.logsave $stem.log
$bin/date >>  $templog
/usr/ucb/mail $email <  $templog

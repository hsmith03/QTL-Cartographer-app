#!INSERTCSHCOMMAND
#           PermuteJZ
#   Copyright (C) 2001 Christopher J. Basten 
# Usage governed by the terms of the GNU  General Public License,  
# version 2 or higher
#  See the file COPYING in this directory
#
#   This file was meant as an example.  You  will need to edit it 
#   to work on your particular system with your data files.
#
#  Start by setting the variables needed.   These were set up for
#  MacOX:  your paths may be different.   
# 
set stem=qtlcart                        #  filename stem
set column=5                            #  LR column to process
set reps=1000                           #  number of permutations
set traits=16                           #  Number of traits in data
set email=basten\@statgen.ncsu.edu      #  email address for notice
set templog=temp.log                    #  temporary log file
set qbin=/Users/basten/bin              #  where are the QTL Cart binaries
set bin=/bin                            #  where are the system programs
#            check this for the commands date, echo, mv, rm and nice
set mail=/usr/bin/mail                  #  the mail command
#                       #  I like /usr/ucb/mail under Solaris. 
# --------------- ------------------- -------------------- --------------
#   Should only need to change what is above.
#
$bin/rm -f $templog                     #  remove old temp log file
$bin/echo "Permutation test started " > $templog   # init. new one
$bin/date >>  $templog                  #  put in the date and params
$bin/echo "Stem: " $stem >> $templog
$bin/echo "Reps: " $reps >> $templog
$bin/echo "Email: " $email >> $templog
$bin/mv $stem.log $stem.logsave         #  save the qtlcart.log file
set j=0                                 # save copies of all the
while ( $j <= $traits )                 # original JZmapqtl output files 
$bin/mv -f $stem.z$j $stem.z$j.save
@ j++
end
set i=1                                 # i counts permutations
$bin/echo "Experimentwise test maxima " > $stem.z0.ewt  # init permutation file
$bin/date >>  $stem.z0.ewt            #  with the date and parameters
$bin/echo "Stem: " $stem >> $stem.z0.ewt
$bin/echo "Reps: " $reps >> $stem.z0.ewt
$bin/echo "Email: " $email >> $stem.z0.ewt
$bin/echo -n "-start" >> $stem.z0.ewt    # This is needed for EWThreshold.pl
while ( $i <= $reps )                    #  main loop.
$qbin/Prune -A -V -i $stem.cro -b 6  >>&  $templog       # Permute
$qbin/JZmapqtl -A -V -t 0   -i $stem.crb  >>&  $templog  # Analyze +traits 
$qbin/GetMaxLR.pl -r $i -C $column < $stem.z0  >> $stem.z0.ewt # Get max. 
set j=0
while ( $j <= $traits )          #  clean up JZmapqtl files
$bin/rm -f  $stem.z$j 
@ j++
end
@ i++
end
set j=0       #  move originals back
while ( $j <= $traits ) 
$bin/mv -f $stem.z$j.save $stem.z$j 
@ j++
end
$bin/echo " " >> $stem.z0.ewt
$bin/echo "-end" >> $stem.z0.ewt
$bin/echo "Now your can run EWThreshold.pl on $stem.z0.ewt" >> $templog
$bin/mv $stem.logsave $stem.log   # move qtlcart.log file back
$bin/date >>  $templog
$mail $email <  $templog  #   notify end of job.

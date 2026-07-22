#!INSERTCSHCOMMAND
#      copyright 2001 Chris Basten
#  Script to test QTL Cart. programs
#  
#   You may need to change some of these to match your system.
#
set bindir=$home/bin
set basedir=$home/qtltest
set delete=/bin/rm
set copy=/bin/cp
set exampledir=$home/qtlcart/QTLCartUnix/example
set log=$basedir/log.txt
#
#  create main directory
#
if ( ! -d $basedir ) mkdir $basedir     
touch $log
#
#   do the tests for the sample.* files
#  
if ( ! -d $basedir/sample ) mkdir $basedir/sample     
cd $basedir/sample
$delete *
$copy $exampledir/sample.* .
$bindir/Rmap -i sample.mps -A     >>  $log
$bindir/Rcross -i sample.raw -A >>  $log
$bindir/Qstats -A >>  $log
$bindir/LRmapqtl -A >>  $log
$bindir/SRmapqtl -A >>  $log
$bindir/Zmapqtl -A >> $log
$bindir/Zmapqtl -A -M 6 >>  $log
$bindir/MImapqtl -A >>  $log
$bindir/Eqtl -A >>  $log
$bindir/Preplot -A >>  $log
#
#   do the tests for the realdat.* files
#  
if ( ! -d $basedir/realdat ) mkdir $basedir/realdat     
cd $basedir/realdat
$delete *
$copy $exampledir/realdat*.* .
$bindir/Rmap -i realdatm.inp -A >>  $log
$bindir/Rcross -i realdatc.inp -A >>  $log
$bindir/Qstats -A >>  $log
$bindir/LRmapqtl -A >> $log
$bindir/SRmapqtl -A >> $log
$bindir/Zmapqtl -A >> $log
$bindir/Zmapqtl -A -M 6 >>  $log
$bindir/MImapqtl -A >> $log
$bindir/Eqtl -A >> $log
$bindir/Preplot -A >> $log
#
#   do the tests for the mletest.* files
#  
if ( ! -d $basedir/mletest ) mkdir $basedir/mletest
cd $basedir/mletest
$delete *
$copy $exampledir/mletest.* .
$bindir/Qstats -X mletest -A >> $log
$bindir/LRmapqtl -A >> $log
$bindir/SRmapqtl -A >> $log
$bindir/Zmapqtl -A >> $log
$bindir/Zmapqtl -A -M 6 >>  $log
$bindir/MImapqtl -A >> $log
$bindir/Eqtl -A >> $log
$bindir/Preplot -A >> $log
#
#   do the tests for the simulations  
#  
if ( ! -d $basedir/sim ) mkdir $basedir/sim 
cd $basedir/sim   
$bindir/Rmap   -A >> $log
$bindir/Rqtl -q 4 -A >> $log
$bindir/Rcross   -A >> $log
$bindir/Qstats -A >> $log
$bindir/LRmapqtl -A >> $log
$bindir/SRmapqtl -A >> $log
$bindir/Zmapqtl -A >> $log
$bindir/Zmapqtl -A -M 6 >>  $log
$bindir/MImapqtl -A >> $log
$bindir/Eqtl -A >> $log
$bindir/Preplot -A >> $log

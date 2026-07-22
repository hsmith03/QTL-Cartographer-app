#!INSERTCSHCOMMAND  
#
#  Run Model 8 iteration
#
#   Usage:
#   Model8  bindir  stem   siglevel  iterations  max_nbp
#     bindir is the binary subdirectory
#     stem is the filename stem
#     siglevel  is the significance level to declare a QTL
#     iterations is the number of iterations
#     max_nbp is the maximal number of background parameters.
#
if ( $1 == '-h' ) then
echo "    Usage:  Model8.csh bindir  stem   siglevel  iterations  max_nbp"
echo "Where"
echo "        bindir  = QTL Cart. binary directory"
echo "          stem  = filename stem"
echo "      siglevel  = Significance level to declare a QTL"
echo "    iterations  = number of iterations"
echo "       max_nbp  = maximal number of background parameters"
echo "    hypothesis  = hypothesis test for Eqtl"
echo " "
echo "Now exiting"
exit
endif
set bindir=$1
set stem=$2
set siglevel=$3
set iterations=$4
set maxnbp=$5
set ihypo=$6
$bindir/Qstats -X $stem -A -V
$bindir/Zmapqtl -A -V -M 3
$bindir/Eqtl -A -V -S $siglevel -I Z -H $ihypo
#
#  Save the original files
#
/usr/bin/mv $stem.eqt $stem.eqt.0
/usr/bin/mv $stem.z $stem.z.0
/usr/bin/cp $stem.sr $stem.sr.0
#
#  Use model 8 iteratively with cofactors from previous run.
#
set i=1
while ( $i < $iterations )
echo "Doing iteration $i"
$bindir/Zmapqtl -A -V -M 8 -n $maxnbp
/usr/bin/rm $stem.sr 
$bindir/Eqtl -A -V -S $siglevel -I Z -H $ihypo
/usr/bin/cp $stem.sr $stem.sr.$i
/usr/bin/mv $stem.eqt $stem.eqt.$i
/usr/bin/mv $stem.z $stem.z.$i
set j=$i
@ j--
$bindir/SRcompare.pl -f $stem.sr.$j < $stem.sr.$i
@ i++
end
/usr/bin/rm $stem.sr 
echo "Finished"

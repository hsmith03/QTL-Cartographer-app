REM
REM Dos batchfile for doing QTL analysis
REM
Rmap -A -V -X realdat -i realdatm.inp -W ..\example
Rcross -A -V -i realdatc.inp
Qstats -A -V
LRmapqtl -A -V
SRmapqtl -A -V -M 2
Zmapqtl -A -V
Zmapqtl -A -V -M 6
Eqtl -A

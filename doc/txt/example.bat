rem  
rem  This is a little batch file that you can run
rem  in a command window to set up an example.
rem
rem  edit this line so that your path to
rem  the QTL Cartographer binaries is correct.
rem
path c:\home\QTLCartWin\bin;%path%
rem
rem  make sure that there is a root directory for
rem  examples.   Here we use c:\home\test
rem
cd c:\home\test
rem
rem  You should have given a command line option
rem  it will create a directory and copy the files to it
rem
mkdir %1
cd %1
rem  
rem  The directory here should be updated to
rem  reflect where the example files reside.
rem
copy c:\home\QTLCartWin\example\%1.* .



# qtlcart
QTL Cartographer Public Version 1.17f



Updated on 19 September 2017

What follows is the old README file. It is obsolete.  Code is now in github and you can 
indicate issues there.  Contact Chris Basten directly at christopher.basten@gmail.com.  




23 March 2004 


                        README file for QTL Cartographer.


QTL Cartographer is a package of programs that will aid in locating the
genes that control quantitative traits using a molecular map of
markers.  It includes some programs to allow simulation studies of
experiments.

Where to get QTL Cartographer:

Via the web, point your browser to ftp://statgen.ncsu.edu/pub/qtlcart
and shift-click on a file to download it.   This probably doesn't work anymore.
If you are reading this, you need to get the source from github at
https://github.com/cbasten/qtlcart 

Otherwise, ftp to statgen.ncsu.edu (152.14.14.17) with

username: ftp
password:  your email address

cd to /pub/qtlcart

ftp> cd /pub/qtlcart


     This is what's in the directory:
----------------------------------------------------------------
----------------------------------------------------------------
Filename                  Explanation
----------------------------------------------------------------
QTLCartWin.zip           32 bit Windows NT/9x version  
QTLCartUnix.tar.Z        Unix distribution
QTLCartMac.sea.hqx       Macintosh distribution
README                   This file
gnuplot-dist/            Distributions of GNUPLOT
   Macgnuplot.zip        Macintosh version 3.7.1
   gnuplot-3.8h.0.tar.gz Unix version 3.8
   gp373w32.zip          Windows version 3.7.3
1.10b/                   Version 1.10b
1.12f/                   Version 1.12f
1.13g/                   Version 1.13g
1.14d/                   Version 1.14d
1.15d/                   Version 1.15d
1.16c/                   Version 1.16c
1.17f/                   Version 1.17f
WQtlSetup.exe            Shengu Wang's GUI front-end for Windows
data/                    Published data sets
----------------------------------------------------------------
----------------------------------------------------------------


****  For MS-Windows 3.1, 3.11 (look in 1.12f/)

Download the files LHARC.EXE and QTLCART.LZH in binary format to your 
DOS machine.  Create a subdirectory for the binaries (say, c:\qtlcart).
Move the files to that subdirectory and uncompress the archive 
QTLCART.LZH with the following  command:

c:\qtlcart\> LHARC x QTLCART

The ten binaries will be uncompressed. 

Now load Windows and start up the File Manager.  Create a program 
group, say "QTL Cartographer", and drag the *.exe files from the File
manager to the new program group.  They will then be ready to use 
with a double click.  You will be presented with a menu of options
to select.

****  For MS-Windows NT/9x:

Download QTLCartWin.zip into a subdirectory and uncompress it by double clicking the
archive.   You may need to download an unzip utility.   


You can now use the UNIX version of QTL Cartographer by first installing the
freeware package Cygwin.   Cygwin is available at http://www.cygwin.com.
If you install the complete distribution, then you will have a gcc compiler
and perl:  This will allow you to run the scripts as well.   Chris highly
recommends that you install and use Cygwin.   The QTL Cartographer binaries
run just as fast within Cygwin as the natively compiled binaries for Windows.  

****  For Macintosh:

PowerMac binaries are in the file QTLCartMac.sea.hqx.  
You will need to unbinhex the download and then double click the
self-extracting binary.   The folder 'bin' has the PowerMac binaries.
The folder 'cbin' has the Carbon binaries for use with MacOS X.
The 'doc' and 'example' folders have documentation and example files.

If you have a Macintosh running MacOS X, you can also download the UNIX
version and install it.   Chris uses the UNIX version on his Macintosh
exclusively.   You can use either GNUPLOT distribution:  The Macintosh 
carbon binaries or the compiled version that uses Aquaterm.   Also, there
is a site that organizes free software for Macintosh OS X users called
fink (http://fink.sourceforge.net).   For OS X users, it is a better way
to get GNUPLOT, R and other software packages.   Mapmaker is also available
through fink, though you have to tweak the fink configuration to allow
for the install of unstable distributions.  Explore the fink site for more
information.  

****  For UNIX:
 
Download the file QTLCartUnix.tar.gz in binary form from statgen.ncsu.edu.
It is in the same directory that this file came from.  Anonymous ftp into
statgen.ncsu.edu and cd to /pub/qtlcart.  On your local machine, create
a subdirectory for the distribution, then move the file QTLCartUnix.tar.gz to
it.  Uncompress and untar the file as follows:

% gunzip QTLCartUnix.tar.gz
% tar xf QTLCartUnix.tar

See the file doc/txt/INSTALL on how to install the programs.

Adobe PDF files containing the documentation are also in the doc/pdf
subdirectory.   


In Chris' opinion, the UNIX version is the de facto standard.   It can be used
under most UNIX and Linux distributions as well as Macintosh OS X and
Windows (with Cygwin installed).   In the future, it is possible that the
binary distributions will be discontinued in favor of a single UNIX distribution.
(The Windows front end produced by Shengchu Wang will still be available as a 
separate program).  Your comments are welcome.

****  Notes:

Remember that the programs all act like Unix programs, and you can get
a list of command line options by invoking the program with the
argument -h.  An example:

brooks:~ % Rmap -h

======================================================================
        QTL Cartographer v. 1.14d, 3/31/2000
        Copyright (C) 1996-1999 C. J. Basten, B. S. Weir and Z.-B. Zeng.
        QTL Cartographer comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions. For details select help.

======================================================================
TIME:    14:19:11 on Wednesday, 21 June 2000
PROGRAM: Rmap
PURPOSE: Create or translate a genetic linkage map
USAGE:   Rmap  [-i  x]  [-o  x]  [-e  x] ...
DEFAULTS:
  [ -i                  ] Input File
  [ -o      qtlcart.map ] Output File
  [ -e      qtlcart.log ] Error File
  [ -s        961611551 ] Random Number Seed
  [ -f                1 ] Map Function [1,8], 1 => Haldane
  [ -p         0.000000 ] Map Function Parameter
  [ -g                1 ] Ouput (1,2,3) => Text, Graphics, Both
  [ -c                4 ] Chromosomes
  [ -m               16 ] Markers per Chromosome
  [-vm         0.000000 ] Standard Deviation of Markers per Chromosome
  [ -d        10.000000 ] Intermarker Distance (cM)
  [-vd         0.000000 ] Standard Deviation of Intermarker Distance (cM)
  [ -t         0.000000 ] Tails (Flanking DNA, in cM)
  [ -M                0 ] Simulation Mode (0,1)

  Also: [-h] for help, [-A] for automatic,  [-V] for non-Verbose
  [-W path] for a working directory, [-R file] to specify a resource
  file and [-X stem] to specify a filename stem.

======================================================================
 Now exiting program without doing any calculations.
  You may need to click on the upper right X to close this screen
======================================================================


The manual  gives a better explanation of the options. A PDF version is in
the manual.pdf file in the doc/pdf subdirectory.  



**Please join the mailing list for QTL Cartographer.  It will be a 
forum for problems you may have in using the programs, and we will
post announcements of updates and bug fixes.  To subscribe, send the
following one line message to MajorDomo@statgen.ncsu.edu:

subscribe qtlcart


**Send any bug reports to qtlcart-bug@statgen.ncsu.edu.  

***************************************************************** 

              Chris Basten        B. S. Weir      Zhao-Bang Zeng

Email:  basten@statgen.ncsu.edu             zeng@statgen.ncsu.edu
Phone:       (919)515-1934                         (919)515-1942
Fax:                           (919)515-7315
                               
Address:
                        Department of Statistics
                     Bioinformatics Research Center
                     North Carolina State University
                         Raleigh, NC 27695-7566
                                  USA

***************************************************************** 


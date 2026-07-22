#!/usr/bin/perl
# 
#  Copyright (C) 2002 Christopher J. Basten 
#  Usage governed by the terms of the GNU  General Public License,  version 2 or higher
#  See the file COPYING in this directory
# 
#  All of the perl scripts will be attached to this file prior to installation.
#  You can make some global changes here.
#  Make sure that your perl command is correct on the first line.  It may be
#  /usr/local/bin/perl or something else.  
#  If you don't have Getopt or it
#  doesn't work on your system, then you have some hacking to do.
#
use Getopt::Std;                     # Do not change this line 
local ($opt_i,$opt_o,$opt_f,$opt_c, $opt_h,$opt_s) = ("u","d","",0,0,0);
$usage = "$0: [-i u||d||m ] [-o u||d||m]   [-h]   < input > output\n";
getopts "i:o:f:chs:" ;
die $usage if ( $opt_h == 1 );

#  First, define the Input_Recod_Separator
$/ = "\r"   if ( $opt_i eq "m"  );    # Mac line-end is \r
$/ = "\n"   if ( $opt_i eq "u"  );    # Unix line-end is \n
$/ = "\r\n" if ( $opt_i eq "d"  );    # Dos line-end is \r\n

# Next, define the Output_Record_Separator  
$\ = "\r"   if ( $opt_o eq "m"  );    # Mac line-end is \r
$\ = "\n"   if ( $opt_o eq "u"  );    # Unix line-end is \n
$\ = "\r\n" if ( $opt_o eq "d"  );    # Dos line-end is \r\n

$skipcode = "XCutXCodeXSkip";
$unskipcode = "XCutXCodeXUnSkip";
$sourcecode = "\t\tCopyright (C) 1994-2005\n\tNorth Carolina State University" if ( $opt_s == 1 ) ;
$sourcecode = "\t\tCopyright (C) 1999-2005\n\tLauren McIntyre and Jun Wu" if ( $opt_s == 2 ) ;
$sourcecode = "\t\tCopyright (C) 2000-2005\n\tPatrick Gaffney and Brian Yandell" if ( $opt_s == 3 ) ;

$header = "/* ------------------------------------------------------ $skipcode
     This file ($opt_f) is part of QTL Cartographer
         
    $sourcecode

  For more information, look at  http://statgen.ncsu.edu/qtlcart or
  email Chris Basten (basten\@statgen.ncsu.edu).   The web site has
  a link to download different versions of the programs.
  
  QTL Cartographer is free software; you can redistribute it
  and/or modify it under the terms of the GNU  General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.
  The GNU General Public License is also available online
  (http://www.gnu.org/licenses/gpl.html).
  
  QTL Cartographer is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied
  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public
  License along with QTL Cartographer; see the file COPYING.  If
  not, write to the Free Software Foundation, Inc., 675 Mass Ave,
  Cambridge, MA 02139, USA.
------------------------------------------------------ $unskipcode */
";

if ($opt_s == 0 ) {


$header = "/* ------------------------------------------------------ $skipcode
     This file ($opt_f) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ $unskipcode */
";


}



$trailer = "
/* ------------------------------------------------------- $skipcode
             End of file $opt_f
------------------------------------------------------------------ */
";

print $header if $opt_c == 1 ;

$skip = 0;
while (<>) {
  chomp;
  $skip = 1 if /$skipcode/ ;
  print if $skip == 0;
  $skip = 0 if /$unskipcode/ ;
}

print $trailer   ;


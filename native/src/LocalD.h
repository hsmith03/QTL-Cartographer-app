/* ------------------------------------------------------ XCutXCodeXSkip
     This file (LocalD.h) is part of QTL Cartographer
         
    		Copyright (C) 1994-2005
	North Carolina State University

  For more information, look at  http://statgen.ncsu.edu/qtlcart or
  email Chris Basten (basten@statgen.ncsu.edu).   The web site has
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
------------------------------------------------------ XCutXCodeXUnSkip */

/* Local definitions for the compiler */

   /*  You can change this to float if you want floating point numbers to be floats rather than doubles.
       This can also be defined on the command line with -DFPN=double or -DFPN=float */
#ifndef FPN
#define FPN double
#endif

    /*Define this if you need the sign() function*/
#ifndef DSIGN
#define DSIGN 1
#endif

    /*Define this if you need div_t, ldiv_t*/
/*
#ifndef DIVT
#define DIVT 1 
#endif
*/
    /*Define this if you need  strlwr()  or strupr()*/
/*
#ifndef ITOA
#define ITOA 1 
#endif
*/
    /* 0 => no debugging output, 1 => memory allocation output */
#ifndef DEBUGGING
#define DEBUGGING 0
#endif

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

/* 
   Define UNIX for the preprocessor if UNIX.
   If this is not UNIX, then it should be Codewarrior.
   
   If using Codewarrior, 
       Define MACWARRIOR for  Macintosh  ( when __MACOS__ is defined)
              WINWARRIOR for  Windows    ( when __MACOS__ is not defined)

*/

#ifdef UNIX  

#define FILESEP '/'

#else

#ifdef  __MACOS__  

#define MACWARRIOR 1 
#define FILESEP ':'

#else

#define WINWARRIOR 1 
#define FILESEP '\\'

#endif
 
#endif

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file LocalD.h
------------------------------------------------------------------ */


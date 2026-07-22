/* ------------------------------------------------------ XCutXCodeXSkip
     This file (mygenome.h) is part of QTL Cartographer
         
    		Copyright (C) 2000-2005
	Patrick Gaffney and Brian Yandell

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

#ifndef MYGENOME_GUARD
#define MYGENOME_GUARD

typedef struct MyGenome { /*Structure to hold a genome defined by the genetic linkage map*/
  int chrom;                 /* Chromosome of the marker interval */
  int markr;                 /* The marker that precedes the interval, can be 0 if borders are allowed */
  double dist;               /* This distance, from marker markr to markr+1 on chromosome chrm, is in M */
  double pos;                /* Position of marker from left telomere in Morgans*/
  int* genotype;             /* marker genotype */
  struct MyGenome *prev;            /* Pointer to previous marker interval */
  struct MyGenome *next;             /* Pointer to next marker interval */
}  mygenome;



#endif


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file mygenome.h
------------------------------------------------------------------ */


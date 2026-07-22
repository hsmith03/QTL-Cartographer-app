/* ------------------------------------------------------ XCutXCodeXSkip
     This file (chrom.h) is part of QTL Cartographer
         
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

#ifndef CHROM_GUARD
#define CHROM_GUARD
/**************************************************************
  File:         chrom.h
  Written by:   Patrick Gaffney
  Modified by:  Amy Jin
  Date:         November 11, 2000/Aug 13,2002
  Version:      0.4

  Purpose:
  -------

    Header file describing chromosome records in MCMC program. 



**************************************************************/


/*#include "Main.h"*/
#include "mygenome.h"


typedef struct chromosomeData CHROMOSOME; 
typedef struct MCQtlData QTL_INFO;


struct chromosomeData
{
   int num;
   int nQtl;
   int nMark;             /* number of markers in the chromosome */
   FPN chromLen;       /* length of the chromosome (in Morgans) */

   mygenome** mark_genome;  /* array of pointers to the genome entries for markers */
   FPN* mark_pos;
  
   QTL_INFO** qtls; 

   FPN** log_prob;

};

#endif





/* ------------------------------------------------------- XCutXCodeXSkip
             End of file chrom.h
------------------------------------------------------------------ */


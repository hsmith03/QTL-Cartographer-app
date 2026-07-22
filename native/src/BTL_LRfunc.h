/* ------------------------------------------------------ XCutXCodeXSkip
     This file (BTL_LRfunc.h) is part of QTL Cartographer
         
    		Copyright (C) 1999-2005
	Lauren McIntyre and Jun Wu

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

/*  
          BTL_LRfunc.h 
          
Created by Jun Wu , Biological Science in Purdue University
This file is part of QTL Cartographer. 
*/




int BTL_calc_lrstats(params *theparams,   individual *individs, int trait, FPN **beta0, FPN **beta1, FPN **fstat, genome *gptr);
void BTL_trans_data(int **dataptr,markermap *themap);
void BTL_print_lrstats(int nn, int trait, markermap *themap, char *minfile, char *iinfile, char *outfile);
int Gauss_elim(FPN **Matrix,  FPN *Ycolumn,int size, FPN *Betta);
int BTL_calc_F2(params *theparams,  individual *individs, int trait, FPN **beta0, FPN **beta1, FPN **beta2,FPN **fstat,FPN **fstatx1,FPN **fstatx2, genome *gptr);
void BTL_print_S(int nn,markermap *themap, FILE *outf, FPN **beta0, FPN **beta1,FPN **beta2,FPN **Rmq, FPN p1,FPN p2,FPN **p3, FPN **fstat, FPN **fstatx1, FPN **fstatx2, int cross);
int calc_Btl(int nn, char *outfile,markermap *themap,int cross,FPN **beta0,FPN **beta1,FPN **beta2,FPN **Rmq,FPN **p3,params *theparams, FPN **fstat, FPN **fstatx1, FPN **fstatx2);
FPN Abs(FPN value);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file BTL_LRfunc.h
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (QSfunc.h) is part of QTL Cartographer
         
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



int move_phenotypes(FPN *trptr, individual *individs, int trait, int nn);
void calc_qstats(FPN *tptr, int nind, FPN *thestats, int kk);

void print_qstats(FPN *stats, params *theparams, int trait, markermap *themap);
void  marker_segregation(individual *individs,markermap *themap,params *theparams);
void miss_mark_summary(individual *individs, markermap *mkptr, params *theparams, int trait);
void ind_data_summary(individual *individs, markermap *mkptr, params *theparams);
void CalcMarkProbs(params *theparams);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file QSfunc.h
------------------------------------------------------------------ */


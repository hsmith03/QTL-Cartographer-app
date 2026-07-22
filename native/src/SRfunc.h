/* ------------------------------------------------------ XCutXCodeXSkip
     This file (SRfunc.h) is part of QTL Cartographer
         
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




FPN init_xmyv(FPN **xm, FPN *yv, genome *gptr, genome *lgptr, int n, individual *individs, markermap *themap, params *theparams);
genome *forward_swr(individual *individs, params *theparams, markermap *themap,
genome *first, linpakws *lnpk, int nbp, int explan, int steps);
genome *backward_swr(individual *individs, params *theparams, markermap *themap, genome *first, linpakws *lnpk, int nbp, int explan, int steps);
genome *for_back_swr(individual *individs, params *theparams, markermap *themap, genome *first, linpakws *lnpk, int nbp, int explan, int steps);

void write_sr_results(genome *agptr, params *theparams, char *srfile, int n, int explan);
int calc_fact(params *theparams);
int calc_samplesize(params *theparams, individual *individs);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file SRfunc.h
------------------------------------------------------------------ */


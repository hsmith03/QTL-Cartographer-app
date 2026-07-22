/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RMfunc.h) is part of QTL Cartographer
         
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




/* These functions create a map of molecular markers. */
markermap *initailize_markermap(int kk, params *theparams);
int mkr_num(int ll,FPN ssigl);
void mrkr_rec(markermap *themap);
void gnuplot_map(params *theparams, markermap *themap);
void mrkr_rec2(markermap *themap);




/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RMfunc.h
------------------------------------------------------------------ */


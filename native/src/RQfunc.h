/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RQfunc.h) is part of QTL Cartographer
         
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






void print_aqtl(params *theparams, char *outfile, int boot);
genome *create_a_qtl(genome *first, aqtl *qtlptr, params *theparams);
FPN qtl_effects(FPN beta, int job);
FPN qtl_dom_effects(FPN b1, FPN b2, int job);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RQfunc.h
------------------------------------------------------------------ */


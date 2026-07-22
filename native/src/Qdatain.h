/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Qdatain.h) is part of QTL Cartographer
         
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


void GetTheModel(params *theparams, char *infile);
void get_knums(char *inputf, markermap *themap);
int get_qtls(aqtl *qtlptr, char *inputf);
int get_knum(char *inputf);
aqtl *get_std_qtls(markermap *themap, char *infile);
aqtl *read_qtls(FILE *fptr, int isu, markermap *themap, int isn);
int movetoqtls(FILE *infile);
void read_interactions(FILE *fptr, markermap *themap, aqtl *qtlptr);


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Qdatain.h
------------------------------------------------------------------ */


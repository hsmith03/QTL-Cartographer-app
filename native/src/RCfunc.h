/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RCfunc.h) is part of QTL Cartographer
         
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




/*   Header files for RCfunc.c:   routines to create a random data set. */
thelines *which_line(int whch, thelines *lines);
thelines *alloc_aline(char *name, individual *iptr, int nn, struct aline *prev, struct aline *next, int which);
void free_aline(thelines *ptr);
void show_lines(thelines *ptr);

void init_pop(individual *p1p2f1, params *theparams);

void calc_phenotypes(individual *indptr, int nn, params *theparams);
int genotype(params *theparams, int mgt, int pgt, FPN mgamete, FPN pgamete);
void create_an_individual(params *theparams, individual *fptr, individual *mptr, individual *optr, genome *genptr);
FPN *calc_genotype(individual *indptr);
FPN calc_genvar(individual *indptr, int nn, int trait);
void assign_indicators(int gt, FPN *xi, FPN *zi);
FPN *calc_genotype2(individual *indptr);

void recombination(genome *genptr,params *theparams );
void do_a_cross(params *theparams, thelines *p1ptr, thelines *p2ptr, thelines *p3ptr, genome *gptr,   int s);
 
void expand_map(params *theparams, genome *gptr, aqtl *qtlptr, markermap *themap);
FPN convert_dist(FPN dist, int flag, int ss);
int sc_genotype(  int mgt, int pgt, FPN mgamete, FPN pgamete);
int first_gam(int gt);
int second_gam(int gt);

void    calc_parental_diffs(markermap *themap,aqtl *theqtls);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RCfunc.h
------------------------------------------------------------------ */


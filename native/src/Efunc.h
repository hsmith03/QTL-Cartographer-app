/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Efunc.h) is part of QTL Cartographer
         
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

void ReEstimateMap(params *theparams  );
void RCDstage1(params *theparams,  char *chptr, char *progname) ;
void RCDstage2(params *theparams,  div_t emethod) ;
void  check_segregation(individual *individs,markermap *themap,params *theparams,genome *gptr);
FPN EstimateRecomb(int l1, int l2,individual *individs,params *theparams,genome **gvect,FPN *lr);

int addlinkgroup(int chrom,int glen,genome **gvect,params *theparams);
int initlinkgroup(int chrom,int glen,genome **gvect,params *theparams);
int dolinkgroup(int chrom,int glen,genome **gvect,params *theparams) ;

genome *reorder_genome(genome *thegenome,params *theparams,int glen,genome **gvecto);
FPN obfuncsar(int chrom,int glen,genome **gvect);
FPN obfuncsal(int chrom,int glen,genome **gvect);
FPN obj_func(params *theparams, int chrom,int glen,genome **gvect);
int FirstLevelMod(params *theparams,markermap *themap,genome **gvect,int glen,individual *thedata,int chrom); 
int SecondLevelMod(params *theparams,markermap *themap,genome **gvect,int glen,individual *thedata,int chrom); 
int ThirdLevelMod(params *theparams,markermap *themap,genome **gvect,int glen,individual *thedata,int chrom); 

void update_recombs(genome *g1ptr, genome *g2ptr,individual *thedata, params *theparams, genome **gvect, int chrom ) ;
void SwapNodeDataUp(individual *thedata,params *theparams, genome *bg1, genome *bg2, markermap *themap); 
void SwapNodesData(individual *thedata,params *theparams, genome *bg1, genome *bg2, markermap *themap);
void SwapNodeDataDown(individual *thedata,params *theparams, genome *bg1, genome *bg2, markermap *themap);
FPN fxlike(int cross,int generation, int l1, int l2, FPN *dcnts, FPN r) ;
FPN bxlike(int generation,FPN *dcnts,  FPN  r);
void m1dotm2(FPN **m1,FPN **m2,FPN **m3 , int ub);
FPN estFXrec(int cross,int generation,int *cnts,FPN *r,int l1, int l2);
void update_chromosome_rec(params *theparams,genome **gvect,int glen,individual *thedata,int chrom);

int GvectNode(genome **gvect,int glen,int chrom, int marker);
FPN Mapobj(params *theparams, int glen,genome **gvect);
void update_interval(genome *lptr, individual *thedata, params *theparams, genome **gvect, int chrom);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Efunc.h
------------------------------------------------------------------ */


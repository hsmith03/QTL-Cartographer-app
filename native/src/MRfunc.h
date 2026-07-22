/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MRfunc.h) is part of QTL Cartographer
         
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


void  AppendAtoB(char *filea,char *fileb);
void EstimateEffects(genome *first,linpakws *lnpk, params *theparams, int trait) ;
long     *GetPositions(char *infile,int ad );
void      GetMRParams(char *infile, FPN *walk, int *otraits, int *traits, int *n);
FPN  **GetTraitData(char *infile, int ntraits, int n,char **names);
int  **GetOTraitData(char *infile, int ntraits, int n,char **names);
void      GetSiteExpecteds(char *infile,int n, genome *gptr,long *positions, int whichsite);
genome   *multiregress(char *infile,  params *theparams, linpakws *lnpk,long *positions,long *dpositions, int trait);
void ShowQTLSinp(genome *first,params *theparams,char *outfile, int i, char *tnames,int where);
void ShowMultiRegress(genome *first,int *pstatus,params *theparams,char *outfile,int i,char *tnames);
void ShowMRParams(params *theparams, char *outfile, linpakws *lnpk,char **otnames );
void     update_pstatus(int *pstatus,int *chromosomes,FPN *locations,long *positions,genome *first,params *theparams) ;
void  GetChromLocales(int *chromosomes,FPN *locations,long *positions,char *infile);
int SetUpMatrix(linpakws *lnpk, params *theparams);
FPN TestThisModel(genome *aptr, genome *dptr, genome *first,linpakws *lnpk, params *theparams,int trait,FPN syy);


void ShowSiteExpecteds(genome *gptr) ;

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MRfunc.h
------------------------------------------------------------------ */


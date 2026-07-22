/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MZfunc.h) is part of QTL Cartographer
         
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



int **jzget_srresults(params *theparams, markermap *themap);
void show_cofactors(FILE *outfile,params *theparams,markermap *themap,genome *first,int **srranks);
int set_cofactors(params *theparams,markermap *themap,genome *first,int **srranks );

void write_position_results(linpakws *lnpk, params *theparams, char *outfile, genome *gptr, FPN abspos, markermap *themap);
void write_jzheader(char *outfile, params *theparams, char *onamae, char *chptr, markermap *themap, int oc,genome *first,int **srranks);
void write_alttraits(char *outfile,params *theparams, markermap *themap,individual *individs, linpakws *lnpk);
long write_altheader(char *outfile,params *theparams,char *onamae, char *chptr,markermap *themap,int oc);
int    do_jzexpecteds(params *theparams,genome *startptr,genome *endptr,markermap *themap, individual *individs,char *outfile, linpakws *lnpk);
FPN expected_int(linpakws *lnpk,params *theparams,markermap *themap,individual *individs,genome *gptr,FILE *outf,int *ipos, FPN excess);
void insert_positions(char *outfile,int positions,long fp);

int add_marker(params *theparams, markermap *themap, individual *individs, linpakws *lnpk, genome *gptr, int row1, int ad, int oti);
int determine_endpoints(params *theparams, markermap *themap, genome **start, genome **end, genome *first, int *do_analysis);
FPN do_jzanalysis(params *theparams, genome *startptr, genome *endptr, markermap *themap, aqtl *theqtls, individual *individs, char *outfile, genome *agptr, linpakws *lnpk );
int     ecm_solve(linpakws *lnpk, int nn, int pp, int ihypo, params *theparams);
FPN   fit_null(linpakws *lnpk, params *theparams, individual *individs, markermap *themap, aqtl *theqtls , genome *first, int *nrows);
int how_many_rows(params *theparams, markermap *themap, individual *individs, int *otr,int marks);
void init_xsave(params *theparams, markermap *themap, individual *individs, linpakws *lnpk);
FPN invdet(FPN **a, FPN **b, int lda, int nt, FPN *z, int *kpvt, int job);
void pick_cofactors(linpakws *lnpk, markermap *themap, aqtl *theqtls, params *theparams, genome *tgptr, genome *agptr );
FPN zmap_int(linpakws *lnpk, params *theparams, markermap *themap, individual *individs, genome *gptr, char *outfile,  int nrows);


int AssignGeneticParameters(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd);
void UnDoTheTrick(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int trait);
void DoTheTrick(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int trait);
void  CreateVarCovar(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd); 
int CalculateTraitlnL(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int trait,int *go_on);
void CalculateJointlnL(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd);
int DoEstep(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int *ktime,FPN s2);
void GxEParameters(linpakws *lnpk,  int cross,   params *theparams,int aa,int dd);



void do_scan(params *theparams,markermap *themap, aqtl *theqtls,individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr,int **srranks);
void do_gbye_same(params *theparams,markermap *themap, aqtl *theqtls,individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr,int **srranks);
void do_gbye_diff(params *theparams,markermap *themap, aqtl *theqtls,individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr,int **srranks);
void do_pleiotropy(params *theparams,linpakws *lnpk,char *progname,char *chptr,int **srranks);
void prep_multiregress(params *theparams,markermap *themap,  individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr);
void write_pleiotropyheader(char *outfile,params *theparams,markermap *themap,genome *first,int **srranks);
FPN pleiotropy_int(linpakws *lnpk,params *theparams,markermap *themap,individual *individs,genome *gptr,char *outfile,int nrows);
FPN PleiotropyCloseLinkage(linpakws *lnpk,params *theparams,markermap *themap,individual *individs,genome *jptr,genome *gptr,char *outfile,int nrows);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MZfunc.h
------------------------------------------------------------------ */


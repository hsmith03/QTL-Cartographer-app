/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Zfunc.h) is part of QTL Cartographer
         
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


int get_srresults(params *theparams, markermap *themap, genome *first);
int get_lratio(params *theparams, FPN *lratio);
int get_CWTR(params *theparams, char *infile, int *pcnts, FPN *lratio);
int get_bootnum(char *infile);
FPN **get_lin_reg(char *lrfile, int chroms, int maxl, int t);

void write_zheader(char *outfile, params *theparams, char *chptr, int eorc, int boots, markermap *themap);
void write_CWTR(params *theparams,genome *startptr, genome *endptr, markermap *themap, FPN deltax, char *outfile, FPN *lratio, int *pcnts, int reps);
void write_maxlr(char *outfile, int jj, FPN *maxlr, int chroms);
void write_zline(params *theparams, char *outfile, genome *gptr, FPN abspos, FPN *estimates);
long write_bootline(params *theparams, char *outfile, char *infile, genome *gptr, FPN abspos, FPN *estimates, long int inset);


int    convert_zfile(params *theparams, char *outfile);
FPN do_zanalysis(params *theparams, genome *startptr, genome *endptr, markermap *themap, aqtl *theqtls, FPN **lin_reg, FPN *yy, individual *individs, char *outfile, FPN *lratio, int *cnts, genome *agptr, linpakws *lnpk);
int    em_solve(linpakws *lnpk, FPN *y, int nn, int pp, int ihypo, params *theparams);
FPN   fit_back_params(params *theparams, individual *individs, markermap *themap, aqtl *theqtls, FPN **xw, FPN *yy, FPN *qraux, FPN *rsd, int chrom, int markr, int **bp, int nbp, genome *gptr);
genome *delete_gnode(genome *gptr);
void   init_phenotypes(params *theparams,individual *individs,FPN **yy,int *samplesize);
void pick_markers(int **bp, markermap *themap, aqtl *theqtls, FPN **lin_reg, params *theparams, genome *tgptr, genome *agptr);

void zmap(params *theparams, linpakws *lnpk, FPN *yy, individual *individs, genome *gptr, int nbp, FPN thetaL, FPN thetaR, FPN ts0, int roffset);




/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Zfunc.h
------------------------------------------------------------------ */


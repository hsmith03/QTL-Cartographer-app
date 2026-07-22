/* ------------------------------------------------------ XCutXCodeXSkip
     This file (EQfunc.h) is part of QTL Cartographer
         
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


typedef struct ZmapqtlResults {
	int chrom;          /* Chromosome of the peak */
	int trait;          /* Trait for this estimate */
	int model;          /* Model for the analysis */
	FPN window;      /* Window size, Models 5-7 */
	int nbp;            /* Number of background parameters in model*/
	long rank;          /* Rank of the node relative to the others */
	long start_offset;  /* File offset...start reading here */
	long end_offset;    /* File offset...end reading here...Not used yet, and not initialized */
	struct ZmapqtlResults *prev;
	struct ZmapqtlResults *next;
}  zresult;

typedef struct JZmapqtlResults {   /*Structure to hold jz results in Eqtl*/
  int chrom;
  int traits;
  FPN position;
  FPN lr10;    /*Likelihood ratios for Hi:Hj  i = 1,2,3,G  j = 0,1,2,E*/
  FPN lr20;
  FPN lr30;
  FPN lr31;
  FPN lr32;
  FPN lrGxE;
  FPN **gparams; /* = dmatrix(1,4,0,traits)    0 is for average, then one per trait. */
                    /* rows are   
                              1  a1
                              2  a3
                              3  d3
                              4  d2   */
  struct JZmapqtlResults *prev;
  struct JZmapqtlResults *next;
  int which;             /* Which qtl are we on? */
}  jzresult;

typedef struct AnalysisResults {
  char *filename;
  long offset;
  int rows;
  int cols;  /*  number of columns for real valued numbers. */
  FPN window;
  int nbp;
  int model;
  int ihypo;
  int trait;   /*  The trait for a single trait analysis, or the number in a multiple trait analysis. */
  int filetype;   /*  integer value for the file type.  */
  int **cm;  /* chromosome, marker matrix.   = imatrix(1,rows,1,2) */
  FPN **results;  /*  results, order is filetype dependant.   = matrix(1,rows,1,19) */
  char **headers;   /*  Column Headers */
  struct AnalysisResults *next;
  struct AnalysisResults *prev;
}  aresults; 



zresult *zresult_rank(zresult *first, int flag);
zresult *zmapqtl_list(params *theparams, int ranker);
zresult *zresult_node(int flag, int c, int t, int m, long int start, long int end, zresult *celle, zresult *prev, zresult *next);
zresult *zresult_elim_nonmodel(zresult *first, params *theparams);
zresult *znode_switch(zresult *tfirst, zresult *biggest);
zresult *zresult_list_abolish(zresult *first);
void print_zresult(zresult *new_node);
int znode_QTLs(zresult *znode, params *theparams);
int znode_QTLs_estimate(zresult *znode, params *theparams, aqtl *theqtls, markermap *themap, int which);
void calc_recomb_dist(aqtl *theqtls, markermap *themap, int numqtls);

int process_bootfile(params *theparams, char *infile, char *outfile, char *progname, char *chptr);
void write_zheader2(char *outfile, params *theparams,   int eorc, int boots);
int process_permEfile(params *theparams, char *efile, char *cfile  );
int process_jzfiles(params *theparams,char *mainfile, genome *first);
void assign_zfilecols(params *theparams, int *lrcol,int *acol,int *dcol, int *r2col, int *tr2col, int *scol);


void    print_jzestimates(params *theparams,int traits,int *whtraits,jzresult *thejzresults,FILE *fptr);
jzresult *jzresult_node(int flag, jzresult *thisnode,int traits);
jzresult *jz_qtlnum(params *theparams,FILE *fptr,int *nqtls,int traits);
void  get_jzestimates(params *theparams,char *mainfile,jzresult *thejzresults,int *whtraits,int wt );
void write_jzmapqtl_ranks(jzresult *thejzresults, params *theparams,genome *first,char *srfile);
genome *NearMark(genome *first, int chrom, FPN cMposition);

void write_zmapqtl_ranks(params *theparams, aqtl *qtlptr,char *srfile);

aresults *AllocAnalysisResults(aresults *first, int rows, FPN window, int nbp, int trait,int filetype, char *filename, int model,int ihypo, long offset);
void UnAllocAnalysisResults(aresults *first);
aresults *AnalysisFirstPass(char *filename);
void AnalysisSecondPass(aresults *first) ;
void ShowAnalysisResults(FILE *out, aresults *cnode,params *theparams);
void SetAnalysisHeaders(aresults *cnode,params *theparams);
void AppendResults(aresults *first, aresults *second);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file EQfunc.h
------------------------------------------------------------------ */


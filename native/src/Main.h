/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Main.h) is part of QTL Cartographer
         
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


#include "LocalD.h"

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifndef UNIX
#include <extras.h>
#endif

/*
#if defined (__WIN32__)
   #include <conio.h>
#endif
*/

#if defined(MACWARRIOR)   
extern void *malloc(size_t);
#include <sioux.h>
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#endif

#if defined(MACWARRIOR)
#include <console.h>
/*char **environ;*/
#endif

/*
#ifdef WIN32
#include <windows.h>
#endif
*/

    /* Define official version */
#ifndef VERSION
#define VERSION "QTL Cartographer v. 1.17j, 28 January 2005"
#endif

/*  These are some default parameters for Simulations...*/
#define MAX_REPS 10001      /*By default, only up to 10,000 repetitions for permutations or bootstraps. Probably shouldn't go higher than 31,000. */
#define REPS 0              /*No reps as starting point Various*/
#define NN 200              /*Sample size   Rcross */
#define NINBRED  1          /*Vestige of something ? */
#define DEFsigl  0.0        /*StDev of markers/chromosome Rmap*/
#define DEFsigs  0.0        /*StDev of intermarker distances Rmap*/
#define DEFm     4          /*Number of chromosomes Rmap*/
#define DEFqtl   9          /*Ave # of QTL per trait Rqtl*/
#define DEFl     16         /*Ave # of markers/chromosome Rmap*/
#define DEFs     10.0       /*Ave intermarker distance   Rmap*/
#define DEFbrdrs 0.0        /*Ave amnt of flanking DNA on each chrom. Rmap*/
#define DEFbeta  2.0        /*Beta parameter for effects Rqtl*/
#define DEFheritability 0.5 /*Heritability Rcross*/
#define SIG_LEVEL 3.84      /*Experimentwise sig. level Eqtl, Preplot*/
#define NUM_SIG 5           /*Number of background markers in CIM. Zmapqtl, JZmapqtl  (nbp)*/
#define WIN_SIZE 10.0       /*Window size in CIM. Zmapqtl, JZmapqtl */
#define DELTAX 2.0          /*Walking speed in cM.  Zmapqtl, JZmapqtl*/
#define MIN_DATA 10         /*Minimum sample size for typing markers. Qstats*/
#define HLINE     70        /*length of horizontal lines in output*/
/*  These are for the ECM algorithm*/
#define STOP_EM  0.00000001
#define STOP_MIM 0.0000000001
#define M_TIME 10000
#define MIN_DIST 0.0001  /*minimum distance from a marker for CIM analysis. */
#define MININTERVAL 0.0001
#define SIGNOFDEVIL 666.0  /* This just indicates that the ECM algorithm failed at a given site. */

#if defined(SIXTYFOURBIT) 
#define MAXQTL 39          /*  Maximum number of QTL to be fitted in MImapqtl */
#define MAXEPISTATICS  190  /* Maximumn number of interaction terms fitted in MImapqtl*/
#else
#define MAXQTL 19          /*  Maximum number of QTL to be fitted in MImapqtl */
#define MAXEPISTATICS  190  /* Maximumn number of interaction terms fitted in MImapqtl*/
#endif
#define MINJOINTFREQ     0.001  /*  Minimum joint frequency for QTL array to be used in MImapqtl*/
#define MAXGT            1000    /*  Maximum number of joint genotypes to consider*/
/*These are for mapfunctions and what missing values are coded as */
#define MAPFUNCTION 1
#define MISS_VAL -1000000.0

/*Some general globals. */
#define MAXNAME 63
#define MAXLINE 511
#define SEED 11221960
#define BIG 2147483647.0   /*  = 2e31 - 1  */
#ifndef PI
#define PI 3.141592654
#endif
#define MAPDELTA 0.00001   /*Used in mapping functions for iterations  Utilities.c */
#define MAXCHROMLEN 10000   /* Maximum chromosome length, in cM.   10000*/
extern char time_buffer[MAXNAME+1];  /*String to keep the time */
extern char gbuffer[MAXLINE+1];      /*Reusable global buffer space  */
extern char gname[MAXNAME+1];        /*Reusable global name space    */
extern char gwarn[MAXNAME+1];        /*Reusable global warning space */

/* Some structures */

/*             QTL                  m3
c1 .......|.....K.....|.............|..........|....K....
   brdrs?              <-----s----->       
c2        |.....K.....|.............|..........|

brdrs are flanking DNA, that is DNA outside the first and last
markers on a chromosome.
*/

typedef struct MapofMarkers {
  int m;         /* number of chromosomes */
  int l;         /* average number of markers per chromosome */
  int maxl;      /* maximum number of markers on any chromosome */
  FPN sigl;   /* variance in number of markers per chromosome.  <= 0 => fixed */
  FPN s;      /* average distance between consecutive markers in cM */
  FPN sigs;   /* variance of s */
  int *mpc;      /* = ivector(1,m); number of markers for each chromosome.  = l forall iff sigl <= 0 */
  FPN **mrf;  /* = dvector(1,m,0,max(mpc)+1); If sigl <= 0, then l is max(mpc), and mpc isn't used.
                      pointer to a matrix of recombination frequencies between markers i and i+1 */ 
  int ml;        /* total number of markers, = m*l iff sigl <= 0 */
  FPN brdrs;  /* Will there be chromosomal material outsite the flanking chromosomal markers? 
                    If this is negative, then there will be no borders.  */
  int *knum;     /* = ivector(1,traits) Number of QTLs on this map for each of the traits */
  int traits;    /* Number of traits to simulate */
  int otraits;   /* Number of other traits */
  int *otypes;   /* Number of classes for each of the other traits */
  char **tnames;  /* Names of the traits */
  FPN *ParentalDiff;   /* = dvector(1,traits).   Difference of means of parental lines */
  char **onames;  /* Names of other traits */
  char **cnames; /* Names of the chromosomes */
  char **names;  /* Names of the markers */
  int  **ttable; /* Table to indicate where in names the marker name is */
  int  **types;  /* Indicate the type of marker:  
                       -1  =>  a-           a dom
                        0  =>  codominant
                        1  =>  A-           A dom           */
  int maxqtl;            /*   max over traits of the knum vector */
} markermap;

typedef struct QTLdata {
  int chrm;        /* this qtl is on chromosome chrm */
  int mrk;         /* it resides after this marker */
  FPN c1;       /* the recombination frequency between marker mrk and the qtl */
  FPN c2;       /* the recombination frequency between marker mrk+1 and the qtl */
  FPN a;        /* the additive effect of the qtl */
  FPN d;        /* the dominance effect of the qtl, 0.0 => no dominace */
  markermap *map;  /* pointer to the map of markers */
  int trait;       /* the trait that this is a qtl for. */
  FPN r2;       /*residuals from the estimation stage */
  FPN tr2;      /*total residuals from the estimation stage*/
  FPN s;        /*test statistic for normality of residuals*/
  int nptrs;       /* How many things are pointing to this node? */
  FPN **episAADD;  /* Pointer to the matrix of epistatic interactions.   = dmatrix(1,themap->knum[trait],1,themap->knum[trait]) */
  FPN **episADDA;  /* Pointer to the matrix of epistatic interactions.   = dmatrix(1,themap->knum[trait],1,themap->knum[trait]) */
}  aqtl;

typedef struct IndividualData {
  int **markers;  /* = gives the values of the markers */
  int **vqtls;    /* = imatrix(1,t,1,k) gives the alleles of the qtls */
  int t;          /* number of quantitative traits measured */
  FPN *y;      /* = dvector(1,t) gives the phenotypes for the t traits, 1 by default */
  char **oy;      /* other traits = cmatrix(1,ot,0,MAXNAME) */
  int *oyt;       /*  oyt = ivector(1,otraits)  */
  FPN *g;      /* gives the genotypic value over all QTLs */
  markermap *map; /* pointer to the map of markers */
  aqtl *qtls;     /* pointer to the array of qtl information */
  char *name;     /* name of this individual */
  char print_flag;/* should we print this individual? Or use it in analysis? */
  int bc;         /*which backcross in Design III.   not used yet.*/
}  individual;

typedef struct TotalGenome { /*Structure to hold a genome defined by the genetic linkage map*/
  int chrom;                 /* Chromosome of the marker interval */
  int markr;                 /* The marker that precedes the interval, can be 0 if borders are allowed */
  FPN dist;               /* This distance, from marker markr to markr+1 on chromosome chrm, is in M */
  FPN pos;                /* Position of marker from left telomere in Morgans*/
  int whichqtl;              /* Indicator of whether there is a qtl on the interval 1 if yes, 0 if no */
  FPN mxo;                /* xo position of maternal    xo = 0.0 if no crossover,   */
  FPN pxo;                /*             or paternal    xo > 0.0 the distance (M) from the marker to the xo */
  char *markername;          /* Name of the left flanking marker */
  FPN *values;            /* Vector of expected values at position pos  or recombination fractions (Emap) */
  int n;                     /*  Sample size or number of markers (Emap)*/
  struct TotalGenome *prev;  /* Pointer to previous marker interval */
  struct TotalGenome *next;  /* Pointer to next marker interval */
}  genome;


typedef struct MIMparam { /* MIM parameters will be held in a linked list. */
  int abd;                 /* abd    use 1 for additive, 2 for dominant and 3 for epistatic*/
  FPN value;                /* current parameter estimate */
  FPN ovalue;          /*previous parameter estimate */
  FPN lnlike;          /*  Log-likelihood if this node is not in the model */
  int qtl1;              /*  which QTL */
  int qtl2;              /* which other QTL if epistatic interaction */
  int gt1;               /*  Genotype index of QTL 1 for use with jpvt */
  int gt2;                /*  Genotype index of QTL 2 for use with jpvt */
  FPN design;        /* design value  for the genotype... updated as needed */
  FPN dbar;          /* mean of designs...used in calculating var-covar matrix*/
  struct QTLdata *qptr1;          /*  Pair of pointers to the QTL  */
  struct QTLdata *qptr2;
  int q1array;                   /* 0 if allocated in MIfunc, and 1 if part of theqtls array */
  int q2array;
  struct MIMparam *prev;  /* Pointer to previous marker interval */
  struct MIMparam *next;  /* Pointer to next marker interval */
}  mimparam;


typedef struct otraittype {   /*Structure to hold information about categorical traits and the categories*/
  char *name;                 /* Categorical (Otrait) trait name */
  int   which;                /* Which Otrait it is for */
  struct otraittype *prev;
  struct otraittype *next;
}  otnode;
  
typedef struct aline {   /*Structure to hold inbred lines that are simulated*/
  char *name;            /*Name of the inbred line*/
  individual *iptr;      /*pointer to the data for that line*/
  int nn;                /*sample size of this data set*/
  struct aline *prev;
  struct aline *next;
  int which;             /* */
}  thelines;

typedef struct Params {  /*Structure to hold parameters*/
  char *resource;           /*file that holds current parameter values */
  char *error;              /* log file*/
  char *map;                /* genetic linkage map*/
  char *mapin;              /* input file for Rmap*/
  char *qtl;                /* genetic model file*/
  char *eqtl;               /* file for estimates of genetic model*/
  char *qtlin;              /* Rqtl input file*/
  char *ifile;              /* data file */
  char *iinfile;            /* input for Rcross*/
  char *qstat;              /* output of Qstats*/
  char *lrfile;             /* output of LRmapqtl*/
  char *srfile;             /* output of SRmapqtl*/
  char *zfile;              /* output of Zmapqtl*/
  char *mimfile;              /* output of MImapqtl*/
  char *stem;               /* stem for filenames*/
  char *mqtfile ;           /* MImapqtl model file (that already exists) */
  char *bayesfile;          /* Bmapqtl.out file */
  char *workdir;            /* working directory*/
  char *term;               /* terminal type (Preplot)*/
  char *tfile;              /* temporary file*/
  char *mimwork;            /* string to encode what MImapqtl should do. */  
  char *mrinput;             /* Output of JZmapqtl and input of MultiRegress*/
  char *mroutput;             /* Output of MultiRegress*/
  
  long seed;       /*random number seed*/
  char *thecross;  /*type of cross*/
  int chrom;       /*chromosomes in map*/
  int wchrom;      /*which chromosome to analyze*/
  int mark;        /*markers per chromosome*/
  FPN vmark;    /*st.dev. of markers/chrom*/
  FPN dist;     /*intermarker distance*/
  FPN vdist;    /*st.dev. of intermarker distances*/
  FPN tail;     /*flanking DNA*/
  int qnum;        /*number of QTL*/
  int dom;         /*dominance flag*/
  FPN beta;     /*parameter for additive effects in simulation*/
  FPN beta1;    /*parameter for dominance effect beta1 in simulation*/
  FPN beta2;    /*parameter for dominance effect beta2 in simulation*/
  int traits;      /*number of traits*/
  int whichtrait;  /*which trait to analyze*/
  int reps;        /*number of repetitions in some simulations*/
  FPN Herit;    /*heritability*/
  FPN Environ;  /*envrinmental variance...overrides heritability*/
  int cross;       /*primary cross type */
  int crosst;      /*primary cross generations   */
  int tcross;      /*type of test cross*/
  int tcrosst;     /*generations of test cross*/
  int nn;           /*sample size*/
  int Model;        /*CIM analysis model*/
  FPN walk;      /*walking speed (cM) for CIM*/
  int Inter;        /*Interactive level flag*/
  int mapfunc;      /*mapping function*/
  int lodflag;      /*whether to convert to LOD scores*/
  FPN maxlr;     /*maximum LR.  Calculated by Eqtl and it ratchets.*/
  FPN maxeffect; /*maximmum additive effect estimate.  Calc. by Eqtl and ratchets.*/

  int gout;         /*output mode flag for Rmap*/
  int Rmode;        /*flag for map simulation type*/
  int verbosity;    /*flag to turn on/off verbose output*/
  int ihypo ;       /*Which hypothesis for Eqtl */
  int boot;         /*Prune flag for what to do*/
  FPN siglevel;  /*Experimentwise sig. level*/
  FPN size;      /*type I error probability*/
  FPN srf1;      /*alpha for adding markers in forward stepwise regression*/
  FPN srb1;      /*alpha for deleting markers in backward elimination regression*/
  int srm;          /*flag for which stepwise regression method to use*/
  int srupper;      /*maximal number of markers that can be added in forward regression*/
  int emethod;      /* Emap method for linkage map creation*/
  int emapobj;      /* Emap objective function*/
  FPN segsize;   /* Emap size of test for segregation*/
  FPN linksize;  /* Emap size of test for linkage */
  
  int nbp;          /* # background markers in CIM*/
  FPN window;    /* window size in CIM*/
  int perms;        /*number of permutations in permutation test*/
  int boots;        /*number of bootstraps to do (actually, just a flag to do one)*/
  int  whoseic;    /*  Which IC for use in MImapqtl and stopping criterion*/
  FPN mimlod;    /*  LOD for adding or dropping parameters*/
  int mimphase;     /* which phase of MImapqtl.  */
  int maxqtl;              /* Maximum number of QTL to be fitted in MImapqtl */
  int maxepistatics;       /* Maximum number of epistatic effects to be fitted in MImapqtl */
  int crosstype;           /* Ultimate cross type, 1, 2, 3, 4, 5 */
  int ngt;                 /* Number of genotypic classes        */
  FPN null_sse;  /*Null hypothesis SSe*/
  FPN total_var; /*phenotypic sample variance of a trait*/
  FPN mapparam;  /*extra parameter needed by some mapping functions*/
  int rwd;          /*rwd secret flag*/

 /* following variables are used in BTL_LRfunc.c and BTL_LRmain.c*/
  FPN p1;        /* distribution of the trait for Parent 1 */
  FPN p2;        /* distribution of the trait for Parent 2 */
  FPN p1_start;  /* start point of distribution range of the trait for Parent1 */ 
  FPN p1_end;    /* end point of distribution range of the trait for Parent1 */ 
  FPN p1_step;   /* step size of distribution range of the trait for Parent1 */ 
  FPN p2_start;  /* start point of distribution range of the trait for Parent2 */ 
  FPN p2_end;    /* end point of distribution range of the trait for Parent2 */ 
  FPN p2_step;   /* step size point of distribution range of the trait for Parent2 */ 
  int    Btl_mode;  /* mode 0 is input specific p1 or p2 value, mode 1 is testing range data*/
  markermap *themap; /* current map */
  genome *thegenome; /* current genome */
  aqtl *theqtls;     /* current set of QTL */
  individual *thedata; /* current data set */
} params;


#include "Utilities.h"
#include "NumRec.h"
#include "params.h"

#include "Blas.h"
#include "Linpak.h"

typedef struct LinpakWorkspace {
  FPN **xx;       /* design matrix...dmatrix(1,p,1,n) */
  FPN **xsave;    /* save design matrix...dmatrix(1,p,1,n) */
  int ldx;           /* leading dimension of xx */
  int n;             /* columns of xx (sample size) */
  int p;             /* rows of xx  (background parameters) */
  int k;             /* should we make this the number of static rows?  mean plus otraits? */
  int t;                    /* number of traits in multitrait analysis */
  FPN **y;               /* trait matrix...dmatrix(0,t,1,n) */
  FPN *wy;              /* working trait vector...dvector(1,n) */
  FPN *qraux;     /* auxiliary vector...dvector(1,p) */
  FPN **rsd;             /* residual vector...dmatrix(0,t,1,n) */
  FPN **wrsd;            /* working residual vector...dmatrix(0,t,1,n) */
  FPN **bb;               /* regression coefficients vector...dmatrix(0,t,1,p) */
  int  **jpvt;             /* pivot matrix...imatrix(0,t,1,p) */
  long  **gtindex;             /* pivot vector...imatrix(0,t,1,p) */
  int **bp;          /* vector to indicate the chromosome and marker of background markers*/
  FPN *pp1;       /* a priori probability of QQ genotypes */
  FPN *pv;        /* a posteriori probability of QQ genotypes */
  FPN *pp2;       /* a priori probability of Qq genotypes */
  FPN *qv;        /* a posteriori probability of Qq genotypes */
  FPN **estimates;       /* parameter estimates for various hypotheses. 
                                                    = dmatrix(0,t,1,9)    */
  FPN **s2;       /* dmatrix(0,t,0,t) the variance-covariance matrix   */  
  FPN **s2i;      /* dmatrix(0,t,0,t) the inverse of the variance-covariance matrix   */  
  int *samplesize;   /* ivector(0,t) sample sizes for each trait.         */ 
  int *kpvt;         /*ivector(0,t) work space for invdet */
  FPN *work;      /*dvector(0,t) work space for invdet */
  int ipos;          /* the number of test positions */
  FPN *lratio;    /* dvector(1, ipos); Keep track of the LR of the sample */
  FPN *ts0;       /* dvector(0, t); null hypothesis Likelihood */
  int *pcnts;        /* ivector(1, ipos); Count number of shuffled sets with LR > SampleLR */
  FPN *maxlr;     /*dvector(0,themap->m) to hold max LR's for perm. test */ 
  FPN *thestats;
} linpakws;

linpakws *cd_jzmapws(linpakws *ilnpk, int n, int p, int t, int cd);
linpakws *cd_mimws(linpakws *ilnpk, int n, int p, int t, int cd);
linpakws *cd_linpakws(linpakws *ilnpk, int n, int p, int cd);
linpakws *cd_multiregressws(linpakws *ilnpk, int n, int p, int t, int cd);
void copy_phenotypes(params *theparams, markermap *themap, individual *individs, linpakws *lnpk);
int how_many_traits(params *theparams, markermap *themap);
int actual_trait(params *theparams,markermap *themap,int trait);
void zplace_qtls(params *theparams, aqtl *qtlptr, genome *first);

#include "MRfunc.h"
#include "MIfunc.h"
#include "MZfunc.h"
#include "Zfunc.h"


#include "RMfunc.h"
#include "RQfunc.h"
#include "RCfunc.h"
#include "EQfunc.h"
#include "QSfunc.h"
#include "LRfunc.h"
#include "SRfunc.h"
#include "BTL_LRfunc.h"

#include "Mdatain.h"
#include "Qdatain.h"
#include "Idatain.h"

#include "Efunc.h"
#include "Genome.h"
#include "MissMark.h"
#include "Otraits.h"

/*  BIM modifications.  */
#include "f2c.h"
#include "revjump.h"
#include "ranlib.h"


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Main.h
------------------------------------------------------------------ */


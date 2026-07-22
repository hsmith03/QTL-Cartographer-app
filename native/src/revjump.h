/* ------------------------------------------------------ XCutXCodeXSkip
     This file (revjump.h) is part of QTL Cartographer
         
    		Copyright (C) 2000-2005
	Patrick Gaffney and Brian Yandell

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

#ifndef REVJUMP_GUARD
#define REVJUMP_GUARD
/**************************************************************
  File:         revjump.h
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:
  -------

    Header file for MCMC program.  See also chrom.h and qtl.h

**************************************************************/


#ifndef VERSION
#include "Main.h"
#endif

#include "chrom.h"
#include "qtl.h"
#include "mygenome.h"
#include "ranlib.h"

#include <assert.h>
int qtl_compare(const void *a, const void *b);

#define MAX_INTERACTION 100
#ifndef MCMCMAXQTL
#define MCMCMAXQTL 30
#define MAX_CHOL (MCMCMAXQTL*2+1)
#endif





/*used to index prior mean,var and log_var (mean, additive and dominance respectively */
#define MU  QTL_NONE
#define ADD  QTL_ADD
#define DOM  QTL_DOM 


#define MAX_PARAM 2
#define MAX_ADD_PARAM 2


/*
#define INFINITY 1E40
*/


/*codes for birth/death moves */
#define C_SWAP 5
#define C_EXCHANGE 4
#define C_UPDATE 3
#define C_BIRTH 1
#define C_DEATH 2

#define NOT_IMPLEMENTED 4
#define BAD_CODE 5




/*--------------------------------------------------------------------------
 *One field in the MCMC_PARAM field contains a code which controls the 
 *analysis prefereence.  Generally it is set to 1, which activates RJ-MCMC,
 *a Poisson prior on the number of QTL, JAYA sampling for new effects (in 
 *RJ-MCMC).  Also additive/dominance effects are sampled as dictated by the
 *design.
 *--------------------------------------------------------------------------
 */

#define MCMC_SELECT_PRIOR_FIELD   0xF      
/*distrib for prior number of QTL  */
#define MCMC_SELECT_SAMPLE_FIELD  0xF0     
/*method in RJ-MCMC of proposing effect */
#define MCMC_SELECT_DOM_FIELD     0xF00    
/*obtain estimates for dominance        */
#define MCMC_SELECT_EFFECT_FIELD  0xF000   
/*type of variance on prior for effects */
#define MCMC_SELECT_PROPOSE_FIELD 0xF0000  
/*proposal from birth death is drawn 
                                    from this distribution*/
#define MCMC_FLAGS                0xF00000 
/*bit flags, described below, to turn
                                     off default features */
#define MCMC_SELECT_MIX           0x3000000
/*change mix from 90-10 birth/death to
								     update (2 bits used), 3rd bit also used as
                                     a flag */
#define MCMC_SELECT_OPTIMIZE      0x30000000  
/*some optimization code bit flags (2 bits used) */
#define MCMC_SELECT_CHROM         0xC0000000                        
                                           


/***************************/
/* MCMC_SELECT_PRIOR_FIELD */
/***************************/
#define SELECT_SIMPLE_MCMC      0         
#define SELECT_NQTL_POISSON     1         
#define SELECT_NQTL_GEOMETRIC   2
#define SELECT_NQTL_UNIFORM     3    
#define SELECT_NQTL_FISCH       4

/****************************/
/* MCMC_SELECT_SAMPLE_FIELD */
/****************************/
#define SELECT_SAMPLE_DECOMP      0x00    /* default, use decomposition */
#define SELECT_SAMPLE_FISCH       0x10 

/**************************/
/* MCMC_SELECT_DOM_FIELD  */
/**************************/
#define SELECT_GET_APPROP         0x000    /* default,get appropriate model */
#define SELECT_GET_ADD_ONLY       0x100  
#define SELECT_GET_DOM_ONLY       0x200    
#define SELECT_RJMCMC             0x300   /* select components to add 
                                             using RJ-MCMC */

/****************************/
/* MCMC_SELECT_EFFECT_FIELD */
/****************************/
#define SELECT_VAR_GAFFNEY        0x0000
#define SELECT_VAR_FISCH          0x1000  

/****************************/
/* MCMC_SELECT_PROPOSE      */
/****************************/
#define SELECT_PROPOSE_PRIOR      0x00      /* default, propose birth-death 
                                               probs from prior */
#define SELECT_PROPOSE_POISSON    0x10000   /* always propose poisson */
#define SELECT_PROPOSE_GEOMETRIC  0x20000   /* always propose geometric */
#define SELECT_PROPOSE_UNIFORM    0x30000   /* always propose uniform */
#define SELECT_PROPOSE_FISCH      0x40000   /* always propose geometric */


/***************************/
/* MCMC_FLAGS              */
/***************************/
/* these flags allow additional options to be set or not (one bit).
   The previous ones allow a range of options (up to 16, x0 to xF,
   4 bits worth)  */

#define MCMC_SELECT_100PCNT_UPDATE 0x100000  /* default update parameters
                                                only if no birth or 
                                                death update */
                                                
#define MCMC_SELECT_BD_NOTEXACT    0x200000  /* here bp+bp <= cval.
                                                Whereas, by default 
                                                we have bp+cp=cval.
                                                cval is set by setting one of
                                                the MCMC_OPTIMIZE flag in 
                                                'revjump' */   
											
#define MCMC_SELECT_SIMPLE_GIBBS   0x400000  /* default is to use a multiple 
                                                parameter update for effects 
                                                (rather than simple Gibbs) */
											          
#define MCMC_SELECT_MOVE_NONE      0x800000    /* default is to propose 
                                                  update on any chrom */

/***************************/
/* MCMC_MIX                */
/***************************/
/* 2 bit field controlling birth-death-update mix */
#define MCMC_OPTIMIZE_DEFAULT     0x0000000    /* standard 90-10 */
#define MCMC_OPTIMIZE_80_20       0x1000000    /* standard 80-20 */
#define MCMC_OPTIMIZE_70_30       0x2000000    /* standard 70-30 */
#define MCMC_OPTIMIZE_50_50       0x3000000    /* standard 50-50 */


/* set this flag to always turn on long-range update (regardless of
  whether a birth/death occured.  Default is to have a long-range
  jump only when updating QTL genotypes */
#define MCMC_OPTIMIZE_BOOST       0x4000000    


/***************************/
/* MCMC_OPTIMIZE           */
/***************************/
/* 2 bit field controlling optimization -- 3 bits*/
#define MCMC_OPTIMIZE_NOPROB     0x10000000  /* turn off prob */
#define MCMC_OPTIMIZE_NOTABLE    0x20000000  /* turn off prob andtable  */
#define MCMC_OPTIMIZE_NONE       0x30000000  /* generate each individual's
                                                prob of QTL geno given
                                                marker */
#define MCMC_FLAG_RANDOM_CHROM   0x40000000  /* default is to select the 
						chromosome with prob proportional
						to its length.  This flag makes
						each chromosome equally likely
						to be chosen (affects acceptance
						ratio) ... lowers it for smaller
						chromosomes.
												*/
#define MCMC_FLAG_NEW_PRIOR   0x50000000  /* selects each chromosome with equal
					     probability, and then position is 
					     uniform over length.  However, this
					     is considered the proior (rather than
					     uniform over the genome).
					  */








typedef struct workData
{
  /* effect and proposed (new) effect */

  /* used in QR decomposition of matrix of genotypes X
	 * with Y = phenotype values, Z = proposed (new) geno,
	 *
	 *  K = Z'RY,   M= 1/(Z'RZ), QZ = Q'Z, QY = Q'Y
	 *  R = I - X{(X'X)^-1}X' = I - QQ', 
	 *
	 *  Note: K = (QZ)'(QY) and M = 1/(QZ'QZ)
	 *        See get_effect(...) in birth.c for full discussion
     */	
  FPN *u;        /*proposed uniform number */
  FPN *resid;    /*residual (y-mu-effect) for all individuals */
  FPN *newResid; /*?????*/
  int*perm_num;   /*use to permute QTL in random_perm_QTL */
  
  /*used in update_params.c */
  FPN *XtY;
  FPN *new_XtY;
  FPN  **XtX;
  FPN  **new_XtX;
  FPN  **chol;
  FPN *p;
  FPN *pvar;
  FPN *pmean;
  FPN *mod_effect;
  FPN *weight;
  FPN *work;
  FPN *oldWeight;   /*needed in proposeFischDeath */
  FPN *oldVar;      /*needed in proposeFischDeath */

  /*all we need to retain from effect update (for birth/death acceptance ratio */
  FPN  bCb;
  FPN  d_invD_d;
  FPN  log_det_Chol;
  FPN  log_det_invD;


  FPN *new_r;
  FPN *log_norm_const;

  QTL_INFO **output_qtls;

  FPN *minProbs;         /*for diagnostics */
  FPN *avgMinProbs; 

} WORK;





typedef struct priorstruct
{
  FPN  mean[3];
  FPN  var[3];
  FPN  alpha[3];          /*prior for variance for add/dom var, 
			       which is Beta(1,beta) */
  FPN  beta[3];          /*prior for variance for add/dom var, 
			      which is Beta(1,beta) */
  FPN  log_var[3];
  int sampleVar[3];


  FPN  log_mu_var;
  FPN  sig_a1;
  FPN  sig_a2;
  FPN  qtl_mean;
  FPN  p;

  int priorDistribution;          /*can be one of  POISSON   Here qtl_mean is the mean of the distrib.UNIFORM   Here qtl_mean holds the maximum number of QTL GEOMETRIC  "   qtl_mean is the mean number of QTL (1/q) Prob of x QTLs = q^x *p (x=0,1,2,...) where q=1-p MIXED     Here qtl_mean is the cutoff before which the prior is uniform.  After there is an geometic decay at rate p. */
} PRIORS;


typedef struct MCMCParam
{
  /*1 => reversible jump, 0 otherwise. */
  int revjump;
  int niter;
  int nby;
  FPN  cval;

  int offset;
  /* int numParam;*/    /*number of parameters to add in birth/death step */
  int addParam;

  FPN  burnIn;
  FPN  preBurnIn;

  FPN  HM;
  FPN  SHM;
  int idx;

  FILE *debug;     /*points to FILE for debug messages */
  FILE *diag;
  FPN  *diag_fY;
  FPN  *maxResid;
  FPN  *pcnt_gt_99;
  FPN  *avg_prob_Q_given_M;
  FPN  *min_prob_Q_given_M;

} MCMC_PARAM;


typedef struct Data {
  /*data for analysis */
  int nn;
  int nChrom;
  int totalMark;
  FPN  totalChromLen;
  FPN *y;
  int *offset;
  FPN  ybar;
  FPN  y_var;
  FPN  y_stdev;
  individual *individs;
  params *theparams;
  markermap *map;

  CHROMOSOME *myChroms;
  FPN *chrom_pos;
   

  /*death-birth stuff */
  FPN *bp;
  FPN *dp;
  FPN *prior_ratio;

  /*some parameters */
  int nQtl;
  FPN  mu;
  FPN  sigmasq;

  /*cross parameters */
  int gmiss;
  int na[3];

  QTL_INFO **myQtls;

} DATA;


typedef struct Missing {
  int lMark;     /*0 => undefined, otherwise there is a valid flanking marker */        
  int rMark;     

  FPN  lPos;   /*position of marker */
  FPN  rPos;

  FPN  scale;  /*In general, we can compute probabilities of a qtl [given
		    flanking markers] to within a constant of proportionality.
		    Only in simple cases (no A- or a- genotypes) can we
		    determine this constant uniquely */

  FPN  **lProb;  /*for each missing marker, we have a vector of 3 values */    
  FPN  **rProb;  /*this give prob of each genotype at the missing marker */    
} MISSING;


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>




#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif


#ifndef MIN
#define MIN(A, B) ((A) > (B) ? (B) : (A))
#endif




#ifndef ABS
#define ABS(A) ( ((A)>0)? (A) : -(A) ) 
#endif

#ifndef EPSILON
#define EPSILON 1E-7
#endif


/*Jaya ususally gives A=genunf() and B=prob[0], then
   returns a marker of genotype 1 wp B and -1 wp (1-B) */
#define random_marker(A,B) ((A) < (B) ? (1) : (-1))
#define ismiss(A) (((A) < -2) ? (1) : (0))


/*mcmc.c */
void mcmc(MCMC_PARAM *myMCMC, DATA *myData, PRIORS *priors,
          CHROMOSOME *chrInfo, WORK *myWork);
void initVars(int nn, FPN  **nmove, int *num_qtl, MCMC_PARAM *myMCMC,
	      FPN  qtl_prior_mn, FPN  cval);
void setupBurnIn (MCMC_PARAM *myMCMC, DATA *myData);
void birth_death(int move_type, DATA *myData, MCMC_PARAM *myMCMC, 
		 WORK *myWork, PRIORS *priors, 
		 CHROMOSOME *chrInfo, FPN   **nmove);
void fixed_locus_update(int nn, int nqtl, QTL_INFO **all_qtls,
			MCMC_PARAM *myMCMC, DATA *myData, PRIORS *priors, 
			WORK *myWork,  
			params *theparams, FPN  **nmove);
void update_effect(int nn, QTL_INFO *thisqtl,  DATA *myData, PRIORS *priors, WORK *myWork);
void update_effects(int nn, int nQtl, int revjump, FPN *y, 
		    FPN *mu, FPN  sigmasq, int *na, 
		    QTL_INFO **qtls, PRIORS *priors,WORK *myWork);
					  
					  

void normalProb(int nn, FPN *resid, FPN  sigmaSq, 
		FPN *diag_fY, FPN *maxResid, FPN *pcnt_gt_99);
void outputResults(FILE *write_res, int nn,
		   int iter, DATA *myData, QTL_INFO **all_qtls, 
		   FPN *resid, QTL_INFO **output_qtls, 
		   PRIORS *priors, MCMC_PARAM *myMCMC, WORK *myWork);
void outputStatistics1(int nn, MCMC_PARAM *myMCMC, params *theparams, 
		       PRIORS *priors, DATA *myData);
void outputStatistics2(int nn, MCMC_PARAM *myMCMC, params *theparams, 
		       FPN   **nmove, int *num_qtl,
		       int countedIter, PRIORS *priors, DATA *myData);
FPN  calc_h2(int nn, FPN *y, FPN  mu, FPN  sigmasq, FPN *resid);
void diagnose(int nn, QTL_INFO **all_qtls, WORK *myWork, 
	      DATA *myData, MCMC_PARAM *myMCMC);
void calcMeanVar(int nn, FPN *vals, FPN *mean, FPN *var);
void updateMean(int i, FPN  val, FPN *mean);


/*lod.c */
FPN  lodnull(int nn, FPN  y_var);
FPN  get_lod(int nn, FPN  sigmasq, FPN  ybar, FPN  y_var,
	       FPN *resid);
/**********************************************************************
    **********************************************************************/



/*birth.c */


int birth(DATA *myData, CHROMOSOME *chrom, PRIORS *priors, 
	  params *theparams, individual *individs, WORK *myWork,
	  MCMC_PARAM *myMCMC);
void setWeights(int nQtl, int *na,
		int revjump, QTL_INFO **all_qtls, 
		FPN *w, QTL_INFO *newQtl);

void get_new_locus(int revjump, 
		   int nChrom, FPN *chrom_pos, FPN  totalChromLen, 
		   CHROMOSOME *chromInfo, CHROMOSOME **chrom, 
		   int *lmark, FPN *pos);
void get_local_locus(CHROMOSOME *chrom, int *lmark, FPN *pos);

void get_new_qtl_genotype(int nn, QTL_INFO *qtl);

FPN  get_effect(long int nn, int nQtl, 
                  FPN *y, FPN  mu, FPN  sigmasq,
		  QTL_INFO *oldQtl, QTL_INFO *newQtl, 
		  QTL_INFO **qtls, int revjump,
		  PRIORS *priors, 
		  int *old_na, FPN *mod_effect,
		  FPN  **XtX, FPN *XtY, FPN  **chol, FPN *p, 
		  FPN *pmean, FPN *pvar, FPN *u, FPN *weight,
		  FPN *bCb, FPN *d_invD_d, 
                  FPN *log_det_Chol, FPN *log_det_invD,
		  WORK *myWork);
				           


QTL_INFO *createQtl(int nn, int qtlNum, QTL_INFO **p_newQtl, 
		    CHROMOSOME *chrom, int lmark, FPN  new_position,
		    int addParam, FPN *a, FPN *w);
void setEffect(int nn, int nQtl, FPN *y, FPN *mod_effect, 
	       QTL_INFO  **all_qtls, FPN *mu, 
	       FPN *w, FPN *resid, int *na);
void setCholParams(WORK *myWork, FPN  bCb, FPN  d_invD_d, 
		   FPN  log_det_Chol, FPN  log_det_invD);
void addQtlToChrom(QTL_INFO *thisqtl);
FPN  calcResidSS(int nn, FPN *resid);
FPN  proposeFischBirth (long int nn, int nQtl, 
			  FPN *y, FPN  mu, FPN  sigmasq,
			  QTL_INFO *newQtl, 
			  QTL_INFO **qtls, int revjump, 
			  PRIORS *priors, int *na,
			  FPN *mod_effect,
			  FPN *pmean, FPN *pvar, FPN *weight,
			  FPN *oldVar, FPN *oldWeight,
			  FPN *resid, FPN *newResid);
FPN  proposeFischDeath (long int nn, int nQtl, 
			  FPN *y, FPN  mu, FPN  sigmasq,
			  QTL_INFO **qtls, int revjump, 
			  PRIORS *priors, int *na,
			  FPN *mod_effect,
			  FPN *pmean, FPN *pvar, FPN *weight,
			  FPN *oldVar, FPN *oldWeight,
			  FPN *resid, FPN *newResid);
int getFischEffect(int nQtl, QTL_INFO **qtls, FPN  mu,
		   PRIORS *priors, FPN *mod_effect, FPN *w, FPN *var);
QTL_INFO *swap_add_dom(int nn, int nQtl, QTL_INFO **qtls,
		       PRIORS *priors, 		              
		       DATA *myData, WORK *myWork, MCMC_PARAM *myMCMC);




/*death.c */
int death(DATA *myData, PRIORS *priors, 
	  params *theparams, WORK *myWork, MCMC_PARAM *myMCMC);
int get_death_interval(int nQtl);

void dropQtl(int nn, int *nQtl, FPN *y, QTL_INFO **all_qtls, 
	     FPN  *mod_effect, FPN *w, FPN *mu, 
	     FPN *resid, int *na);
int removeQtlFromList(int nQtl, QTL_INFO *thisqtl, QTL_INFO **qtls);
void moveQtlToEndofList(int nQtl, QTL_INFO *thisqtl, QTL_INFO **qtls);
QTL_INFO *removeIQtlFromList(int nQtl, int i, QTL_INFO **qtls);


/*accept_birth.c */
int select_move(int nQtl, FPN *bp, FPN *dp);

void calcResid2(int nn, int nQtl, FPN *y, FPN *modified_parms,
		QTL_INFO **qtls, FPN *newResid);
FPN  get_log_proposal_ratio(int nQtl, FPN *bp, FPN *dp, 
			      FPN *priorRatio);
FPN  get_log_position_ratio(int revjump, CHROMOSOME *chrom, 
			      DATA *myData);
FPN  get_prior_position_ratio(FPN  totalChromLen, int nChrom, 
				int nMark, FPN  IMlen);



/* update_parms.c */
FPN  update_mu(int nn, FPN  mu, FPN  sigmasq, FPN *resid,
		 FPN  mu_prior_mn, FPN  mu_prior_var);
FPN  update_add_effect(int nn, int nQtl, FPN  sigmasq,
			 QTL_INFO *modQtl, FPN *resid,
			 FPN  effect_prior_mn, FPN  effect_prior_var);
FPN  update_dom_effect(int nn, int nQtl, FPN  sigmasq,
			 QTL_INFO *modQtl, FPN *resid,
			 FPN  effect_prior_mn, FPN  effect_prior_var);
FPN  Gibbs_update_sigmasq(int nn, FPN *resid,
			    FPN  prior_sigsq_a1,
			    FPN  prior_sigsq_a2);
FPN  MH_update_sigmasq(int nn, FPN *resid,
			 FPN  sigmasq,
			 FPN  deltaS2, FPN  maxS2);
FPN  qtl_prior_effect_ratio (FPN  effect_prop, FPN  w_prop,
			       FPN  effect, FPN  w,
			       FPN  effect_mean, FPN  effect_var, 
			       FPN  log_effect_var);
FPN  calc_prior_effect_term(FPN  effect, FPN  w, FPN  effect_mean,
			      FPN  log_effect_var, FPN  effect_var);

FPN  update_effect_prior_var(int type, int nQtl,  QTL_INFO **all_qtls,
			       PRIORS *priors, FPN  y_var);

void setUpGeno(mygenome *qptr, int **prevGeno, int **nextGeno, int **geno);
void swapRowAndCol(int nn, FPN  **a, int row1, int row2);
void moveQtlToEndOfXtX(int nQtl, QTL_INFO **qtls, int changeQtl, 
		       FPN  **XtX, FPN *XtY, int nterm);
void removeRowCol(int n, FPN  **a, int row);
FPN  log_determinant(int nterm, FPN *p, FPN  sigmasq);
QTL_INFO *long_range_update(int nn, int nQtl, QTL_INFO **qtls, 
			    FPN *chrom_pos, FPN  totalChromLen,
			    int nChrom, CHROMOSOME *chromInfo, PRIORS *priors, 		              
			    DATA *myData, WORK *myWork, MCMC_PARAM *myMCMC);



int cholesky(int n, FPN  **a, FPN *p);
int incremental_cholesky(int n, FPN  **chol, FPN *p, FPN  *newcol);
void choleskySolve(int n, FPN  **a, FPN *p, FPN *b, FPN *x);
void forwardSubs(int n, FPN  **a, FPN *p, FPN *b);
void backwardSubs(int n, FPN  **a, FPN *p, FPN *b);
void mycopy(int n, int m, FPN  **to, FPN  **from);
int setAddDomCovMatrix(int nn, int nQtl, QTL_INFO **selectedQtls, 
		       FPN *y, FPN  **XtX, FPN *XtY, int *na);
				       
void addColToAddDom(int nn, int old_nQtl, int *old_na, 
		    QTL_INFO **qtls, QTL_INFO *thisqtl,
		    FPN *y, FPN *XtY, FPN  **XtX);

void setAddDomDiag_Row1(int nn, int idx, QTL_INFO *thisqtl, 
			FPN *y, FPN  **a, FPN *wr);
int generate_effects(int nterm, FPN  **XtX, FPN *XtY, 
		     FPN *pmean, FPN *pvar,FPN  sigma2, 
		     FPN  **chol, FPN *p, 
		     FPN *u, FPN *mod_effect,
		     FPN *bCb, FPN *d_invD_d,
		     FPN *log_det_Chol, FPN *log_det_invD, int genBeta);
void calc_interaction(int idx, int idx2, FPN  **a, 
		      QTL_INFO *thisqtl, QTL_INFO *bqtl, int nn);
void setPriorMeanVar(int nQtl, int revjump, 
		     QTL_INFO **qtls, QTL_INFO *newQtl, PRIORS *priors, 
		     FPN *w, int *na, FPN *k, FPN *v);



/*update_qtl.c */
void setValidFlag(QTL_INFO *thisqtl, int revjump);
FPN  update_lambda_qtl(QTL_INFO *thisqtl, QTL_INFO **all_qtls, 
			 DATA *myData, WORK *myWork, int revjump);
FPN  get_new_lambda(int nMark, FPN  lambda);
int setNewFlank(mygenome *qptr, int *prevRenew, int *nextRenew,
		FPN *prevDist, FPN *nextDist, int revjump);

void SwapGenome(mygenome **g1, mygenome **g2);
void Swap3DTable(FPN  ****g1, FPN  ****g2);
void Swap2DTable(FPN  ***g1, FPN  ***g2);
void SwapIVec(int **g1, int **g2);
void SwapDVec(FPN  **g1, FPN  **g2);
void SwapInt(int *g1, int *g2);
void SwapDble(FPN *g1, FPN *g2);
void swapQtl(QTL_INFO **g1, QTL_INFO **g2);
void swapQtlData(QTL_INFO *Qtl1, QTL_INFO *Qtl2);

FPN  calc_flank_dist(mygenome *qptr);
void get_new_pos(CHROMOSOME *chrom, FPN  oldPos, FPN  delta, int *newMark, FPN *newPos);
int accept_new_lambda(QTL_INFO *thisqtl, QTL_INFO **all_qtls, int nn, int *offset, int newMark, FPN  newPos, int gmiss,params *theparams, individual *individs, int revjump);
void Gibbs_update_geno(QTL_INFO *thisqtl, int nn, int *offset, FPN  **resid, FPN  sigmasq,  params *theparams, individual *individs, int gmiss,
		       FPN *new_r, int **genotype, FPN  **newResid, FPN *log_norm_const, int revjump);
void propose_qtl_geno(QTL_INFO *thisqtl, int nn, int *offset,
		      FPN *resid, FPN  sigmasq,  
		      params *theparams, individual *individs, int gmiss,
		      FPN *new_r, int *genotype, FPN *newResid,
		      FPN *log_norm_const, int revjump);
FPN  gen_qprob (FPN *a, int geno, FPN  sigmasq, FPN  resid,  FPN  *prob, FPN  *qprob, FPN *new_r, int gmiss);
void setupTable(int gmiss, params *theparams, int renewPrev, int renewNext, QTL_INFO *thisqtl);
void calc_trans_prob(params *theparams, FPN  rec,FPN  **mrec);
						
void normalizeProb(int gmiss, FPN  ***condProb);
void setLambda(QTL_INFO *thisqtl, int mark, FPN  pos, int qtlNum);
void copySwapGenoTo(QTL_INFO *newQtl, QTL_INFO *thisqtl);
void genProbs(int nn, int *offset, params *theparams, individual *individs, int gmiss, QTL_INFO *thisqtl, int revjump);




/*initvals.c */
void GenGenotype(int gmiss, FPN *pr, int *genotype);
void initQtl(int nn, int *offset, QTL_INFO *thisqtl, params *theparams, individual *individs, int gmiss, int revjump);
FPN  IMdist(QTL_INFO *thisqtl);
int *genotype1(QTL_INFO *thisqtl);
int *igenotype(QTL_INFO *thisqtl);
int  **p_igenotype(QTL_INFO *thisqtl);
void dgenotype(QTL_INFO *thisqtl, int nn, FPN *dz);
int qtlnum(QTL_INFO *thisqtl);
FPN  haldane( FPN  dist_val );
FPN  gammln2(FPN  xx);

/*read_data.c */
void getdata( MCMC_PARAM *myMCMC, DATA *myData, PRIORS *priors, CHROMOSOME *chromInfo,  FPN  **new_lambda, int **new_chrom);

/*qtl_move.c */
FPN  calc_intermarker_dist(CHROMOSOME *chrom, int m);
mygenome *findInsertPos(FPN  pos, mygenome *lmark, int chrom);
void insertQtl(mygenome *qptr, int lmark, FPN  pos, CHROMOSOME *chrom, int qtlNum);
void removeQtl(mygenome *qptr);
int checkIntegrity(int nQtl, CHROMOSOME *chrom);
void replaceQtl(mygenome *old, mygenome *new);
void restoreQtl(mygenome *qptr);
int EQUALS(FPN  x, FPN  y);
int LE(FPN  x, FPN  y);
int LT(FPN  x, FPN  y);
FPN  getPos(CHROMOSOME *chrom, FPN  lambda, mygenome **lptr);

/*numcmp.c */
int binSearch(int nval, FPN *vals, FPN  searchVal);
int numcmp(FPN  *v1, FPN  *v2);
void checkResid(int nn, int nQtl, FPN  mu, FPN *y, QTL_INFO **allQtls, FPN *resid, int gmiss);
void calc_cond_prob2(params *theparams, int bc, mygenome *gptr,int kk, FPN  *pAA,FPN  *pAa,FPN  *paa);
void selfed_f_tpm2(FPN  **sftpm, FPN  r, int t);
void bselfed_f_tpm2(FPN  **sftpm, FPN  r, int t, int bc);
void AinR2(FPN  (*rr)[3], FPN  **aa);
void AdotB2(FPN  **rr, FPN  (*aa)[3], FPN  (*bb)[3]);

/*experiment.c */
void Cholesky(int n, FPN **z, FPN  **L);
void Inverse(int n, FPN  **z, FPN  **inverse_z);
void Multiply(int n, int m,  FPN  **a, FPN *b, FPN *c);
FPN  Determinant(int n, FPN  **z);

/*main.c */
void setupMCMC(MCMC_PARAM *myMCMC, DATA *myData, WORK *myWork, PRIORS *priors, params *theparams,   CHROMOSOME **pChromInfo);
void copyGenomeTo(mygenome *mystartptr, genome *startptr, mygenome *prev);
void setupWork(int nn, WORK *myWork);
void setupBirthDeathProbs(int revjump, FPN  qtl_param,  FPN  cval, FPN  **bp, FPN  **dp, FPN  **priorRatio);
void noramlizeBirthDeath(FPN  *bp, FPN  *dp, FPN  cval);
void setupDiagnostics(int nn, MCMC_PARAM  *myMCMC, params *theparams);
void setupQtl(DATA *myData, MCMC_PARAM *myMCMC, WORK *myWork, CHROMOSOME *chromInfo, PRIORS *priors, int *new_chrom, FPN  *new_lambda, params *theparams);			  
void setupChromosomes(DATA *myData, MCMC_PARAM *myMCMC, CHROMOSOME**pChromInfo,params *theparams);
void setupTraitData(DATA *myData, params *theparams );
void random_perm_QTL(int nQtl, QTL_INFO **qtls, QTL_INFO **buff, int *perm_num);
void printX(int nn, int nQtl, QTL_INFO **qtls);
void printXtX(int nterm, FPN  **XtX);

/*diagnostics */
void checkCholesky(int nQtl,int revjump, QTL_INFO **qtls, PRIORS *priors, FPN  *weight, int *na, FPN  *pmean,FPN  *pvar,
		   FPN  **XtX, FPN  *XtY, FPN  sigmasq, FPN  **chol, FPN  *p, FPN  *u, FPN  *mod_effect,					    
           FPN  *bCb, FPN  *d_invD_d, FPN  *log_det_Chol, FPN  *log_det_invD, WORK *myWork);


#endif









/* ------------------------------------------------------- XCutXCodeXSkip
             End of file revjump.h
------------------------------------------------------------------ */


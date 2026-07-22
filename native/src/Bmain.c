/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Bmain.c) is part of QTL Cartographer
         
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

#include "Main.h"


/**************************************************************
  File:         main.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:
  -------

    This is the MAIN program for reversible jump mcmc
    computations.



  Modified by Chris Basten so as to integrate with QTL Cartographer
  Distribution.  December 2004 - January 2005.  
  
**************************************************************/


int nniter, nnby;
 

void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams,  int flag)
{
  
  int ii, jj;
  if (flag == 1) {
    strcpy(opt[1], "-i");
    strcpy(opt[2], "-o");
    strcpy(opt[3], "-e");
    strcpy(opt[4], "-m");
    strcpy(opt[5], "-s");
    strcpy(opt[6], "-c");
    strcpy(opt[7], "-t");
    strcpy(opt[8], "-x");
    strcpy(opt[9], "-d");

    strcpy(opt_e[1], "Input File");
    strcpy(opt_e[2], "Output File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Genetic Linkage Map File");
    strcpy(opt_e[5], "Random Number Seed");
    strcpy(opt_e[6], "Trait to analyze");
    strcpy(opt_e[7], "Chromosome to analyze (>0)"); 
    strcpy(opt_e[8], "Number of iterations");
    strcpy(opt_e[9], "Record increment");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      *(*(opt_v + ii) + jj) = '\0';
  strcpy(opt_v[1], theparams->ifile);
  strcpy(opt_v[2], theparams->bayesfile);
  strcpy(opt_v[3], theparams->error);
  strcpy(opt_v[4], theparams->map);
  sprintf(opt_v[5], "%ld", theparams->seed);
  sprintf(opt_v[6], "%d", theparams->whichtrait);
  sprintf(opt_v[7], "%d", theparams->wchrom);
  
  sprintf(opt_v[8], "%d", nniter);  
  sprintf(opt_v[9], "%d", nnby);
  /*These were commented:  why? */
}

void update_params(char **opt_v,  params *theparams)
{
  strcpy(theparams->ifile, opt_v[1]);
  
  strcpy(theparams->bayesfile, opt_v[2]);
  strcpy(theparams->error, opt_v[3]);
  strcpy(theparams->map, opt_v[4]);
  theparams->seed = atol(opt_v[5]);
  theparams->whichtrait = atoi(opt_v[6]);
  theparams->wchrom = atoi(opt_v[7]);

  nniter = atoi(opt_v[8]);
  nnby = atoi(opt_v[9]);
}



int main(int argc, char *argv[])
{
  CHROMOSOME *chromInfo;
  DATA myData;
  WORK myWork;
  MCMC_PARAM myMCMC;
  PRIORS priors;
  

  /* ============== inserted QtlCart Code  ================== */

  char *chptr, *purpose;
  int nopts, automatic;
  params *theparams;
#if defined(MACWARRIOR)  
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  whichprogram = 11;    
  nniter = 1000;/*  What should the defaults be?   None had been set*/
  nnby = 100;   /*  Also, code to set them had been commented out.  Why? */

  purpose = cvector(0,MAXNAME);  
  strcpy(purpose, "Do Bayesian MCMC QTL Analysis");
  theparams = NULL;
  nopts = 9;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2(); 
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);  
/*  theparams->perms = theparams->boots = 0;   Are these vestiges of Zmapqtl? */
  
/* Initailize the data structures...Updated in 2003. */
  theparams->theqtls = NULL;
  GetTheMap(theparams, theparams->map);
  GetTheData(theparams, theparams->ifile, 0 );
  theparams->thegenome = create_genome(theparams->themap);   
  
  /* ============== end QtlCart Code  ================== */
  
  /* setup number generator */
  while (theparams->seed < 0) theparams->seed += 2147483563;
  if (theparams->seed == 0) theparams->seed = time(NULL);
  
  if ((theparams->seed & 16)==0)
    setall(theparams->seed/4+13, theparams->seed/3 + 211);
  else
    setall(theparams->seed/3 + 211, theparams->seed/4+13);
  
  advnst((theparams->seed & 31)+3);    
  /* here  */
  setupMCMC(&myMCMC, &myData, &myWork, &priors,theparams,   &chromInfo);
  mcmc (&myMCMC, &myData, &priors, chromInfo, &myWork);

/*  Clean up qtlcart memory structures. */
  write_trailer(theparams,chptr,1);
  free_cvector( purpose,0,MAXNAME);
  theparams = create_params(theparams, 0, NULL);
  
  return(1);
}


void setupMCMC(MCMC_PARAM* myMCMC, DATA* myData, WORK* myWork, PRIORS* priors, params* theparams,  CHROMOSOME** pChromInfo)
{
  FPN *new_lambda;
  int *new_chrom;
  
  /* a little help optimizing code for 2-genotype crosses, like BC, DH */
    switch (theparams->cross) {
    case 1:  myData->gmiss = -1; break;
    case 2:  myData->gmiss =  1; break;
    case 5:  myData->gmiss =  0;  break;
    default: myData->gmiss = -2;
    };

   setupTraitData(myData, theparams );
   setupChromosomes(myData, myMCMC, pChromInfo, theparams);
   setupWork(myData->nn, myWork);
   
   /* miscellaneous */
   myData->theparams = theparams;
   myData->individs = theparams->thedata;
   myData->bp=myData->dp=myData->prior_ratio=(FPN *)NULL;  /* initialize */
   
  
   getdata(myMCMC, myData, priors, *pChromInfo,   &new_lambda, &new_chrom ); 
   
   setupQtl(myData, myMCMC,  myWork, *pChromInfo, priors, new_chrom, new_lambda,
 		    theparams );
   
   
#ifndef MCMC_OMIT_DEBUG
   checkResid(myData->nn, myData->nQtl, myData->mu, myData->y, myData->myQtls, 
	      myWork->resid,myData->gmiss);
#endif
   setupDiagnostics(myData->nn, myMCMC, theparams);
   priors->priorDistribution= myMCMC->revjump & 0xF;
   
   /* setup ... */
   if ((myMCMC->revjump & MCMC_SELECT_MIX) == MCMC_OPTIMIZE_80_20) 
     myMCMC->cval = -0.8;
   else if ((myMCMC->revjump & MCMC_SELECT_MIX) == MCMC_OPTIMIZE_70_30) 
     myMCMC->cval = -0.7;
   else if ((myMCMC->revjump & MCMC_SELECT_MIX) == MCMC_OPTIMIZE_50_50) 
     myMCMC->cval = -0.5;
   else 
     myMCMC->cval = -0.9;  /* default */

   /* decide whether to make the birth/death (BD) probs, bp+dp, exactly
      equal to myMCMC->cval (default) or to make it an inequality */
   if (myMCMC->revjump & MCMC_SELECT_BD_NOTEXACT) myMCMC->cval = -myMCMC->cval;
   
   setupBirthDeathProbs(myMCMC->revjump, priors->qtl_mean, 
			myMCMC->cval, &myData->bp, &myData->dp, &myData->prior_ratio);/* here*/
   
   myMCMC->HM = myMCMC->SHM = 0;
   myMCMC->idx = 0;
   

   if (myMCMC->revjump & MCMC_SELECT_SIMPLE_GIBBS) 
     {
       myMCMC->revjump |= MCMC_SELECT_MOVE_NONE;
       myMCMC->revjump |= SELECT_SAMPLE_FISCH;
   }
}





void setupWork(int nn, WORK* myWork)
{
  int i;
  
  
   /* residuals */
  myWork->resid = dvector(1,nn);
  myWork->newResid = dvector(1,nn);
  
  
  myWork->u = dvector(1,MAX_CHOL);
  myWork->pmean = dvector(1,MAX_CHOL);
  myWork->pvar = dvector(1,MAX_CHOL);
  myWork->p = dvector(1,MAX_CHOL);
  myWork->XtY = dvector(1,MAX_CHOL);
  myWork->new_XtY = dvector(1,MAX_CHOL);
  myWork->XtX = dmatrix(1,MAX_CHOL, 1, MAX_CHOL);
  myWork->new_XtX = dmatrix(1,MAX_CHOL, 1, MAX_CHOL);
  myWork->chol = dmatrix(1,MAX_CHOL, 1, MAX_CHOL);
  myWork->perm_num = ivector(1,MCMCMAXQTL);
  
  /* diagnostics */
  myWork->minProbs = dvector(1,nn);
  myWork->avgMinProbs = dvector(1,nn);
  
  /* for qtl updating */
  myWork->new_r = dvector(-1,1);
  myWork->log_norm_const = dvector(1,nn);
  myWork->output_qtls = (QTL_INFO**)malloc((MCMCMAXQTL) * sizeof(QTL_INFO*));
  myWork->output_qtls--;
  for (i=1; i<= MCMCMAXQTL; i++) myWork->output_qtls[i]=(QTL_INFO*)NULL;
  
  myWork->mod_effect = dvector(0,MAX_CHOL);
  myWork->weight = dvector(0,MAX_CHOL);
  myWork->work = dvector(0,MAX_CHOL);
  myWork->oldWeight = dvector(0,MAX_CHOL);
  myWork->oldVar = dvector(0,MAX_CHOL);
  myWork->weight[0] = 1.0;          /* the weight for the prior mean variance */
}




void noramlizeBirthDeath(FPN *bp, FPN *dp, FPN cval)
  /* if cval is positive, then we ensure that they are no more
   * than cval.   The latter ensures that the prior ratio cancels
   * out the ratio of death to birth proposal probabilities in the
   * acceptance ratio for RJ-MCMC (and is standard Green 1995).
   *
   * If cval is negative, we ensure that death and birth probs 
   * add up to exactly cval 
   */
{
    int i;
    FPN mval;
    
    if (cval > 0 && cval <= 1)
      {
	mval = 0;
	for (i=0; i<=MCMCMAXQTL; i++)
	  mval = MAX(mval, (bp[i] + dp[i]));
	mval = cval/mval;
	
	for (i=0; i<=MCMCMAXQTL; i++)
	  {
	    bp[i] *= mval; 
	    dp[i] *= mval;
	  }
      }
    else {
      cval = -cval;
      if (cval <= 0) cval = 0.9;
      for (i=0; i<=MCMCMAXQTL; i++)
	{
	  bp[i] = cval * bp[i]/(bp[i] + dp[i]);
	  dp[i] = cval - bp[i];
	}
    }
}




void setupBirthDeathProbs(int revjump, FPN qtl_param, 
			  FPN cval,
			  FPN **bp, FPN **dp, FPN **priorRatio)
{
  int i;
  int priorType = (revjump & MCMC_SELECT_PRIOR_FIELD);
  int proposalType = (revjump & MCMC_SELECT_PROPOSE_FIELD);
  
  /* birth/death proposal probabilities */
   if (!*bp) *bp = dvector(0,MCMCMAXQTL);
   if (!*dp) *dp = dvector(0,MCMCMAXQTL);
   if (!*priorRatio) *priorRatio = dvector(1,MCMCMAXQTL);
   (*dp)[0] = (*bp)[MCMCMAXQTL] = 0.0;
   
   if (priorType==0)
     {
       for (i=0; i<=MCMCMAXQTL; i++) (*bp)[i] = (*dp)[i] = (*priorRatio)[i] = 0.0;
       return;
     }
   
   for (i=0; i<MCMCMAXQTL; i++)
     {
       switch (priorType) 
	 {
	 case SELECT_NQTL_POISSON:
		  (*priorRatio)[i+1] = qtl_param/(i+1);
		  break;
		  
	 case SELECT_NQTL_GEOMETRIC:
	   (*priorRatio)[i+1] = 1.0/qtl_param;
  	      break;
	      
	 case SELECT_NQTL_FISCH:
	   (*priorRatio)[i+1] = qtl_param/(i+1)/(i+1);
	   break;
	   
	 case SELECT_NQTL_UNIFORM:
	   (*priorRatio)[i+1] = 1;
	   break;
	 default: printf("Should have 1,2,3 or 4 (Poisson,Geometric,"
			 "Fisch or Uniform for rightmost (hex) "
			 "digit of revjump.  Choosing poisson as prior.\n");
	   (*priorRatio)[i+1] = qtl_param/(i+1);
	 };
       
       switch (proposalType)		  
	 {
	 case SELECT_PROPOSE_PRIOR:
	   (*bp)[i] = MIN(1.0, (*priorRatio)[i+1]);
	   (*dp)[i+1] = MIN(1.0, 1.0/(*priorRatio)[i+1]);
	   break;
	 case SELECT_PROPOSE_POISSON:
	   (*bp)[i] = MIN(1.0, qtl_param/(i+1));
	   (*dp)[i+1] = MIN(1.0, (i+1)/qtl_param);
	   break;
	 case SELECT_PROPOSE_GEOMETRIC:
	   (*bp)[i] = MIN(1.0, 1.0/qtl_param);
	   (*dp)[i+1] = MIN(1.0, qtl_param);
	   break;
	 case SELECT_PROPOSE_FISCH:
	   (*bp)[i] = MIN(1.0, qtl_param/(i+1)/(i+1));
	   (*dp)[i+1] = MIN(1.0, (i+1)*(i+1)/qtl_param);
	   break;
	 case SELECT_PROPOSE_UNIFORM:
	   (*bp)[i] = 1;
	   (*dp)[i+1] = 1;
	   break;
	 default: printf("Should have 0,1,2,3,4 (Default,Poisson,Geometric,"
			 "Fisch or Uniform for 5th (hex) "
			 "digit of revjump.  Choosing default.\n");
	   (*bp)[i] = MIN(1.0, (*priorRatio)[i+1]);
	   (*dp)[i+1] = MIN(1.0, 1.0/(*priorRatio)[i+1]);
	  };
       
       
       
     }
   noramlizeBirthDeath(*bp, *dp, cval);
}





void  copyGenomeTo(mygenome* mystartptr, genome* startptr, mygenome* prev)
{
  mystartptr->chrom = startptr->chrom;
  mystartptr->dist = startptr->dist;
  mystartptr->markr = startptr->markr;
  mystartptr->pos = startptr->pos; 
  mystartptr->genotype = (int*)NULL;
  mystartptr->next = (mygenome*)NULL;
  mystartptr->prev = prev;
  if (prev) prev->next = mystartptr;
}





void setupDiagnostics(int nn, MCMC_PARAM* myMCMC, params* theparams)
{
  int i,j,k;
  char debugStr[3] = "deb";
  char diagStr[4] = "diag";
  
  k=1;
  if (theparams->error)	
    {
      /* determine if we want full debug info */
      for (i=0; i<(int)strlen(theparams->error) && theparams->error[i] != '.'; i++) k+=1 ;
      for (j=0, i++; i<(int)strlen(theparams->error) && j<3 && 
	     debugStr[j] == theparams->error[i]; i++, j++) k+=1 ;
      if (j==3) myMCMC->debug=fopen(theparams->error,"w");
      else myMCMC->debug = (FILE*)NULL;
      
      /* determine if we want full diagnostic info */
      for (i=0; i<(int)strlen(theparams->error) && theparams->error[i] != '.'; i++) k+=1 ;
      for (j=0, i++; i<(int)strlen(theparams->error) && j<3 && 
	     diagStr[j] == theparams->error[i]; i++, j++) k+=1 ;
      if (j==3) myMCMC->diag=fopen(theparams->error,"w");
      else myMCMC->diag = (FILE*)NULL;
    }
  else 
    {
      myMCMC->debug = (FILE*)NULL;
      myMCMC->diag = (FILE*)NULL;
    }
  myMCMC->diag_fY = dvector(1,nn);
  myMCMC->min_prob_Q_given_M = dvector(1,nn);
  myMCMC->avg_prob_Q_given_M = dvector(1,nn);
  myMCMC->pcnt_gt_99 = dvector(1,nn);
  myMCMC->maxResid = dvector(1,nn);
  
}




void setupQtl(DATA* myData, MCMC_PARAM* myMCMC, WORK* myWork, 
	      CHROMOSOME* chromInfo, PRIORS* priors, 
	      int* new_chrom, FPN *new_lambda, 
	      params* theparams )			  
{
  QTL_INFO *newQtl;
  CHROMOSOME* chrom;
  int nn = myData->nn;
  int i, c, lmark;
  
  /* qtl info .. we allocate two more than normal for work variables */
  myData->myQtls = (QTL_INFO**)malloc((MCMCMAXQTL+2) * sizeof(QTL_INFO*));
  for (i=0; i< MCMCMAXQTL+2; i++) myData->myQtls[i]=(QTL_INFO*)NULL;
  createQtl(theparams->nn, 0, &myData->myQtls[0], NULL, 0, 0.0, QTL_NONE, NULL, NULL);
  
  
  myData->na[ADD]=myData->na[DOM]=0;
  
  for( i=1; i<=myData->nQtl; i++ ) 
    {
      /* allow for errors, and special case 
	 new_chrom[i] == 0 => use qtl index i as chromosome number
	 new_lambda[i] == 0 => use middle of chromosome as starting position
      */
      c = new_chrom[i] - (myMCMC->offset - 1);
      if (c < 1 || c > myData->nChrom)
	c = ignuin(1,myData->nChrom);
      chrom = &chromInfo[c];
      
      if (new_lambda[i] <= 0.0 || new_lambda[i] >= chrom->chromLen)
	get_local_locus(chrom, &lmark, &new_lambda[i]);
      else
	lmark = binSearch(chrom->nMark, chrom->mark_pos, new_lambda[i]);
      
      newQtl = createQtl(nn, i, &myData->myQtls[i], chrom, 
			 lmark, new_lambda[i], myMCMC->addParam, NULL, NULL);
      if (newQtl->flag & QTL_ADD) myData->na[QTL_ADD]++;
      if (newQtl->flag & QTL_DOM) myData->na[QTL_DOM]++;
      
      
      if (myMCMC->addParam == QTL_NONE) {newQtl->flag = QTL_ADD_DOM;}  
      /* force it to fully fit*/
      
      initQtl(myData->nn, myData->offset, newQtl, 
	      theparams, theparams->thedata, myData->gmiss, myMCMC->revjump);
      newQtl->chrom->nQtl++;
    }
  
  if (myData->nQtl > 0) free_dvector(new_lambda,1,myData->nQtl);
  
  update_effects(nn, myData->nQtl, myMCMC->revjump, myData->y, 
		 &myData->mu, myData->sigmasq, myData->na, 
		 myData->myQtls, priors, myWork);
  
}





void setupChromosomes(DATA* myData, MCMC_PARAM* myMCMC, CHROMOSOME** pChromInfo, params* theparams)
{
  int i,c, idx, k;
  int *geno;
  genome *startptr;
  mygenome *mystartptr, *prev;
  CHROMOSOME *chromInfo;
  CHROMOSOME *chrom;
  int nn = myData->nn;
  int *offset = myData->offset;
  
  /* setup the number of chromosomes to analyze */
  if (theparams->wchrom == 0) {
    /* analyze all chromosomes */
    myData->nChrom = theparams->chrom;
    myMCMC->offset = 1;   /* first chromosome in our data structure
					       is number 1 */
  }
  else {
    myData->nChrom = 1;
    myMCMC->offset = theparams->wchrom;
  }
  
  
  /* set the pointer into the genome (marker) structure */
  startptr = theparams->thegenome;
  while (startptr && startptr->chrom != myMCMC->offset)
    startptr = startptr->next;
  
  /* allocate space for chromosome information */
  *pChromInfo = (CHROMOSOME*)malloc(myData->nChrom * sizeof(CHROMOSOME)); 
  chromInfo = --(*pChromInfo);   /* so indexing is 1,2,.... */
  myData->myChroms = chromInfo;
  myData->totalMark = 0;
  myData->totalChromLen = 0;
  myData->chrom_pos = dvector(1,myData->nChrom+1);
  myData->chrom_pos[1] = 0.0;
  
  
  
  for (c=1, idx=myMCMC->offset; c<= myData->nChrom; c++, idx++) 
    {
      chrom = &chromInfo[c];
      chrom->num = idx;
      chrom->nMark = theparams->themap->mpc[idx];
      chrom->nQtl=0;
      
      chrom->mark_genome = (mygenome**)malloc(sizeof(mygenome*) * theparams->themap->mpc[idx]);
      chrom->mark_pos = dvector(1,chrom->nMark);
      chrom->mark_genome--;
      prev = (mygenome*)NULL;
      
      /* make a reference array for marker info */
      for (k=1; k<= chrom->nMark; k++) {
        mystartptr = chrom->mark_genome[k] = (mygenome*)malloc(sizeof(mygenome));
	copyGenomeTo(mystartptr, startptr, prev);
	chrom->mark_pos[k] = mystartptr->pos;
	prev = mystartptr;
	
        if (!startptr || startptr->chrom != idx) {
	  printf("Too few markers for chromosome %ld.  Expected %d, got %d",
		 idx, chrom->nMark, k);
	  exit(1);
        }
	
	/*--------------------------------------------------------
	 * may seem a waste since we already have this marker
	 * information in a structure. we again copy the entire 
	 * marker set.  However, we have the concept of flanking 
	 * markers, and now two ivectors contain ALL the relevant 
	 * marker information MOST (we have coded specially 
	 * for missing data!) of the time
	 *--------------------------------------------------------
	 */
	geno = chrom->mark_genome[k]->genotype = ivector(1,nn);
	for (i=1; i<= nn; i++) 
	  geno[i] = theparams->thedata[offset[i]].markers[idx][k];
	
	startptr = startptr->next;
      }
      
      
      /* capture some overall statistics (why not!) */
      /* keep some chromosome length statistics (needed for position proposal) */
      chrom->chromLen = chrom->mark_pos[chrom->nMark] - chrom->mark_pos[1];
      myData->chrom_pos[c+1] = myData->chrom_pos[c] + chrom->chromLen;
      
      myData->totalMark += chrom->nMark;
      myData->totalChromLen += chrom->chromLen;
      
      /* set up record of chromosomes on the QTL */
      chrom->qtls = (QTL_INFO**)malloc(MCMCMAXQTL * sizeof(QTL_INFO*));	 
      chrom->qtls--;
      for (i=1; i<= MCMCMAXQTL; i++) chrom->qtls[i]=(QTL_INFO*)NULL;
    }
}

void setupTraitData(DATA* myData, params* theparams )
{
  int nn, i, j;
  
  nn = theparams->nn;
  
  myData->y = dvector(1, theparams->nn);
  myData->offset = ivector(1,theparams->nn);
  
  for (i=1,j=1; i<=theparams->nn; i++)
    {
      /* copy y data */
      myData->y[j] = theparams->thedata[i].y[theparams->whichtrait];   
      myData->offset[j] = i;
      
      if (myData->y[j] >= MISS_VAL) 
	j++;
      else 
	nn--;
      
    }
  
  if (nn<0) {
    printf("Exiting, all trait data is missing\n"); exit(2);
  }
  
  myData->nn = j-1;
  assert(nn == myData->nn);
  
  calcMeanVar(myData->nn, myData->y, &myData->ybar, &myData->y_var);
  myData->y_stdev = sqrt(myData->y_var);
}  


















/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Bmain.c
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (mcmc.c) is part of QTL Cartographer
         
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

/************************************************************
    File:       mcmc.c
    Written by: Patrick Gaffney
    Date:       November 11, 2000
    Version:    0.4

    Purpose:

      Function MCMC gets the states of the required 
      Markov chain from its equilibrium distribution.
************************************************************/
#define PRE_BURN_IN .15   /* accelerated percent of burn-in     */
#define BURN_IN .10       /* burn-in is BURN_IN percent of total*/

#include "revjump.h"
#include "ranlib.h"
#include <assert.h>


void mcmc(MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, CHROMOSOME* chrInfo, WORK* myWork) 
{
  int move_type, num_qtl[MCMCMAXQTL+1];
  FPN LOD;

  FILE *write_res;
  int iter, i;
  int nn = myData->nn;
  FPN** nmove = dmatrix(1,5,1,2);
  QTL_INFO** all_qtls = myData->myQtls;
  int countedIter=0;
  int always_update = 1;                   /* for pre-burn-in */  
  /* a FPN negative, forces update of peramaters
     regardless of whether birth/death ocurred */
  int always_update_long = (myMCMC->revjump & MCMC_OPTIMIZE_BOOST);
  int rj_mcmc = (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD);
  QTL_INFO* lqtl;
  int k;

  write_res = fopen(myData->theparams->bayesfile, "w");
  fprintf(write_res, "niter nqtl iqtl chrom LOD mu sigmasq addvar domvar add dom locus esth");
  if (!rj_mcmc) fprintf(write_res, " HM SHM\n");
  else fprintf(write_res, "\n");


  /* adjust to get posterior */
  priors->sig_a1 += nn / 2.0;
  initVars(nn, nmove, num_qtl, myMCMC, priors->qtl_mean, 1.0);
  setupBurnIn (myMCMC, myData);

  
  outputStatistics1(nn, myMCMC, myData->theparams, priors, myData);
	               
  
  
  setupBirthDeathProbs(SELECT_NQTL_UNIFORM, priors->qtl_mean, 
		       0.9, &myData->bp, &myData->dp, &myData->prior_ratio);




  for( iter=(int)-myMCMC->burnIn; iter < myMCMC->niter; iter++)
    {

      /* restore the user's choice */
      if (iter == (int)(-myMCMC->burnIn + myMCMC->preBurnIn))
	{
	  if (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD)
	    setupBirthDeathProbs(myMCMC->revjump, priors->qtl_mean,   myMCMC->cval, &myData->bp, &myData->dp,  &myData->prior_ratio);
	  always_update = (myMCMC->revjump & MCMC_SELECT_100PCNT_UPDATE);
	}


      if (iter % 5000 == 0)  /* perform integrity checks */
	{
	  printf ("%d.",iter/1000); fflush(stdout);
	  checkResid(myData->nn, myData->nQtl, myData->mu, myData->y, myData->myQtls, myWork->resid, myData->gmiss);
	  for (i=1; i<= myData->nChrom; i++)
	    checkIntegrity(myData->nQtl, &chrInfo[i]);
	} 

      /****************************************
					       determine whether birth or death step
					       or whether to update parameters only
      ****************************************/

      /*********************************************************************
      **********************************************************************
      BIRTH AND DEATH STEPS ARE REQUIRED ONLY WHEN DOING REVERSIBLE
      JUMP MCMC COMPUTATIONS, AS THE FOLLOWING IF STATEMENT INDICATES.

      FOR REVERSIBLE JUMP MCMC, CHOOSE MOVE TYPE (BIRTH, DEATH OR
      STAY) AND THEN PROCEED. 
      **********************************************************************
      *********************************************************************/


      if (rj_mcmc)  {
	    move_type = select_move( myData->nQtl, myData->bp, myData->dp);
	    if (move_type != C_UPDATE)
	      birth_death(move_type, myData, myMCMC, myWork, priors, chrInfo, nmove);
	  }
      else
	    move_type = C_UPDATE;                         /* always update for fixed MCMC */


      if (always_update) move_type = C_UPDATE;        /* force 100% updates */
		

	

      if (move_type == C_UPDATE)
	    fixed_locus_update(nn, myData->nQtl, all_qtls, myMCMC, myData, priors, myWork, myData->theparams, nmove);
      else if (myData->nQtl >=1 && always_update_long) {
	    nmove[C_EXCHANGE][1]++;
          
	    lqtl = long_range_update(nn, myData->nQtl, all_qtls, myData->chrom_pos, myData->totalChromLen,  myData->nChrom, myData->myChroms, priors, myData, myWork, myMCMC);	
	    if (lqtl) {
	      nmove[C_EXCHANGE][2]++;
	      setValidFlag(lqtl, myMCMC->revjump);         
	    }
		  
	  } 



      /* record frequency of the number of QTL after burn-in (every entry) */
      if (iter >= 0) num_qtl[myData->nQtl]++;

      /* output results and possibly perform diagnostics */
      if(iter % myMCMC->nby == 0  && iter >= 0)
	{
	  countedIter++;
	  outputResults(write_res, nn, iter, myData, all_qtls, myWork->resid,  myWork->output_qtls, priors, myMCMC, myWork);
	  if (myMCMC->diag) 
	    diagnose(nn, all_qtls, myWork, myData, myMCMC);
	}


#ifndef MCMC_OMIT_DEBUG
      for (k=1; k<= myData->nQtl; k++)
	{
	  if ((!(EQUALS(all_qtls[k]->qptr->dist, all_qtls[k]->nextDist)) ||
	       !(EQUALS(all_qtls[k]->qptr->prev->dist, all_qtls[k]->prevDist))) 
	      && (all_qtls[k]->probValid && all_qtls[k]->chrom->nQtl == 1))
	    printf("miss\n");
	}
      checkResid(nn, myData->nQtl, myData->mu, myData->y, 
		 all_qtls, myWork->resid, myData->gmiss);
#endif


      if (myMCMC->debug)
	{
	  if(myData->nQtl > 0) 
	    LOD = get_lod(nn, myData->sigmasq, myData->ybar, myData->y_var, myWork->resid);
	  else
	    LOD = 0.0;
	  fprintf(myMCMC->debug, " .. LOD = %6.2lf\n", LOD);
	}

    }  /* end iter */

  outputStatistics2(nn, myMCMC, myData->theparams,  nmove, num_qtl, countedIter,  priors, myData);


  /* close results file */
  fflush(write_res);
  fclose(write_res);
}







void birth_death(int move_type, DATA* myData, MCMC_PARAM* myMCMC, WORK* myWork, PRIORS* priors, 
		 CHROMOSOME* chrInfo, FPN** nmove)
{
  int i;
  static char* moves[] = {"-","B","D","U"};


  if (myMCMC->debug)
    {
      fprintf (myMCMC->debug, "%d%s  [", myData->nQtl, moves[move_type]);
      for (i=1; i< myData->nQtl; i++)
	fprintf(myMCMC->debug, "{%d}%5.2f,",myData->myQtls[i]->chrom->num,
		myData->myQtls[i]->qptr->pos);
      fprintf(myMCMC->debug, "{%d}%5.2f] ",myData->myQtls[myData->nQtl]->chrom->num,
	      myData->myQtls[myData->nQtl]->qptr->pos);
    }


  nmove[move_type][1] += 1;

  if (move_type==C_BIRTH)
    {
      if (birth(myData, chrInfo, priors, myData->theparams, myData->individs, myWork, myMCMC))
	nmove[C_BIRTH][2] += 1;
    }
  else if (move_type==C_DEATH) 
    {
      if (death(myData, priors, myData->theparams, myWork, myMCMC))
	nmove[C_DEATH][2] += 1;
    }

  if (myData->nQtl >= 1 && 
      ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_RJMCMC))
    {
      nmove[C_SWAP][1]+=1;
      if (swap_add_dom(myData->nn, myData->nQtl, 
		       myData->myQtls, priors, myData, myWork, myMCMC))
	nmove[C_SWAP][2] += 1;
    }
				   

}




void fixed_locus_update(int nn, int nqtl, QTL_INFO** all_qtls,MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, WORK* myWork, params* theparams, FPN** nmove)
{
  int i; 
  QTL_INFO* lqtl;
  FPN sum;
  int simpleGibbs = myMCMC->revjump & MCMC_SELECT_SIMPLE_GIBBS;

  if ( theparams->nn <=0 )
    i =1;  

  if(nqtl > 1) 
    random_perm_QTL(nqtl, myData->myQtls, myWork->output_qtls, myWork->perm_num);

  /* we must update the sigmasq and effect_prior_var BEFORE update_effects 
     to retain the validity of the Cholesky decomposition 
  */
  /*  myData->sigmasq = MH_update_sigmasq(nn, myWork->resid,
      myData->sigmasq, 
      myData->y_var/20.0, 1.5 * myData->y_var);
  */
  myData->sigmasq = Gibbs_update_sigmasq(nn, myWork->resid, priors->sig_a1, priors->sig_a2);

  if (priors->sampleVar[ADD] || priors->sampleVar[DOM])
    {
      if (myMCMC->addParam == QTL_NONE || myMCMC->addParam == QTL_ADD_DOM)
	{
	  if (priors->sampleVar[ADD]) 
	    update_effect_prior_var(ADD, nqtl, all_qtls, priors, myData->y_var);
	  if (priors->sampleVar[DOM]) 
	    update_effect_prior_var(DOM, nqtl, all_qtls, priors, myData->y_var);
	}
      else if (priors->sampleVar[myMCMC->addParam])
	update_effect_prior_var(myMCMC->addParam, nqtl, all_qtls, priors, myData->y_var);
    }


  for (i=1, sum=0; i<=nqtl; i++) 
    {
      lqtl = all_qtls[i];	
      sum += update_lambda_qtl(lqtl, all_qtls, myData, myWork, myMCMC->revjump);

      if (!(EQUALS(lqtl->qptr->dist, lqtl->nextDist)) ||
	  !(EQUALS(lqtl->qptr->prev->dist, lqtl->prevDist)))
	printf("miss\n");


      if (simpleGibbs)
	update_effect(nn, lqtl, myData, priors, myWork);	
      setValidFlag(lqtl, myMCMC->revjump);         
	                        /* invalidate the cache of log_prob of QTL given flanking
				   markers, for safety if warranted (i,e. more than one
				   QTL on this chromosome, or if feature disabled) */

    }


  if (!simpleGibbs)
    {
      update_effects(nn, nqtl, myMCMC->revjump, 
		     myData->y, &myData->mu, myData->sigmasq, myData->na,
		     all_qtls, priors, myWork);

    }
  else
    myData->mu = update_mu(nn, myData->mu, myData->sigmasq, myWork->resid,
			   priors->mean[MU], priors->var[MU]);



  /* now propose long-range updates (if appropriate) */
  if (nqtl >= 1 && 
      !(myMCMC->revjump & MCMC_SELECT_MOVE_NONE)) 
    {	
      nmove[C_EXCHANGE][1]++;
      lqtl = long_range_update(nn, nqtl, all_qtls, 
			       myData->chrom_pos, myData->totalChromLen, 
			       myData->nChrom, myData->myChroms,
			       priors, myData, myWork, myMCMC);					  
      if (lqtl)
	{
	  nmove[C_EXCHANGE][2]++;
	  setValidFlag(lqtl, myMCMC->revjump);         

#ifndef MCMC_OMIT_DEBUG
	  if (!(EQUALS(lqtl->qptr->dist, lqtl->nextDist)) ||
	      !(EQUALS(lqtl->qptr->prev->dist, lqtl->prevDist)))
	    printf("miss\n");
#endif
	}


    }

		  
  /* keep tally */
  if (nqtl >= 1) 
    { 
      nmove[C_UPDATE][1]++;
      nmove[C_UPDATE][2] += sum/myData->nQtl;
    }
	
  /* verify */  
#ifndef MCMC_OMIT_DEBUG
  checkResid(nn, nqtl, myData->mu, myData->y, 
	     all_qtls, myWork->resid, myData->gmiss);
#endif

}






void update_effect(int nn,  QTL_INFO* lqtl, DATA* myData, PRIORS* priors, WORK* myWork)
{
  if (lqtl->flag & QTL_ADD)
    lqtl->a[QTL_ADD] = update_add_effect(nn, myData->nQtl, myData->sigmasq, lqtl, 
					 myWork->resid, priors->mean[QTL_ADD], 
					 lqtl->w[QTL_ADD] * priors->var[QTL_ADD]);

  if (lqtl->flag & QTL_DOM)
    lqtl->a[QTL_DOM] = update_dom_effect(nn, myData->nQtl, myData->sigmasq, lqtl, 
					 myWork->resid, priors->mean[QTL_DOM], 
					 lqtl->w[QTL_DOM] * priors->var[QTL_DOM]);
}


int numcmp(FPN *v1, FPN *v2)
{
  if(*v1<*v2)
    return -1;
  else if(*v1>*v2)
    return 1;
  else
    return 0;
}


int binSearch(int nval, FPN* vals, FPN searchVal)
     /* assumptions:
      searchVal lies between vals[1] and vals[nval]

    inputs:
	    nvals    :  number of intervals
		vals     :  array containing nval interval bounds
		searchVal:  value whose interval we're seeking
    returns:
	  the interval in which searchVal lies, i.e. 
	      i : vals[i]<=searchVal<vals[i+1]
*/
{
  int low, high, mid;
  low = 1;
  high = nval;
  do {
    mid = (low+high)/2;
    if (searchVal < vals[mid]) {
      if (high==mid) 
	break; 
      high=mid;
    }
    else {
      if (low==mid) break; 
      low=mid;
    }
  } while (1);

  return low;
}


FPN lodnull(int nn, FPN y_var)
{
  return -(nn/2)*log(y_var) - nn/2;
}



FPN get_lod(int nn, FPN sigmasq, FPN ybar, FPN y_var,
	       FPN* resid)
{
  FPN lod, null_lod, work, sum;
  static FPN logten = 2.302585;
  int i;

  sum = 0.0;
  for (i=1; i<= nn; i++) {
    work = resid[i];
    work *= work;
    sum += work;
  }

  lod = -(nn/2)* log(sigmasq) - 0.5*sum/sigmasq;
  null_lod = lodnull(nn, y_var);

  return (lod - null_lod)/logten;
}


void normalProb(int nn, FPN* resid, FPN sigmaSq, FPN* diag_fY,
		FPN* maxResid, FPN* pcnt_gt_99)
{
  /* increment the f(y) diagnostic by */

  static FPN logTwoPi = 1.837877066;
  FPN sqrtSigmaSq = sqrt(sigmaSq);
  int i;
  FPN z;
  FPN zlim;
		
  if (40*nn < 4000)	zlim = 1.412979 * pow(40.0*nn,0.1122);   
	                        /* approx to find z such that P(z<-zlim or z>zlim) = p = 1/(20*nn)
	                           reasonable for 0.01 < p < 0.0005 */
  else zlim = 2.221161 * pow(40.0*nn, 0.05646);
	                        /* approx to find z such that P(z<-zlim or z>zlim) = p = 1/(20*nn)
	                           reasonable for 0.0005 < p < 0.000001 */
  if (zlim < 2) zlim=2;

  for (i=1; i<=nn; i++)
    {
      z = resid[i]/sqrtSigmaSq;
      if (fabs(z) >  fabs(maxResid[i]))  maxResid[i] = z;

		
      if (fabs(z) > zlim) pcnt_gt_99[i] += 1.0;

      diag_fY[i] += exp(-logTwoPi -resid[i]*resid[i]/(2*sigmaSq))/sqrtSigmaSq;
    }
}




void diagnose(int nn, QTL_INFO** all_qtls, WORK* myWork, DATA* myData, MCMC_PARAM* myMCMC)
{
  int i,j;
  QTL_INFO* lqtl;
  int* geno;
  FPN pr;

  normalProb(nn, myWork->resid, myData->sigmasq, 
	     myMCMC->diag_fY, myMCMC->maxResid, myMCMC->pcnt_gt_99);
  for (i=1; i<=nn; i++)
    {
      myWork->avgMinProbs[i]=0.0;
      myWork->minProbs[i] = 0;
    }

  /* now compute a few statistics re. QTL prob given markers, across QTLs */
  for (j=1; j<=myData->nQtl; j++)
    {
      lqtl = all_qtls[j];
      geno = igenotype(lqtl);

      /* calculate prob of QTL genos given flanking markers */
      genProbs(nn, myData->offset, myData->theparams, myData->individs, 
	       myData->gmiss, lqtl, myMCMC->revjump);
					       

      lqtl->probValid = 0;

      for (i=1; i<=nn; i++)
	{		  
	  pr = lqtl->log_prob[j][geno[i]];
	  if (pr < myWork->minProbs[j]) 
	    myWork->minProbs[j] = pr;
	  myWork->avgMinProbs[j] += exp(pr);				      
	}
    }
		 
  if (myData->nQtl > 0) 
    {
      for (j=1; j<=nn; j++)
	{
	  myMCMC->min_prob_Q_given_M[j] += exp(myWork->minProbs[j]);
	  myMCMC->avg_prob_Q_given_M[j] += myWork->avgMinProbs[j]/myData->nQtl;
	}
    }
}	




void outputResults(FILE* write_res, int nn,
		   int iter, DATA* myData, QTL_INFO** all_qtls, FPN* resid, 
		   QTL_INFO** output_qtls, PRIORS* priors, MCMC_PARAM* myMCMC,
		   WORK* myWork)
{
  int iqtl,j;
  QTL_INFO* lqtl;
  FPN H, LOD;
  FPN residSS, val;
  static FPN logTwoPi = 1.837877066;
  int rj_mcmc = (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD);

  /* order QTLs for output */
  for(iqtl=1; iqtl<=myData->nQtl; iqtl++)
    {
      /* look for place to insert QTL in ordered sequence*/ 
      for(j=iqtl; ; j--)
	{
	  if (j==1 ||
	      all_qtls[iqtl]->chrom->num > output_qtls[j-1]->chrom->num ||
	      (all_qtls[iqtl]->chrom->num == output_qtls[j-1]->chrom->num &&
	       all_qtls[iqtl]->qptr->pos > output_qtls[j-1]->qptr->pos))
	    {
              output_qtls[j] = all_qtls[iqtl];
	      break;
	    }
	  else
	    output_qtls[j] = output_qtls[j-1];  /* open a space */
	}
    }

  if(myData->nQtl > 0) 
    {
      LOD = get_lod(nn, myData->sigmasq, myData->ybar, myData->y_var, resid);
      H = calc_h2(nn, myData->y, myData->mu, myData->sigmasq, resid);
    }
  else
    LOD = 0.0;

  if (!rj_mcmc)
    {
      myMCMC->idx++;
      residSS = calcResidSS(nn, myWork->resid);
      val = ((nn/2.0) * (logTwoPi + log(myData->sigmasq))) + 
	residSS/(2.0*myData->sigmasq);
      updateMean(myMCMC->idx, exp(val), &myMCMC->HM);

      val = priors->sig_a1 * log(1 + residSS/(2*priors->sig_a2));
      val += gammln2(priors->sig_a1 - (nn/2.0));
      val += ((nn/2.0)*(logTwoPi + log(priors->sig_a2)));
      val -= gammln2(priors->sig_a1);
      updateMean(myMCMC->idx, exp(val), &myMCMC->SHM);		
    }

  /* rewrite by Yandell in new format 30 mar 1999, PGA 11-18-2000 */
  if(myData->nQtl==0)
    {
      fprintf(write_res, "%d 0 1 . %15.4f %12.4f %12.4f ",
	      iter, LOD, myData->mu, myData->sigmasq);
      if (myMCMC->addParam == QTL_ADD_DOM || myMCMC->addParam == QTL_NONE)
	fprintf(write_res,"%12.4f %12.4f ", priors->var[QTL_ADD]/myData->y_var, 
		priors->var[QTL_DOM]/myData->y_var);
      else if (myMCMC->addParam == QTL_ADD)
	fprintf(write_res,"%12.4f . ", priors->var[QTL_ADD]/myData->y_var);
      else
	fprintf(write_res,". %12.4f ", priors->var[QTL_DOM]/myData->y_var);
      if (!rj_mcmc) fprintf(write_res,". . . 0 %12.4f %12.4f\n",myMCMC->HM, myMCMC->SHM);
      else fprintf(write_res,". . . 0\n");
    }
  else
    for(iqtl=1; iqtl<=myData->nQtl; iqtl++)
      {
	lqtl = output_qtls[iqtl];

	/* output preamble ... including number of chromosome, LOD, mu, sigmasq */
	if (iqtl == 1)
	  {
	    fprintf(write_res, "%d %d %d %d %15.4f %12.4f %12.4f", 
		    iter, myData->nQtl, iqtl, lqtl->chrom->num,
		    LOD, myData->mu, myData->sigmasq);
     	    if (myMCMC->addParam == QTL_ADD_DOM || myMCMC->addParam == QTL_NONE)
	      fprintf(write_res,"%12.4f %12.4f ", priors->var[ADD]/myData->y_var, 
		      priors->var[DOM]/myData->y_var);
	    else if (myMCMC->addParam == QTL_ADD)
	      fprintf(write_res,"%12.4f . ", priors->var[QTL_ADD]/myData->y_var);
	    else
	      fprintf(write_res,". %12.4f ", priors->var[QTL_DOM]/myData->y_var);
	  }
	else
	  fprintf(write_res, "%d %d %d %d . . . . . ", iter, myData->nQtl, 
		  iqtl, lqtl->chrom->num);

	/* output effects */
	if (lqtl->flag & QTL_ADD)
	  fprintf(write_res, " %12.4f",lqtl->a[QTL_ADD]);
	else 
	  fprintf(write_res, ". ");

	if (lqtl->flag & QTL_DOM)
	  fprintf(write_res, " %12.4f",lqtl->a[QTL_DOM]);
	else 
	  fprintf(write_res, " .");

	/* for first QTL (with nQtl>0) output QTL pos and heriability */
	if (iqtl!=1)
	  {
	    if (!rj_mcmc) fprintf(write_res, " %10.4f . . .\n", lqtl->qptr->pos*100);
	    else fprintf(write_res, " %10.4f .  \n", lqtl->qptr->pos*100);
	  }
	else if (!rj_mcmc) 
	  fprintf(write_res, " %10.4f %8.4f %12.4f %12.4f\n",lqtl->qptr->pos*100,  H,
		  myMCMC->HM, myMCMC->SHM);
	else
	  fprintf(write_res, " %10.4f %8.4f\n",lqtl->qptr->pos*100,  H);
      }
}



void outputStatistics1(int nn, MCMC_PARAM* myMCMC, params* theparams, 
		       PRIORS* priors, DATA* myData)
{
  char oFile[256];
  FILE* write_nmove;

  strcpy(oFile, theparams->bayesfile);
  strcat(oFile, ".mov");

  write_nmove = fopen(oFile, "w");

  fprintf(write_nmove,"QtlNew Version 6.01  4-4-2001\n"
	  "===========================\n\n");
  fprintf(write_nmove,"MCMC Parameters:\n\tBurnin:    \t%d\n\t"
	  "PreBurnIn: \t%d\n\t"
	  "Iters:     \t%d\n\t"
	  "StepSize:  \t%d\n\t"
	  "cval:      \t%6.4f\n\t"
	  "seed:      \t%d\n\t"
	  "ifile:     \t%s\n\t"
	  "bfile:     \t%s\n\t"
	  "mfile:     \t%s\n",
	  (int)myMCMC->burnIn, (int)myMCMC->preBurnIn, 
	  myMCMC->niter, myMCMC->nby, myMCMC->cval, 
	  theparams->seed, theparams->ifile, 
	  theparams->bayesfile, theparams->map);

  fprintf(write_nmove,"Prior Parameters:\n\t"
	  "qtl mean:  \t%f\n\t"
	  "Distrib:   \t%d\n\t"
	  "revjump:   \t%X\n\t"
	  "a1:        \t%f\n\t"
	  "a2:        \t%f\n\t"
	  "Mu_mean:   \t%f\n\t"
	  "Mu_stdev:  \t%g\n\t"
	  "A_Eff_mu:  \t%f\n\t"
	  "Eff_stdev: \t%g\n\t"
	  "D_Eff_mu:  \t%f\n\t"
	  "Eff_stdev: \t%g\n\t"
	  "A_alpha:   \t%g\n\t"
	  "A_beta:    \t%g\n\t"
	  "D_alpha:   \t%g\n\t"
	  "D_beta:    \t%g\n",
	  priors->qtl_mean,priors->priorDistribution, 
	  myMCMC->revjump,
	  priors->sig_a1, priors->sig_a2,
	  priors->mean[MU], sqrt(priors->var[MU]),	                  
	  priors->mean[ADD], sqrt(priors->var[ADD]),
	  priors->mean[DOM], sqrt(priors->var[DOM]),
	  priors->alpha[ADD], priors->beta[ADD],
	  priors->alpha[DOM], priors->beta[DOM]);

  fprintf(write_nmove,"Data Parameters:\n\t"
	  "Y_mean:    \t%g\n\t"
	  "Y_stdev:   \t%g\n\t"
	  "N:         \t%d\n\t"
	  "nChrom:    \t%d\n\t"
	  "TotalMark: \t%d\n\t"
	  "Length:    \t%f\n\t",
	  myData->ybar, myData->y_stdev, 
	  myData->nn,myData->nChrom,
	  myData->totalMark, myData->totalChromLen);
  
  fclose(write_nmove);
}



void outputStatistics2(int nn, MCMC_PARAM* myMCMC, params* theparams, 
		       FPN** nmove, int* num_qtl,
		       int countedIter, PRIORS* priors, DATA* myData)
{
  int i;
  char oFile[256];
  FILE* write_nmove;
  FPN hm = 1/myMCMC->HM;
  FPN shm = 1/myMCMC->SHM;

  strcpy(oFile, theparams->bayesfile);
  strcat(oFile, ".mov");

  if (myMCMC->diag && countedIter)
    {
      fprintf(myMCMC->diag, "fY   pcnt99   maxZ   avgQ_M  minQ_M\n");
      for (i=1; i<=nn; i++)
	fprintf(myMCMC->diag, "%7.4f %7.4f %7.4f %7.4f %7.4f\n", 
		myMCMC->diag_fY[i]/countedIter,
		myMCMC->pcnt_gt_99[i]/countedIter, myMCMC->maxResid[i],
		myMCMC->avg_prob_Q_given_M[i]/countedIter, 
		myMCMC->min_prob_Q_given_M[i]/countedIter);
    }


  write_nmove = fopen(oFile, "a");
  fprintf(write_nmove,"HM:        \t%g\n\t"
	  "SHM:       \t%g\n",
	  hm, shm);

  fprintf(write_nmove,"\n\n");
  fprintf(write_nmove,"          birth   death   locus   exchange\nselect ");
  for(i=1; i<=5; i++)
    fprintf(write_nmove, " %9.2f", nmove[i][1]);
  fprintf(write_nmove,"\naccept ");
  for(i=1; i<=5; i++)
    fprintf(write_nmove, " %9.2f", nmove[i][2]);
  fprintf(write_nmove,"\n\nnqtl frequency\n");
  for(i=0; i<=MCMCMAXQTL; i++)
    fprintf(write_nmove, "%4d %d\n", i, num_qtl[i]);

  fclose(write_nmove);
}




FPN calc_h2(int nn, FPN* y, FPN mu, FPN sigmasq, FPN* resid)
{
  int i;
  FPN sum, mean, add;

  if (nn <=0) return 0;

  /* calculate the genetic variance, i.e. var(gene effects) */
  for (i=1, sum=0.0; i<= nn; i++) 
    sum += y[i] - mu - resid[i] ; 
  mean = sum/ nn;

  for (i=1, sum=0.0; i<= nn; i++) 
    {
      add = y[i] - mu - resid[i] - mean; 
      sum += add*add;
    }
  sum /= nn;       /* estimate for genetic variance */

  return sum / (sum + sigmasq);
	
}


void initVars(int nn, FPN** nmove, int* num_qtl, MCMC_PARAM* myMCMC,
	      FPN qtl_prior_mn, FPN cval)
{
  int i;

  for(i=1; i<=2; i++)
    nmove[C_BIRTH][i] = nmove[C_DEATH][i] = nmove[C_UPDATE][i] =
      nmove[C_EXCHANGE][i] = nmove[C_SWAP][i] = 0.0;

  for (i=0; i<=MCMCMAXQTL; i++)
    num_qtl[i] = 0;

  if (myMCMC->diag)
    {
      for (i=1; i<=nn; i++)
	{
	  myMCMC->diag_fY[i]=0.0;
	  myMCMC->maxResid[i]=0.0;
	  myMCMC->avg_prob_Q_given_M[i]=0.0;
	  myMCMC->min_prob_Q_given_M[i]=0.0;
	  myMCMC->pcnt_gt_99[i]=0.0;
	}
    } 
}





void setupBurnIn (MCMC_PARAM* myMCMC, DATA* myData)
{
  FPN ratio, burnIn, preBurnIn;;

  
  /* BurnIn: a ratio of the number of iterations (10%) */
  if (myMCMC->burnIn <0) burnIn = BURN_IN * myMCMC->niter;
  else burnIn = myMCMC->burnIn * myMCMC->niter;

  if (myMCMC->preBurnIn <0) preBurnIn = myData->nChrom * 100;   
  else preBurnIn = myMCMC->preBurnIn * burnIn;

  if (myMCMC->burnIn < 0 || myMCMC->preBurnIn < 0)
    {
      /* PreBurnIn: get 50 birth proposals per chromosome on average */

      /* now check if PreBurnIn not long enough */
      ratio = preBurnIn/(burnIn * PRE_BURN_IN);
      if (ratio > 1) 
	{
	  burnIn *= ratio; 
	  myMCMC->niter = (int)(myMCMC->niter * ratio); 
	}
      else if (ratio < 1)
	{
	  preBurnIn = PRE_BURN_IN * burnIn;
	}
    }


  /*  printf("%lf %lf\n",burnIn, preBurnIn);
   */  myMCMC->burnIn = floor(burnIn);
  myMCMC->preBurnIn = floor(preBurnIn);
  /*  printf("%lf %lf\n",myMCMC->burnIn, myMCMC->preBurnIn);*/

  printf("preBurnIn: %d\nBurnIn:    %d\n",(int)myMCMC->preBurnIn,(int)myMCMC->burnIn);
}

/*
            h2 = LOD * log(10);
            h2 = 2*(h2-1)/(h2+myData->nn-1);
*/

void calcMeanVar(int nn, FPN* vals, FPN* mean, FPN* var)
{
  int i;
  FPN sum,t;

  for (i=1, sum=0.0; i<=nn; i++)
    sum += vals[i];
  *mean = sum/nn;

  for (i=1, sum=0.0; i<=nn; i++)
    {
      t = vals[i] - *mean;
      sum += t*t;
    }
  *var = sum/nn;
}


void updateMean(int i, FPN val, FPN* mean)
     /* i=1,2,.... 
   user must initialize mean=0, so when called with i=1, mean= val */
{
  FPN delta = (val - *mean) / i;
  *mean += delta;
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file mcmc.c
------------------------------------------------------------------ */


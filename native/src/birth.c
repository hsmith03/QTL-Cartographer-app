/* ------------------------------------------------------ XCutXCodeXSkip
     This file (birth.c) is part of QTL Cartographer
         
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

/****************************************************
  File:         birth.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:

    Function BIRTH performs the "birth step" for
    adding a new QTL to the model.
****************************************************/

#include "revjump.h"
#include "ranlib.h"

int birth(DATA* myData, CHROMOSOME* chromInfo, PRIORS* priors, 
		   params* theparams, individual* individs, WORK* myWork, 
		   MCMC_PARAM* myMCMC)
{
  FPN new_position;
  int lmark;
  FPN uni_ran; 
  long int nn = myData->nn;
  int nQtl = myData->nQtl;
  FPN sigmasq = myData->sigmasq;
  QTL_INFO** all_qtls = myData->myQtls; 
  QTL_INFO* newQtl;
  CHROMOSOME* chrom;
  FPN log_accept_birth_prob, log_effect_ratio, log_proposal, log_position;
  FPN bCb, d_invD_d, log_det_Chol, log_det_invD;

  /* work structures we'll use */
  FPN* mod_effect = myWork->mod_effect;    /* new grand mean and effects */
  FPN** XtX = myWork->XtX;                 /* XtX matrix */
  FPN* XtY = myWork->XtY;                 /* XtX vector */
  FPN** chol = myWork->chol;               /* new Cholesky decomposition */
  FPN* p = myWork->p;                      /* new diagonal for the Cholesky */
  FPN* pmean = myWork->pmean;              /* mean of grand mean/effects */
  FPN* pvar = myWork->pvar;                /* variance of grand mean/effects */
  FPN* u = myWork->u;                      /* iid normal values */
  FPN* weight = myWork->weight;            /* new weights */


  /* GET NEW BIRTH INTERVAL */
  get_new_locus(myMCMC->revjump, myData->nChrom, myData->chrom_pos, 
	            myData->totalChromLen, chromInfo, &chrom, &lmark, &new_position);

  newQtl = createQtl(nn, nQtl+1, &all_qtls[nQtl+1], chrom, 
	                 lmark, new_position, myMCMC->addParam, (FPN*)NULL, (FPN*)NULL);
  
                                  

  /* GET GENOTYPES AT NEW POSITION */
  initQtl(myData->nn, myData->offset, newQtl, 
	      theparams, individs, myData->gmiss, myMCMC->revjump);


  log_effect_ratio = get_effect(nn, nQtl, myData->y, myData->mu, sigmasq,
                                NULL, newQtl, all_qtls, myMCMC->revjump,priors, myData->na, mod_effect,
  			                    XtX, XtY, chol, p, pmean, pvar, u, weight,
			                    &bCb, &d_invD_d, &log_det_Chol, &log_det_invD, myWork);

                     
  
/************************************************
  If jacobian = 0, then acceptance probability 
  is 0.
************************************************/
  log_proposal = get_log_proposal_ratio(nQtl, myData->bp, myData->dp, myData->prior_ratio);
  log_position = get_log_position_ratio(myMCMC->revjump, newQtl->chrom, myData);

  log_accept_birth_prob = log_proposal + log_position + log_effect_ratio;
										 

    uni_ran = genunf(0.0,1.0);
    uni_ran = log(uni_ran);

    if(uni_ran < log_accept_birth_prob){
	  setEffect(nn, nQtl+1, myData->y, myWork->mod_effect, all_qtls, &myData->mu, 
		        myWork->weight, myWork->resid, myData->na);
      setCholParams(myWork, bCb, d_invD_d, log_det_Chol, log_det_invD);
		        
      myData->nQtl++;
	  newQtl->chrom->nQtl++;
      setValidFlag(newQtl, myMCMC->revjump);         

#ifndef MCMC_OMIT_DEBUG
     checkResid(myData->nn, myData->nQtl, myData->mu, myData->y, 
		         myData->myQtls, myWork->resid, myData->gmiss);
	 checkIntegrity(myData->nQtl, newQtl->chrom);
#endif


	  return 1;
    }
	else {

	  removeQtl(newQtl->qptr);
#ifndef MCMC_OMIT_DEBUG
	 checkIntegrity(myData->nQtl, newQtl->chrom);
#endif
	  return 0;
	}
}




void get_new_locus(int revjump, int nChrom, FPN* chrom_pos, FPN totalChromLen, 
				   CHROMOSOME* chromInfo, CHROMOSOME** chrom, 
				   int* lmark, FPN* pos)
{

  /* alternative coding */
  int c;

  if (revjump & MCMC_FLAG_RANDOM_CHROM)
  {
	  c = ignuin(1,nChrom);
      *chrom = &chromInfo[c];
	  *pos = genunf(0,(*chrom)->chromLen);
  }
  else
  {
     *pos = genunf(0.0, totalChromLen);
     c = binSearch(nChrom+1, chrom_pos, *pos);
     *pos -= chrom_pos[c];
     *chrom = &chromInfo[c];
  }


  *lmark = binSearch((*chrom)->nMark, (*chrom)->mark_pos, *pos);
  
}


void get_local_locus(CHROMOSOME* chrom, int* lmark, FPN* pos)
{
  *pos = genunf(0.0, chrom->chromLen);
  *lmark = binSearch(chrom->nMark, chrom->mark_pos, *pos);
}














void setEffect(int nn, int nQtl, FPN* y, FPN* mod_effect, QTL_INFO** all_qtls, 
			   FPN* mu, FPN* w, FPN* resid, int* na)
{
  /* assign all qtls the mod_effect, sets the weights using 'w' and evaluates the
	 number of additive and domiance effects, na[ADD] and na[DOM] respectively */

	int j, idx;
	QTL_INFO* thisqtl;
	int type;

	*mu = mod_effect[1];
	na[ADD]=na[DOM]=0;

	for (j=1, idx=2; j<= nQtl; j++) 
	{
		thisqtl = all_qtls[j];
		for (type=QTL_ADD; type <= QTL_DOM; type++)
		   if (thisqtl->flag & type) 
		   {
			  na[type]++;
		  	  thisqtl->a[type] = mod_effect[idx];
			  thisqtl->w[type] = w[idx];
			  idx++;
		   }
	}

	/* finally calculate the residuals */
	calcResid2(nn, nQtl, y, mod_effect, all_qtls, resid);
}




QTL_INFO* createQtl(int nn, int qtlNum, QTL_INFO** p_newQtl, CHROMOSOME* chrom, 
					int lmark, FPN new_position,
					int addParam, FPN* a, FPN* w)
{

  QTL_INFO* newQtl;
  int i;
  FPN u;
   

  if (!*p_newQtl) {
	  newQtl = (*p_newQtl) = (QTL_INFO*)malloc(sizeof(QTL_INFO));
	  newQtl->qptr = (mygenome*)malloc(sizeof(mygenome));
	  newQtl->qptr->genotype = ivector(1,nn);
      newQtl->log_condProb = (FPN***)malloc(3*sizeof(FPN**));
      newQtl->log_condProb++;
      for (i=-1;i<=1;i++) newQtl->log_condProb[i] = dmatrix(-1,1,-1,1);
      newQtl->transProbR = dmatrix(0,2,0,2);
      newQtl->transProbL = dmatrix(0,2,0,2);
	  newQtl->log_prob = (FPN**)malloc(nn * sizeof(FPN*));
	  newQtl->log_prob--;
	  for (i=1; i<=nn; i++) newQtl->log_prob[i] = (FPN*)NULL;
	  newQtl->missing_prob = dmatrix(1,nn,-1,1);
      newQtl->a = dvector(QTL_ADD, QTL_DOM);
      newQtl->w = dvector(QTL_ADD, QTL_DOM);
  }
  else
	  newQtl = *p_newQtl;

  newQtl->chrom = chrom;

  /* initialize values to that required by setBirthWeights */
  newQtl->a[QTL_ADD] = (a==(FPN*)NULL) ? 0: a[QTL_ADD];
  newQtl->a[QTL_DOM] = (a==(FPN*)NULL) ? 0: a[QTL_DOM];
  newQtl->w[QTL_ADD] = (w==(FPN*)NULL) ? 1: w[QTL_ADD];
  newQtl->w[QTL_DOM] = (w==(FPN*)NULL) ? 1: w[QTL_DOM];

  if (addParam == QTL_NONE)        /* reversible jump over additive/dominance effects */
  {
	  u = RANF();
	  if (u < 3.0/7) 
		  newQtl->flag = QTL_ADD;
      else if (u < 4.0/7)
		  newQtl->flag = QTL_ADD_DOM;
	  else
		  newQtl->flag = QTL_DOM;
  }
  else
	  newQtl->flag = addParam;

  newQtl->nParam = (newQtl->flag == QTL_ADD_DOM)? 2:1;
  newQtl->nInt = 0;


  setLambda(newQtl, lmark, new_position, qtlNum); 

  /* set to flag that tables need be calculated */
  newQtl->prevDist = -1.0;
  newQtl->nextDist = -1.0;
  newQtl->probValid = 0;


  return newQtl;
}














FPN  get_effect(long int nn, int nQtl, 
				           FPN* y, FPN mu, FPN sigmasq,
				           QTL_INFO* oldQtl, QTL_INFO* newQtl, 
						   QTL_INFO** qtls, int revjump,
				           PRIORS* priors, 
				           int* old_na, FPN* mod_effect,
						   FPN** XtX, FPN* XtY, FPN** chol, FPN* p, 
						   FPN* pmean, FPN* pvar, FPN* u, FPN* weight,
						   FPN* bCb, FPN* d_invD_d, 
                      	   FPN* log_det_Chol, FPN* log_det_invD,
						   WORK* myWork)
/*=====================================================================

   Description
   -----------

   Inputs
   ------
         nn      ... the number of inidividuals
         nQtl    ... number of QTLs in original model (before birth, death, etc.)
         y       ... their trait values
         sigmasq ... model environmental variance
		 oldQtl  ... pointer to old QTL (for death or update ... NULL for birth)
		 newQtl  ... pointer to new QTL (for birth or update ... NULL for death)
         qtls    ... an array containing pointers to the nQtl QTLs
         myWork  ... structure containing all relevant          
					                and dominance components 
		 priors  ... prior data
         na      ... number of effects in nQtl model

  
   Outputs
   -------
    	 myWork->w          .. weights for new model            
	     myWork->mod_effect .. new effects

   Returns 
   -------
		 logarithm of effect ratio (part of acceptance ratio) correctly adjusted for birth
		 or death.

         
      
    
  =====================================================================*/
{
	int i, info;
	int na[3];
	int getFischProposal = (revjump & SELECT_SAMPLE_FISCH);

	if (oldQtl && newQtl && getFischProposal) return -1E30;


	for (i=1; i<=2; i++) na[i] = old_na[i];  /* protect the na vector */

   
    #ifndef MCMC_OMIT_DEBUG 
	if (!getFischProposal)
	     checkCholesky(nQtl,revjump, qtls, priors, weight,na, pmean,pvar,
	                   XtX, XtY, sigmasq,chol, p, u, mod_effect,					    
                       bCb, d_invD_d, log_det_Chol, log_det_invD, myWork);
	#endif


   if (oldQtl) 
   {
	/* set the y[k] times genotype values for existing QTL, and the
       interaction terms between existing QTL and new one */
   	    if (oldQtl->flag & ADD) na[ADD]--;
        if (oldQtl->flag & DOM) na[DOM]--;
		nQtl--;

        if (getFischProposal)
           return proposeFischDeath(nn, nQtl, y, mu, sigmasq, qtls,
		  				            revjump, priors, na, mod_effect, 
							 	    pmean, pvar, weight, 
									myWork->oldVar, myWork->oldWeight,
									myWork->resid, myWork->newResid);
   }



   if (newQtl)
   {
  	  /* set the sum of y[k] times the genotype values for new QTL, and
	     also the new covariance matrix values diagonal terms.
	     IMPORTANT.. nQtl now is the number of valid QTL in qtls array.
	     newQtl is not qtls[nQtl+1], in general.
	   */
 	   addColToAddDom(nn, nQtl, na, qtls, newQtl,y, XtY, XtX);					       		

	   if (newQtl->flag & ADD) na[ADD]++;
	   if (newQtl->flag & DOM) na[DOM]++;
	   nQtl++;

       if (getFischProposal)
          return proposeFischBirth (nn, nQtl, y, mu, sigmasq, newQtl, qtls,
						  revjump, priors, na, mod_effect, pmean,
						  pvar, weight, 
						  myWork->oldVar, myWork->oldWeight,
						  myWork->resid, myWork->newResid);
   }

   /* update the cholesky with new genotype value(s) */
    setPriorMeanVar(nQtl,revjump, qtls, newQtl, priors, weight, na, pmean,pvar);
    info = generate_effects(na[ADD]+na[DOM]+1, XtX, XtY, 
 	                 pmean, pvar, sigmasq,
					 chol, p, u, mod_effect,					    
                     bCb, d_invD_d, log_det_Chol, log_det_invD, 1);

  if (info)
  {
	printX(nn,nQtl,qtls);
	printf("It's been bad, laddie...\n");
  }


    return 0.5 * (na[ADD]+na[DOM]-old_na[ADD]-old_na[DOM]) * log(sigmasq) 
		   + 0.5 * (*bCb - myWork->bCb) 
		   - 0.5*(*d_invD_d - myWork->d_invD_d)
		   - (*log_det_Chol - myWork->log_det_Chol)   /* recall these are inverses */
		   + 0.5 * (*log_det_invD - myWork->log_det_invD);

}








QTL_INFO* long_range_update(int nn, int nQtl, QTL_INFO** qtls,
					  FPN* chrom_pos, FPN totalChromLen,
					  int nChrom, CHROMOSOME* chromInfo, PRIORS* priors, 		              
					  DATA* myData, WORK* myWork, MCMC_PARAM* myMCMC)
/*=======================================================================================
 *
 *
 *
 *
 *  Exchange with any QTL in genome:
 *		 long_range_update(nn, nqtl, thisqtl, all_qtls,
 *					   myData->chrom_pos, myData->totalChromLen, 
 *					   myData->nChrom, myData->myChroms,
 *		               priors, theparams, individs, myData, myWork, myMCMC);					  
 *
 *
 *  Update on one chromosome only:
 *		 long_range_update(nn, nqtl, thisqtl, all_qtls,
 *					   NULL, 0, 0, NULL,
 *		               priors, theparams, individs, myData, myWork, myMCMC);					  
 *
 *=========================================================================================
 */
		              
{
  FPN log_accept_update_prob, new_position, log_effect_ratio ;
  FPN uni_ran; 
  int lmark, i;

  /* things for the new QTL */
  CHROMOSOME* chrom;
  QTL_INFO* oldQtl;
  QTL_INFO* thisqtl;

  /* work structures we'll use */
  FPN* mod_effect = myWork->mod_effect;    /* new grand mean and effects */
  FPN** new_XtX = myWork->new_XtX;         /* new XtX matrix (copy of XtX, then adjust) */
  FPN** chol = myWork->chol;               /* new Cholesky decomposition */
  FPN* new_XtY = myWork->new_XtY;          /* new XtX matrix (copy of XtX, then adjust) */
  FPN* p = myWork->p;                      /* new diagonal for the Cholesky */
  FPN* pmean = myWork->pmean;              /* mean of grand mean/effects */
  FPN* pvar = myWork->pvar;                /* variance of grand mean/effects */
  FPN* u = myWork->u;                      /* iid normal values */
  FPN* weight = myWork->weight;            /* new weights */
  params* theparams = myData->theparams;
  individual* individs = myData->individs; 

  /* some variables we need for the acceptance ratio */
  FPN bCb, d_invD_d, log_det_Chol, log_det_invD;
  int nterm = myData->na[ADD] + myData->na[DOM] + 1;
  int changeQtl;

  changeQtl = ignuin(1,nQtl);
  if (changeQtl != nQtl)
     moveQtlToEndOfXtX(nQtl, qtls, changeQtl, myWork->XtX, myWork->XtY, nterm);					 

  /* remove selected QTL from record structure */
  thisqtl = qtls[nQtl];
  removeQtl(thisqtl->qptr);  
  thisqtl->chrom->nQtl--;
    
  /* GET NEW BIRTH INTERVAL */
  if (chromInfo && chrom_pos && nChrom)
     get_new_locus(myMCMC->revjump, nChrom, chrom_pos, totalChromLen,
	               chromInfo, &chrom, &lmark, &new_position);
  else
  {
	 chrom = thisqtl->chrom;
	 get_local_locus(chrom, &lmark, &new_position);
  }

 

  mycopy(nterm, nterm,new_XtX, myWork->XtX);
  for (i=1; i<= nterm; i++) new_XtY[i] = myWork->XtY[i];
                                  /* don't want to corrupt the XtX unless we make the move */


  /* create and initialize the new QTL, we call it oldQtl since we will swap data
     with that of the oldQtl */
  oldQtl = createQtl(nn, nQtl, &qtls[0], chrom, lmark, new_position, 
	               thisqtl->flag, thisqtl->a, thisqtl->w);
  oldQtl->chrom->nQtl++;

  swapQtlData(oldQtl, thisqtl);  /* swap the records ... so keep the old data safe and place
                                 new record at the end of the list */

  
  initQtl(nn, myData->offset, thisqtl, 
	      theparams, individs, myData->gmiss, myMCMC->revjump);

  /* now get the effect + params required for the acceptance ratio */
  log_effect_ratio = get_effect(nn, nQtl, myData->y, myData->mu,
	                            myData->sigmasq, oldQtl, thisqtl, 
                                qtls, myMCMC->revjump, priors, myData->na,mod_effect,
 			                    new_XtX, new_XtY,chol, p, pmean, pvar, u, weight,
			                    &bCb, &d_invD_d, &log_det_Chol, &log_det_invD, myWork);
							    

		
    log_accept_update_prob =  log_effect_ratio;
    uni_ran = genunf(0.0,1.0);
    uni_ran = log(uni_ran);

    if(uni_ran < log_accept_update_prob)
	{
	  Swap2DTable(&myWork->XtX, &myWork->new_XtX);
	  SwapDVec(&myWork->XtY, &myWork->new_XtY);
      setCholParams(myWork, bCb, d_invD_d, log_det_Chol, log_det_invD);

	  

	  /* important to have swapped data before setting effect, since now
	     the first nQtl elements of 'qtls' hold the real QTLs */
	  setEffect(nn, nQtl, myData->y, myWork->mod_effect, qtls, &myData->mu, myWork->weight, 
		        myWork->resid, myData->na);

#ifndef MCMC_OMIT_DEBUG
	  checkIntegrity(nQtl, thisqtl->chrom);
	  if (oldQtl->chrom != thisqtl->chrom) checkIntegrity(nQtl, oldQtl->chrom);
      checkResid(nn, nQtl, myData->mu, myData->y, 
		         qtls, myWork->resid, myData->gmiss);
#endif
	  return thisqtl;
    }
	else 
	{	  
      removeQtl(thisqtl->qptr);           /* return our genome data structures to normal */
      restoreQtl(oldQtl->qptr);            
	  thisqtl->chrom->nQtl--;
	  oldQtl->chrom->nQtl++;
      swapQtlData(oldQtl, thisqtl);  /* swap the records ... so keep the old data safe */
   
#ifndef MCMC_OMIT_DEBUG
	  checkIntegrity(nQtl, thisqtl->chrom);
	  if (oldQtl->chrom != thisqtl->chrom) checkIntegrity(nQtl, oldQtl->chrom);
#endif

	  return NULL;
	}
}



QTL_INFO* swap_add_dom(int nn, int nQtl, QTL_INFO** qtls,
					   PRIORS* priors, 		              
					   DATA* myData, WORK* myWork, MCMC_PARAM* myMCMC)
/*=======================================================================================
 *
 *
 *
 *
 *  Exchange with any QTL in genome:
 *		 long_range_update(nn, nqtl, thisqtl, all_qtls,
 *					   myData->chrom_pos, myData->totalChromLen, 
 *					   myData->nChrom, myData->myChroms,
 *		               priors, theparams, individs, myData, myWork, myMCMC);					  
 *
 *
 *  Update on one chromosome only:
 *		 long_range_update(nn, nqtl, thisqtl, all_qtls,
 *					   NULL, 0, 0, NULL,
 *		               priors, theparams, individs, myData, myWork, myMCMC);					  
 *
 *=========================================================================================
 */
		              
{
  FPN log_accept_update_prob, log_effect_ratio ;
  FPN uni_ran; 
  int lmark, i;

  /* things for the new QTL */
  CHROMOSOME* chrom;
  QTL_INFO* oldQtl;
  QTL_INFO* thisqtl;

  /* work structures we'll use */
  FPN* mod_effect = myWork->mod_effect;    /* new grand mean and effects */
  FPN** new_XtX = myWork->new_XtX;         /* new XtX matrix (copy of XtX, then adjust) */
  FPN** chol = myWork->chol;               /* new Cholesky decomposition */
  FPN* new_XtY = myWork->new_XtY;          /* new XtX matrix (copy of XtX, then adjust) */
  FPN* p = myWork->p;                      /* new diagonal for the Cholesky */
  FPN* pmean = myWork->pmean;              /* mean of grand mean/effects */
  FPN* pvar = myWork->pvar;                /* variance of grand mean/effects */
  FPN* u = myWork->u;                      /* iid normal values */
  FPN* weight = myWork->weight;            /* new weights */
  params* theparams = myData->theparams;
  individual* individs = myData->individs; 
  FPN log_position;

  /* some variables we need for the acceptance ratio */
  FPN bCb, d_invD_d, log_det_Chol, log_det_invD;
  int nterm = myData->na[ADD] + myData->na[DOM] + 1;
  int changeQtl;

  changeQtl = ignuin(1,nQtl);
  if (changeQtl != nQtl)
     moveQtlToEndOfXtX(nQtl, qtls, changeQtl, myWork->XtX, myWork->XtY, nterm);					 

  oldQtl = qtls[nQtl];
  chrom = oldQtl->chrom;
  removeQtl(oldQtl->qptr);  
 

  mycopy(nterm, nterm,new_XtX, myWork->XtX);
  for (i=1; i<= nterm; i++) new_XtY[i] = myWork->XtY[i];
                                  /* don't want to corrupt the XtX unless we make the move */

  /* create and initialize the new QTL */
  lmark = binSearch(chrom->nMark, chrom->mark_pos, 
	                oldQtl->qptr->pos);
  thisqtl = createQtl(nn, nQtl, &qtls[0], chrom, lmark, 
	               oldQtl->qptr->pos, myMCMC->addParam, 
				   (FPN*)NULL, (FPN*)NULL);

  /* copy the genotype */
  for (i=1; i<=nn; i++) thisqtl->qptr->genotype[i] = oldQtl->qptr->genotype[i];
	  	               
  log_position = get_log_position_ratio(myMCMC->revjump, thisqtl->chrom, myData) -
                 get_log_position_ratio(myMCMC->revjump, oldQtl->chrom, myData);


  /* now get the effect + params required for the acceptance ratio */
  log_effect_ratio = get_effect(nn, nQtl, myData->y, myData->mu,
	                            myData->sigmasq, oldQtl, thisqtl, 
                                qtls, myMCMC->revjump, priors, myData->na,mod_effect,
 			                    new_XtX, new_XtY,chol, p, pmean, pvar, u, weight,
			                    &bCb, &d_invD_d, &log_det_Chol, &log_det_invD, myWork);
							    
		
    log_accept_update_prob =  log_effect_ratio + log_position;
    uni_ran = genunf(0.0,1.0);
    uni_ran = log(uni_ran);

    if(uni_ran < log_accept_update_prob)
	{
	  Swap2DTable(&myWork->XtX, &myWork->new_XtX);
	  SwapDVec(&myWork->XtY, &myWork->new_XtY);
      setCholParams(myWork, bCb, d_invD_d, log_det_Chol, log_det_invD);

      swapQtlData(oldQtl,thisqtl);

	  /* important to have swapped data before setting effect, since now
	     the first nQtl elements of 'qtls' hold the real QTLs */
	  setEffect(nn, nQtl, myData->y, myWork->mod_effect, qtls, &myData->mu, myWork->weight, 
		        myWork->resid, myData->na);


#ifndef MCMC_OMIT_DEBUG
	  checkIntegrity(nQtl, oldQtl->chrom);
      checkResid(nn, nQtl, myData->mu, myData->y, 
		         qtls, myWork->resid, myData->gmiss);
#endif

	  return oldQtl;
    }
	else 
	{
      removeQtl(thisqtl->qptr);           /* return our genome data structures to normal */
      restoreQtl(oldQtl->qptr);            

#ifndef MCMC_OMIT_DEBUG
	  checkIntegrity(nQtl, oldQtl->chrom);
#endif
	  return NULL;
	}
}




void setCholParams(WORK* myWork, FPN bCb, FPN d_invD_d, 
			 FPN log_det_Chol, FPN log_det_invD)
{
      myWork->bCb = bCb; 
      myWork->d_invD_d = d_invD_d;
	  myWork->log_det_Chol = log_det_Chol; 
	  myWork->log_det_invD = log_det_invD;
}





void setWeights(int nQtl, int* na,
				  int revjump, QTL_INFO** all_qtls, FPN* w,
				  QTL_INFO* newQtl)
			 	  
/*
  input:
    nQtl    ... number of QTLs in model after jump or transition
	            (if newQtl is set, then nQtl-1 are located in 
				qtls, and newQtl holds the nQtl QTL)
    n_a     ... number of additive terms in model with nQtl QTLs
	n_d     ... number of dominance terms in model with nQtl QTLs
	all_qtls .. pointer to array containing pointers to QTL records
	            Note all_qtls[nQtl+1] points to proposed (birth) QTL

  output:
	w       ... new weights.  

*/
{
   QTL_INFO* thisqtl;
   int i,j;
   int weightType = (revjump & MCMC_SELECT_EFFECT_FIELD);
   

   


   switch (weightType)
   {
   case SELECT_VAR_FISCH: 
       for (i=2; i<= na[ADD]+na[DOM]+1; i++) w[i]=1.0;
	   break; 

   case SELECT_VAR_GAFFNEY:
   


	   for (i=1, j=2; i<=nQtl-1; i++)
	   {
		   thisqtl = all_qtls[i];
		   if (thisqtl->flag & QTL_ADD) w[j++] = 1.0/na[ADD];
		   if (thisqtl->flag & QTL_DOM) w[j++] = 1.0/na[DOM];
	   }	   
       if (!newQtl && nQtl >0) newQtl = all_qtls[nQtl];
	   if (newQtl)
	   {
	      if (newQtl->flag & QTL_ADD) w[j++] = 1.0/na[ADD];
	      if (newQtl->flag & QTL_DOM) w[j++] = 1.0/na[DOM];
	   }


	   if (na[ADD]+na[DOM]+2 != j)
		   printf("Error in setWeights: %d + %d != %d\n", na[ADD], na[DOM], j-2);

	   break; 

   default:
	   printf("weightType assumed to be SELECT_VAR_FISCH\n");
       for (i=2; i<= na[ADD]+na[DOM]+1; i++) w[i]=1.0;
   };
}


int getFischEffect(int nQtl, QTL_INFO** qtls, FPN mu, 
			       PRIORS* priors, FPN* mod_effect, FPN* w, FPN* var)
{
	int i, idx;

	mod_effect[1] = mu;
	for (i=1, idx=2; i<= nQtl; i++)
	{
		if (qtls[i]->flag & QTL_ADD) 
		{
			mod_effect[idx]= qtls[i]->a[ADD];
			w[idx] = qtls[i]->w[ADD];
			var[idx] = w[idx] * priors->var[ADD];
			idx++;
		}
		if (qtls[i]->flag & DOM) 
		{
			mod_effect[idx]= qtls[i]->a[DOM];
			w[idx] = qtls[i]->w[DOM];
			var[idx] = w[idx] * priors->var[DOM];
			idx++;
		}
	}
	return idx-1;
}


FPN proposeFischBirth (long int nn, int nQtl, 
				           FPN* y, FPN mu, FPN sigmasq,
				           QTL_INFO* newQtl, 
						   QTL_INFO** qtls, int revjump, 
				           PRIORS* priors, int* na,
				           FPN* mod_effect,
						   FPN* pmean, FPN* pvar, FPN* weight,
						   FPN* oldVar, FPN* oldWeight,
						   FPN* resid, FPN* newResid)
/*
    nQtl .. number of QTL in model we're jumping to
    na  ... number of effects of type in model with nQtl QTL
  */

{
	int i;
	FPN val, temp;
	int nterm = na[ADD]+na[DOM]+1 - newQtl->nParam;  
	         /* number of parameters in model we're jumping from */
	
    getFischEffect(nQtl-1, qtls, mu, priors, mod_effect, oldWeight, oldVar);	
    setPriorMeanVar(nQtl, revjump, qtls, NULL, priors, 
					 weight, na, pmean, pvar);

	/* calculate the prior ratio */
	for (i=2, val=0; i<= nterm; i++) 
	{
		temp = (mod_effect[i] - pmean[i]);
		val += 0.5 * temp*temp*(1/oldVar[i] - 1/pvar[i]);
		val += 0.5*log(oldWeight[i]/weight[i]);
	}
	
	/* generate the birth effects, if a birth step */
	if (newQtl->flag & QTL_ADD) 
	{
		nterm++;
	    mod_effect[nterm] = gennor(pmean[nterm], sqrt(pvar[nterm])); 
	}
	if (newQtl->flag & QTL_DOM) 
	{
		nterm++;
		mod_effect[nterm] = gennor(pmean[nterm], sqrt(pvar[nterm])); 
	}

	calcResid2(nn,nQtl,y,mod_effect,qtls,newResid);
	return val - (calcResidSS(nn,newResid) - calcResidSS(nn, resid)) / (2.0 * sigmasq);
}



FPN proposeFischDeath (long int nn, int nQtl, 
				           FPN* y, FPN mu, FPN sigmasq,
						   QTL_INFO** qtls, int revjump, 
				           PRIORS* priors, int* na,
				           FPN* mod_effect,
						   FPN* pmean, FPN* pvar, FPN* weight,
						   FPN* oldVar, FPN* oldWeight,
						   FPN* resid, FPN* newResid)
/*
    nQtl .. number of QTL in model we're jumping to
    na  ... number of effects of type in model with nQtl QTL
  */

{
	int i;
	FPN val, temp;
	int nterm = na[ADD]+na[DOM]+1;  
	         /* number of parameters in model we're jumping from */
	
    getFischEffect(nQtl+1, qtls, mu, priors, mod_effect, oldWeight, oldVar);	
    setPriorMeanVar(nQtl, revjump, qtls, NULL, priors, 
					 weight, na, pmean, pvar);

	/* calculate the prior ratio */
	for (i=2, val=0; i<= nterm; i++) 
	{
		temp = (mod_effect[i] - pmean[i]);
		val += 0.5 * temp*temp*(1/oldVar[i] - 1/pvar[i]);
		val += 0.5*log(oldWeight[i]/weight[i]);
	}
	
	calcResid2(nn,nQtl,y,mod_effect,qtls,newResid);
	return val - (calcResidSS(nn, newResid) - calcResidSS(nn, resid)) / (2.0 * sigmasq);
}



FPN calcResidSS(int nn, FPN* resid)
{
	FPN residSS = 0.0;
	int i;

	for (i=1; i<=nn; i++)
		residSS += resid[i]* resid[i];

	return residSS;
}


void checkCholesky(int nQtl,int revjump, QTL_INFO** qtls, PRIORS* priors, 
				   FPN* weight, int* na, FPN* pmean,FPN* pvar,
	               FPN** XtX, FPN* XtY, FPN sigmasq,
				   FPN** chol, FPN* p, FPN* u, FPN* mod_effect,					    
                   FPN* bCb, FPN* d_invD_d, FPN* log_det_Chol, 
				   FPN* log_det_invD, WORK* myWork)
{
	int getFischProposal = (revjump & SELECT_SAMPLE_FISCH);

	if (!getFischProposal)
	{
	   setPriorMeanVar(nQtl,revjump, qtls, NULL, priors, weight,na, pmean,pvar);
       generate_effects(na[ADD]+na[DOM]+1, XtX, XtY, pmean, pvar, sigmasq,
 	                    chol, p, u, mod_effect,					    
                        bCb, d_invD_d, log_det_Chol, log_det_invD, 0);
	   if (!EQUALS(*d_invD_d, myWork->d_invD_d) ||
           !EQUALS(*log_det_Chol,myWork->log_det_Chol) ||	                                       
           !EQUALS(*log_det_invD,myWork->log_det_invD) ||
   		   !EQUALS(*bCb,myWork->bCb))
	   {
		   printf("Cholesky not correct on entry to gen_effects ... \n"
		          "recall you must recompute Cholesky (i.e. call\n"
			      "setPriorMeanVar(),generate_effects() and setCholParams()\n"
			      " or update genotypes with update_lambda_qtl(), \n"
			      "after updating sigmasq\n");
	   }
	}
}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file birth.c
------------------------------------------------------------------ */


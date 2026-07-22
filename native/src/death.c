/* ------------------------------------------------------ XCutXCodeXSkip
     This file (death.c) is part of QTL Cartographer
         
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
  File:       death.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:

    Function DEATH performs the "death step" for
    removing an existing QTL.
****************************************************/

#include "revjump.h"
#include "ranlib.h"




int death(DATA* myData, PRIORS* priors, 
		   params* theparams, WORK* myWork, MCMC_PARAM* myMCMC)
{
  FPN log_accept_death_prob, log_effect_ratio, 
	     log_proposal, log_position;
  FPN uni_ran;
  int i;
  /* get elements */
  int nQtl = myData->nQtl;
  int nn = myData->nn;
  FPN mu = myData->mu;
  FPN sigmasq = myData->sigmasq;
  QTL_INFO** all_qtls = myData->myQtls;
  QTL_INFO* deathQtl;
  CHROMOSOME* chrom;
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
  int nterm = myData->na[ADD] + myData->na[DOM] + 1;
    int changeQtl;

  if ( theparams->nn > 0 )
    i=1;

  changeQtl = ignuin(1,nQtl);
  if (changeQtl != nQtl)
     moveQtlToEndOfXtX(nQtl, all_qtls, changeQtl, XtX, XtY, nterm);					 

  deathQtl = all_qtls[nQtl];   /* recall we've randomized order */
  chrom = deathQtl->chrom;

  log_position = get_log_position_ratio(myMCMC->revjump, deathQtl->chrom, myData);

  log_effect_ratio = get_effect(nn, nQtl, myData->y, mu, sigmasq,
	                            deathQtl, NULL, all_qtls, myMCMC->revjump,
			                    priors, myData->na, mod_effect,
  			                    XtX, XtY, chol, p, pmean, pvar, u, weight,
								&bCb, &d_invD_d, &log_det_Chol, &log_det_invD, myWork);


                     
  log_proposal = get_log_proposal_ratio(nQtl-1, myData->bp, myData->dp, myData->prior_ratio);

  log_accept_death_prob = -log_proposal - log_position + log_effect_ratio;


/* DECIDE WHETHER OR NOT TO ACCEPT NEW QTL */

/**************************************************
proceed with accept/reject death only if jacobian
is not 0. if jacobian is 0, reject death.
**************************************************/


    uni_ran = genunf(0.0,1.0);
    uni_ran = log(uni_ran);

    if(uni_ran < log_accept_death_prob) 
	{
      dropQtl(nn, &myData->nQtl, myData->y, all_qtls, mod_effect, weight, 
		      &myData->mu, myWork->resid, myData->na);

#ifndef MCMC_OMIT_DEBUG
	 checkIntegrity(myData->nQtl, deathQtl->chrom);
#endif

	 setCholParams(myWork, bCb, d_invD_d, log_det_Chol, log_det_invD);
	  return 1;
	}
	return 0;
}




void dropQtl(int nn, int* nQtl, FPN* y, QTL_INFO** all_qtls, 
	         FPN *mod_effect, FPN* w, FPN* mu, FPN* resid, int* na)
{
  /*  recall  mu = effect[0]  and modified_mu = modified_parms[0]
      also modified_parms[*nqtl+1] contains the effect of the new qtl */

  QTL_INFO* dropQtl = all_qtls[*nQtl];
   
  /* remove its record from genome list */
  removeQtl(dropQtl->qptr);
  if (*nQtl == 0)
	  printf ("cannot delete non-existant QTL\n");
  else
     (*nQtl)--;

 
  /* set values */
  setEffect(nn, *nQtl, y, mod_effect, all_qtls, mu,
	        w, resid, na);
  dropQtl->chrom->nQtl--;
}




/***************************************
      Function GET_DEATH_INTERVAL
    randomly draws a new QTL interval 
***************************************/

int get_death_interval( int nQtl )
{
  return (int)((genunf(0.0, 1.0) * nQtl) + 1);    /* recall qtls are numbered 1,2,... */
}





int removeQtlFromList(int nQtl, QTL_INFO* thisqtl, QTL_INFO** qtls)
{
	int i,j;

	for (i=1; i<= nQtl; i++)
	{
		if (qtls[i] == thisqtl)
		{
			for (j=i+1; j<=nQtl; j++) 
				qtls[j-1] = qtls[j];
			break;
		}
	}
	return i;
}






/* ------------------------------------------------------- XCutXCodeXSkip
             End of file death.c
------------------------------------------------------------------ */


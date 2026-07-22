/* ------------------------------------------------------ XCutXCodeXSkip
     This file (read_data.c) is part of QTL Cartographer
         
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

/**************************************************************
  File:         read_data.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:
  -------
    getdata of type VOID.This reads 
    trait, marker, distance and initial values from different 
    files. 



**************************************************************/

#include "revjump.h"
#include "ranlib.h"


/*************************************************
             FUNCTION GETDATA
*************************************************/

void getdata( MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, 
              CHROMOSOME* chromInfo,   FPN** new_lambda, int** new_chrom )
{

  /* open files to read and write */
  int end_of_file1;
  FILE *read_n;

  int i;
  int type;
  FPN burn_in, preburn_in;
  FPN prior_sd[3];



  /***************************************** 
               READ NTRAIT-NMARK DATA
  *****************************************/

  read_n = fopen("nval.dat", "r");
  end_of_file1 = fscanf( read_n, " %X %lf %lf %d %d %d", 
		  &myMCMC->revjump, &burn_in, &preburn_in, &myMCMC->niter, &myMCMC->nby, 
		  &myData->nQtl);


  /*------------------------------------------------------------------ 
   * determine how many parameters to add each time in birth/death step
   * if revjump field is set for RJ_MCMC, we only have birth/death 
   * of one parameter at a time.  So addParam = 1.  If gmiss is not -2,
   * then we have only 2 genotypes, so addParam = 1.  Otherwise, we
   * will favor adding an additive and dominance parameter (addParam=2)
   *------------------------------------------------------------------ 
   */	  
  myMCMC->burnIn = (burn_in >= 0 && burn_in < 1)? burn_in: -1;
  myMCMC->preBurnIn = (preburn_in >= 0 && preburn_in < 1)? preburn_in: -1;


  if (myData->gmiss == -2 && 
	  (myMCMC->revjump & MCMC_SELECT_DOM_FIELD) != SELECT_RJMCMC &&
	  ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_GET_APPROP))
	   myMCMC->addParam = QTL_ADD_DOM;
  else
  {
	  if ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_GET_DOM_ONLY)
	     myMCMC->addParam = QTL_DOM;
	  else if ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_RJMCMC && myData->gmiss == -2)
	     myMCMC->addParam = QTL_NONE;
	  else
	     myMCMC->addParam = QTL_ADD;
  }

  
/*************************************************************************
   the parameter revjump is binary. revjump=1 means do 
   reversible jump mcmc computations. if revjump=0, no reversible
   jump mcmc, do only regular mcmc as in satagopan et al (1996).
*************************************************************************/

 
  /* initialize our data structures */
  if (myData->nQtl > 0) 
    {  
      *new_lambda = dvector(1,myData->nQtl);
      *new_chrom = ivector(1,myData->nQtl);
      
      end_of_file1 = fscanf( read_n, "%d", &(*new_chrom)[1]); 
      if ((*new_chrom)[1] <=0)
	{
	  for( i=2; i<= myData->nQtl; i++ ) (*new_chrom)[i] = (*new_chrom)[1];
	}
      else
	for( i=2; i<= myData->nQtl; i++ )
	  end_of_file1 = fscanf( read_n, "%d", &(*new_chrom)[i] );
      
      end_of_file1 = fscanf( read_n, "%lf", &(*new_lambda)[1] );
      if ((*new_lambda)[1] <=0)
	for( i=2; i<= myData->nQtl; i++ ) (*new_lambda)[i] = (*new_lambda)[1];
      else
	{
	  (*new_lambda)[1] /= 100.0;
	  for( i=2; i<= myData->nQtl; i++ ) 
	    {
	      end_of_file1 = fscanf( read_n, "%lf", &(*new_lambda)[i] );
	      (*new_lambda)[i] /= 100.0;
	    }
	}
    }
  
  
  end_of_file1 = fscanf( read_n, "%lf %lf" , &myData->mu, &myData->sigmasq ); 
  end_of_file1 = fscanf( read_n, "%lf %lf %lf %lf", 
			 &priors->mean[MU], &prior_sd[MU],
			 &priors->sig_a1, &priors->sig_a2 );
  
  if (priors->sig_a2 < 0) priors->sig_a2 = -priors->sig_a2 * myData->y_var;
  
  end_of_file1 = fscanf( read_n, "%lf %lf %lf %lf", 
			 &priors->mean[ADD], &prior_sd[ADD], 
			 &priors->mean[DOM], &prior_sd[DOM] );
  
  
  /* qtl prior mean is required only if doing reversible jump mcmc
     computations as the following "if" statement indicates */
  
  if( (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD) >= 1) {
    end_of_file1 = fscanf( read_n, "%lf", &priors->qtl_mean );
  }
  else{
    priors->qtl_mean = myData->nQtl;
  }
  
  
  /* check for defaults (-1 in standard deviation field => 
     use the mean field (0<mean<1 typically) times the data
     variance to give the variance for the prior variance */
  if (myData->sigmasq<0) 
    {
      myData->sigmasq = myData->mu * myData->y_var;
      myData->mu = myData->ybar;
    }
  
  priors->sampleVar[MU] = 0;
  
  if (prior_sd[MU] < 0) 
    {
    priors->var[MU] = priors->mean[MU] * myData->y_var;
    priors->mean[MU] = myData->ybar;
  }
  else
    priors->var[MU] = prior_sd[MU]*prior_sd[MU];
  
  priors->log_var[MU] = log(priors->var[MU]);
  
  
  
  for (type = ADD; type <= DOM; type++)
    {
      priors->sampleVar[type] = 0;
      priors->alpha[type]= priors->beta[type]=0.0;
    
    if (LE(prior_sd[type],0)) 
      {
	if (!LE(priors->mean[type],0)) 
	  priors->var[type] = priors->mean[type] * myData->y_var;
	else
	  {
	    priors->sampleVar[type] = 1;
	    
	    if (EQUALS(prior_sd[type],0) || !LT(priors->mean[type], 0))
	      {
		priors->alpha[type] = 2;
		priors->beta[type] = 10;   /* this gives Beta(1,5) as a prior with mean 0.25 */
	      }
	    else 
	      {
		priors->alpha[type] = -priors->mean[type];
		priors->beta[type] = -prior_sd[type];
	      }
	    
	    priors->var[type] = genbet(priors->alpha[type], priors->beta[type]) * 2 * myData->y_var;	
	  }
	priors->mean[type] = 0;
      }
    else
      priors->var[type] = prior_sd[type]*prior_sd[type];
    
    priors->log_var[type] = log(priors->var[type]);
    }
  
  fclose(read_n);
}






/* ------------------------------------------------------- XCutXCodeXSkip
             End of file read_data.c
------------------------------------------------------------------ */


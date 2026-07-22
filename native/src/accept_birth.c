/* ------------------------------------------------------ XCutXCodeXSkip
     This file (accept_birth.c) is part of QTL Cartographer
         
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
  File:       accept_birth.c
  Written by: Patrick Gaffney
  Date:       November 11, 2000
  Version:    0.4

  Purpose:

    Function ACCEPT_BIRTH computes the acceptance
    probability for the birth of a new QTL.
****************************************************/

#include "revjump.h"
#include "ranlib.h"

void checkResid(int nn, int nQtl, FPN mu, FPN* y, QTL_INFO** allQtls, FPN* resid,int gmiss)
{
  FPN sum;
  int i,j;
  int* geno;
  QTL_INFO* thisqtl;
  if ( gmiss == 0 )
    i = 1;
  for(i=1;i<=nn;i++){
    sum = mu;
    for(j=1; j<=nQtl; j++) 
      {
	thisqtl = allQtls[j];
	geno = igenotype(thisqtl);
	sum += EFFECT_VALUE(geno[i], thisqtl->a);
      }
    sum = y[i] - sum;
    if (!EQUALS(resid[i], sum)) 
      {
	printf("==>%d .. %f %f",i,sum, resid[i]);
      }
  }
}



int select_move(int nQtl, FPN* bp, FPN* dp)
     /*---------------------------------------------------------------
 * Description
 *    Selects the move type (1=>birth, 2=>death, 3=>none)
 *
 * chrom           pointer to chromosome where death/birth is    
 *                 occuring                                      
 * qtl_prior_mean  prior mean for the number of qtl (for POISSON!)
 * cval            Multiplier for birth and death probabilities.
 *                 Must be chosen so that the product cval*
 *                 (birth_prob+death_prob)<1. Recommended cval=0.4
 *                 If cval is in error, the maximum value is chosen
 *                                                               
 *--------------------------------------------------------------*/								    
{
  FPN uni_ran;
  int move_type;
  


  uni_ran = genunf(0.0, 1.0);
  if( uni_ran < bp[nQtl])
    move_type = C_BIRTH;
  else if( uni_ran < (bp[nQtl] + dp[nQtl]) )
    move_type = C_DEATH;
  else
    move_type=C_UPDATE;


  return move_type;
}





FPN get_log_proposal_ratio(int nQtl, FPN* bp, FPN* dp, FPN* priorRatio)
{
  /* nQtl ... number of QTL in reduced model */
  return log( (dp[nQtl+1]/bp[nQtl]) * priorRatio[nQtl+1]);
}



FPN get_log_position_ratio(int revjump, CHROMOSOME* chrom, DATA* myData)
{
  if (revjump & MCMC_FLAG_RANDOM_CHROM)
    return log(chrom->chromLen * myData->nChrom / myData->totalChromLen);
  else
    return 0;   /* we sample the position from the prior */
}

void calcResid2(int nn, int nQtl, FPN* y, FPN* modified_parms,
		QTL_INFO** qtls, FPN* newResid)
{
  int i,j,idx;
  QTL_INFO* thisqtl;
  FPN a[QTL_ADD_DOM];
  int* geno;

  for(i=1;i<=nn;i++) newResid[i] = y[i] - modified_parms[1];

  for(j=1, idx=2; j<=nQtl; j++) 
    {
      thisqtl = qtls[j];
      geno = igenotype(thisqtl);

      /* find which terms are applicable for the QTL */
      a[QTL_ADD] = (thisqtl->flag & QTL_ADD)? modified_parms[idx++]: 0;
      a[QTL_DOM] = (thisqtl->flag & QTL_DOM)? modified_parms[idx++]: 0; 

      for (i=1;i<=nn;i++) 
	newResid[i] -= EFFECT_VALUE(geno[i], a);
    }
}








/* ------------------------------------------------------- XCutXCodeXSkip
             End of file accept_birth.c
------------------------------------------------------------------ */


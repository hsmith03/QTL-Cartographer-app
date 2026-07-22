/* ------------------------------------------------------ XCutXCodeXSkip
     This file (INITVALS.c) is part of QTL Cartographer
         
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
  File:       initgeno.c
  Written by: Patrick Gaffney
  Date:       November 11, 2000
  Version:    0.4

  Purpose:
  -------
    These functions find the index of missing markers and traits.

**************************************************************/

#include "revjump.h"
#include "ranlib.h"

/********************************************
         Function INITQTL_GENO

  initializes qtl genotypes given their 
  location and flanking marker genotypes
  for the given chromosome
********************************************/


void initQtl(int nn, int* offset, QTL_INFO* thisqtl, params* theparams, individual* individs, int gmiss, int revjump)
{
  /* here we initialize the QTL on the basis of its flanking markers 
     only (and don't use trait values as in propose_qtl_geno() in
	file update_qtl.c) 
*/
	int i;
	int* geno;
	mygenome* qptr = thisqtl->qptr;

    genProbs(nn, offset, theparams, individs, gmiss, thisqtl, revjump);   
    thisqtl->prevDist = qptr->prev->dist;
    thisqtl->nextDist = qptr->dist;
	
	
	geno = igenotype(thisqtl);

	for (i=1; i<=nn; i++)
	   GenGenotype(gmiss, thisqtl->log_prob[i], &geno[i]);
}




void GenGenotype (int gmiss, FPN* log_pr, int* geno)				 
     /* assumes pr indexed by -1,0 or 1 */
     /* gmiss is set only for crosses with only 2-genotype progeny, to
   indicate the missing genotype (-1,0 or 1) */
{
  FPN u;
  
  u = genunf(0,1);	
  
  if (log(u) < log_pr[1] && gmiss !=1) *geno = 1;
  else if (log(1-u) <= log_pr[-1] && gmiss != -1) *geno = -1;
  else *geno = 0;
}




int* igenotype(QTL_INFO* thisqtl)
/* the first element (zero index) of array contains the first genotype value */
{
	return thisqtl->qptr->genotype;
}


int** p_igenotype(QTL_INFO* thisqtl)
/* the first element (zero index) of array contains the first genotype value */
{
	return &thisqtl->qptr->genotype;
}


int* genotype1(QTL_INFO* thisqtl)
/* the first element (zero index) of array contains the first genotype value */
{
	return (&thisqtl->qptr->genotype)[1];
}



void dgenotype(QTL_INFO* thisqtl, int nn, FPN* dz)
/* the first element (zero index) of array dz contains the first genotype value */
{
	int i;
	int* geno = igenotype(thisqtl);

	for (i=0; i<nn; i++)
		dz[i] = (FPN)(geno[i]);
}



int qtlnum(QTL_INFO* thisqtl)
{
	return -thisqtl->qptr->markr;
}



FPN IMdist(QTL_INFO* thisqtl)
{
	CHROMOSOME* chrom = thisqtl->chrom;
	int mark = binSearch(chrom->nMark, chrom->mark_pos, thisqtl->qptr->pos);
	return chrom->mark_pos[mark+1] - chrom->mark_pos[mark];
}






void random_perm_QTL(int nQtl, QTL_INFO** qtls, QTL_INFO** buff, int* perm_num)
/* randomly permute the QTLs (for updating) */
{
	int i,j, val;

	for (i=1; i<=nQtl; i++) perm_num[i] = i;
    for (i=nQtl; i>=1; i--) 
      {
        j = ignuin(1,i);
	val = perm_num[j];

	SwapInt(&perm_num[i], &perm_num[j]);  /* now move entry containing val to end */
	buff[i] = qtls[val];
	buff[i]->qptr->markr = -i;            /* make sure marker numbers consistant */
	/* though that now doesn't mean much   */
      }
    for (i=1; i<=nQtl; i++) 
      qtls[i] = buff[i];                       /* restore our list */
}



void printX(int nn, int nQtl, QTL_INFO** qtls)
{
  FILE* f;
  QTL_INFO* thisqtl;
  int i,j;

  f = fopen("x.txt","w");
  for (i=1; i<=nn; i++)
    {
      fprintf(f,"1  ");
      for (j=1; j<=nQtl; j++)
	{
	  thisqtl = qtls[j];
	  if (thisqtl->flag & QTL_ADD) fprintf(f,"%d  ",thisqtl->qptr->genotype[i]);
	  if (thisqtl->flag & QTL_DOM) 
	    fprintf(f,"%5.2lf  ",(double) DOM_GENO_VALUE(thisqtl->qptr->genotype[i]));
	}
      fprintf(f,"\n");
    }
  fclose(f);
}

void printXtX(int nterm, FPN** XtX)
{
  int i,j;
  FILE* f;
	  
  f = fopen("x.xtx","w");
	  
  for (i=1; i<=nterm; i++)
    {
      for (j=1; j<=nterm; j++)
	fprintf(f,"%lf ",XtX[i][j]);
      fprintf(f,"\n");
    }
  fclose(f);
}



FPN gammln2(FPN xx)
{
  FPN x,y,tmp,ser;
  static FPN cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file INITVALS.c
------------------------------------------------------------------ */


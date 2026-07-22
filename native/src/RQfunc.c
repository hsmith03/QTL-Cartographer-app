/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RQfunc.c) is part of QTL Cartographer
         
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

#include "Main.h"
/* 

Functions to simulate a genetic model for a quantitative trait.


Printing routines are now in Qdatain.c.   This file is only required by Rqtl.

*/


/*
 * This calculates the allele effect of a QTL. job determines the
 * distribution. 
 * 
 
   if (beta < -10.0)
     job = 1;
   else if (beta <= 0.0)
     job = 2;
   else
     job = 3;

 
    job = 
       1 => effect = 0.25   
       2 => effect = 0.25 * N(0,1) 
       3 => effect = effect = gamma(alpha,beta)
 *
 * where alpha is scaled so that the distribution has a mean of 1.
 */

FPN qtl_effects(FPN beta,int job)
{
  FPN aa, bb, pp, effect;
  if (job == 1)
    effect = (FPN) 0.25;
  else if (job == 2)
    effect = (FPN) 0.25 * (FPN) fabs(gasdev(&job));
  else {
    if (beta > (FPN) 0.25)
      effect = (FPN) 0.25 * gamgbh(beta, job);
    else {
      aa = (FPN) 1.0 / beta;
      bb = (FPN) -1.0 / ((FPN) 1.0 - beta);
      pp = gammp(beta, (FPN) 1.0);
      effect = (FPN) 0.25 * gamnl1(beta, aa, bb, pp, job);
    }
  }
  return (effect);
}



/*
 * This calculates the dominance effect of a QTL. job determines the
 * distribution.
         1 => effect =  0.0  and there is no dominance
 job =   2 => effect =  1.0  Complete dominance of A over a
         3 => effect = -1.0  Complete dominance of a over A
         4 => effect is beta distributed with parameter b1 and b2.
 */

FPN qtl_dom_effects(FPN b1, FPN b2, int job)
{
  FPN aa, bb, pp, effect, dom1, dom2;
  if (job == 1)
    effect = (FPN) 0.0;
  else if (job == 2)
    effect = (FPN) 1.0;
  else if (job == 3)
    effect = (FPN) -1.0;
  else {
    if (b1 > (FPN) 0.25)
      dom1 = gamgbh(b1, 10);
    else {
      aa = (FPN) 1.0 / b1;
      bb = (FPN) -1.0 / ((FPN) 1.0 - b1);
      pp = gammp(b1, (FPN) 1.0);
      dom1 = gamnl1(b1, aa, bb, pp, 10);
    }
    if (b2 > 0.25)
      dom2 = gamgbh(b2, 10);
    else {
      aa = (FPN) 1.0 / b1;
      bb = (FPN) -1.0 / ((FPN) 1.0 - b2);
      pp = gammp(b2,(FPN)  1.0);
      dom1 = gamnl1(b2, aa, bb, pp, 10);
    }
    effect = (FPN) 2.0 * dom1 / (dom1 + dom2) - (FPN) 1.0;
  }
  return (effect);
}



/*
 * This places a QTL on the genome at a random place. first points to the
 * genome.  The random spot is uniform over the genome, with the caveat that
 * no QTL will be in the same interval as an existing QTL, nor in an interval
 * adjioning an interval that contains a QTL. QTLs can exist in the tails
 * outside the flanking markers of a chromosome.  Once the QTL is places, the
 * recombination frequencies between it an the flanking markers are
 * calculated, and a random allelic effect is calculated by qtl_effects.
 */

genome *create_a_qtl(genome *first,aqtl *qtlptr,params *theparams)
{
  genome *gptr;
  FPN position, total, tgenome;
  int qtljob,k;
  k=0;
 /* determine what kind of allelic effect the QTL will
 
 have
 
 */
  if (theparams->beta < -10.0)
    qtljob = 1;
  else if (theparams->beta <= 0.0)
    qtljob = 2;
  else
    qtljob = 3;
 /* Calculate total length of genome in Morgans */
  gptr = first;
  tgenome = (FPN) 0.0;
  while (gptr != NULL) {
    tgenome = tgenome + gptr->dist;
    gptr = gptr->next;
  }
  gptr = first;
  position = tgenome * ranf(1);
  total = (FPN) 0.0;
  while (total < position && gptr != NULL) {	/* This determines the marker following the QTL */
    total = total + gptr->dist;
    gptr = gptr->next;
  }
  if (gptr != NULL)	/* if qtpr is NULL, then the QTL is on the tail of the last chromosome */
    gptr = gptr->prev;
  else
    for (gptr = first; gptr->next != NULL; gptr = gptr->next) k+=1;
  qtlptr->a = qtl_effects(theparams->beta, qtljob);
  qtlptr->d = qtl_dom_effects(theparams->beta1, theparams->beta2, theparams->dom);
  qtlptr->chrm = gptr->chrom;
  qtlptr->mrk = gptr->markr;
  qtlptr->c2 = mapfunc((total - position), -1);
  qtlptr->c1 = mapfunc((gptr->dist - total + position), -1);
  if (qtlptr->map->brdrs > 0.0) {
    if (qtlptr->mrk == 0)	/* Are we on a left border? */
      qtlptr->c1 = (FPN) 0.5;

    if (qtlptr->mrk == qtlptr->map->mpc[qtlptr->chrm])  /* Are we on a right border? */
	  qtlptr->c2 = (FPN) 0.5;
  }

 /* Now eliminate the correct nodes.  We don't want QTLs in
    the same interval.
    
    
    If theparams->Rmode == 0, then we don't wnat QTL in
    the same interval as a pre-existing QTL or in an adjacent interval.
    We need to eliminate at least one, and possibly 3 nodes. 
    
    If theparams->Rmode == 1 then we don't wnat QTL in
    the same interval.
    
    We eliminate nodes only if they are on the same chromosome
    and the markers are consecutive. 
    
    
    
    */
  if (gptr->prev == NULL) {	/* if we are on the leftmost tail of the chromosomes, eliminate 2				
				intervals */
    if (theparams->Rmode == 0 && gptr->chrom == gptr->next->chrom && gptr->markr + 1 == gptr->next->markr)
      gptr = elim_gnode(gptr);  /*  Eliminate the node */
    gptr = elim_gnode(gptr);  

  }
  else if (gptr->next == NULL) {	/* if we are on the rightmost tail of the chromosomes,
					   eliminate 2 intervals */
    if (theparams->Rmode == 0 && gptr->chrom == gptr->prev->chrom && gptr->markr - 1 == gptr->prev->markr)
      gptr = elim_gnode(gptr);
    gptr = elim_gnode(gptr);
    gptr = first;
  }
  else {	/* else we eliminate up to 3 intervals, if they are all consecutive and on the same
		   chromosome */
    if (theparams->Rmode == 0 && gptr->chrom == gptr->prev->chrom && gptr->markr - 1 == gptr->prev->markr)
      gptr = elim_gnode(gptr->prev);

    if (theparams->Rmode == 0 && gptr->chrom == gptr->next->chrom && gptr->markr + 1 == gptr->next->markr)
      gptr = elim_gnode(gptr);

    gptr = elim_gnode(gptr);
    gptr = first;

  }
  return (gptr);	/* make sure to return a pointer to the beginning of the genome */
}




/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RQfunc.c
------------------------------------------------------------------ */


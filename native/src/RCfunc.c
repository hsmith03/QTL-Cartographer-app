/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RCfunc.c) is part of QTL Cartographer
         
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

   Subroutines for simulating   data sets 

   All output routines are now in Idatain.c.   This means that Rcross is the
   only module that requires this file.   
*/

/*
  Expand the map for RI lines.  theparams->crosst should be 0, 1, or 2.
  if 0, it means RI lines created random F2, followed by purification (Trudy MacKay) = Doubled Haploid
  if 1, it means RI lines created by selfing.
  if 2, it means RI lines created by sib mating.
*/
void expand_map(params *theparams,genome *gptr,aqtl *qtlptr,markermap *themap)
{
  int kk,ii;
  genome *tgptr;
  kk = 0;
  for (ii = 1; ii <= themap->traits; ii++)
    kk = kk + themap->knum[ii];
  for (ii = 1; ii <= kk; ii++) {
    qtlptr[ii].c1 = convert_dist( qtlptr[ii].c1, 1, theparams->crosst);
    qtlptr[ii].c2 = convert_dist( qtlptr[ii].c2, 1, theparams->crosst);    
  }  
  for ( tgptr = gptr ; tgptr != NULL ; tgptr = tgptr->next ) {    
    tgptr->dist = convert_dist(tgptr->dist, 0, theparams->crosst);        
    tgptr->pos  = convert_dist(tgptr->pos, 0, theparams->crosst);     
    tgptr->mxo  = convert_dist(tgptr->mxo, 0, theparams->crosst);    
    tgptr->pxo  = convert_dist(tgptr->pxo, 0, theparams->crosst);     
  } 
}

/*
  if flag = 0, dist is in morgans
  if flag = 1, dist is in rec
  ss indicates selfed RI (1) or sibmated (2).
  return units given
*/
FPN convert_dist(FPN dist,int flag,int ss)
{
  FPN rec;
    if ( flag == 0 ) /*convert to recombination freq.*/
      rec =    mapfunc(dist,-1);  
    else 
      rec = dist; 
    if (ss == 2 ) /*sib mating ri lines*/
      rec = (FPN)4.0*rec / ((FPN)1.0+(FPN)6.0*rec);
    else if ( ss == 1 ) /*selfed ri lines*/
      rec = (FPN)2.0*rec / ((FPN)1.0+(FPN)2.0*rec);
    else 
      rec = dist;
    if ( flag == 0 ) /*convert to morgans.*/
      rec = mapfunc(rec,1);

  return(rec);
}

/*  Return a pointer to the line with index whch */
thelines *which_line(int whch, thelines *lines)
{
  thelines *lptr;
  for (lptr = lines; lptr != NULL; lptr = lptr->next)
    if (lptr->which == whch)
      return (lptr);
  return (NULL);
}

/*  allocate space for a line */
thelines *alloc_aline(char *name,individual *iptr,int nn,struct aline *prev,struct aline *next,int which)
{
  thelines *ptr;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  ptr = (thelines *) malloc((size_t) sizeof(thelines));
#else
  ptr = (thelines *) malloc((unsigned) sizeof(thelines));
#endif
    if ( debugging > 2 ) {
        sprintf(gwarn,"In alloc_aline(), allocated 1 line  at %p\n",ptr);
        MemoryAccount(gwarn);
    }
  if (!ptr)
    nrerror("allocation failure in alloc_aline()");
  ptr->name = cvector(0, MAXNAME);
  strcpy(ptr->name, name);
  ptr->iptr = iptr;
  ptr->nn = nn;
  ptr->prev = prev;
  ptr->next = next;
  ptr->which = which;
  return (ptr);
}

/*  Deallocate space for a line */
void free_aline(thelines *ptr)
{
  if (ptr->name != NULL)
    free_cvector(ptr->name, 0, MAXNAME);
  if (ptr->iptr != NULL)
    free_indvector(ptr->iptr - 1, ptr->nn);
    if ( debugging > 2 ) {
        sprintf(gwarn,"In free_aline(), deallocated 1 line  at %p\n",ptr);
        MemoryAccount(gwarn);
    }
  if (ptr != NULL)
    free((char *) ptr);
}

/* Show the lines currently existing */
void show_lines(thelines *ptr)
{
  thelines *tptr;
  printf("\n\n\n\tThese are the lines we have so far\n\n");
  for (tptr = ptr; tptr != NULL; tptr = tptr->next)
    printf("\n\t%3d.  %s", tptr->which, tptr->name);
  printf("\n\n\n");
}



/*
   Do a cross.  pg1 is a pointer to parental population 1.
                pg2 is a pointer to parental population 2.
                off is the offspring population pointer.
                there are pg1n, pg2n and offn individuals in these three populations.
                pg1 and pg2 can be the same population.
                
                
                This now includes selfing.  cross ==3 implies a selfed F2, while
                cross==4 implies a random F2.  
                
                the s is a flag to indicate whether it is selfing or not.  test crosses
                to a selfed F line require this flag != 1.
*/
void do_a_cross(params *theparams,thelines *p1ptr,thelines *p2ptr,thelines *p3ptr,genome *gptr, int s)
{
  FPN malefemale;
  int p1i, p2i, ii;
  individual *pg1, *pg2, *off;
  int pg1n, pg2n, offn,start;
/*  thelines *tp1ptr;*/
    start = 0;
    pg1 = p1ptr->iptr;
    pg1n = p1ptr->nn;
    pg2 = p2ptr->iptr;
    pg2n = p2ptr->nn;
    off = p3ptr->iptr;
    offn = p3ptr->nn;
    if ( s == -1 ) 
      offn = p3ptr->nn/2;
    else if ( s == -2 )
      start =  p3ptr->nn/2 ;
    if (pg1n == 1)
      p1i = 0;
    if (pg2n == 1)
      p2i = 0;
    for (ii = start; ii < offn; ii++) {
      if ( theparams->cross == 6 )
        malefemale = (FPN)0.1;
      else
        malefemale = ranf(2);
      if (pg1n > 1)
        p1i = (int) ((FPN) pg1n * ranf(p1i));
      if ( s==1  ) 
        p2i = p1i;
      else if (pg2n > 1)
        p2i = (int) ((FPN) pg2n * ranf(p2i));
      recombination(gptr, theparams );
      if (malefemale < 0.5)
        create_an_individual(theparams, &pg1[p1i], &pg2[p2i], &off[ii], gptr);
      else
        create_an_individual(theparams,&pg2[p2i], &pg1[p1i], &off[ii], gptr);
      if ( s == -1 || s == -2 )
        off[ii].bc = -s;
    }

}

/*
 * Perform recombination for a pair of individuals.  If L is the length of
 * a chromosome, then there will be Poisson(L,vL) crossovers in that chromosome
 * in each parent.
 *
 *  genptr is just a linked list representing the genome.
 *  For each chromosome, we calculate xo1 and xo2, which are poisson distributed
 *  random variables determining the number of crossovers on that chromosome.
 *  xo1 is the number of the paternal gamete, while xo2 is that for the maternal
 *  gamete.  These crossovers are placed uniformly on the genome.
 */
void recombination(genome *genptr,params *theparams )
{
  FPN tgenome, *rec1, *rec2, howfar;
  int ii, xo1, xo2, r1, r2, chromosome;
  genome *gptr, *thischromptr, *nextchromptr;

  chromosome = 1;
  gptr = genptr;
  while (chromosome <= theparams->themap->m) {
/* Do crossing over for each chromosome in turn.  */
    thischromptr = gptr;
    tgenome = (FPN)0.0;
    while (gptr != NULL && gptr->chrom == chromosome) {
      gptr->mxo = (FPN) 0.0;
      gptr->pxo = (FPN) 0.0;
      tgenome = tgenome + gptr->dist;
      gptr = gptr->next;
    }
    nextchromptr = gptr;
    ii = 10;
    xo1 = (int) poidev(tgenome, &ii);	/* number of crossovers for the paternal gamete */
    xo2 = (int) poidev(tgenome, &ii);	/* number of crossovers for the maternal gamete */
    if (xo1 > 0) {
      rec1 = ran_arry(xo1);	/* array of x01 uniform random numbers.  Will determine */
      for (ii = 1; ii <= xo1; ii++)	/* the placement of the xo1 crossovers on this chromosome */
	    rec1[ii] = tgenome * rec1[ii];
      if (xo1 > 1)
	    sort(xo1, rec1);
      rec1[0] = tgenome * (FPN) 2.0;
    }
    if (xo2 > 0) {	/* same as above with maternal chromosome */
      rec2 = ran_arry(xo2);
      for (ii = 1; ii <= xo2; ii++)
	    rec2[ii] = tgenome * rec2[ii];
      if (xo2 > 1)
	    sort(xo2, rec2);
      rec2[0] = tgenome * (FPN)2.0;
    }

    r1 = r2 = 1;
    howfar = (FPN) 0.0;
    gptr = thischromptr;
    while (gptr != NULL && gptr->chrom == chromosome) {
/*
   Go along the chromosome.  Determine if a crossover occurred in some interval.
   Place it's position in M in the gptr->pxo field.  Position of crossovers for
   paternal and maternal gametes.
*/
      howfar = howfar + gptr->dist;
      if (xo1 > 0)
	    if (howfar > rec1[r1]) {
	      gptr->pxo = rec1[r1] - howfar + gptr->dist;
	      r1 = r1 + 1;
	      if ( r1 > xo1)
	        r1 = 0;
	    }
      if (xo2 > 0)
	    if (howfar > rec2[r2]) {
	      gptr->mxo =  rec2[r2] - howfar + gptr->dist;
	      r2 = r2 + 1;
	      if ( r2 > xo2)
	        r2 = 0;
	    }
      if (r1 == 0 && r2 == 0)
	    gptr = NULL;
      else
	    gptr = gptr->next;
    }


    if (xo1 > 0)
      free_dvector(rec1, 0, xo1);
    if (xo2 > 0)
      free_dvector(rec2, 0, xo2);
    chromosome = chromosome + 1;
    gptr = nextchromptr;
  }
}

/*
 * The following takes two individuals and 'mates' them to produce an
 * offspring.  Inputs are pointers to the parents and to the offspring, as
 * well as a genome.
 *
 *
 *
 *
 * The following are assumed: Value  Genotype  Explanation
 *
 * 3     AA       homozygote, Allele A from both parents
   2     Aa       heterozygote, with the A having come from the paternal parent
   1     aA       heterozygote, with the A having come from the maternal parent
   0     aa       homozygote, Allele a from both parents
 *
 * For a pair of loci, 22 or 11 mean coupling, while 21 or 12 mean repulsion
 * gametes.
 */
void create_an_individual(params *theparams,individual *fptr,individual *mptr,individual *optr,genome *genptr)
{
  int ii, jj, i, j, pgt, mgt, wqtl,onchrom;
  FPN pgamete, mgamete, tpgamete, tmgamete;
  genome *gptr;
  optr->map = fptr->map;
  optr->qtls = mptr->qtls;
  gptr = genptr;
  onchrom = 0;
  while (gptr != NULL) {
    ii = gptr->chrom;
    jj = gptr->markr;
/*    if ((fptr->map->brdrs > 0.0 && jj == 0) || (fptr->map->brdrs <= 0.0 && jj == 1)) {*/
    if ( ii != onchrom ) {   /* reset the gamete at each chromosome boundary */
      mgamete = ranf(10);
      pgamete = ranf(10);
      onchrom = ii;
    }	/* if gamete <= 0.5, then a 1 contributes an 'a' while a 2 contributes an 'A'. if gamete >
	   0.5, then a 1 contributes an 'A' while a 2 contributes an 'a'. gamete should be renewed
	   for each chromosome randomly, and should be switched for each XO. */
    mgt = mptr->markers[ii][jj];
    pgt = fptr->markers[ii][jj];
    if ( theparams->cross == 6 )
      optr->markers[ii][jj] = sc_genotype( mgt, pgt, mgamete, pgamete);
    else
      optr->markers[ii][jj] = genotype(theparams,mgt, pgt, mgamete, pgamete);
    if (gptr->whichqtl > 0) {	/* if there is at least one QTL in the interval... */
      wqtl = 0;
      for (i = 1; i <= optr->map->traits; i++)
	for (j = 1; j <= optr->map->knum[i]; j++) {
	  wqtl = wqtl + 1;
	  if ( optr->qtls[wqtl].chrm == gptr->chrom && optr->qtls[wqtl].mrk == gptr->markr) {
	    mgt = mptr->vqtls[i][j];
	    pgt = fptr->vqtls[i][j];
	    tpgamete = pgamete;
	    tmgamete = mgamete;
	    if (gptr->markr != 0) {
	      if (gptr->pxo > (FPN) 0.0 && mapfunc( fptr->qtls[wqtl].c1, 1) > gptr->pxo)
		    tpgamete = (FPN) 1.0 - pgamete;
	      if (gptr->mxo > (FPN) 0.0 && mapfunc( mptr->qtls[wqtl].c1, 1) > gptr->mxo)
		    tmgamete = (FPN) 1.0 - mgamete;
	    }
	    else {
	      if (gptr->pxo > (FPN) 0.0 && gptr->dist - mapfunc(fptr->qtls[wqtl].c2, 1) > gptr->pxo)
		    tpgamete = (FPN) 1.0 - pgamete;
	      if (gptr->mxo > (FPN) 0.0 && gptr->dist - mapfunc(mptr->qtls[wqtl].c2, 1) > gptr->mxo)
		    tmgamete = (FPN) 1.0 - mgamete;
	    }

	    optr->vqtls[i][j] = genotype(theparams, mgt, pgt, tmgamete, tpgamete);
	  }
	}
    }
    if (gptr->pxo > (FPN) 0.0)
      pgamete = (FPN) 1.0 - pgamete;
    if (gptr->mxo > (FPN) 0.0)
      mgamete = (FPN) 1.0 - mgamete;
    gptr = gptr->next;
  }
  optr->g = calc_genotype(optr);
}


/*
 *
 * The following are assumed: Value  Genotype  Explanation
 *
 *
   3     AA       homozygote, Allele A from both parents
   2     Aa       heterozygote, with the A having come from the paternal parent
   1     aA       heterozygote, with the A having come from the maternal parent
   0     aa       homozygote, Allele a from both parents
 *
 *
 * mgt is the maternal genotype, mgamete tells which gamete is used pgt is the
 * paternal genotype, pgamete tells which gamete is used
 */

int genotype(params *theparams,int mgt,int pgt,FPN mgamete,FPN pgamete)
{
  int ogt, mgam, pgam;
  switch (mgt) {
   case 3: mgam = 1; break;
   case 2:
    if (mgamete <= 0.5)  mgam = 1;
    else                 mgam = 0;
    break;
   case 1:
    if (mgamete <= 0.5)  mgam = 0;
    else                 mgam = 1;
    break;
   case 0:               mgam = 0; break;
   default:              mgam = -1; break;
  }
  switch (pgt) {
   case 3:               pgam = 1; break;
   case 2:
    if (pgamete <= 0.5)  pgam = 1;
    else                 pgam = 0;
    break;
   case 1:
    if (pgamete <= 0.5)  pgam = 0;
    else                 pgam = 1;
    break;
   case 0:               pgam = 0; break;
   default:              pgam = -1; break;
  }
  ogt = -1;
  if (pgam == 1 && mgam == 1)       ogt = 3;
  else if (pgam == 1 && mgam == 0)  ogt = 2;
  else if (pgam == 0 && mgam == 1)  ogt = 1;
  else if (pgam == 0 && mgam == 0)  ogt = 0;
  if ( theparams->cross == 5 ) {
    if ( pgam == 1 )                ogt = 3;
    else                            ogt = 0;
  }
  return (ogt);
}




/*
 * This produces the initial Parental lines and the F1 generation
 *
 * P1 is the high line.  All of its molecular markers have a value of 3, which
 * means AA.  QTLs have the same alleles.
 *
 * P2 is the low line.  All of its modlecular markers have a value of 0, which
 * means aa.  QTLs have the same alleles.
 *
 * F1 is the F1 hybrid line.  All of its molecular markers have a value of 1 or
 * 2, which means Aa.  QTLs have the same alleles. There will be two F1
 * types, those with all 1's and those with all 2's. Those with all 1's had
 * fathers from P1, while those with all 2's had fathers from P2
 *
 *  For Sue Carson:  Change the initial states of the P1, P2 and F1 populations
 *  so that the QTL alleles are configured randomly...P1, P2 & F1 will have a combination
 *  of 0, 1, 2, 3 at the QTL locus.
 Ultimately, there will be up to four alleles at each marker locus, A1, A2, A3 and A4.
 Abbreviate A1.A1 as 11, A1.A2 as 12, etc.  Order is important internally, so 
 12 and 21 are different.  There are 16 genotypes at each marker:
 
Genotype code:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
Genotype:      11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44

The first number indicates the paternally derivied allele, and the second the maternally
derived allele. 


Crosses:

  SC1:  Assume that the paternal line has genotypes 12 or 21 at all marker loci
  (experimental design), and that the maternal line has genotype 22 at all loci.
  Offspring will be 12 or 22 depending upon whether the paternal gamete was 1 or 2.
  1 is dominant to 2.  The phase of the markers are random over loci.  All QTL loci
  are QQ, Qq, qQ or qq, chosen randomly. (Codes 3, 2, 1, 0, respectively).  The ouput
  will be phenotypes with the paternal allele (1 or 2).
 
  SC2:  Assume that the maternal line has genotypes 34 or 43 at all marker loci
  (experimental design), and that the paternal line has genotype 33 at all loci.
  Offspring will be 33 or 34 depending upon whether the maternal gamete was 3 or 4.
  4 is dominant to 3.  The phase of the markers are random over loci.  All QTL loci
  are QQ, Qq, qQ or qq, chosen randomly. (Codes 3, 2, 1, 0, respectively).  The ouput
  will be phenotypes with the maternal allele (3 or 4 recoded to 0 or 1, respectively).
 
 
 */
void init_pop(individual *p1p2f1, params *theparams)
{
  int ii, jj, kk, mrkrs;
  FPN dgamete;
  if ( theparams->cross == 6 ) {  /* Do Sue Carson's lines   */
  
	  for (ii = 1; ii <= (p1p2f1 + 1)->map->traits; ii++)
	    for (kk = 1; kk <= *((p1p2f1 + 1)->map->knum + ii); kk++) { /*QTLs are 0 1 2 or 3*/
	      p1p2f1[1].vqtls[ii][kk] = (int) ( 4.0 * ranf(ii));
	      p1p2f1[2].vqtls[ii][kk] = (int) ( 4.0 * ranf(ii));
	      p1p2f1[3].vqtls[ii][kk] = (int) ( 4.0 * ranf(ii));
	      p1p2f1[4].vqtls[ii][kk] = (int) ( 4.0 * ranf(ii));
	    }
	  for (ii = 1; ii <= p1p2f1[1].map->m; ii++) {
	    mrkrs = p1p2f1[1].map->mpc[ii];
	    for (jj = 1; jj <= mrkrs; jj++) {
	      dgamete = ranf(ii);
	      if (theparams->crosst == 1 && dgamete <= 0.5 ) {
	        p1p2f1[1].markers[ii][jj] = 2; /* SC1: 50% 12 x 22 */
	        p1p2f1[2].markers[ii][jj] = 6;
	      }
	      else if (theparams->crosst == 1 && dgamete > 0.5 ) {
	        p1p2f1[1].markers[ii][jj] = 5; /* SC1: 50%  21 x 22 */
	        p1p2f1[2].markers[ii][jj] = 6;
	      }
	      else if (theparams->crosst == 2 && dgamete <= 0.5 ) {
	        p1p2f1[1].markers[ii][jj] = 12; /* SC2: 50% 33 x 34 */
	        p1p2f1[2].markers[ii][jj] = 11;
	      }
	      else if (theparams->crosst == 2 && dgamete > 0.5 ) {
	        p1p2f1[2].markers[ii][jj] = 15; /* SC2: 50% 33 x 43 */
	        p1p2f1[1].markers[ii][jj] = 11;
	      }
	    }
	    p1p2f1[3].markers[ii][jj] = p1p2f1[4].markers[ii][jj] =  1;
	  }
	  for ( ii=1; ii<=4; ii++)
	    p1p2f1[ii].g = calc_genotype(&p1p2f1[ii]);
  }
  else {  /* do the default inbred lines. */
	  for (ii = 1; ii <=  p1p2f1[1].map->traits; ii++)
	    for (kk = 1; kk <=  p1p2f1[1].map->knum[ii]; kk++) {
	      p1p2f1[1].vqtls[ii][kk] = 3;
	      p1p2f1[2].vqtls[ii][kk] = 0;
	      p1p2f1[3].vqtls[ii][kk] = 1;
	      p1p2f1[4].vqtls[ii][kk] = 2;
	    }
	  for (ii = 1; ii <= p1p2f1[1].map->m; ii++) {
	      mrkrs = p1p2f1[1].map->mpc[ii];
	    for (jj = 1; jj <= mrkrs; jj++) {
	      p1p2f1[1].markers[ii][jj] = 3;
	      p1p2f1[2].markers[ii][jj] = 0;
	      p1p2f1[3].markers[ii][jj] = 1;
	      p1p2f1[4].markers[ii][jj] = 2;
	    }
	  }
	  for ( ii=1; ii<=4; ii++)
	    p1p2f1[ii].g = calc_genotype(&p1p2f1[ii]);
  }
}


/*
 * The genonotypes are just the sum of the additive effects over all QTLs,
 *
  Genotype    AA    Aa    aa
  Value        a    ad    -a 
  
  This routine used in version 1.14 and earlier.  
 */
FPN *calc_genotype2(individual *indptr)
{
  FPN ka;
  int kk, ii, wqtl;
  if (indptr->g == NULL)
    indptr->g = dvector(1, indptr->map->traits);

  wqtl = 0;
  for (ii = 1; ii <= indptr->map->traits; ii++) {
    ka = (FPN) 0.0;
    for (kk = 1; kk <= indptr->map->knum[ii]; kk++) {
      wqtl = wqtl + 1;
      if ( indptr->vqtls[ii][kk] == 3)	/* 3 => AA, so add  a */
	    ka = ka + indptr->qtls[wqtl].a;
      else if (indptr->vqtls[ii][kk] == 1 || indptr->vqtls[ii][kk] == 2)
	    ka = ka + indptr->qtls[wqtl].a * indptr->qtls[wqtl].d ;	/* 1 or 2 => Aa, so add a*d */
      else if (indptr->vqtls[ii][ kk] == 0 )
	    ka = ka - indptr->qtls[wqtl].a  ;	/* 0 => aa, so subtract a */
    }
    indptr->g[ii] = ka;
  }
  return (indptr->g);
}

/*
   Fore each trait, sum the additive effects and multiply by two.  
   
   This is the expected difference in the phenotype between the two parental lines.  

*/
void    calc_parental_diffs(markermap *themap,aqtl *theqtls) {
  int wo,kk,ii; 
  
  themap->ParentalDiff = dvector(1,themap->traits);
  
    wo = 0;
    for (kk = 1; kk <= themap->traits; kk++) {
      themap->ParentalDiff[kk] = (FPN) 0.0;
      for (ii = 1; ii <= themap->knum[kk] ; ii++) {
	    wo = wo + 1;
        themap->ParentalDiff[kk] = themap->ParentalDiff[kk] + theqtls[wo].a ;
      }
      themap->ParentalDiff[kk] = (FPN) 2.0 * themap->ParentalDiff[kk] ;
  
   }
  
}
/*

  The genotypes  follow the model
 
 G = mu + sum_i a_i x_i   + sum_i d_i z_i  
        + sum_{i=1,k-1} sum_{j=i+1,k} ( bAA_ij x_i x_j + bAD_ij x_i z_j + bDA_ij z_i x_j  bDD_ij z_i z_j )
        
       QQ   Qq   qq
      ---------------
 x_i   1     0   -1
 z_i  -1/2  1/2 -1/2
 
a_i is additive effect at locus i
d_i is dominance effect at locus i
bAA_ij is additive by additive effect between loci i and j 
bAD_ij is additive by dominance effect between loci i and j 
bDA_ij is dominance by additive effect between loci i and j 
bDD_ij is dominance by dominance effect between loci i and j 


The mean mu is zero...should it be an option?

This used in version 1.15 and later.
*/
FPN *calc_genotype(individual *indptr)
{
  FPN ka,xi,xj,zi,zj;
  int i,j, ii,  qi,base;
  if (indptr->g == NULL)
    indptr->g = dvector(1, indptr->map->traits);
  base = 0;
  for (ii = 1; ii <= indptr->map->traits; ii++) {
    ka = (FPN) 0.0;
    for (i = 1; i <= indptr->map->knum[ii]; i++) {
      qi = base+i;
      assign_indicators(indptr->vqtls[ii][i], &xi, &zi);
      ka = ka + xi * indptr->qtls[qi].a + zi * indptr->qtls[qi].d ;
      for (j=i+1; j<=indptr->map->knum[ii]; j++ ) {
        assign_indicators(indptr->vqtls[ii][j], &xj, &zj);
        ka = ka + xi * xj *  indptr->qtls[qi].episAADD[j][i] ;
        ka = ka + xi * zj *  indptr->qtls[qi].episADDA[j][i] ;  
        ka = ka + zi * xj *  indptr->qtls[qi].episADDA[i][j] ;
        ka = ka + zi * zj *  indptr->qtls[qi].episAADD[i][j] ;  
      }
    }
    indptr->g[ii] = ka;
    base = base + indptr->map->knum[ii];
  }
  return (indptr->g);
}


/*
Indicator variable for locus with genotype gt

gt     3   2,1    0
       QQ   Qq   qq
      ---------------
 x_i   1     0   -1
 z_i  -1/2  1/2 -1/2
*/

void assign_indicators(int gt, FPN *xi, FPN *zi) {
  if ( gt == 3) {	/* 3 => QQ   */
    *xi = (FPN) 1.0;  
    *zi = (FPN) -0.5;
  }
  else if (gt == 1 || gt == 2) { /* 1,2 => Qq or qQ  */
    *xi = (FPN) 0.0;  
    *zi = (FPN) 0.5;
  }
  else if (gt == 0 ) { /* 3 => qq   */
    *xi = (FPN) -1.0;  
    *zi = (FPN) -0.5;
  }
}


/*
 * The phenotypes are just the sum of the additive effects over all QTLs,
 * plus some normal environmental variance, with mean zero and variance
 * determined by the genetic variance and heritability.
 */
void calc_phenotypes(individual *indptr,int nn,params *theparams)
{
  FPN genvar, envar;
  int ii, tt;
  envar = (FPN) 1.0;
  for (tt = 1; tt <= indptr[1].map->traits; tt++) {
    if ( theparams->Environ > 0.0) 
      envar = theparams->Environ;
    else if ( theparams->Herit > 0.0) {
      genvar = calc_genvar(indptr, nn, tt);
      envar = genvar * ((FPN) 1.0 / theparams->Herit - (FPN) 1.0);
    }
	if ( envar > 0.0 ) 
	  envar = (FPN) sqrt(envar);
    for (ii = 1; ii <= nn; ii++)
      indptr[ii].y[tt] = indptr[ii].g[tt] + envar * gasdev(&ii);
  }
}

/*
 * This calculates the genetic sample variance for the QTLs
 */
FPN calc_genvar(individual *indptr,int nn, int trait)
{
  FPN sum, ssum, genvar;
  int ii;
  sum = ssum = (FPN) 0.0;
  for (ii = 1; ii <= nn; ii++) {
    sum = sum + indptr[ii].g[trait];
    ssum = ssum + indptr[ii].g[trait] * indptr[ii].g[trait];
  }
  genvar = (ssum - sum * sum / (FPN) nn) / (FPN) (nn - 1);
  return (genvar);
}




/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RCfunc.c
------------------------------------------------------------------ */


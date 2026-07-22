/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MIfunc.c) is part of QTL Cartographer
         
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


/* functions to implement multiple interval mapping. */


#include "Main.h"



/*  
put phenotypes in lnpk->y.   
Assume that individuals with missing phenotypes have been  moved to the end of the array, 
after theparams->nn.
*/
void copy_traits(params *theparams,markermap *themap,individual *individs,linpakws *lnpk)
{
  int ii;
  ii = themap->traits;
  for ( ii=1; ii<=theparams->nn; ii++ ) 
    lnpk->y[theparams->whichtrait][ii]  = individs[ii].y[theparams->whichtrait] ;
  lnpk->samplesize[theparams->whichtrait] = theparams->nn;
}


/*
This subroutine converts an integer code into a multilocus
genotype.    If there are m QTL in the model, then there
are  2^m multilocus QTL genotypes for crosses with two genotypic
classes (BC and RI lines), or 3^m multilocus QTL genotypes for 
Fx lines.  Note that this means that we can handle only 31 QTL
for BC and RI lines, and only 19 for Fx lines.   We would need to go
with long FPNs to accomadate more QTL.  

You input the following:

genotype = ivector(1,loci);  will hold the multilocus genotype
ngt      = 2 or 3:  a 2 for BC or RI, and a 3 for Fx's
j          is an integer code for the mulilocus genotype
loci       is the number of loci

Each position in the genotype vector will hold a 1, 0 or -1 indicating
the genotype at that locus.   These codes will stand for the
following:

-------------------------
           Cross
      -------------------
Code  BC1   BC2  RI    Fx
-------------------------
 1    QQ    Qq   QQ    QQ
 0    Qa    qq   qq    Qa
-1                     qq
-------------------------
ngt    2     2    2     3
-------------------------
-------------------------

This can be used to enumerate all the possibilities in 
a branching tree diagram in which each node has either two
or three branches.

There are limits to the number of loci that can be coded
and decoded.


Machine        BC/RI    Fx
32bit          31       19
64bit          63       39


*/
void WhichMultilocusGT(int loci, long j, int *genotype, int ngt) {
  long total,row;
  int l;
  row=j;
  total = 1;
  if ( ngt == 2 ) {
	  for ( l=1; l<loci; l++ )
	    total *= 2;                /* One-half the number of multilocus genotypes */
	  for ( l=1; l<=loci; l++ ) {
	    if ( row <= total ) 
	      genotype[l] = 0;
	    else {
	      row = row - total;
	      genotype[l] = 1;   
	    }
	    total = total/2;
	  }
  }
  else if ( ngt == 3 ) {
	  for ( l=1; l<loci; l++ )
	    total *= 3;                /* One-third the number of multilocus genotypes */
	  for ( l=1; l<=loci; l++ ) {
	    if ( row <= total ) 
	      genotype[l] = -1;
        else if ( row <= 2*total ) {
          row = row - total;
          genotype[l] = 0;       
        }
	    else {
	      row = row - 2*total;
	      genotype[l] = 1;   
	    }
	    total = total/3;
	  } 
  }
}

/*
   This takes a vector representing a multilocus genotype and returns
   the integer code for it.   
*/
long LongMultilocusGT(int loci, int *genotype, int ngt) {
  long j,factor;
  int i;

	if ( ngt == 2 ) {  
	  factor = 1;
	  j=1;
	  for ( i=loci; i>0; i-- ) {
	    if ( genotype[i] == 1 )
	      j = j+factor;
	    factor = factor*ngt; 
	  }
	}
	else if (ngt == 3) {  
	  factor = 1;
	  j=1;
	  for ( i=loci; i>0; i-- ) {
	    if ( genotype[i] == 1 ) 
	      j = j+2*factor;
	    else if ( genotype[i] == 0 ) 
	      j = j+factor;
	    factor = factor*ngt; 
	  }
	}  
  return(j);
}


/*
  Indicate whether a QTL for this trait resides in the interval 
 */
void mimplace_qtls(params *theparams,mimparam *mimptr,genome *first)
{
 
  genome *gptr;
  mimparam *lmimptr;
  int ii;
  ii = theparams->traits;
  for ( gptr = first; gptr != NULL; gptr = gptr->next )
    gptr->whichqtl = 0;
  for ( lmimptr=mimptr->next; lmimptr != NULL ; lmimptr = lmimptr->next )   {
      gptr = first;
      while (gptr != NULL)
        if ( gptr->chrom == lmimptr->qptr1->chrm && gptr->markr == lmimptr->qptr1->mrk ) {
          gptr->whichqtl = lmimptr->qtl1;
          gptr = NULL;
        }
        else
          gptr = gptr->next;          
  }  
}


/*  Write a header for the mim output file */
void write_mimheader(char *outfile,params *theparams,char *onamae, char *chptr,markermap *themap,int oc)
{
  int ii;
  ii = themap->traits;
  if (oc == 1 )
    print_head(onamae,outfile,chptr,1,100,theparams);  
}


/* 
    This is a driver for actually doing the MI mapping.  
    
*/
void do_mimanalysis(params *theparams,markermap *themap,aqtl *theqtls,individual *individs,genome *agptr,linpakws *lnpk)
{
  FPN initlogln,genvar;
  mimparam *mimptr;
  int k,jj, nparams,m,maxparams, savemaxqtl;
  FILE *outf,*moutf;
  k=0;
  outf = fileopen(theparams->mimfile, "a");
  if (outf == NULL)
      outf = fileopen(theparams->mimfile, "w");
  if ( theparams->verbosity == 1 ) {
        putline(stdout,'-',79);
	    printf("\n\n\tThis is an analysis for trait number %d named %s.\n",theparams->whichtrait,themap->tnames[theparams->whichtrait]);
        putline(stdout,'-',79);
        printf("\n");
  }
/*    Should we limit the number of parameters?  Jansen reccomends  nparams <= 2 sqrt(n)   */  
  savemaxqtl = theparams->maxqtl;
  maxparams =   (int)  floor(  2.0 * sqrt((FPN) lnpk->samplesize[theparams->whichtrait]) ) ;
  if ( theparams->maxqtl*(theparams->ngt-1) > maxparams ) {
    savemaxqtl = theparams->maxqtl;
    theparams->maxqtl = maxparams / ( theparams->ngt-1) ;
    sprintf(gwarn,"\n\n WARNING: Readusting maximum QTL downward from %d to %d \n",savemaxqtl, theparams->maxqtl);         
    LogTheError(theparams->error, gwarn);
    if (theparams->verbosity == 1 ) 
      printf("\n\n WARNING: Readusting maximum QTL downward from %d to %d \n",savemaxqtl, theparams->maxqtl);         
  }


/* Convert theqtls into a model */
  mimptr = ConstructMIMparams(theparams,themap,theqtls);
  nparams = HowManyQTL(mimptr,0);
  if ( theparams->mimwork[1] == 'M' ) {
	  if ( theparams->verbosity == 1 && theparams->mimwork[1] == 'M' ) 
	    printf("\n\n 1. Loading initial model for analysis... ");
	  if ( nparams > 0 ) 
	    PrintFinalModel(outf,theparams,themap,mimptr,"The Initial Model is",lnpk->samplesize[theparams->whichtrait]);
	  if ( theparams->verbosity == 1 ) 
	    printf(" DONE.");
  }
  else if ( theparams->verbosity == 1 )
    printf("\n\n 1. Skipped loading initial model for analysis... ");

/*  Rrefine the estimates of the given model */
  if ( theparams->mimwork[2] == 'P' ) {
	  if ( theparams->verbosity == 1 ) 
	    printf("\n 2. Refining estimates of all parameters in this model... ");
	  if ( nparams > 0  ) {
	    initlogln = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	    PrintFinalModel(outf, theparams,themap,mimptr,"After PARAMETER refinement,",lnpk->samplesize[theparams->whichtrait]);
	  }
	  if ( theparams->verbosity == 1  ) 
	    printf(" DONE.");
  }
  else if ( theparams->verbosity == 1 )
    printf("\n\n 2. Skipped refining estimates of parameters... ");

/*  Next, refine the estimates for the position of the QTLs in this model. */
  if ( theparams->mimwork[3] == 'R' ) {
	  if ( theparams->verbosity == 1    )  
	    printf("\n 3. Refining QTL position estimates within current intervals... ");
	  if ( nparams > 0  ) {
	    JitterPositions(theparams,themap, individs,mimptr,agptr,lnpk );
	    initlogln = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	    PrintFinalModel(outf, theparams,themap,mimptr,"After within-interval POSITION refinement, ",lnpk->samplesize[theparams->whichtrait]);
	  }
	  if ( theparams->verbosity == 1   ) 
	    printf(" DONE."); 
  }
  else if ( theparams->mimwork[3] == 'A' ) {
	  if ( theparams->verbosity == 1    )  
	    printf("\n 3. Refining QTL position estimates by testing adjacent intervals... ");
	  if ( nparams > 0  ) {
        TestAdjacentIntervals(theparams,themap,individs,mimptr,agptr,lnpk) ;
	    PrintFinalModel(outf, theparams,themap,mimptr,"After adjacent-interval POSITION refinement, ",lnpk->samplesize[theparams->whichtrait]);
	  }
	  if ( theparams->verbosity == 1   ) 
	    printf(" DONE."); 
  }
  else if ( theparams->verbosity == 1 )
    printf("\n\n 3. Skipped refining QTL position estimates... ");

/*  Recheck the significance of each QTL. */
  if ( theparams->mimwork[4] == 'T' ||   theparams->mimwork[4] == 'D' ||   theparams->mimwork[4] == 'E' ) {
	  if ( theparams->verbosity == 1 ) 
	    printf("\n 4. Checking each parameter for significance... ");
	  ConfirmCurrentModel(outf,theparams,themap,individs,mimptr,agptr,lnpk);
	  nparams = HowManyQTL(mimptr,0);
	  if ( nparams > 0 ) {
	    initlogln = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	    PrintFinalModel(outf, theparams,themap,mimptr,"After CONFIRMATION of parameters,",lnpk->samplesize[theparams->whichtrait]);  
	  }
	  if ( theparams->verbosity == 1 ) 
	    printf(" DONE.");
  }
  else if ( theparams->verbosity == 1 )
    printf("\n\n 4. Skipped checking each parameter for significance... ");
  
/*  Search for more QTL.  */
  if ( theparams->mimwork[5] == 'S' || theparams->mimwork[5] == 'A' || theparams->mimwork[5] == 'D' ) {
	  if ( theparams->verbosity == 1 ) 
	    printf("\n 5. Searching for more QTL... ");
	  FindMoreQTL(outf, theparams,themap,individs,mimptr,agptr,lnpk,nparams );
	  initlogln = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	  if ( theparams->mimwork[0] != 'S' )
	    PrintFinalModel(outf, theparams,themap,mimptr,"After SEARCHING for more QTL,",lnpk->samplesize[theparams->whichtrait]);  
	  if ( theparams->verbosity == 1 ) 
	    printf(" DONE.");
  }
  else if ( theparams->verbosity == 1 )
    printf("\n\n 5. Skipped searching for more QTL... ");

/* Search for epistatic terms*/
  if ( theparams->mimwork[6] == 'E' || theparams->mimwork[6] == 'B'  || theparams->mimwork[6] == 'U' ) {
	  if ( theparams->verbosity == 1 ) 
	    printf("\n 6. Searching for Epistatic interactions via %c...",theparams->mimwork[6]);
	  m = HowManyQTL(mimptr,1);  jj=0;
	  if ( m > 1 ) {
	    nparams = HowManyQTL(mimptr,0);
	    if ( theparams->mimwork[6] == 'B' || theparams->mimwork[6] == 'U' ) 
          jj = BackElimEpistasis( theparams,themap,individs,mimptr,agptr,lnpk,nparams,initlogln);
        
	    if ( theparams->verbosity == 1 && jj == -2 ) 
	      printf("\nToo many epistatic terms to do backward elimination...trying forward search.  " );
	    if ( theparams->mimwork[6] == 'E' || jj == -2 ) 
	      SearchForEpistasis(theparams,themap,individs,mimptr,agptr,lnpk,nparams,initlogln);

	    initlogln = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	    PrintFinalModel(outf, theparams,themap,mimptr,"After SEARCHING for epistatic terms, ",lnpk->samplesize[theparams->whichtrait]);  
	  }  
	  else if ( theparams->verbosity == 1 ) 
	      printf(" NO interactions to find!");  

	  if ( theparams->verbosity == 1 ) 
	    printf(" DONE.");
  }
  else if ( theparams->verbosity == 1 )
    printf("\n\n 6. Skipped searching for Epistatic interactions... ");

/* Calculate the breeding values and the Variance-Covariance matrix */
  if ( theparams->mimwork[7] == 'C'  ) {
	  if ( theparams->verbosity == 1 ) 
	    printf("\n 7. Calculating Variance-Covariance matrix and Breeding Values...");
      UpdateGTindex(mimptr);
	  m = HowManyQTL(mimptr,0);
	  QTLpostions(theparams,themap,mimptr,agptr);
	  if ( m > 0  ) {
	    CalcGTvalues(theparams,mimptr,lnpk);
	    ShowGTvalues(outf,theparams,lnpk,individs);
	    genvar = CalcVarCovar(theparams,mimptr,lnpk);
	    ShowVarCovar(outf,theparams,themap,mimptr,lnpk,genvar);  
        moutf = fileopen(theparams->tfile, "a");  
	    PrintRqtlout(moutf,theparams,themap,mimptr,2);
        fileclose(theparams->tfile, moutf);
	  }
	  else {
	    if ( theparams->verbosity == 1 ) 
	      printf(" NO QTLs left!");  
	  }
	  if ( theparams->verbosity == 1 ) 
	    printf(" DONE.");
  }
  else if ( theparams->mimwork[7] == 'R'  ) {
	  if ( theparams->verbosity == 1 ) 
	    printf("\n 7. Calculating residuals for the current model...   ");
      moutf = fileopen(theparams->tfile, "a");  
	  PrintRqtlout(moutf,theparams,themap,mimptr,2);
      fileclose(theparams->tfile, moutf);
      CalculateResiduals(theparams,themap,individs,mimptr,agptr,lnpk);
	  if ( theparams->verbosity == 1 ) 
	    printf(" DONE.");
  }
  else if ( theparams->verbosity == 1 )
    printf("\n\n 7. Skipped calculating Variance-Covariance matrix/Breeding Values or residuals... ");

/* Say goodbye and clean up. */    
  if ( theparams->verbosity == 1 ) {
      printf("\n");
      putline(stdout,'-',79);
      printf("\n\n\tMultiple Interval Mapping for trait %d named %s is complete.\n",theparams->whichtrait,themap->tnames[theparams->whichtrait] );
      putline(stdout,'-',79);
      printf("\n");
  }
  fileclose(theparams->mimfile, outf);
  mimptr->abd = 10; /*  delete_mimparam won't delete a mimptr if abd == 0 (an anchor), so we unprotect it.*/
  while ( (mimptr = delete_mimparam(mimptr)) != NULL ) k+=1;
  theparams->maxqtl = savemaxqtl;
}

/* Find terminal node of mim chain */
mimparam *TerminalMIMnode(mimparam *mimptr) {
  mimparam *tmimptr;
  int k;
  k=0;
  for ( tmimptr=mimptr; tmimptr->next != NULL ; tmimptr = tmimptr->next) k+=1;
  return(tmimptr);
}

/*
    Determine the position of each QTL from the left telomere in cM
*/
void QTLpostions(params *theparams,markermap *themap,mimparam *mimptr,genome *agptr) {
  mimparam *tmimptr;
  genome *tgptr;
  int i;
  i=themap->traits;
  i=theparams->traits;
  for ( tmimptr=mimptr->next; tmimptr != NULL ; tmimptr = tmimptr->next ) 
    if ( tmimptr->abd == 1 ) 
      for ( tgptr=agptr; tgptr != NULL ; tgptr = tgptr->next ) 
        if ( tgptr->chrom == tmimptr->qptr1->chrm && tgptr->markr == tmimptr->qptr1->mrk ) 
          tmimptr->qptr1->s = (FPN) 100.0 * tgptr->pos + mapfunc(tmimptr->qptr1->c1 ,2);	  
}	  
	  

/*
  Search all intervals for new QTL.
*/
void FindMoreQTL(FILE *outf, params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk,int nparams ) {
  int go_on,whqtl;
  time_t itptr,ctptr;
  mimplace_qtls(theparams,mimptr,agptr);
  go_on = nparams;
  go_on = 1;
  whqtl = HowManyQTL(mimptr,1);
  if ( theparams->verbosity == 1 ) 
    time(&itptr);
  if ( theparams->mimwork[5] == 'S' || theparams->mimwork[5] == 'A' )
    while ( go_on == 1 )   {
      if ( theparams->verbosity == 1 )
        printf("\n\tAt %2d QTL: ",whqtl);
      go_on = ScanGenomeQTL(outf,theparams,themap,individs,mimptr,agptr,lnpk);
      if ( theparams->verbosity == 1 ) {
        time(&ctptr);
        itptr = ctptr - itptr;
        printf("(%ld sec)",(long) itptr);
        itptr = ctptr;
        whqtl +=1;
      }
    }
  else if ( theparams->mimwork[5] == 'D' )
    while ( go_on == 1 )   
      go_on = SearchForDominance(theparams,themap,individs,mimptr,agptr,lnpk);
}

/*
  Search all intervals for new QTL.
  
  Do this based on theparams->mimwork[5]
  
    if S, then search for QTL as additive-dominance pairs
    if A, then search for additve effects
*/
int ScanGenomeQTL(FILE *outf,params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) {
  mimparam *amimptr,*dmimptr,*tmimptr;
  genome *gptr,*tgptr;
  int whqtl,lparams,foundone;
  FPN localic,d1,dist,pos,bestic,bestd,newic;

  whqtl = HowManyQTL(mimptr,1);
  if ( whqtl >= theparams->maxqtl ) {
    if (theparams->verbosity == 1 ) 
          printf(" Maxed Out on QTL search ");
    return(0);
  }
  
    foundone = 0;
  bestic = localic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;

  if ( theparams->mimwork[0] == 'S' ) {
    lparams = HowManyQTL(mimptr,0);
    if ( themap->tnames != NULL )
      fprintf(outf, "\n-trait       %5d      Analyzed trait [%s]", theparams->whichtrait,themap->tnames[theparams->whichtrait]);
    else
      fprintf(outf, "\n-trait       %5d      Trait analyzed", theparams->whichtrait);
    fprintf(outf, "\n-maxqtl      %5d    Maximum number of allowed QTL in the model",theparams->maxqtl);
    fprintf(outf, "\n-maxepis     %5d    Maximum number of epistatic terms allowed",theparams->maxepistatics);
    fprintf(outf, "\n-xic         %5d    Code for the IC criterion",theparams->whoseic);
    fprintf(outf, "\n-walk        %5.1f    Walking speed for position refinement and QTL search, in cM",theparams->walk);
    fprintf(outf, "\n-LRthresh    %5.1f    Likelihood ratio threshold for adding/deleting a QTL",theparams->mimlod);
    fprintf(outf, "\n-workcode    %s       Code indicating what to do",theparams->mimwork);
    fprintf(outf, "\n-params      %5d    Number of parameters in initial model",lparams);
    fprintf(outf, "\n-modelfile   %s  \n#   current LR %f",theparams->mqtfile,localic);
    fprintf(outf,"\nchrom  mark  position         IC          LR\n-s");
  }
  

  tmimptr = TerminalMIMnode(mimptr);
  whqtl = HowManyQTL(mimptr,1);
  whqtl = -1 * (whqtl + 1);
  amimptr = insert_mimparam(tmimptr,1, (FPN) 0.0,  whqtl, 0) ;
  amimptr->qptr1->map = themap;
  amimptr->qtl1 = -1 * whqtl ; 

  if ( theparams->mimwork[5] != 'A' ) {
	  if ( theparams->crosstype == 3 || theparams->crosstype == 4 ) {
	    dmimptr = insert_mimparam(amimptr,2, (FPN) 0.0,  whqtl, 0) ;
	    dmimptr->qtl1 = -1 * whqtl ; 
	    dmimptr->qptr1 = amimptr->qptr1;
	  }
  }
  bestd = (FPN) 0.0;
  
  for ( gptr=agptr ; gptr  != NULL ; gptr  = gptr->next ) {  /*Scan genome for best new QTL*/
    if ( gptr->whichqtl == 0 ) {
	  amimptr->qptr1->chrm = gptr->chrom;
	  amimptr->qptr1->mrk = gptr->markr;
	  amimptr->qptr1->a = (FPN) 0.0;
	  dist = 100 * gptr->dist;
	  d1 = (FPN) MIN_DIST ;
	  for ( pos = (FPN) MIN_DIST ; pos < dist; pos = pos + theparams->walk ) {
	    amimptr->qptr1->c1 = mapfunc( pos, -2 );
	    amimptr->qptr1->c2 = mapfunc( dist-pos, -2);
	    newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
        if ( theparams->mimwork[0] == 'S' )         
          fprintf(outf,"\n%5d %4d %10.5f %12.2f  %10.4f",gptr->chrom,gptr->markr, 100*gptr->pos + pos , newic, localic - newic);
	    if ( newic < bestic ) {
	      tgptr = gptr;
	      bestd = pos;
	      bestic = newic; 
	    }
	  } 
	}     
  }
  
  if ( (localic-bestic) > theparams->mimlod ) {/* We found another QTL */
	  amimptr->qptr1->chrm = tgptr->chrom;
	  amimptr->qptr1->mrk = tgptr->markr;
	  amimptr->qptr1->a = (FPN) 0.0;
	  amimptr->qptr1->c1 = mapfunc( bestd, -2 );
	  amimptr->qptr1->c2 = mapfunc( ((FPN) 100.0 * tgptr->dist - bestd), -2);
	  amimptr->qptr1->r2 = lrtolod(localic-bestic);
	  newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
	  tgptr->whichqtl = amimptr->qtl1;
	  if ( theparams->verbosity == 1 ) 
	    printf("+");
	  foundone = 1;
  }
  else { /* No new QTL */
    if ( theparams->mimwork[5] != 'A' ) {
      if (theparams->crosstype == 3 || theparams->crosstype == 4 ) 
        tmimptr = delete_mimparam(dmimptr);
    }
    tmimptr = delete_mimparam(amimptr);
  }
  if ( theparams->mimwork[0] == 'S' ) {
    foundone = 0;  /*  Just do one pass if S in position 0 */
    fprintf(outf,"\n-e\n");
  }
  whqtl = HowManyQTL(mimptr,1);
  return(foundone);
}

/*
  Check all QTL for dominance effects
  
*/
int SearchForDominance(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) {
  mimparam *amimptr,*dmimptr,*tmimptr;
  int foundone,m,d,i;
  FPN localic,bestic,newic;
  
  foundone = 0;
  bestic = localic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
  tmimptr = mimptr;
  m = HowManyQTL(mimptr,1);
  d = 0;
  for ( i=1; i<=m ; i++ ) {  /* how many QTL don't yet have dominance effects? */
    dmimptr = WhichMIMnode(mimptr,i,0,2,1);
    if ( dmimptr->abd == 0 && (theparams->crosstype == 3 || theparams->crosstype == 4) )
      d = d+1;
  }
  if ( d > 0 ) {	/* if more than 0, then test them. */	  
	  dmimptr = insert_mimparam(mimptr,2, (FPN) 0.0, d, 0) ;
	  for ( i=1; i<=m ; i++ ) {
	    amimptr = WhichMIMnode(mimptr,i,0,1,1);
	    dmimptr = WhichMIMnode(mimptr,i,0,2,1);
	    if ( dmimptr->abd == 0 && (theparams->crosstype == 3 || theparams->crosstype == 4) ) {
		  dmimptr->qtl1 = i ; 
		  dmimptr->qptr1 = amimptr->qptr1;
	      newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
          amimptr->lnlike = localic - newic;
          if ( newic < bestic ) {
            bestic = newic;
            tmimptr = amimptr;
          }
		}
		else
		  amimptr->lnlike = (FPN) 0.0;		
	  }  
  }
  else /* if they all have such effects, we are done. */
    return(foundone);
  
  if ( tmimptr->abd != 0 && tmimptr->lnlike > theparams->mimlod ) {/* We found a new dominance effect */
    dmimptr->qtl1 = tmimptr->qtl1;
    dmimptr->qptr1 = tmimptr->qptr1;
	dmimptr->qptr1->tr2 = lrtolod(tmimptr->lnlike);
    foundone = 1;
	newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
	if ( theparams->verbosity == 1 ) 
	    printf("d");
  }  
  else  /* No new dom. effect */
    tmimptr = delete_mimparam(dmimptr);
  return(foundone);
}
/* 
  Given a model, evaluate it and return the IC. 
*/
FPN EvalModel(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk  ) {
    int i,jj,lparams;
    FPN newlogln,ysum,ysum2,en,normfact,mean,variance;
    en = (FPN) 1.0/(FPN) lnpk->samplesize[theparams->whichtrait] ;
    lparams = HowManyQTL(mimptr,0);

    if ( lparams > 0 ) {
      UpdateGTindex(mimptr);
      CalculatePriors(theparams,themap,individs,mimptr,agptr,lnpk);
      CalcMeanVar(theparams,mimptr,lnpk); 
      jj =  RefineEstimates( theparams,mimptr,lnpk);
      lparams = HowManyQTL(mimptr,0);
      newlogln = CalcIC(lnpk->s2[theparams->whichtrait][0] , lparams, lnpk->samplesize[theparams->whichtrait], theparams->whoseic);
    }
    else {  /*  For a null model (no parameters aside from the mean)*/
      ysum = ysum2 = (FPN) 0.0;
      for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ ) {
        ysum = ysum + lnpk->y[theparams->whichtrait][i];
        ysum2 = ysum2 + (FPN) pow(lnpk->y[theparams->whichtrait][i] , 2.0);
      }
      mean = ysum * en ;
      variance = (ysum2 -  ysum * ysum * en  ) * en ;
      normfact = (FPN) 1.0 /(FPN) sqrt( 2.0 * (FPN) PI * variance ) ;

      newlogln = (FPN) 0.0;
      for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ ) 
	    newlogln =  newlogln + (FPN) log( normfact *  exp( -0.5 * pow(lnpk->y[theparams->whichtrait][i] - mean,2.0)/variance) );	 
      newlogln = CalcIC(newlogln , lparams, lnpk->samplesize[theparams->whichtrait], theparams->whoseic);
    }
    return(newlogln);
}

/*  
  Search for epistatic terms 
  Need to cut out when number of epistatics exceeds theparams->maxepistatics
*/
void SearchForEpistasis(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk,int nparams,FPN initlogln) {
  mimparam *rmimptr,*smimptr,*tmimptr,*lmimptr,*ismimptr;
  int lparams,nepis,maxepis,maxparams,nqtl,i,j,gt1,gt2,abd;
  FPN localic,minic,previc,**episAADD,**episADDA;
  localic = initlogln;
  i=nparams;
  nqtl = HowManyQTL(mimptr,1);
  episAADD = dmatrix(1,nqtl,1,nqtl);
  episADDA = dmatrix(1,nqtl,1,nqtl);
  
  maxparams =   (int)  floor(  2.0 * sqrt((FPN) lnpk->samplesize[theparams->whichtrait]) ) ; /* 2 sqrt(n) is maximum allowed parameters */
  lparams = HowManyQTL(mimptr,0);
  maxepis = maxparams - lparams;   /* 2 sqrt(n) - lparams is what is left over for fitting epistatic interactions. */
  if ( theparams->maxepistatics <= maxepis )  /*  Downratchet maxepis if it is bigger than requested */
    maxepis = theparams->maxepistatics;
  else { /*  Indicate whether we had use fewer than requested. */
    sprintf(gwarn,"\n\n WARNING: Readusting maximum epistatic interactions from  %d to %d \n",theparams->maxepistatics, maxepis);         
    LogTheError(theparams->error, gwarn);
    if (theparams->verbosity == 1 ) 
      printf("\n\n WARNING: Readusting maximum epistatic interactions from  %d to %d \n",theparams->maxepistatics, maxepis);     
  }    
  
  
  nepis = HowManyQTL(mimptr, -3);
  while  ( nepis < maxepis ) {
    minic = previc = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
    for ( i=1; i<=nqtl; i++ )
      for ( j=1; j<=nqtl; j++ )
        episAADD[i][j] = episADDA[i][j] = (FPN) 0.0;
    tmimptr = TerminalMIMnode(mimptr);
    lmimptr = insert_mimparam(tmimptr, 0,(FPN) 0.0,0,0 ) ;   
	for ( rmimptr = tmimptr ;  rmimptr->abd != 0 ; rmimptr = rmimptr->prev ) 	
	  for ( smimptr = rmimptr->prev ; smimptr->abd != 0 ; smimptr = smimptr->prev ) 
	    if ( rmimptr->abd == 1 && smimptr->abd == 1 ) {  /*  We have a pair of QTL, check all epistatic interactions that aren't part of the model. */
	      lmimptr->qtl1 = rmimptr->qtl1;
	      lmimptr->qtl2 = smimptr->qtl1;
          lmimptr->qptr1 = rmimptr->qptr1;
          lmimptr->qptr2 = smimptr->qptr1;
          lmimptr->qptr1->nptrs = lmimptr->qptr1->nptrs + 1;
          lmimptr->qptr2->nptrs = lmimptr->qptr2->nptrs + 1;  
          ismimptr = WhichMIMnode(mimptr,rmimptr->qtl1,smimptr->qtl1,3,2);
          if ( ismimptr->abd == 0 ) 
              episAADD[rmimptr->gt1][smimptr->gt1] = CheckEpistasis(theparams,themap,individs,mimptr,rmimptr,smimptr,lmimptr,agptr,lnpk,3);             
		  if (  theparams->crosstype == 3 || theparams->crosstype == 4 ) { /* Search for ad, da and dd effects if three genotypes */
             ismimptr = WhichMIMnode(mimptr,rmimptr->qtl1,smimptr->qtl1,6,2);
             if ( ismimptr->abd == 0 )
                episAADD[smimptr->gt1][rmimptr->gt1] = CheckEpistasis(theparams,themap,individs,mimptr,rmimptr,smimptr,lmimptr,agptr,lnpk,6);
             ismimptr = WhichMIMnode(mimptr,rmimptr->qtl1,smimptr->qtl1,6,2);
             if ( ismimptr->abd == 0 )
                episADDA[smimptr->gt1][rmimptr->gt1] = CheckEpistasis(theparams,themap,individs,mimptr,rmimptr,smimptr,lmimptr,agptr,lnpk,5);
             ismimptr = WhichMIMnode(mimptr,rmimptr->qtl1,smimptr->qtl1,4,2);
             if ( ismimptr->abd == 0 )
                episADDA[rmimptr->gt1][smimptr->gt1] = CheckEpistasis(theparams,themap,individs,mimptr,rmimptr,smimptr,lmimptr,agptr,lnpk,4);
		   }
		   lmimptr->qptr1->nptrs = lmimptr->qptr1->nptrs - 1;
           lmimptr->qptr2->nptrs = lmimptr->qptr2->nptrs - 1; 
	     }      
/*  Now, look at all IC's that are not 0.0... find the largest.   If */
    gt1 = gt2 = abd = 0;
    for ( i=1; i<nqtl; i++ )
      for ( j=i+1; j<=nqtl; j++ ) {
          if ( episAADD[j][i] > (FPN) 0.0 && episAADD[j][i] < minic ) {
            minic = episAADD[j][i];
            gt1 = i;
            gt2 = j;
            abd = 3;
          }
          if ( episAADD[i][j] > (FPN) 0.0 && episAADD[i][j] < minic  ) {
            minic = episAADD[i][j];
            gt1 = i;
            gt2 = j;
            abd = 6;
          }
          if ( episADDA[i][j] > (FPN) 0.0 && episADDA[i][j] < minic  ) {
            minic = episADDA[i][j];
            gt1 = i;
            gt2 = j;
            abd = 4;
          }
          if ( episADDA[j][i] > (FPN) 0.0 && episADDA[j][i] < minic ) {
            minic = episADDA[j][i];
            gt1 = i;
            gt2 = j;
            abd = 5;
          }        
      }
	  if ( previc - minic  > theparams->mimlod ) {
	    nepis = nepis + 1;
        rmimptr = WhichMIMnode(mimptr,gt1,0,1,-1);
        smimptr = WhichMIMnode(mimptr,gt2,0,1,-1);
	    lmimptr->qtl1 = rmimptr->qtl1;
	    lmimptr->qtl2 = smimptr->qtl1;
        lmimptr->qptr1 = rmimptr->qptr1;
        lmimptr->qptr2 = smimptr->qptr1;
        lmimptr->qptr1->nptrs = lmimptr->qptr1->nptrs + 1;
        lmimptr->qptr2->nptrs = lmimptr->qptr2->nptrs + 1;  
        lmimptr->abd = abd;	 
        lmimptr->lnlike = previc - minic;
        if (theparams->verbosity == 1 ) 
          printf("+");
	  }
	  else {  /*  Didn't find another epistatic interaction. */
	    lmimptr->qptr1 = NULL;
	    lmimptr->qptr2 = NULL;
	    lmimptr = delete_mimparam(lmimptr);
	    nepis = maxepis + 1;
      }
  } 
  free_dmatrix(episAADD,1,nqtl,1,nqtl);
  free_dmatrix(episADDA,1,nqtl,1,nqtl);

}

/*
   Check the epistatic term
*/
FPN CheckEpistasis(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,mimparam *rmimptr, mimparam *smimptr,mimparam *lmimptr,genome *agptr,linpakws *lnpk,int epitype) {
  mimparam *ismimptr;
  FPN newic;
  lmimptr->abd = 0;
  ismimptr = WhichMIMnode(mimptr,rmimptr->gt1,smimptr->gt1,epitype,-2);
  if ( ismimptr->abd == 0 ) {
    lmimptr->abd = epitype;
    newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
  }
  else
    newic = (FPN) 0.0;
  return(newic);
}


/*  
  Search for epistatic terms 
  Need to cut out when number of epistatics exceeds theparams->maxepistatics
*/
int BackElimEpistasis(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk,int nparams,FPN initlogln) {
  mimparam *rmimptr,*smimptr,*tmimptr,*lmimptr,*minmimptr;
  
  int lparams,nepis,maxepis,maxparams,nqtl,go_on,i;
  FPN localic,minlr,previc,newic,lrratio;
  i=nparams;
  localic = initlogln;
  nepis = HowManyQTL(mimptr, -3);
  if ( nepis > 0 )  /* Assume that there are no epistatic interactions yet.  */
    return(-1);
    
  nqtl = HowManyQTL(mimptr,1);
  if ( theparams->ngt == 2 )  /* How many could there be? */
     nepis = nqtl * ( nqtl - 1 ) / 2 ;
  else 
     nepis = nqtl * ( nqtl - 1 ) * 2 ;
  if ( theparams->mimwork[6] == 'B' )
    maxparams =   (int)  floor(  2.0 * sqrt((FPN) lnpk->samplesize[theparams->whichtrait]) ) ; /* 2 sqrt(n) is maximum allowed parameters */
  else
    maxparams =  lnpk->samplesize[theparams->whichtrait] - 1;
  lparams = HowManyQTL(mimptr,0);
  maxepis = maxparams - lparams;   /* 2 sqrt(n) - lparams is what is left over for fitting epistatic interactions. */
  
  if ( nepis > maxepis )  /*  There are too many for a backward elimination*/
    return(-2);

/*
   First, we put in a node for all possible epistatic interactions.   We will only
   get here if there aren't too many.   
*/
  lmimptr = TerminalMIMnode(mimptr);

  UpdateGTindex(mimptr);
  for ( rmimptr = lmimptr ;  rmimptr->abd != 0 ; rmimptr = rmimptr->prev ) 	
	for ( smimptr = rmimptr->prev ; smimptr->abd != 0 ; smimptr = smimptr->prev ) 
	  if ( rmimptr->abd == 1 && smimptr->abd == 1 ) {   /* Put in a node for each interaction term. */
        tmimptr = TerminalMIMnode(mimptr);
        insert_epistatic(tmimptr, rmimptr, smimptr, 3 );
        if ( theparams->ngt ==	3 ) {  /* For 3 genotypic classes, put in AD, DA and DD terms */
          insert_epistatic(tmimptr, rmimptr, smimptr, 4 );
          insert_epistatic(tmimptr, rmimptr, smimptr, 5 );
          insert_epistatic(tmimptr, rmimptr, smimptr, 6 );
        }  	  
	  }
/*
    Now we need to do the backward elimination.   
    
    1. Check the IC for the full model
    2. For each epistatic term, pull it, calculate the IC, and put it back
    3. Find the minimum IC.   If below the threshold, delete it permanently, and go to 1.
                              Else, quit loop
           
*/	  
  go_on = nepis;
  while  ( go_on >0 ) {
  /*  This is step 1. */
    UpdateGTindex(mimptr);
    previc = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
    minlr = (FPN) BIG;
    minmimptr = NULL;	  
  /*  This is step 2. */
    tmimptr = TerminalMIMnode(mimptr);
    for ( rmimptr = tmimptr ;  rmimptr->abd != 0 ; rmimptr = rmimptr->prev ) 
	  if ( rmimptr->abd > 2 ) {
	      PullMIMnode(rmimptr);
	      newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	      lrratio = newic - previc ;
	      PushMIMnode(rmimptr);
	      rmimptr->lnlike = lrratio;
	      if ( lrratio < minlr ) {
	        minlr = lrratio;
	        minmimptr = rmimptr;
	      }
          UpdateGTindex(mimptr);                
      }
  /*  This is step 3. */
    if ( minmimptr != NULL ) {
      if ( minlr < theparams->mimlod ) {
        go_on = go_on - 1;
        minmimptr = delete_mimparam(minmimptr);
        if (theparams->verbosity == 1 ) 
          printf("-");        
      }
      else
        go_on = 0;
    }
  } 
  return(0);
}

/*
Insert an epistatic term.   Uses insert_mimparam, but also sets the pointers to the
additive effects and increments their counters.   

*/
void insert_epistatic(mimparam *tmimptr, mimparam *rmimptr, mimparam *smimptr, int abdtype ) {
	mimparam *lmimptr;
	
    lmimptr = insert_mimparam(tmimptr, abdtype,(FPN) 0.0,rmimptr->qtl1, smimptr->qtl1 ) ;   
    lmimptr->qptr1 = rmimptr->qptr1;
    lmimptr->qptr2 = smimptr->qptr1;
    lmimptr->qptr1->nptrs = lmimptr->qptr1->nptrs + 1;
    lmimptr->qptr2->nptrs = lmimptr->qptr2->nptrs + 1;  
}

/*
  Calculate the IC
  
  We are currently using IC = -2(log(L(k) - k c(n)/2 )
    where L(k) is the likelihood with a k parameter model and sample size n
    
    c(n) is a penalty function
  We may modify this for new options on comparing models.  
  
*/
FPN CalcIC( FPN logln, int nparams, int nn, int whoseic) {
  FPN ic;
  if ( whoseic == 1 )
    ic = -(FPN) 2.0 * (logln  -  (FPN) 0.5 * (FPN) nparams * (FPN) log((FPN) nn) );       /*  c(n) = log(n) */
  else if ( whoseic == 2 )
    ic = -(FPN) 2.0 * (logln  -    (FPN) nparams   );                            /*  c(n) = 2 */
  else if ( whoseic == 3 )
    ic = -(FPN) 2.0 * (logln  -    (FPN) nparams * (FPN) log(log((FPN) nn)) );      /*  c(n) = 2 log(log(n)) */
  else if ( whoseic == 4 )
    ic = -(FPN) 2.0 * (logln  -    (FPN) nparams * (FPN) log((FPN) nn) );       /*  c(n) = 2 log(n) */
  else if ( whoseic == 5 )
    ic = -(FPN) 2.0 * (logln  -  (FPN)1.5 * (FPN) nparams * (FPN) log((FPN) nn) );       /*  c(n) = 3 log(n) */
  else if ( whoseic == 6 )
    ic = -(FPN) 2.0 * (logln    );                                                  /*  c(n) = 0 */
  
  return(ic);
}

/*
  Find the mimnode of type abd corresponding to the QTL.
  
  if 
     indic      what to check  qtl against
     
      -2          tptr->gt1, tptr->gt2
      -1          tptr->gt1
       1          tptr->qtl1
       2          tptr->qtl1, tptr->qtl2
*/
mimparam *WhichMIMnode(mimparam *mimptr, int qtl1, int qtl2, int abd, int indic ) {

  mimparam *tptr;

  for ( tptr=mimptr->next; tptr != NULL ; tptr = tptr->next )
    if ( indic == 1 ) {
      if ( tptr->abd == abd && tptr->qtl1 == qtl1 )
        return(tptr); /*  Return the correct node where qtl1 matches qtl */
    }
    else if ( indic == -1 )  {
      if ( tptr->abd == abd && tptr->gt1 == qtl1 )
        return(tptr); /*  Return the correct node where gt1 matches qtl  */
    }
    else if ( indic == -2 )  {
      if ( tptr->abd == abd  && tptr->gt1 == qtl1 && tptr->gt2 == qtl2 )
        return(tptr); /*  Return the correct node  where gt2 matches qtl */
    }
    else if ( indic == 2 )  {
      if ( tptr->abd == abd && tptr->qtl1 == qtl1 && tptr->qtl2 == qtl2 )
        return(tptr); /*  Return the correct node  where qtl2 matches qtl */
    }
      
  return(mimptr);  /*  Return the anchor if none found. */
}



/*
   This is meant to check each parameter for significance.   
   If the parameter is not significant, then remove it from the chain.   
*/
void ConfirmCurrentModel(FILE *outf,params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) {
  mimparam  *lptr,*aptr,*dptr;
  int go_on,nparams,comparewith ;
  FPN minlod;
  comparewith = 1;
  if ( theparams->mimwork[4] == 'D' )
    comparewith = 2;
  else if ( theparams->mimwork[4] == 'E' )
    comparewith = 3; 
  go_on = 1;
  while ( go_on == 1 ) {
    nparams = HowManyQTL(mimptr,0);
    if ( nparams > 0 ) {
      CheckThisModel(outf,theparams,themap,individs,mimptr,agptr,lnpk);
      minlod = theparams->mimlod;
      aptr = mimptr;
      for ( lptr=mimptr->next; lptr != NULL ; lptr = lptr->next ) 
        if (   ((lptr->abd == comparewith) || ( comparewith > 2 && lptr->abd > 2 )) && lptr->lnlike < minlod ) { 
        /*  find minimum that is less than threshold */
          minlod = lptr->lnlike;
          aptr = lptr;
        }
      if ( aptr->abd != 0  ) { /*  Found one.   */
        if ( theparams->mimwork[4] == 'D'  ) {
          aptr = delete_mimparam(aptr);
          if ( theparams->verbosity == 1 ) 
            printf("d");   
        }  
        else if (   theparams->mimwork[4] == 'E' ) {
          aptr = delete_mimparam(aptr);
          if ( theparams->verbosity == 1 ) 
            printf("e");   
        } else {
          dptr = WhichMIMnode(mimptr, aptr->gt1,0, 2, -1);
          if ( dptr != NULL )
            dptr = delete_mimparam(dptr);
          aptr = delete_mimparam(aptr);
          if ( theparams->verbosity == 1 ) 
            printf("-");   
        }  
      }
      else 
          go_on = 0;     
     } 
     else {
       go_on = 0;
       if ( theparams->verbosity == 1 ) 
          printf(" NO QTLs left!");     
     }
   }
}

/*
   Calculate a likelihood for each parameter dropped.

   At present, this assumes that all parameters are additive or dominant.   
   Need to generalize to epistatics as well.  

   For dominance, both the additive and dominance parameters are dropped, and the likelihood is placed in the
   additive spot.  Thus, we check only the additives in the calling routine ConfirmCurrentModel.
   
To Do  (Does it make sense?):   
   Do this based on the value of theparams->mimwork[4]
   
   If T, additive and dominant effects are tested together
   If D, only dominance effects are tested
   If E, only epistatic effects are tested
   
*/
void CheckThisModel(FILE *outf,params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) {
  mimparam  *aptr,*dptr;
  int i,m , nepis;
  FPN lrratio,lod, newic,initic;

  putline(outf, '*', 85);
  if ( theparams->mimwork[4] == 'D' )
    fprintf(outf,"\n\t  This Section checks each dominance effect for significance");
  else
    fprintf(outf,"\n\t  This Section checks each QTL for significance");
  putline(outf, '=', 85);
  fprintf(outf,"\n  QTL Type   Est      C   M    c1    c2    IC(k)   IC(k-1)    LR       LOD");
  putline(outf, '-', 85);
  initic =  EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
/*  */
  if ( theparams->mimwork[4] == 'T' ) { /* Check the additive/dominance effects,   */
	  m = HowManyQTL(mimptr,1);
	  for ( i=1; i<=m; i++) {
	    aptr = WhichMIMnode(mimptr,i,0,1,-1);
	    nepis = CheckForEpistatics( mimptr, aptr);
	    if (aptr->abd == 1 && nepis == 0 ) {
	      PullMIMnode(aptr);
	      dptr = WhichMIMnode(mimptr,i,0,2,-1);
	      if (  dptr->abd == 2 )
	        PullMIMnode(dptr);
	      newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	      PrintMIMnode(outf,aptr); 
	      lrratio = newic - initic ;
	      lod = lrratio * (FPN) 0.217;
	      fprintf(outf, "  %g   %g   %g  %g",initic,newic,lrratio,lod);
	      if ( dptr->abd == 2 )
	        PushMIMnode(dptr);
	      PushMIMnode(aptr);
	      aptr->lnlike = lrratio;
	      aptr->qptr1->tr2 = lod;
          UpdateGTindex(mimptr);
	    }  
	    else if ( aptr->abd == 1 )
	      aptr->lnlike = (FPN) 10.0 * theparams->mimlod ; 
	  }
  }
  else if (  theparams->mimwork[4] == 'D' ) { /* Check the dominance effects only */
	  m = HowManyQTL(mimptr,2);
	  for ( i=1; i<=m; i++) {
	    dptr = WhichMIMnode(mimptr,i,0,2,-1);
	    if (dptr->abd == 2 ) {
	      PullMIMnode(dptr);
	      newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	      PrintMIMnode(outf,dptr); 
	      dptr->lnlike = lrratio = newic - initic ;
	      lod = lrratio * (FPN) 0.217;
	      fprintf(outf, "  %g   %g   %g  %g",initic,newic,lrratio,lod);
	      PushMIMnode(dptr);
          UpdateGTindex(mimptr);
	    }  
	  }
  }
  else if (  theparams->mimwork[4] == 'E' ) { /* Check the epistatic effects only */
	  for ( dptr = mimptr->next ; dptr != NULL ; dptr = dptr->next ) {

	    if (dptr->abd > 2 ) {
	      PullMIMnode(dptr);
	      newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
	      PrintMIMnode(outf,dptr); 
	      dptr->lnlike = lrratio = newic - initic ;
	      lod = lrratio * (FPN) 0.217;
	      fprintf(outf, "  %g   %g   %g  %g",initic,newic,lrratio,lod);
	      PushMIMnode(dptr);
          UpdateGTindex(mimptr);
	    }  
	  }
  }

  putline(outf, '=', 85);
  fprintf(outf,"\n\nCompleted confirmation of parameters\n\n");
  putline(outf, '*', 85);
  fprintf(outf,"\n\n");

}


/*
  Assume that aptr is a node for an additive effect.
  Check all epistatic effects for whether they involve this
  putative gene.  Return the number of epistatic effects that
  involve this gene.  

*/
int CheckForEpistatics(mimparam *mimptr, mimparam *aptr) {
  mimparam *lptr;
  int nepis ;
  nepis = 0;

  for ( lptr = mimptr->next; lptr != NULL ; lptr = lptr->next ) 
    if ( lptr->abd > 2 )  
      if ( lptr->gt1 == aptr->gt1 || lptr->gt2 == aptr->gt1 )
        nepis = nepis+1;

  return(nepis);
}


/* For each QTL, map in entire interval and determine the best position.  */
void  JitterPositions(params *theparams,markermap *themap ,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk ) {
  mimparam  *lptr ;
  FPN minln;
  for ( lptr = mimptr->next; lptr != NULL ; lptr=lptr->next)     /*  For each QTL, test different positions in the interval. */
    if ( lptr->abd == 1 ) /*  This is an additive effect.  Do a jitter here.  */
      minln =  DoAJitter(theparams,themap,individs,mimptr,lptr,agptr,lnpk);
}


/* For this QTL, map in entire interval and determine the best position.  */
FPN  DoAJitter(params *theparams,markermap *themap, individual *individs,mimparam *mimptr,mimparam *lptr,genome *agptr,linpakws *lnpk) {
  FPN d1,d2,dist,pos,minln,bestd1,newlogln,newic,minic;
  mimparam *dptr;
/*   First, pull the node corresponding to this QTL (additive and dominance nodes) and get an IC for the model. */
  PullMIMnode(lptr);
  dptr = WhichMIMnode(mimptr,lptr->gt1,0,2,-1);
  if (  dptr->abd == 2 )
	        PullMIMnode(dptr);
  minic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk);
  if ( dptr->abd == 2 )
	        PushMIMnode(dptr);
  PushMIMnode(lptr);
  minln = (FPN) 0.0;
/*  minic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
  printf("\n   %f  ",minln);*/
  bestd1 = d1 = mapfunc( lptr->qptr1->c1, 2);
  d2 = mapfunc(lptr->qptr1->c2, 2);
  dist = d1 + d2;
/*   Now, check the interval for the best position of this QTL */
  for ( pos = (FPN) MIN_DIST ; pos < dist; pos = pos + theparams->walk ) {
    lptr->qptr1->c1 = mapfunc( pos, -2 );
    lptr->qptr1->c2 = mapfunc( dist-pos, -2);
    newic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
    newlogln = minic - newic;
/*    printf(" %f",newlogln);*/
    if ( newlogln > minln ) {
      bestd1 = pos;
      minln = newlogln;
    }
  }
/*  mypause();*/
  d2 = dist - bestd1;  /*  Time to update the best position. */
  lptr->qptr1->c1 = mapfunc(bestd1, -2);
  lptr->qptr1->c2  = mapfunc(d2, -2);  
  minln = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
  if ( theparams->verbosity == 1 ) {
        if ( bestd1 < d1 )
          printf("<");
        else if ( bestd1 > d1 )
          printf(">");
        else
          printf("=");
  }

  return(minln);       
}


/* 

For each QTL, determine whether an adjacent interval is a better place for the QTL.  

        prioric = posteric = thisic = Eval current model
        
      1.  Is there an open interval prior to the QTL interval on the same chromosome?   Use DoAJitter to get prioric
      2.  Is there an open interval after the QTL interval on the same chromosome?      Use DoAJitter to get posteric


   Do we drop and delete the lptr?  or just change the qtl data?

*/
void  TestAdjacentIntervals(params *theparams,markermap *themap, individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) {
  mimparam  *lptr ;
  int mark,isprior,isposter;
  FPN thisic, prioric, posteric,c1,c2;
  for ( lptr = mimptr->next; lptr != NULL ; lptr=lptr->next)     /*  For each QTL, test the adjacent intervals if they are open. */
    if ( lptr->abd == 1 ) { /*  This is an additive effect.  Do a test here.  */
      mark = lptr->qptr1->mrk;
      c1 = lptr->qptr1->c1;
      c2 = lptr->qptr1->c2;      
      prioric = posteric = thisic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
      lptr->qptr1->mrk = mark -1; /* change info in  lprt->qtl1 */
      isprior = isclearint(theparams,mimptr,lptr,agptr);
      if (  isprior == 1 )   
        prioric  =  DoAJitter(theparams,themap, individs,mimptr,lptr,agptr,lnpk);
      else
        prioric = prioric + (FPN) 1.0;
      lptr->qptr1->mrk = mark + 1;
      isposter = isclearint(theparams,mimptr,lptr,agptr);     
      if ( isposter == 1 ) 
        posteric =  DoAJitter(theparams,themap,individs,mimptr,lptr,agptr,lnpk);
      else
        posteric = posteric + (FPN) 1.0;
      
      if ( prioric < thisic && prioric <= posteric ) {
        lptr->qptr1->mrk = mark -1; 
        isprior = isclearint(theparams,mimptr,lptr,agptr);
        if ( theparams->verbosity == 1 )
          printf("<-");
      }
      else if ( posteric < thisic && posteric < prioric ) {
        lptr->qptr1->mrk = mark + 1;
        isposter = isclearint(theparams,mimptr,lptr,agptr);     
        if ( theparams->verbosity == 1 )
          printf("->");
      }
      else {
        lptr->qptr1->mrk = mark;            
        isposter = isclearint(theparams,mimptr,lptr,agptr);     
        if ( theparams->verbosity == 1 )
          printf(".=.");
      }
  }
  thisic = EvalModel(theparams,themap,individs,mimptr,agptr,lnpk) ;
}

/*
    Want to know whether the QTL in lptr has been put in an interval that
    already has a QTL, or if it has been put in an interval that doesn't exist.
    
    Return a 1 if it is in an interval without a pre existing QTL,
    else return a 0.

*/
int isclearint(params *theparams, mimparam *mimptr,mimparam *lptr,genome *agptr) {
  genome *tgptr;
  mimparam *tmim;
  int intexists,qtlexists,i;
  FPN intlength;
  i=theparams->traits;
  intexists =0;
  
  for ( tgptr=agptr; tgptr != NULL ; tgptr = tgptr->next )
    if (  tgptr->chrom ==lptr->qptr1->chrm && tgptr->markr ==  lptr->qptr1->mrk ) {
      intexists = 1;
      intlength = tgptr->dist;
    }

  qtlexists = 0;
  if ( intexists == 0 ) /* This interval doesn't exist */
    return(0);    
  else {
    for ( tmim = mimptr->next; tmim != NULL ; tmim = tmim->next ) 
      if ( tmim->abd == 1 && tmim != lptr ) { /*  This interval already has a QTL in it. */
        if ( tmim->qptr1->chrm == lptr->qptr1->chrm &&   tmim->qptr1->mrk ==   lptr->qptr1->mrk )
          qtlexists = 1;     
      }  
  }
  if ( qtlexists == 1 )
    return(0);
  lptr->qptr1->c1 = mapfunc((FPN) MIN_DIST , -2);
  lptr->qptr1->c2 =  mapfunc(   intlength - (FPN) 0.01 * (FPN) MIN_DIST  , -1 );
  return(1);
}


/* Iterate on E and M steps to refine the estimates of a model.  Set debugging to 2 to see the likelihood change
    over the iterations. */
int RefineEstimates(params *theparams,mimparam *mimptr, linpakws *lnpk) {
  int jj,ii,thestep;

  jj = 1;
  if ( theparams->verbosity == 1  && debugging < 2) 
    thestep = Rotator(0);
  
  if (  debugging > 1 )
    printf("\n  Here is how the likelihood is changing\n  Step      NewIC      Change   ");
  do {  /*  refine estimates*/    
    if ( theparams->verbosity == 1 && debugging < 2) 
      thestep = Rotator(thestep);
    mimptr->ovalue = lnpk->s2[theparams->whichtrait][0];  /*  Previous likelihood */
    EStepPi(theparams,mimptr,lnpk);
    MStepE(theparams,mimptr,lnpk);
    CalcMeanVar(theparams,mimptr,lnpk);  
    mimptr->value = lnpk->s2[theparams->whichtrait][0];  /*  New likelihood */
    if (  debugging > 1 )
      printf("\n  %d    %20.12f  %20.12f", jj, mimptr->value, mimptr->value - mimptr->ovalue );
    ii = jj;
  } while (  (jj = ShouldWeContinue(mimptr,jj) ) > 1 );
  if (  debugging > 1 )
    printf("\n");
  
  if ( theparams->verbosity == 1  && debugging < 2)  
    thestep = Rotator(1);
  return(ii);
}

/*  
      Print out the model 
*/
void PrintFinalModel(FILE *outf, params *theparams,markermap *themap, mimparam *mimptr,char *outstring, int nn) {

  mimparam *gptr;
  int count,k;
  putline(outf, '#', 80);
  fprintf(outf,"\n\t%s\n#",outstring );
  
  count = HowManyQTL(mimptr, 0);
  if ( count > 0 ) {
    if ( themap->tnames != NULL )
      fprintf(outf, "\n-trait       %5d      Analyzed trait [%s]", theparams->whichtrait,themap->tnames[theparams->whichtrait]);
    else
      fprintf(outf, "\n-trait       %5d      Trait analyzed", theparams->whichtrait);
    fprintf(outf, "\n-maxqtl      %5d    Maximum number of allowed QTL in the model",theparams->maxqtl);
    fprintf(outf, "\n-maxepis     %5d    Maximum number of epistatic terms allowed",theparams->maxepistatics);
    fprintf(outf, "\n-xic         %5d    Code for the IC criterion",theparams->whoseic);
    fprintf(outf, "\n-walk        %5.1f    Walking speed for position refinement and QTL search, in cM",theparams->walk);
    fprintf(outf, "\n-LRthresh    %5.1f    Likelihood ratio threshold for adding/deleting a QTL",theparams->mimlod);
    fprintf(outf, "\n-workcode    %s       Code indicating what to do",theparams->mimwork);
    fprintf(outf, "\n-modelfile   %s  \n#",theparams->mqtfile);
	  count = HowManyQTL(mimptr, 1);
	  fprintf(outf, "\n-Aqtl     %3d    Number of QTL with additive effects",count);
	  count = HowManyQTL(mimptr, 2);
	  fprintf(outf, "\n-Dqtl     %3d    Number of QTL with dominance effects",count);
	  count = HowManyQTL(mimptr, 3);
	  fprintf(outf, "\n-AAqtl    %3d    Number of QTL with additive by additive effects",count);
	  count = HowManyQTL(mimptr, 4);
	  fprintf(outf, "\n-ADqtl    %3d    Number of QTL with additive by dominance effects",count);
	  count = HowManyQTL(mimptr, 5);
	  fprintf(outf, "\n-DAqtl    %3d    Number of QTL with dominance by additive effects",count);
	  count = HowManyQTL(mimptr, 6);
	  fprintf(outf, "\n-DDqtl    %3d    Number of QTL with dominance by dominance effects",count);
	  fprintf(outf, "\n________________");
	  count = HowManyQTL(mimptr, 0);
	  fprintf(outf, "\n-Total    %3d    Number of parameters in this model",count);
      count =   (int)  floor(  2.0 * sqrt((FPN) nn) ) ;
	  fprintf(outf, "\n-Maximum  %3d    Number of parameters allowed in this model (2sqrt(n))",count);
	  PrintModelHeader(outf);
	  for ( gptr=mimptr->next; gptr != NULL ; gptr = gptr->next )  
	    PrintMIMnode(outf, gptr); 
	  fprintf(outf, "\n-e");
  }
  else
      fprintf(outf, "\n\n\nNo QTL found at this stage...\n\n");

  putline(outf, '#', 80);
  
  	    if ( themap->ParentalDiff != NULL ) {
	      strcat(gname, "y\0");
          k = OttoJones(theparams ,outf, themap->ParentalDiff, gname,mimptr );
        }

  putline(outf, '#', 80);
  
  fflush(outf);
} 

/*  
      Print out the model in Rqtl.out format.  This generally goes at the end
      and can be used as a new initial model. 
*/
void PrintRqtlout(FILE *outf, params *theparams,markermap *themap,mimparam *mimptr,int position) {
  mimparam *aptr,*dptr;
  int m,i,r,s;
  if ( position == 1 ) {
    fprintf(outf,"\n#\n#");
    fprintf(outf,"\n-start model\n#");
    fprintf(outf,"\n#  The model defined here is in Rqtl.out format and can be used as");
    fprintf(outf,"\n#  an input model for Rcross or for another run of MImapqtl.  Preplot will process it as well.");
    fprintf(outf,"\n#\n-t     %d   is the number of traits\n#\n#",theparams->traits);
  }
  else if (position == 2) {
  	  m = HowManyQTL(mimptr,1);
	  fprintf(outf,"\n-k    %d",m); 
	  fprintf(outf,"\n# for trait -named %s which is  -number     %d",themap->tnames[theparams->whichtrait],theparams->whichtrait);
	  fprintf(outf,"\n#\n#      #  ..Chrom..Markr.    .RecombiL.    .RecombiR.    .Additive.    .Dominance");
	  for ( i=1; i<=m; i++ ) {
	    aptr = WhichMIMnode(mimptr, i,0,1,-1);
	    dptr = WhichMIMnode(mimptr, i,0,2,-1);
	    if ( dptr->abd == 0 )
		  fprintf(outf, "\n-l %5d  %5d  %5d   %12.7f  %12.7f  %12.7f        0.0000", i, aptr->qptr1->chrm, aptr->qptr1->mrk, aptr->qptr1->c1, aptr->qptr1->c2, aptr->value);
	    else
		  fprintf(outf, "\n-l %5d  %5d  %5d   %12.7f  %12.7f  %12.7f  %12.7f", i, aptr->qptr1->chrm, aptr->qptr1->mrk, aptr->qptr1->c1, aptr->qptr1->c2, aptr->value, dptr->value);
	      
	  }
	  fprintf(outf,"\n#\n#   QTL1   QTL2   Type       Value  \n#");
	  for ( r=1; r<m; r++ )
	    for ( s=r+1; s<=m ; s++ ) {
	      aptr = WhichMIMnode(mimptr, r,s,3,-2);
	      if ( aptr->abd != 0 )
		    fprintf(outf, "\n-i %5d  %5d    AA   %12.4f ",  aptr->gt1, aptr->gt2, aptr->value );
		  else {
	        aptr = WhichMIMnode(mimptr, s,r,3,-2);
	        if ( aptr->abd != 0 )
		      fprintf(outf, "\n-i %5d  %5d    AA   %12.4f ",  aptr->gt2, aptr->gt1, aptr->value );	  
		  }
	      aptr = WhichMIMnode(mimptr, r,s,4,-2);
	      if ( aptr->abd != 0 )
		    fprintf(outf, "\n-i %5d  %5d    AD   %12.4f ",  aptr->gt1, aptr->gt2, aptr->value );
		  else {
	        aptr = WhichMIMnode(mimptr, s,r,4,-2);
	        if ( aptr->abd != 0 )
		      fprintf(outf, "\n-i %5d  %5d    AD   %12.4f ",  aptr->gt2, aptr->gt1, aptr->value );	  
		  }
	      aptr = WhichMIMnode(mimptr, r,s,5,-2);
	      if ( aptr->abd != 0 )
		    fprintf(outf, "\n-i %5d  %5d    DA   %12.4f ",  aptr->gt1, aptr->gt2, aptr->value );
		  else {
	        aptr = WhichMIMnode(mimptr, s,r,5,-2);
	        if ( aptr->abd != 0 )
		      fprintf(outf, "\n-i %5d  %5d    DA   %12.4f ",  aptr->gt2, aptr->gt1, aptr->value );	  
		  }
	      aptr = WhichMIMnode(mimptr, r,s,6,-2);
	      if ( aptr->abd != 0 )
		    fprintf(outf, "\n-i %5d  %5d    DD   %12.4f ",  aptr->gt1, aptr->gt2, aptr->value );    
		  else {
	        aptr = WhichMIMnode(mimptr, s,r,6,-2);
	        if ( aptr->abd != 0 )
		      fprintf(outf, "\n-i %5d  %5d    DD   %12.4f ",  aptr->gt2, aptr->gt1, aptr->value );	  
		  }
	    }
	  fprintf(outf,"\n#\n#********* ********* ********* ********* \n#");
  }
  else {
    fprintf(outf,"\n#\n#End of this block of information\n#\n-end model\n#\n");
  }
  fflush(outf);
} 

/*  
         For a model, print a header for the table of parameters
*/
void PrintModelHeader(FILE *outf) {
  fprintf(outf, "\n#\n#   Here is a summary of the QTL ");
  fprintf(outf, "\n--------------------------------------------------------------------------------------");
  fprintf(outf, "\n                   Position (Main effect)                 Epistatic effect Position   ");
  fprintf(outf, "\n                   -----------------------------   -----------------------------------");
  fprintf(outf, "\n QTL Effect Value     C   M      c1         c2     (QTL) (C)   (M)     (c1)       (c2) \n-s");
}

/*  
     Print a node in the mimparam chain 
*/
void PrintMIMnode(FILE *outf,mimparam *gptr) {

  if ( gptr->abd != 0 ) {   
    fprintf(outf,"\n %3d  ",gptr->qtl1);           
    if ( gptr->abd == 1 ) 
      fprintf(outf," A  ");
    else if ( gptr->abd == 2 )
      fprintf(outf," D  ");
    else if ( gptr->abd == 3 )
      fprintf(outf," AA ");
    else if ( gptr->abd == 4 )
      fprintf(outf," AD ");
    else if ( gptr->abd == 5 )
      fprintf(outf," DA ");
    else if ( gptr->abd == 6 )
      fprintf(outf," DD ");
    fprintf(outf," %6.3f   %3d %3d %10.8f %10.8f", gptr->value,gptr->qptr1->chrm,gptr->qptr1->mrk,gptr->qptr1->c1,gptr->qptr1->c2);
    if (    gptr->abd > 2 )    
      fprintf(outf,"  %3d  %3d   %3d  %10.8f %10.8f", gptr->qtl2,gptr->qptr2->chrm,gptr->qptr2->mrk,gptr->qptr2->c1,gptr->qptr2->c2);
  }  
}

/*
   Count the number of nodes of the specified parameter type abd, 
   or count all nodes.
               0    all nodes
    abd   =    1    additive   a nodes
               2    dominance  d nodes
               3    a x a nodes
               4    a x d nodes          
               5    d x a nodes
               6    d x d nodes

Assume that mimptr is the anchor...start at the next one.
*/
int HowManyQTL(mimparam *mimptr,int abd) {
  mimparam *lptr;
  int count, mabd;
  
  count = 0;
  if ( abd < 0 ) {
    mabd = -abd;
    for ( lptr=mimptr->next; lptr != NULL ; lptr = lptr->next ) 
      if ( lptr->abd >= mabd ) 
        count = count+1;
  }
  else  {
    for ( lptr=mimptr->next; lptr != NULL ; lptr = lptr->next ) 
      if ( lptr->abd == abd || abd == 0) 
        count = count+1;
  }
  return(count);

}

/* 
   Convert theqtls to an mimptr for the trait of interest.    
   
   Need to put in epistatic effects.  
*/
mimparam *ConstructMIMparams(params *theparams, markermap *themap, aqtl *theqtls) {
  mimparam *mimptr;
  int ii, k,kk,r,s;
  
  mimptr = insert_mimparam(NULL , 0, (FPN) 1.0 , 0, 0);  /*  abdtype = 0 means that this is an anchor.*/
  if ( theqtls == NULL )
    return(mimptr);  /*   If there is no initial model, return an anchor to build upon. */
  k = 0;  
  for (kk = 1; kk <= themap->traits; kk++) {
    for (ii = 1; ii <= themap->knum[kk]; ii++) {
	  k = k + 1;
	  if ( kk == theparams->whichtrait ) {
	    mimptr = insert_mimparam(mimptr, 1, theqtls[k].a, ii, 0);   /* Put in a node for the additive effect. */
	    mimptr->qptr1 = &theqtls[k]; 
	    mimptr->q1array = 1; 
	    mimptr->qptr1->nptrs  = mimptr->qptr1->nptrs + 1;
	    if ( ( theparams->crosstype == 3 || theparams->crosstype == 4 ) && theqtls[k].d != (FPN) 0.0) { /* Put in a node for the dominance effect if Fx cross */
	      mimptr = insert_mimparam(mimptr, 2, theqtls[k].d, ii, 0);
	      mimptr->qptr1 = &theqtls[k]; 
	      mimptr->q1array = 1; 
	      mimptr->qptr1->nptrs  = mimptr->qptr1->nptrs + 1;
	      if ( ii == themap->knum[kk] ) 
	        for ( r=2; r<=themap->knum[kk]; r++ ) 
	          for (s=1; s<r ; s++ ) {
	            if ( theqtls[k].episAADD[s][r] != (FPN) 0.0 ) {
	              mimptr = insert_mimparam(mimptr, 6, theqtls[k].episAADD[s][r], s, r);   /* Put in a node for the dominance by dominance effect. */
   	              mimptr->qptr1 = &theqtls[k-themap->knum[kk]+s]; mimptr->qptr1->nptrs = mimptr->qptr1->nptrs + 1;
   	              mimptr->qptr2 = &theqtls[k-themap->knum[kk]+r]; mimptr->qptr2->nptrs = mimptr->qptr2->nptrs + 1; 
	              mimptr->q1array = 1; mimptr->q2array = 1; 
	            }
	            if ( theqtls[k].episADDA[s][r] != (FPN) 0.0 ) {
	              mimptr = insert_mimparam(mimptr, 5, theqtls[k].episADDA[s][r] , s, r);   /* Put in a node for the dominance by additive effect. */
   	              mimptr->qptr1 = &theqtls[k-themap->knum[kk]+s]; mimptr->qptr1->nptrs = mimptr->qptr1->nptrs + 1;
   	              mimptr->qptr2 = &theqtls[k-themap->knum[kk]+r]; mimptr->qptr2->nptrs = mimptr->qptr2->nptrs + 1; 
	              mimptr->q1array = 1; mimptr->q2array = 1; 
	            }
	            if ( theqtls[k].episADDA[r][s] != (FPN) 0.0 ){
	              mimptr = insert_mimparam(mimptr, 4, theqtls[k].episADDA[r][s] , s, r);   /* Put in a node for the additive by dominance effect. */
   	              mimptr->qptr1 = &theqtls[k-themap->knum[kk]+s]; mimptr->qptr1->nptrs = mimptr->qptr1->nptrs + 1;
   	              mimptr->qptr2 = &theqtls[k-themap->knum[kk]+r]; mimptr->qptr2->nptrs = mimptr->qptr2->nptrs + 1; 
	              mimptr->q1array = 1; mimptr->q2array = 1; 
	            }
	          }
	    }
	    if ( ii == themap->knum[kk] ) 
	      for ( r=2; r<=themap->knum[kk]; r++ ) 
	        for (s=1; s<r ; s++ ) 
	          if ( theqtls[k].episAADD[r][s] != (FPN) 0.0 )  {
	            mimptr = insert_mimparam(mimptr, 3, theqtls[k].episAADD[r][s], s, r);   /* Put in a node for the additive by additive effect. */
   	            mimptr->qptr1 = &theqtls[k-themap->knum[kk]+s]; mimptr->qptr1->nptrs = mimptr->qptr1->nptrs + 1;
   	            mimptr->qptr2 = &theqtls[k-themap->knum[kk]+r]; mimptr->qptr2->nptrs = mimptr->qptr2->nptrs + 1; 
	            mimptr->q1array = 1; mimptr->q2array = 1; 
	          }
	  }
	}
  }
  while (   mimptr->prev != NULL )
     mimptr = mimptr->prev ;
  UpdateGTindex(mimptr);
  return(mimptr);
}

/*
     Need to have an indexing system for the genotypes.  This updates it.
     
     qtl1, qtl2 reference the original numbers of the QTL.
     gt1, gt2 are a new numbering system that gets updated.
*/
void UpdateGTindex(mimparam *mimptr) {
  mimparam *tptr,*aptr; 
  int qtlnum;
  qtlnum = 0;
  for ( tptr = mimptr->next; tptr != NULL ; tptr = tptr->next ) 
    if (tptr->abd == 1 ) {  /* additive effects are redone */
      qtlnum = qtlnum + 1;
      tptr->gt1 = qtlnum;  
      tptr->gt2 = 0;
    }
  for ( tptr = mimptr->next; tptr != NULL ; tptr = tptr->next ) 
    if (tptr->abd == 2 ) {  /* dominance effects are redone */
      aptr =  WhichMIMnode(mimptr, tptr->qtl1,0, 1,1);
      tptr->gt1 = aptr->gt1;
      tptr->gt2 = 0;
    }  
  for ( tptr = mimptr->next; tptr != NULL ; tptr = tptr->next ) 
    if (tptr->abd > 2 ) {  /* epistatic effects are redone */
      aptr =  WhichMIMnode(mimptr, tptr->qtl1,0, 1,1);
      tptr->gt1 = aptr->gt1;
      aptr =  WhichMIMnode(mimptr, tptr->qtl2,0, 1,1);
      tptr->gt2 = aptr->gt1;      
    }
}
 
/* 
     Allocate space for a genome node and set elements to specified values  
*/
mimparam *insert_mimparam(mimparam *inptr, int abdtype, FPN initvalue, int qtl1, int qtl2) 
{
  mimparam *gptr;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  gptr = (mimparam *) malloc( (size_t) sizeof(mimparam));
#else
  gptr = (mimparam *) malloc( (unsigned) sizeof(mimparam));
#endif
 if ( debugging > 2) {
        sprintf(gwarn,"In alloc_mimparam(), allocated 1 mimparam node at %p\n",gptr);
        MemoryAccount(gwarn);
 }

  gptr->abd = abdtype;     /* 1, 2, 3, 4, 5, 6 for a, d, aa, ad, da, dd */
  gptr->ovalue = gptr->value = initvalue;    /* initial parameter value*/
  gptr->dbar = gptr->lnlike = (FPN) 0.0;
  gptr->gt1 = gptr->qtl1 = qtl1;  /* which QTL is this for */
  gptr->gt2 = gptr->qtl2 = qtl2;  /* if epistatic, which other QTL */
  gptr->design = 0;  /*  Design value for this parameter */
  gptr->qptr1 = NULL; gptr->q1array = 1;
  gptr->qptr2 = NULL; gptr->q2array = 1;
  if ( qtl1 < 0 && abdtype == 1) {
#if defined(MACWARRIOR) || defined(WINWARRIOR)  
    gptr->qptr1 = (aqtl *) malloc((size_t)  sizeof(aqtl) );
#else
    gptr->qptr1 = (aqtl *) malloc((unsigned)  sizeof(aqtl) );
#endif
    gptr->q1array = 0;
    CleanQTL( gptr->qptr1 , 1, NULL);
    gptr->qptr1->nptrs = 1;
  }
  gptr->prev = NULL;   /* Pointer to previous node */
  gptr->next = NULL;   /* Pointer to next node */
  if ( inptr != NULL ) {
    gptr->prev = inptr;
    gptr->next = inptr->next;
    inptr->next = gptr;
    if ( gptr->next != NULL )
      gptr->next->prev = gptr;
  }
/*  PrintMIMnode( stdout, gptr, NULL );*/
  return(gptr);
}

/*
    Delete an mimparam node.
 
    If inptr an interior or terminal node, then return previous node.
    if inptr is an anchor node, then return same.    
 */
mimparam *delete_mimparam(mimparam *inptr)
{
  mimparam *gptr;
  gptr = inptr;
  if ( gptr != NULL && gptr->abd == 0 )
    return(gptr);
  if ( inptr != NULL ) {
    if ( inptr->prev != NULL ) {
      gptr = inptr->prev;
      gptr->next = inptr->next;
      if ( gptr->next != NULL )
        gptr->next->prev = gptr;    
    }
    else if ( inptr->next != NULL ) {
      gptr = inptr->next;
      gptr->prev = NULL;    
    }
    else  {
      gptr = NULL;
    }  
    if ( debugging > 2 ) {
          sprintf(gwarn,"In delete_mimparam(), deallocated 1 mimparam node at %p\n",inptr);
          MemoryAccount(gwarn);
    }
    if ( inptr->qptr1 != NULL && inptr->qptr1->nptrs == 0 && inptr->abd == 1 ) 
        free((char *) inptr->qptr1 );  /*  We delete the aqtl structure if this is the last thing to point to it. */
    else if ( inptr->qptr1 != NULL )
      inptr->qptr1->nptrs = inptr->qptr1->nptrs - 1; /*  Otherwise, we decrement the indicator of how many things point to it. */
    if (  inptr->qptr2 != NULL && inptr->qptr2->nptrs == 0  && inptr->abd == 1 ) 
        free((char *) inptr->qptr2 );
    else if ( inptr->qptr2 != NULL )
      inptr->qptr2->nptrs = inptr->qptr2->nptrs - 1;

    free((char *) inptr);
  }
  return(gptr);
}

/*  Pull a  mimparam pointer from  a chain. */
void PullMIMnode(mimparam *lptr)
{
  if ( lptr->prev != NULL )  
    lptr->prev->next = lptr->next;
  if ( lptr->next != NULL )
    lptr->next->prev = lptr->prev; 
}

/*  Push a  mimparam pointer into a chain.   */
void PushMIMnode(mimparam *lptr)
{
  if ( lptr->prev != NULL )  
    lptr->prev->next = lptr ;
  if ( lptr->next != NULL )
    lptr->next->prev = lptr; 
}

/*    Design parameters xi, zi, betaij

For BC, RI lines,   we have a pair of genotypes   {(QQ, Qq), (Qq, qq), (QQ,qq)}.  Call them Q and q

       Q        q
 x    1/2    -1/2
 
For loci i,j    Use    xi xj  as the epistatic term  aa
 
For Fx lines, we have three genotypes (QQ, Qq, qq).   

        QQ       Qq      qq
  x      1        0      -1
  z     -1/2     1/2    -1/2

Parameterization  
--------------------------------------------------       
--------------------------------------------------       
                       a for crosstype
              ------------------------------------
              BC1   BC2     SFx     RFx      RIx
       d      1      2       3       4        5 
--------------------------------------------------       
qq   -1/2     0    -1/2     -1      -1       -1 
Qq    1/2   -1/2    1/2      0       0        0
QQ   -1/2    1/2     0       1       1        1 
--------------------------------------------------       
--------------------------------------------------       

For loci i,j, use    xi xj for aa,   xi zj for ad, zi xj for da and zi zj for dd

Note that for BC2, qq and Qq are encoded by 0 and 1, respectively, thus the column becomes  

-1/2        0
 1/2  =>  -1/2
  0        1/2   
  
And for RIx, qq and QQ are encoded by 0 and 1, respectively, thus the column becomes  

-1        0
 0  =>   -1        
 1        1   
                                                  
                                                           
*/

void init_design(mimparam *gptr, int *genotype, int crosstype) {
  mimparam *lgptr;
  FPN xz[3][6] = { { -(FPN) 0.5,  (FPN) 0.0,  (FPN) 0.0, -(FPN) 1.0, -(FPN) 1.0,  (FPN) 0.0 },
                      {  (FPN)0.5, -(FPN)0.5, -(FPN)0.5,  (FPN) 0.0,  (FPN) 0.0, -(FPN) 1.0 },
                      { -(FPN)0.5,  (FPN)0.5,  (FPN)0.5,  (FPN) 1.0,  (FPN) 1.0,  (FPN) 1.0 }   };  
  for ( lgptr=gptr; lgptr != NULL ; lgptr = lgptr->next) {
    if ( lgptr->abd == 1 )    /* 1   for an additive effect for QTL1    a  */
      lgptr->design = xz[genotype[lgptr->gt1]+1][crosstype] ;
    else if ( lgptr->abd == 2 ) /* 2   for a dominance effect for QTL1  d  */
      lgptr->design = xz[genotype[lgptr->gt1]+1][0] ;
    else if ( lgptr->abd == 3 )  /* 3   for additive by additive effect of QTL 1 and QTL 2     aa   */
      lgptr->design = xz[genotype[lgptr->gt1]+1][crosstype] * xz[genotype[lgptr->gt2]+1][crosstype] ;
    else if ( lgptr->abd == 4 ) /* 4   for additive by dominance effect of QTL 1 and QTL 2     ad  */
      lgptr->design = xz[genotype[lgptr->gt1]+1][crosstype] * xz[genotype[lgptr->gt2]+1][0] ;
    else if ( lgptr->abd == 5 ) /* 4   for dominance by additive effect of QTL 1 and QTL 2     da   */
      lgptr->design = xz[genotype[lgptr->gt1]+1][0] * xz[genotype[lgptr->gt2]+1][crosstype] ;
    else if ( lgptr->abd == 6 ) /* 3   for dominance by dominance effect of QTL 1 and QTL 2    dd   */
      lgptr->design = xz[genotype[lgptr->gt1]+1][0] * xz[genotype[lgptr->gt2]+1][0] ;
  
  }
}


/* Calculate the prior probabilities of the QTL given the flanking markers and
   map distance to the flanking markers. 
   
Genotypes look like this:
-------------------------
           Cross
      -------------------
Code  BC1   BC2  RI    Fx
-------------------------
 1    QQ    Qq   QQ    QQ
 0    Qa    qq   qq    Qa
-1                     qq
-------------------------
ngt    2     2    2     3
-------------------------
-------------------------

we return probablilites in lnpk->rsd   as [3]->pQQ, [2]->pQq and [1]->pqq
-------------------------
           Cross
      -------------------
Row   BC1   BC2  RI    Fx
-------------------------
 3    QQ    Qq   QQ    QQ
 2    Qa    qq   qq    Qa
 1          qq   qq    qq
-------------------------
-------------------------
   
  Why?   The WhichMultiLocusGT doesn't know about crosses.  If two classes, only 0 or 1 are returned.  
  We take into account the crosses here.  This also allows init_design to work without dealing with crosses 
*/
void CalcPriorMatrix(params *theparams,genome *gptr,individual *individs,int jj, linpakws *lnpk, mimparam *mptr) {
  int ii;
  mimparam *tmptr;
  genome   *tgptr,lg,rg;
  ii = 0;
  for ( tmptr=mptr; tmptr !=NULL ; tmptr = tmptr->next ) {
    if ( tmptr->abd == 1 ) {  
      ii = ii+1;
      tgptr = aqtl_genomenode(gptr, tmptr->qptr1);      
      if ( tgptr == NULL )
        exit(2);
      add_virtual_node(tgptr, &lg , &rg, tmptr->qptr1->c1, tmptr->qptr1->c2);  
      calc_cond_prob(theparams,&rg,individs,jj, &lnpk->rsd[3][ii], &lnpk->rsd[2][ii], &lnpk->rsd[1][ii]);
      del_virtual_node(tgptr);
      if (theparams->crosstype == 2 ) {
        lnpk->rsd[3][ii] = lnpk->rsd[2][ii];   
        lnpk->rsd[2][ii] = lnpk->rsd[1][ii];
      }
      else if ( theparams->crosstype == 5 ) 
        lnpk->rsd[2][ii] = lnpk->rsd[1][ii];
    }
  }
}

/* 
  Calculate the prior probabilities of the QTL given the flanking markers and
   map distance to the flanking markers.  
   
   This routine is no longer used...it was mainly for debugging purposes.   
*/
void ShowPriorMatrix(params *theparams,  linpakws *lnpk, mimparam *mptr) {
  int ii;
  FPN largest,total;
  mimparam *tmptr;
  FILE *outptr;
  
  outptr = fileopen( theparams->error, "a" );
  
  fprintf(outptr, "\n#  Prior Probabilities\n#QTL    pAA       pAa       paa");
  ii = 0;
  total = (FPN) 1.0;
  for ( tmptr=mptr; tmptr !=NULL ; tmptr = tmptr->next ) {
    if ( tmptr->abd == 1 ) {
        ii = ii+1;
        if ( lnpk->rsd[3][ii] > lnpk->rsd[2][ii] )
          largest = lnpk->rsd[3][ii];
        else
          largest = lnpk->rsd[2][ii] ;
        if ( lnpk->rsd[1][ii] > largest )
          largest = lnpk->rsd[1][ii];
        total = total *  largest;
        fprintf(outptr, "\n%3d  %8.5f  %8.5f  %8.5f",ii,lnpk->rsd[3][ii],lnpk->rsd[2][ii],lnpk->rsd[1][ii]);
    }
  }
  fprintf(outptr,"\n#   The most common joint genotype has frequency %g",total);
  fileclose( theparams->error, outptr);
}

/*

   For all QTL in the model, calculate the most probable joint genotypes for each  individual.   
   There will be some threshold for a joint genotype to be considered in the mapping.  It is
   0.01 at the time of coding.    Look for the MINJOINTFREQ variable in Main.h to change it.
   
   The lnpk->xsave will have the priors of the genotypes in lnpk->gtindex (they are indexed). 
*/
void CalculatePriors(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) {
  int  jj,m,k,cntr,traits;
  long ltypes, lii,lcnt;
  FPN freq,freqsum;
  m = HowManyQTL(mimptr,1);
  ltypes = (long)  pow( (FPN) theparams->ngt, (FPN) m );
  traits = themap->traits;
  for ( jj=1; jj<=theparams->nn; jj++) {
    if ( lnpk->xsave[jj] != NULL )   /*  Releasing and reallocating memory might be a problem in Windows....*/
      free_dvector(lnpk->xsave[jj],1,lnpk->gtindex[jj][0]) ;   
    if ( lnpk->xx[jj] != NULL )   
      free_dvector(lnpk->xx[jj],1,lnpk->gtindex[jj][0]) ;   
    if ( lnpk->gtindex[jj] != NULL )  
      free_lvector(lnpk->gtindex[jj],0,lnpk->gtindex[jj][0]) ;   
    CalcPriorMatrix(theparams,agptr,individs,jj,lnpk,mimptr);
    lcnt = 0L; cntr=0; freqsum=(FPN) 0.0;
    for (  lii = 1L ; lii <= ltypes; lii++ ) {
      WhichMultilocusGT(m, lii, lnpk->jpvt[theparams->whichtrait] ,theparams->ngt);
      freq = (FPN) 1.0;    
      for ( k=1; k<=m; k++ ) 
        freq = freq * lnpk->rsd[ lnpk->jpvt[theparams->whichtrait][k] + 2 ][k] ;
      if ( freq > (FPN) MINJOINTFREQ ) {
        cntr = cntr+1;
        if ( cntr <= (int) MAXGT ) {
          lnpk->gtindex[0][cntr] = lii;
          lnpk->xsave[0][cntr] = freq;  /*Make initial posteriors equal to priors */
          freqsum = freqsum + freq;
          lcnt = lcnt + 1L;
        }
        else {
          sprintf(gwarn,"\n\n WARNING:  %d QTL model exceeded max number of genotypes for individual %d\n",m,jj);         
          LogTheError(theparams->error, gwarn);
          if (theparams->verbosity == 1 ) 
            printf("\n WARNING:  %d QTL model exceeded max number of genotypes for individual %d\n",m,jj);
          lnpk->gtindex[0][0] = (long) MAXGT ; /*  This is something that may need to be fixed... */
        }
      }
    }    
    lnpk->gtindex[jj] = lvector(0,cntr);   
    if ( cntr <= (int) MAXGT ) 
      lnpk->gtindex[jj][0] = cntr;
    else
      lnpk->gtindex[jj][0] =  (long) MAXGT ;
    lnpk->xsave[jj] = dvector(1,cntr);   
    lnpk->xx[jj] = dvector(1,cntr);   
    if ( freqsum < 0.95 ) {      
      sprintf(gwarn,"\n\n WARNING:  %d QTL model priors for individual %d sum to %f\n",m,jj,freqsum);         
      LogTheError(theparams->error, gwarn);
    }
    for ( k=1; k<= (int) lnpk->gtindex[jj][0] ; k++ )   {/*  normalize the priors. */
      lnpk->xx[jj][k] = lnpk->xsave[jj][k] = lnpk->xsave[0][k] / freqsum ;
      lnpk->gtindex[jj][k] = lnpk->gtindex[0][k];
    }
    if ( debugging > 1 )
    	ShowPriorMatrix(  theparams,  lnpk, mimptr); 
  }
      if ( debugging > 1 )
        ShowSigGenotypes(  theparams,lnpk);  
}


/*  
    Once we have figured out the priors on all joint genotypes, write out the 
    integer codes for those with frequencies greater than MINFREQ 
    
    Nothing uses this now...it was mainly for debugging purposes.  
    
*/
void ShowSigGenotypes( params *theparams, linpakws *lnpk) {
  int cntr,jj,k;  
  FILE *outf;
  
  outf = fileopen( theparams->error, "a" );
  cntr = 0;
  fprintf(outf, "\n#\n#  Summary of which joint genotypes have frequencies greater than %f",(FPN) MINJOINTFREQ);
  fprintf(outf, "\n#\n# Individual   # GT  \n#------------------------->");
  for ( jj=1; jj<=theparams->nn; jj++) {
    if ( lnpk->gtindex[jj][0] > cntr ) 
      cntr = lnpk->gtindex[jj][0];
    fprintf(outf,"\n#  %5d   %4ld : ",jj, lnpk->gtindex[jj][0] );
    for ( k=1; k<= lnpk->gtindex[jj][0] ; k++ )
      fprintf( outf, " %ld", lnpk->gtindex[jj][k] );
  }
  fprintf(outf, "\n#\n#  The maximum number of sig. joint genotypes was %d \n#",cntr);
  fileclose( theparams->error, outf );
}


/*   
   This does equations 6 and 7 from Zeng, Kao and Basten, Genet. Res, 74:281
   
   Assume that mimptr is the anchor...start at next one.

*/
void  CalcMeanVar(params *theparams,mimparam *mimptr,linpakws *lnpk) {
  int m,i,j;
  FPN  rsum, jrsum, ijrsum, rssum, jrssum, ijrssum, ysum, ysum2, yijrsum;
  FPN mean, variance;
  mimparam *rmimptr,*smimptr;
  m = HowManyQTL(mimptr,1);
  
  yijrsum = ysum = ysum2 = ijrssum = ijrsum = (FPN) 0.0;
  for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ ) {
    ysum = ysum + lnpk->y[theparams->whichtrait][i];
    ysum2 = ysum2 + lnpk->y[theparams->whichtrait][i] * lnpk->y[theparams->whichtrait][i];
    jrsum = jrssum  = (FPN) 0.0;
    for ( j=1; j<= (int) lnpk->gtindex[i][0] ; j++ ) {      
      WhichMultilocusGT(m, lnpk->gtindex[i][j], lnpk->jpvt[theparams->whichtrait], theparams->ngt);      
      init_design(mimptr, lnpk->jpvt[theparams->whichtrait], theparams->crosstype);      
      rsum = rssum = (FPN) 0.0;
      for ( rmimptr=mimptr->next ; rmimptr != NULL ; rmimptr = rmimptr->next ) {
        rsum = rsum + rmimptr->value * rmimptr->design ;
        for ( smimptr=mimptr->next; smimptr != NULL ; smimptr = smimptr->next ) 
          rssum = rssum + rmimptr->value * rmimptr->design * smimptr->value * smimptr->design ;
      }
      jrsum = jrsum + rsum * lnpk->xx[i][j];
      jrssum = jrssum + rssum *  lnpk->xx[i][j];      
    }
    ijrsum  = ijrsum +   jrsum;
    ijrssum = ijrssum + jrssum;
    yijrsum = yijrsum + jrsum * lnpk->y[theparams->whichtrait][i] ;  
  }
  lnpk->s2[0][theparams->whichtrait] = mean = ( ysum - ijrsum ) / (FPN) lnpk->samplesize[theparams->whichtrait]  ;
  lnpk->s2[theparams->whichtrait][theparams->whichtrait] = variance = (ysum2 - (FPN) 2.0*mean*ysum + (FPN) lnpk->samplesize[theparams->whichtrait] *mean*mean -  (FPN) 2.0*yijrsum + (FPN) 2.0*mean*ijrsum + ijrssum) /  (FPN) lnpk->samplesize[theparams->whichtrait] ;  
}


/*   
   This does equations 18 and 19 from Zeng, Kao and Basten, Genet. Res, 74:281
   
   The Variance-Covariance matrix will be in lnpk->s2i
*/
FPN  CalcVarCovar(params *theparams,mimparam *mimptr,linpakws *lnpk) {
  int m,mt,i,j,r,s;
  FPN genvar,phenvar;
  mimparam *rmimptr , *smimptr;
  m = HowManyQTL(mimptr,1);
  mt = HowManyQTL(mimptr,0);
  lnpk->s2i = dmatrix(0,mt,0,mt);
  lnpk->ldx = mt;
  for ( r=0; r<=mt ; r++ ) 
    for ( s=0; s <= mt; s++ )  
      lnpk->s2i[r][s] = (FPN) 0.0 ;

  
  phenvar = (FPN) 0.0;
  for ( i=1; i<= lnpk->samplesize[theparams->whichtrait] ; i++ ) {
    phenvar = phenvar + (lnpk->y[theparams->whichtrait][i] - lnpk->pp1[0]) *  (lnpk->y[theparams->whichtrait][i] - lnpk->pp1[0]) ;
    for ( j=1; j <= (int) lnpk->gtindex[i][0] ; j++ ) {      
      WhichMultilocusGT(m, lnpk->gtindex[i][j], lnpk->jpvt[theparams->whichtrait], theparams->ngt);      
      init_design(mimptr, lnpk->jpvt[theparams->whichtrait], theparams->crosstype);   
      r=0;
      for ( rmimptr=mimptr->next ; rmimptr != NULL ; rmimptr = rmimptr->next )  {
        r = r + 1;
        lnpk->s2i[0][r] = lnpk->s2i[0][r] + lnpk->xx[i][j] * (FPN) pow( (rmimptr->design - rmimptr->dbar) * rmimptr->value , 2.0 ) ;      
        s = r;
        for ( smimptr=rmimptr->next ; smimptr != NULL; smimptr = smimptr->next ) {
          s = s + 1;
          lnpk->s2i[r][s] = lnpk->s2i[r][s] + lnpk->xx[i][j] * ((rmimptr->design - rmimptr->dbar)  * rmimptr->value) *((smimptr->design - smimptr->dbar)  * smimptr->value) ;
        }
      }
    }  
  }
/*

structure of s2i:   

row    0   1  2  ...   mt
col 0      variance of each param in row 0, 1-mt
    1       covariances of parameters 
    2             in the upper triangle
    .
    .
    .   R2 analog of 
        of covariances in
    mt   lower triangle
       ^ R2 of each param in col 1, rows 1-mt


[0][0] contains the phentypic variance
[1][1] contains the genetic variance

*/  
  genvar = (FPN) 0.0;
  for ( r=1; r<=mt ; r++ ) {/*  var/covar in upper triangle*/
    lnpk->s2i[0][r] = lnpk->s2i[0][r] / (FPN) lnpk->samplesize[theparams->whichtrait] ;
    genvar = genvar + lnpk->s2i[0][r];
    for ( s=r+1; s <= mt; s++ ) {
      lnpk->s2i[r][s] = (FPN) 2.0 * lnpk->s2i[r][s] / (FPN) lnpk->samplesize[theparams->whichtrait] ;
      genvar = genvar + lnpk->s2i[r][s];
    }
  }
  lnpk->s2i[1][1] = genvar;
  phenvar = phenvar /  (FPN) lnpk->samplesize[theparams->whichtrait];
  for ( r=1; r<=mt ; r++ ) {  /*  R2 values in lower triangle*/
    lnpk->s2i[r][0] = lnpk->s2i[0][r] / phenvar ;
    for ( s=r+1; s <= mt; s++ ) 
      lnpk->s2i[s][r] =  lnpk->s2i[r][s] / phenvar;
  }
  lnpk->s2i[0][0] = phenvar;
  return(genvar);
}

/*   
   This shows the  the Variance-Covariance matrix 
*/
void  ShowVarCovar(FILE *outf, params *theparams,markermap *themap,mimparam *mimptr,linpakws *lnpk,FPN genvar) {
  FPN *rowsums,doofus;
  int m,mt,r,s;
  mimparam *rmimptr , *smimptr;
  UpdateGTindex(mimptr);
  m = HowManyQTL(mimptr,1);
  mt = HowManyQTL(mimptr,0);
  rowsums = dvector(0,mt);
  doofus = genvar;
  PrintFinalModel(outf, theparams,themap, mimptr,"We use this model for the Variance-Covariance matrix.",lnpk->samplesize[theparams->whichtrait]);

  fprintf(outf,"\n\n\tThis is the Variance-Covariance Matrix.\n");
  fprintf(outf,"\nPhenotypic Variance: %12.4g\nGenetic Variance:    %12.4g\nResidual Variance:   %12.4g\n",lnpk->s2i[0][0] ,lnpk->s2i[1][1], lnpk->s2i[0][0] -lnpk->s2i[1][1]);

  r =   0;
  putline(outf, '-', 12+mt*13);
  putline(outf, '-', 12+mt*13);
  fprintf(outf,"\nQTL(s) Type");
  for ( rmimptr = mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next   ) 
    if ( rmimptr->abd < 3 ) 
      fprintf(outf,"     %3d     ",rmimptr->qtl1);  
    else         
      fprintf(outf,"   %3dx%3d   ",rmimptr->qtl1,rmimptr->qtl2);  
  putline(outf, '-', 12+mt*13);
  
  
  r =   0;
  
  for ( rmimptr = mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next   ) {
    r = r+1;
    if ( rmimptr->abd < 3 ) 
      fprintf(outf,"\n%3d    .",rmimptr->qtl1);  
    else         
      fprintf(outf,"\n%3dx%3d.",rmimptr->qtl1,rmimptr->qtl2);  
    if ( rmimptr->abd == 1 ) 
      fprintf(outf," A");
    else if ( rmimptr->abd == 2 )
      fprintf(outf," D");
    else if ( rmimptr->abd == 3 )
      fprintf(outf,"AA");
    else if ( rmimptr->abd == 4 )
      fprintf(outf,"AD");
    else if ( rmimptr->abd == 5 )
      fprintf(outf,"DA");
    else if ( rmimptr->abd == 6 )
      fprintf(outf,"DD");
    s = 0;
    for ( smimptr = mimptr->next; smimptr != NULL ; smimptr = smimptr->next   ) {
      s = s+1;
      if ( s < r ) 
        fprintf(outf,"             " );
      else if ( s == r ) 
        fprintf(outf," %12.4g",lnpk->s2i[0][r] );
      else 
        fprintf(outf," %12.4g",lnpk->s2i[r][s]);
      if ( s > r ) 
        rowsums[s] = rowsums[s] + (FPN)0.5* lnpk->s2i[r][s] ;
      else if ( s == r )
        rowsums[s] = rowsums[s] + lnpk->s2i[0][r] ;
      else
        rowsums[s] = rowsums[s] + (FPN)0.5* lnpk->s2i[s][r] ;
    }
  }
  putline(outf, '-', 12+mt*13);
  fprintf(outf,"\n  Sum     ");
  for ( r=1; r<=mt; r++ ) {
    rowsums[0] = rowsums[0] + rowsums[r];
    fprintf(outf," %12.4g",rowsums[r]);
  }
  fprintf(outf,"\nTotal      %12.4g",rowsums[0]);

  putline(outf, '*', 12+mt*13);
  
  for ( r=0; r<=mt; r++ )
    rowsums[r] = (FPN) 0.0;  
  fprintf(outf,"\n\n\tHere are the R2 values \n");
  fprintf(outf,"\nGenetic:  %12.4f\nResidual: %12.4f\n",lnpk->s2i[1][1]/lnpk->s2i[0][0],   (lnpk->s2i[0][0] -lnpk->s2i[1][1])/lnpk->s2i[0][0] );
  r =   0;
  putline(outf, '-', 12+mt*13);
  putline(outf, '-', 12+mt*13);
  fprintf(outf,"\nQTL(s) Type");
  for ( rmimptr = mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next   ) 
    if ( rmimptr->abd < 3 ) 
      fprintf(outf,"     %3d     ",rmimptr->qtl1);  
    else         
      fprintf(outf,"   %3dx%3d   ",rmimptr->qtl1,rmimptr->qtl2);  
  putline(outf, '-', 12+mt*13);
  r = 0;
  for ( rmimptr = mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next   ) {
    r = r+1;
    if ( rmimptr->abd < 3 ) 
      fprintf(outf,"\n%3d    .",rmimptr->qtl1);  
    else         
      fprintf(outf,"\n%3dx%3d.",rmimptr->qtl1,rmimptr->qtl2);  
    if ( rmimptr->abd == 1 ) 
      fprintf(outf," A");
    else if ( rmimptr->abd == 2 )
      fprintf(outf," D");
    else if ( rmimptr->abd == 3 )
      fprintf(outf,"AA");
    else if ( rmimptr->abd == 4 )
      fprintf(outf,"AD");
    else if ( rmimptr->abd == 5 )
      fprintf(outf,"DA");
    else if ( rmimptr->abd == 6 )
      fprintf(outf,"DD");
    s = 0;
    for ( smimptr = mimptr->next; smimptr != NULL ; smimptr = smimptr->next   ) {
      s = s+1;
      if ( s < r ) 
        fprintf(outf,"             " );
      else if ( s == r ) 
        fprintf(outf," %12.4f",lnpk->s2i[r][0] );
      else 
        fprintf(outf," %12.4f",lnpk->s2i[s][r]);
      if ( s > r ) 
        rowsums[s] = rowsums[s] + (FPN)0.5* lnpk->s2i[s][r] ;
      else if ( s == r ) 
        rowsums[s] = rowsums[s] + lnpk->s2i[r][0] ;
      else
        rowsums[s] = rowsums[s] + (FPN)0.5* lnpk->s2i[r][s] ;
    }
  }
  putline(outf, '-', 12+mt*13);
  fprintf(outf,"\n  Sum     ");
  for ( r=1; r<=mt; r++ ) {
    rowsums[0] = rowsums[0] + rowsums[r];
    fprintf(outf," %12.4f",rowsums[r]);
  }
  fprintf(outf,"\nTotal      %12.4f",rowsums[0]);
  fprintf(outf,"\n\n");
  putline(outf, '*', 12+mt*13);
  
  fprintf(outf,"\n\n\t        Estimates of QTL positions, effects and interactions ");
  r =   0;
  putline(outf, '-', 100);
  putline(outf, '-', 100);
  fprintf(outf,"\nQTL(pair)     Type     Chrom.  Marker   Position     dIC         Effect            Effect (%%)");
  putline(outf, '-', 100);
  r = 0;
  for ( rmimptr = mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next   ) {
    r = r+1;
    if ( rmimptr->abd < 3 ) 
      fprintf(outf,"\n%3d             ",rmimptr->qtl1);  
    else         
      fprintf(outf,"\n%3dx%-3d         ",rmimptr->qtl1,rmimptr->qtl2);  
    if ( rmimptr->abd == 1 ) 
      fprintf(outf," A");
    else if ( rmimptr->abd == 2 )
      fprintf(outf," D");
    else if ( rmimptr->abd == 3 )
      fprintf(outf,"AA");
    else if ( rmimptr->abd == 4 )
      fprintf(outf,"AD");
    else if ( rmimptr->abd == 5 )
      fprintf(outf,"DA");
    else if ( rmimptr->abd == 6 )
      fprintf(outf,"DD");
    if ( rmimptr->abd < 3 ) 
      fprintf(outf,"      %4d  %4d  %12.4f %10.4f  %12.4f     %12.1f", rmimptr->qptr1->chrm,rmimptr->qptr1->mrk, rmimptr->qptr1->s, rmimptr->qptr1->r2,rmimptr->value, (FPN) 100.0*rowsums[r]  );
    else /*(FPN) 100.0*lnpk->s2i[r][0]  or rowsums[r] */
      fprintf(outf,"                               %10.4f  %12.4f     %12.1f", rmimptr->lnlike, rmimptr->value, (FPN) 100.0*rowsums[r] );
  }
  putline(outf, '-', 100);
  fprintf(outf,"\n\n");

  free_dvector(rowsums,0,mt);
}

/*   
   This shows the breeding values for the individuals.   
   These refer to the equations in Zeng, Kao and Baston
(1999) see
http://statgen.ncsu.edu/zeng/GeneticalResearch-99-ZKB.pdf

*/
void  ShowGTvalues(FILE *outf,params *theparams,linpakws *lnpk,individual *individs) {
  int i,j,linelen;
  linelen = 55;
  putline(outf, '-', linelen);
  putline(outf, '-', linelen);
  fprintf(outf,"\n\nThe breeding values of the individuals are calculated");
  fprintf(outf,"\naccording to Zeng, Kao and Basten (1999) Estimating the"); 
  fprintf(outf,"\ngenetic architecture of quantitative traits. Genetical");
  fprintf(outf,"\nResearch, 74:279-290.   See that paper for equations ");
  fprintf(outf,"\nnumbered 14 and 15.");
  putline(outf, '-', linelen);
  fprintf(outf,"\nIndividual     Obs. Phenotype Equation 14  Equation 15" );
  putline(outf, '-', linelen);
  j=0;
  for ( i=1; i<= theparams->nn ; i++ ) 
    if (individs[i].print_flag == 'y') {
      j = j+1;
      if ( individs[i].name != NULL ) 
        fprintf(outf,"\n%12s %12.4f %12.4f   %12.4f",individs[i].name,individs[i].y[theparams->whichtrait],lnpk->pp1[j],lnpk->pp2[j]);
      else
        fprintf(outf,"\n  %8d   %12.4f %12.4f   %12.4f",i,individs[i].y[theparams->whichtrait],lnpk->pp1[j],lnpk->pp2[j]);
    }
  putline(outf, '-', linelen);
  fprintf(outf,"\n\n");
}

/*   
   This does equations 14 and 15 from Zeng, Kao and Basten, Genet. Res, 74:281
   In addition, calculates the yi and Dr needed for equations 18 and 19.

lnpk->pp1[0] holds ybar
*/
void  CalcGTvalues(params *theparams,mimparam *mimptr,linpakws *lnpk) {
  int m,i,j;
  mimparam *rmimptr ;
  m = HowManyQTL(mimptr,1);
  for ( rmimptr=mimptr; rmimptr != NULL ; rmimptr = rmimptr->next )
    rmimptr->dbar = (FPN) 0.0;
  lnpk->pp1[0] =  (FPN) 0.0;
  for ( i=1; i<= lnpk->samplesize[theparams->whichtrait] ; i++ ) {
    lnpk->pp1[0] = lnpk->pp1[0] + lnpk->y[theparams->whichtrait][i];
    lnpk->pp2[i] = lnpk->pp1[i] =  lnpk->s2[0][theparams->whichtrait] ;
    for ( j = 1; j <= (int) lnpk->gtindex[i][0] ; j++ ) {      
      WhichMultilocusGT(m, lnpk->gtindex[i][j], lnpk->jpvt[theparams->whichtrait], theparams->ngt);      
      init_design(mimptr, lnpk->jpvt[theparams->whichtrait], theparams->crosstype);   
      for ( rmimptr=mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next ) {
         rmimptr->dbar = rmimptr->dbar + lnpk->xx[i][j] * rmimptr->design;  /* sum_i sum_j pi_ij D_jr / n  needed in 18 and 19 */
         lnpk->pp1[i] = lnpk->pp1[i] + lnpk->xx[i][j] *  rmimptr->design * rmimptr->value ;   /*  Equation 14 */   
         lnpk->pp2[i] = lnpk->pp2[i] + lnpk->xsave[i][j] *  rmimptr->design * rmimptr->value ;   /*  Equation 15 */     
      }
    }  
  }
  for ( rmimptr=mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next )
    rmimptr->dbar = rmimptr->dbar / (FPN) lnpk->samplesize[theparams->whichtrait] ;  /* These are the Dbar's needed for 18 and 19*/
  lnpk->pp1[0] = lnpk->pp1[0] / (FPN) lnpk->samplesize[theparams->whichtrait] ;        /*  Mean values of the traits needed in 18 and 19 */
}


/*  
  This is equation 5 in Zeng, Kao and Basten
   It does the M step in the EM algorithm.  
   It essentially updates the parameter estimates  
*/
void MStepE(params *theparams,mimparam *mimptr,linpakws *lnpk) {

  FPN ijDsum, ijrsum, ssum;
  int i, j,m;
  mimparam *rmimptr, *smimptr;
  m = HowManyQTL(mimptr,1);
  
  for ( rmimptr=mimptr->next; rmimptr != NULL ; rmimptr = rmimptr->next ) {
    ijDsum = ijrsum = (FPN) 0.0;
	for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ ) {
	  for ( j=1; j<= (int) lnpk->gtindex[i][0] ; j++ ) {      
	    WhichMultilocusGT(m,   lnpk->gtindex[i][j], lnpk->jpvt[theparams->whichtrait], theparams->ngt);      
	    init_design(mimptr, lnpk->jpvt[theparams->whichtrait], theparams->crosstype);      
	    ssum =   (FPN) 0.0;
	    for ( smimptr=mimptr->next; smimptr != NULL ; smimptr = smimptr->next ) 
	      if ( smimptr != rmimptr )
	        ssum = ssum + smimptr->value * smimptr->design;
	    ijrsum = ijrsum +   lnpk->xx[i][j] * rmimptr->design * ( lnpk->y[theparams->whichtrait][i] - lnpk->s2[0][theparams->whichtrait] - ssum ) ;
	    ijDsum = ijDsum + lnpk->xx[i][j] *  rmimptr->design *  rmimptr->design ;	
	  }
    }
    rmimptr->ovalue = rmimptr->value;
    if ( ijDsum != (FPN) 0.0 )
       rmimptr->value = ijrsum / ijDsum ;
    else {
      if (theparams->verbosity == 1 ) 
        printf("!");
        
      sprintf(gwarn,"\n\n WARNING: ijDsum is zero in M step of EM algorithm. \n");         
      LogTheError(theparams->error, gwarn);
    }
  }
}


/*  
  Equation 4 from Zeng, Kao and Basten
  
  This is the E step in the EM algorithm and updates the posteriors in lnpk->xx (called \pi_{ij}^{[t+1]} on the
  left hand side of equation 4).
*/
void  EStepPi(params *theparams,mimparam *mimptr,linpakws *lnpk) {

  FPN ijsum,ssum,mean,variance,normfact;
  int i, j,m;
  mimparam  *smimptr;
  m = HowManyQTL(mimptr,1);  
  variance = lnpk->s2[theparams->whichtrait][theparams->whichtrait];
  normfact = (FPN) 1.0 / (FPN)sqrt( 2.0 * (FPN) PI * variance ) ;
  lnpk->s2[theparams->whichtrait][0] = (FPN) 0.0;
  for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ ) {
    ijsum = (FPN) 0.0;
	for ( j=1; j<= (int) lnpk->gtindex[i][0] ; j++ ) {      
	    WhichMultilocusGT(m, lnpk->gtindex[i][j], lnpk->jpvt[theparams->whichtrait], theparams->ngt);      
	    init_design(mimptr, lnpk->jpvt[theparams->whichtrait], theparams->crosstype);      
	    ssum =   (FPN) 0.0;
	    for ( smimptr=mimptr->next; smimptr != NULL ; smimptr = smimptr->next ) 
	        ssum = ssum + smimptr->value * smimptr->design;
	    mean = lnpk->s2[0][theparams->whichtrait] +  ssum;
	     
	    lnpk->xx[i][j] =  lnpk->xsave[i][j] * normfact *  (FPN) exp( -0.5 * pow( lnpk->y[theparams->whichtrait][i] - mean, 2.0)/variance);	 
	    ijsum = ijsum +  lnpk->xx[i][j];    
	}
	if (ijsum > (FPN) 0.0 )
	  lnpk->s2[theparams->whichtrait][0] = lnpk->s2[theparams->whichtrait][0] + (FPN) log( ijsum );
	for ( j=1; j<=lnpk->gtindex[i][0] ; j++ )       /*  Normalize, and keep track of log of the likelihood */
        lnpk->xx[i][j] = lnpk->xx[i][j] / ijsum ;	
  }
}

/*
    Need to decide on a stopping criterion for the EM algorithm.  


 1.   If  |lnL(t)| > 1.0 and | lnL(t+1) - lnL(t) | / |lnL(t)|   < STOP_MIM , then stop. 
      If  |lnL(t)| <= 1.0 and | lnL(t+1) - lnL(t) |    < STOP_MIM , then stop. 
      The likelihood of the model at step t [L(t)] is kept in the ovalue
      variable of the anchor node (mimptr->ovalue).  The updated likelihood [L(t+1)] 
      is stored in the value variable of the anchor node (mimptr->value).
    
 2.   Do at least 2 and no more than M_TIME iterations. 
*/
int ShouldWeContinue(mimparam *mimptr,int jj) {
  FPN delta;
  if ( jj < 3 ) 
    return( jj+1 );
  else if ( jj > (int) M_TIME ) 
    return(0);
  delta = (FPN) fabs(mimptr->value - mimptr->ovalue );
  if ( (FPN) fabs(mimptr->ovalue) > (FPN) 1.0 )
    delta = delta / (FPN)  fabs(mimptr->ovalue) ;
  
  if ( delta < (FPN) STOP_MIM )
    return(0);
  else
    return(jj+1);  
}


/* 
  Given a model, calculate the residuals
  
  Should this be the residuals plus the expected value of the current model
  for each individual?
*/
void CalculateResiduals(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) {
    int i,j,m,jj;
    FPN  rsum;
    mimparam *rmimptr ;
    m = HowManyQTL(mimptr,0);

    if ( m > 0 ) {     /*  Do this only if there are some QTLs */
      UpdateGTindex(mimptr);
      CalculatePriors(theparams,themap,individs,mimptr,agptr,lnpk);  
      for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ ) {
        rsum = (FPN) 0.0;
        for ( j=1; j<= (int) lnpk->gtindex[i][0] ; j++ ) {      
          WhichMultilocusGT(m, lnpk->gtindex[i][j], lnpk->jpvt[theparams->whichtrait], theparams->ngt);      
          init_design(mimptr, lnpk->jpvt[theparams->whichtrait], theparams->crosstype);      
          for ( rmimptr=mimptr->next ; rmimptr != NULL ; rmimptr = rmimptr->next ) 
            rsum = rsum + rmimptr->value * rmimptr->design  * lnpk->xx[i][j] ;          
         }
         lnpk->y[theparams->whichtrait][i] = lnpk->y[theparams->whichtrait][i] - rsum;  /* Replace trait values with residuals */
      }
  }
  else {     /*  Otherwise, subtract the mean */
    rsum = (FPN) 0.0;
    for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ )
      rsum = rsum + lnpk->y[theparams->whichtrait][i] ;
    rsum = rsum / (FPN) lnpk->samplesize[theparams->whichtrait]  ;
    for (   i = 1; i <= lnpk->samplesize[theparams->whichtrait] ; i++ )
      lnpk->y[theparams->whichtrait][i] =  lnpk->y[theparams->whichtrait][i] - rsum;
  }
  
  jj = 0;
  for ( i=1; i<= theparams->nn; i++ ) 
    if ( individs[i].print_flag == 'y' ) {
      jj = jj+1; /* Replace trait values with residuals in individs structure */
      individs[i].y[theparams->whichtrait] = lnpk->y[theparams->whichtrait][jj];    
    }      

}




/* constants */
#define CHISQR1 3.841
/* for 95% confidence limit */
#define ITERATE 20


/* ***********************

OttoJones was modified by Chris Basten from main.  It is now a function that can be
inserted into QTL Cartographer.   Probably will go into Eqtl.   

main essentially provides a minimal user interface, some error checking, and a test of
the functions provided. Each function has a reference to the equation in paper, Otto and Jones 2000
Please note that there are two essential definitions above.
These are 'magic numbers' required in several of the sub routines.

 User friendly version of O J estimator 24/7/02 version 0.2  jones and otto.
cojo@ucdavis.edu


As inputs,

   outfile         a string to indicate the output file...use "-" for stdout.  
   fParentalDiff   difference in parental lines, no units
   fNumQtlFound    number of QTL in model
   fMinEffect      smallest (absolute value?) QTL additive effect
   fMeanEffect     mean of QTL additive effects...what if some are negative?
   scUseEstThress  indicator of whether to estimate threshold or use fMinEffect.

************************ */

int OttoJones(params *theparams,FILE *fileptr, FPN *fParentalDiff, char *scUseEstThress, mimparam *mimptr )
{
	FPN fNumQtlFound,  fMinEffect,  fMeanEffect,theeffect;
    FPN fTotNumQTl, fLower95, fUpper95;
    FPN dTau, dMeanEffMissed, dSegVarMissed,  dMinEffectUsed,  dEstMinEffect;
    int wo, kk, k, ii, lb, ub;
    mimparam *aptr;
    if ( fileptr == NULL )
      return(-1);
    fprintf(fileptr,"\n");
    fprintf(fileptr,"Welcome to the OJ estimator!\n");
    fprintf(fileptr,"See Otto and Jones, 2000, Genetics 156:2093-2107.\n\n");
    fprintf(fileptr,"This is a pared down version of what is available in the Mathamatica package\n");
    fprintf(fileptr,"(see http://www.zoology.ubc.ca/~otto/Research).\n\n");
    fprintf(fileptr,"This implimentation assumes an exponential distribution of QTL effects (CV=1).\n");
    fprintf(fileptr,"For more complicated distributions please consult us (cojo@ucdavis.edu).\n\n");
    fprintf(fileptr,"ASSUMPTIONS:\n");
    fprintf(fileptr,"1.  Fixed allelic differences in parents\n");
    fprintf(fileptr,"2.  Little to no transgressive segregation\n");
    fprintf(fileptr,"3.  Exponential distibution of QTL effects\n");
    fprintf(fileptr,"4.  The number of QTLs is substantially less than the number of intervals\n");
    fprintf(fileptr,"5.  The probablilty of detecting a QTL rises from near zero to one around the threshold of detection\n");
    fprintf(fileptr,"\n");
    fprintf(fileptr,"Several of these assumptions were relaxed in Otto and Jones (see text for details).\n");
    fprintf(fileptr,"\n");	 

    lb = 1;
    ub = theparams->themap->traits;
    if ( mimptr != NULL )
      lb = ub = theparams->whichtrait;
    wo = 0;
    for (kk = lb; kk <= ub; kk++) {
      if ( mimptr == NULL )
        k = theparams->themap->knum[kk];
      else
        k = HowManyQTL(mimptr, 0);
      if ( k > 0 ) {
      
	        fprintf(fileptr," \n\n" );
	      if ( theparams->themap->tnames != NULL )
	        fprintf(fileptr,"    The following is for trait number %d named %s\n",kk,theparams->themap->tnames[kk]);
	      else
	        fprintf(fileptr,"    The following is for trait number %d\n",kk);

	      fNumQtlFound = (FPN) k; 
	      fMeanEffect = (FPN) 0.0;
	      fMinEffect = (FPN) FLT_MAX ;
	      for (ii = 1; ii <= k; ii++) {
	        if ( mimptr == NULL ) {
		      wo = wo + 1;	
		      theeffect = theparams->theqtls[wo].a;
		    }
		    else {
	          aptr = WhichMIMnode(mimptr, ii,0,1,-1);
	          theeffect = aptr->value;
	        }
	          if ( fMinEffect > (FPN) fabs(theeffect)   )
	            fMinEffect = (FPN) fabs(theeffect) ;
	          fMeanEffect = fMeanEffect + theeffect ;
	      }

	        
	        
	      fMeanEffect = fMeanEffect / fNumQtlFound; 
	      fprintf(fileptr,"    You found %d QTL for the current trait with a parental difference of %f.\n",k, fParentalDiff[kk]);
	      fprintf(fileptr,"    The average effect is %f and the minimum effect is %f\n  ",fMeanEffect,fMinEffect);

		  if (fNumQtlFound < 3) 
			    fprintf(fileptr,"No estimates because the estimator performs poorly when less than 3 QTLs are detected.\n  ");
	      else if (fMeanEffect < fMinEffect) 
			    fprintf(fileptr,"No estimates because the minimum estimated effect MUST be below the average detected effect.\n");
		  else {
		    if ( fMeanEffect < (1.25 * fMinEffect) )
			  fprintf(fileptr,"The estimator may not perform well when the threshold is less than 25%% of mean detected.\n");
		    fprintf(fileptr,"\n\n");
	/*  Have user determine how they wish to estimate the threshold of sensitivity for the analysis: default should be yes*/
		    if (strcmp(scUseEstThress, "n") == 0)
			  dMinEffectUsed = fMinEffect;
		    else if (strcmp(scUseEstThress, "N") == 0)
			  dMinEffectUsed = fMinEffect;
		    else {
			  dEstMinEffect = estimatethreshold(fNumQtlFound, fMeanEffect, fMinEffect);
			  fprintf(fileptr,"Your estimated QTL detection threshold is:			 %g\n\n", dEstMinEffect);
			  dMinEffectUsed = dEstMinEffect;
		    }	/* find tau */
		    dTau = calctau(fMeanEffect, dMinEffectUsed);	
		/* estimate number of QtL  */
		    fTotNumQTl = esttotalqtl(fParentalDiff[kk], fMeanEffect, dMinEffectUsed);
		    fprintf(fileptr,"Estimated  number of QTLs missed:					 %3.3f\n", (fTotNumQTl - fNumQtlFound));
		    fprintf(fileptr,"Estimated total number of QTLs:						 %3.3f\n\n", fTotNumQTl);
		/* estimate confidence intervals */
		    fprintf(fileptr,"These are the 95%% confidence intervals for the estimated true number of QTLs\n");
		    fLower95 = estimatelowerbound(fNumQtlFound, fMeanEffect, dMinEffectUsed, fParentalDiff[kk], fTotNumQTl);
		    fprintf(fileptr,"lower 95%% CI: 										  %3.3f\n", fLower95);
		    fUpper95 = estimateupperbound(fNumQtlFound, fMeanEffect, dMinEffectUsed, fParentalDiff[kk]);
		    if (fUpper95 < 0)
			  fprintf(fileptr," Upper 95%% CI exceeded 1000 QTLs.  The number cannot be effectively estimated (check for error)\n");
		    else
			  fprintf(fileptr,"Upper 95%% CI: 				   						 %3.3f\n\n", fUpper95);
		    dMeanEffMissed = meaneffectundetected(fMeanEffect, dTau);
		    fprintf(fileptr,"Mean effect of factors missed:						  %g", dMeanEffMissed);
		    dSegVarMissed = segvarundetected(dTau);
		    fprintf(fileptr," \n");
		    fprintf(fileptr,"Proportion of segregation variance explained by undetected factors %g\n", dSegVarMissed);
		  }
      }
	}
    return(0);
}


/*  eqn 6 p 2097 */
FPN esttotalqtl(FPN fParentalDiff, FPN fMeanEffect, FPN dMinEffectUsed)
{
    FPN fTotNumQTl;
    fTotNumQTl = -666;
    fTotNumQTl = fParentalDiff / (fMeanEffect - dMinEffectUsed);
    return fTotNumQTl;
}

/*  eqn 12 p 2098 */
FPN estimatethreshold(FPN fNumQtlFound, FPN fMeanEffect, FPN fMinEffect)
{
    FPN dSoughtMin;
    FPN dEstimateThresh = 666;
    FPN dMinMin, dMaxMin,  dOldMin;
    FPN dNumer,  dDenom;
    int iTest;
    iTest = 0;
    dOldMin = -666;
    dSoughtMin = 0;
    dSoughtMin = (FPN) fMinEffect / 2;
    dMinMin = 0;
    dMaxMin = fMinEffect;
    while (iTest < ITERATE)    {
/* ITERATE is defined as the number of iterations of where the estimate does not change
    required before exiting the loop.  see the #define macro at the top */
	  dSoughtMin = ((dMaxMin - dMinMin) / 2) + dMinMin;
	  if (dSoughtMin == dOldMin)
	    iTest++;
	  dNumer = (fMeanEffect - (FPN) exp((fMinEffect * fNumQtlFound) / (fMeanEffect - dSoughtMin)) * (fMeanEffect - dSoughtMin) - dSoughtMin + fMinEffect * fNumQtlFound);
	  dDenom = (fNumQtlFound * ((FPN)1 - (FPN)exp((fMinEffect * fNumQtlFound) / (fMeanEffect - dSoughtMin))));
	  dEstimateThresh = (dNumer / dDenom) - fMinEffect + dSoughtMin;	/* ok */
	  dOldMin = dSoughtMin;
	  if (dEstimateThresh > 0)
	    dMaxMin = dSoughtMin;
	  else
	    dMinMin = dSoughtMin;
    }
    return dSoughtMin;
}

/* eqn 8/9 ; p2097 */
FPN estimatelowerbound(FPN fNumQtlFound, FPN fMeanEffect, FPN dMinEffectUsed, FPN fParentalDiff, FPN fTotNumQTl)
{
    FPN fEstimateLwr = 666;
    FPN fSoughtLow, fMinMin,  fMaxMin, fOldMin;
/*FPN fNumer, fDenom;*/
    int iTest;
    iTest = 0;
    fOldMin = -666;
    fMinMin = 0;
    fMaxMin = fTotNumQTl;
    while (iTest < ITERATE)  {	/* ITERATE is defined as the number of
	   iterations of where the estimate does not change  equired before
	   exiting the loop.  see the #define macro  at the top */
	fSoughtLow = ((fMaxMin - fMinMin) / 2) + fMinMin;
	if (fSoughtLow == fOldMin)
	    iTest++;	/* CHISQR1  is defined above  it is simply  the chi square value
			   with 1 degree of  freedom */
	fEstimateLwr = ((FPN) CHISQR1 * fParentalDiff + (FPN) 2 * fParentalDiff * fNumQtlFound - (FPN) 2 * fMeanEffect * fSoughtLow * fNumQtlFound + (FPN) 2 * dMinEffectUsed * fSoughtLow * fNumQtlFound - (FPN) 2 * fParentalDiff * fNumQtlFound * (FPN) log(fParentalDiff) + (FPN) 2 * fParentalDiff * fNumQtlFound * (FPN) log((fMeanEffect - dMinEffectUsed) * fSoughtLow));
	fOldMin = fSoughtLow;
	if (fEstimateLwr > 0)
	{
	    fMaxMin = fSoughtLow;
	}
	else
	{
	    fMinMin = fSoughtLow;
	}
    }
    return fSoughtLow;
}
/* eqn 8/9 ; p2097 */
FPN estimateupperbound(FPN fNumQtlFound, FPN fMeanEffect, FPN dMinEffectUsed, FPN fParentalDiff)
{
    FPN fEstimateHi = 666;
    FPN fSoughtHi, fMinMin,  fMaxMin,  fOldMin;
/*FPN fNumer, fDenom;*/
    int iTest, iTest2 = 1;
    iTest = 0;
    fOldMin = -666;
    fMinMin = 0;
    fMaxMin = 1000;	/* this is arbitrary,  but highly unlikely  to be exceeded */
    while (iTest < ITERATE) {	
 /* ITERATE is defined as the number of iterations of where
	   the estimate does not change equired before exiting the loop.
	   see the #define macro  at the top */
	  fSoughtHi = ((fMaxMin - fMinMin) / 2) + fMinMin;
	  if (fSoughtHi == fOldMin)
	    iTest++;
	  fEstimateHi = ((FPN) CHISQR1 * fParentalDiff + (FPN) 2 * fParentalDiff * fNumQtlFound - (FPN) 2 * fMeanEffect * fSoughtHi * fNumQtlFound + 2 * dMinEffectUsed * fSoughtHi * fNumQtlFound - 2 * fParentalDiff * fNumQtlFound * (FPN) log(fParentalDiff) + (FPN) 2 * fParentalDiff * fNumQtlFound * (FPN) log((fMeanEffect - dMinEffectUsed) * fSoughtHi));
	  fOldMin = fSoughtHi;
	  if (fEstimateHi < 0)
	    fMaxMin = fSoughtHi;
	  else
	    fMinMin = fSoughtHi;
    }
    if (fSoughtHi > 999)  {
	  fSoughtHi = -1 * fSoughtHi;	/* returns a negative
					   number if arbitrary
					   limit is exceeded */
    }
    return fSoughtHi;
}

/* eqn 10; p2097 */
FPN meaneffectundetected(FPN fMeanEffect, FPN dTau)
{
    FPN dAnswer;
    dAnswer = (fMeanEffect * ((FPN)1 - (((FPN)1 - dTau) / ((FPN)1 - (FPN) exp(-1 * ((1 - dTau) / dTau))))));
    return dAnswer;
} 

/* bottom of p2097 */ 
FPN segvarundetected(FPN dTau)
{
    FPN dAnswer;
    dAnswer = (FPN) 1 - ((FPN) 0.5 * ( (FPN) exp((-1 * ((1 - dTau) / dTau))) * (1 + (1 / (dTau * dTau)))));
    return dAnswer;
}



/* eqn 10 ; p2097 */ 
FPN calctau(FPN fMeanEffect, FPN dMinEffectUsed)
{
    FPN dTau;
    dTau = (FPN) (fMeanEffect - dMinEffectUsed) / fMeanEffect;
    return dTau;
}


/* constants*/
#undef CHISQR1  
#undef ITERATE 

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MIfunc.c
------------------------------------------------------------------ */


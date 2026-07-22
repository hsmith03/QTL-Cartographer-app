/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MLnpkws.c) is part of QTL Cartographer
         
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



/*  allocate workspace for analyses. */
#include "Main.h"


/*

This was written by Chris Basten.  It creates some workspace for Linpak.
*/
linpakws *cd_linpakws(ilnpk,n,p,cd)
  linpakws *ilnpk;
  int n,p,cd;
{
  linpakws *olnpk;

  olnpk = NULL;
  if ( cd == 0 && ilnpk != NULL ) {
    if ( ilnpk->xx != NULL )
      free_dmatrix(ilnpk->xx,1,ilnpk->p,1,ilnpk->n);
    if ( ilnpk->xsave != NULL )
      free_dmatrix(ilnpk->xsave,1,ilnpk->p,1,ilnpk->n);
    if ( ilnpk->y != NULL )
      free_dmatrix(ilnpk->y,1,1,1,ilnpk->n);
    if ( ilnpk->rsd != NULL )
      free_dmatrix(ilnpk->rsd,1,1,1,ilnpk->n);
    if ( ilnpk->bb != NULL )
      free_dmatrix(ilnpk->bb,1,1,1,ilnpk->p);
    if ( ilnpk->qraux != NULL )
      free_dvector(ilnpk->qraux,1,ilnpk->p);
    if ( ilnpk->jpvt != NULL )
      free_imatrix(ilnpk->jpvt,1,1,1,ilnpk->p);
    if ( ilnpk->pp1 != NULL)
      free_dvector(ilnpk->pp1,1,ilnpk->n);
    if ( ilnpk->pp2 != NULL)
      free_dvector(ilnpk->pp2,1,ilnpk->n);
    if ( ilnpk->pv != NULL)
      free_dvector(ilnpk->pv,1,ilnpk->n);
    if ( ilnpk->qv != NULL)
      free_dvector(ilnpk->qv,1,ilnpk->n);
    if ( ilnpk->estimates != NULL)
      free_dmatrix(ilnpk->estimates,1,2,0,18);
    if ( ilnpk->bp != NULL)
      free_imatrix(ilnpk->bp, 1, 2, 1, ilnpk->k );;
    if ( ilnpk->wrsd != NULL)
      free_dmatrix(ilnpk->wrsd,1,1,1,ilnpk->n);
    if ( ilnpk->wy != NULL)
      free_dvector(ilnpk->wy,1,ilnpk->n);
    if ( debugging > 2 ) {
/*    char buffer[MAXNAME];*/
        sprintf(gwarn,"In cdlinpakws(), deallocated 1 linpakworkspace at %p\n", ilnpk);
        MemoryAccount(gwarn);
    }
    free( (char *) ilnpk );
  }
  else if ( cd == 1 ) {
#if defined(MACWARRIOR) || defined(WINWARRIOR)  
    olnpk = (linpakws *) malloc((size_t) sizeof(linpakws));
#else
    olnpk = (linpakws *) malloc((unsigned)   sizeof(linpakws));
#endif
    if ( debugging > 2 ) {
/*    char buffer[MAXNAME];*/
        sprintf(gwarn,"In cdlinpakws(), allocated 1 linpakworkspace at %p\n", olnpk);
        MemoryAccount(gwarn);
    }
    olnpk->xx = dmatrix(1,p,1,n);
    olnpk->xsave = dmatrix(1,p,1,n);
    olnpk->y = dmatrix(1,1,1,n);
    olnpk->rsd = dmatrix(1,1,1,n);
    olnpk->bb = dmatrix(1,1,1,p);
    olnpk->qraux = dvector(1,p);
    olnpk->jpvt = imatrix(1,1,1,p);
    olnpk->n = olnpk->ldx = n;
    olnpk->p = olnpk->k = p;
    olnpk->pp1 = NULL;
    olnpk->pp2 = NULL;
    olnpk->pv = NULL;
    olnpk->qv = NULL;
    olnpk->estimates = dmatrix(1,2,0,18);
    olnpk->thestats = olnpk->estimates[2] ;
    olnpk->bp = NULL;
    olnpk->wrsd = NULL;
    olnpk->wy = NULL;
  }
  return(olnpk);
}

/*
This was written by Chris Basten.  It creates some workspace for Linpak.
There is also a lot of extra space that is required for some of the 
subroutines in the QTL Cartographer system.  Windows doesn't like it
when you allocate and deallocate large chunks of memory.  It fragements the
memory and causes general protection faults.  This subroutine will create
workspace if cd==1 and dellocate it if cd==0.  You need to give an
input pointer (which can be NULL if cd==1) and it returns an output
pointer (NULL if cd==0).

The use of this is to allocate the space in the main section of a program, and
to pass the pointer to subroutions that can then access the memory.

ilnpk  input pointer
n      sample size
p      rows in design matrix
t      traits in multitrait analysis
olnpk  output pointer
*/
linpakws *cd_jzmapws(linpakws *ilnpk, int n, int p, int t, int cd)
{
  linpakws *olnpk;

  olnpk = NULL;
  if ( cd == 0 && ilnpk != NULL ) {
    if ( ilnpk->xx != NULL )
      free_dmatrix(ilnpk->xx,1,ilnpk->p,1,ilnpk->n);
    if ( ilnpk->xsave != NULL )
      free_dmatrix(ilnpk->xsave,1,ilnpk->p,1,ilnpk->n);
    if ( ilnpk->y != NULL )
      free_dmatrix(ilnpk->y,0,ilnpk->t,1,ilnpk->n);
    if ( ilnpk->rsd != NULL )
      free_dmatrix(ilnpk->rsd,0,ilnpk->t,1,ilnpk->n);
    if ( ilnpk->bb != NULL )
      free_dmatrix(ilnpk->bb,0,ilnpk->t,1,ilnpk->p);
    if ( ilnpk->qraux != NULL )
      free_dvector(ilnpk->qraux,1,ilnpk->p);
    if ( ilnpk->jpvt != NULL )
      free_imatrix(ilnpk->jpvt,0,ilnpk->t,1,ilnpk->p);
    if ( ilnpk->pp1 != NULL)
      free_dvector(ilnpk->pp1,1,ilnpk->n);
    if ( ilnpk->pp2 != NULL)
      free_dvector(ilnpk->pp2,1,ilnpk->n);
    if ( ilnpk->pv != NULL)
      free_dvector(ilnpk->pv,1,ilnpk->n);
    if ( ilnpk->qv != NULL)
      free_dvector(ilnpk->qv,1,ilnpk->n);
    if ( ilnpk->estimates != NULL)
      free_dmatrix(ilnpk->estimates,0,ilnpk->t,1,12);
    if ( ilnpk->bp != NULL)
      free_imatrix(ilnpk->bp, 1, 2, 1, ilnpk->k );;
    if ( ilnpk->wrsd != NULL)
      free_dmatrix(ilnpk->wrsd,0,ilnpk->t,1,ilnpk->n);
    if ( ilnpk->wy != NULL)
      free_dvector(ilnpk->wy,1,ilnpk->n);
    if ( ilnpk->s2 != NULL)
      free_dmatrix(ilnpk->s2,0,ilnpk->t,0,ilnpk->t);
    if ( ilnpk->s2i != NULL)
      free_dmatrix(ilnpk->s2i,0,ilnpk->t,0,ilnpk->t);
    if ( ilnpk->samplesize != NULL)
      free_ivector(ilnpk->samplesize,0,ilnpk->t);
    if ( ilnpk->pcnts != NULL )
      free_ivector(ilnpk->pcnts,1,ilnpk->ipos);
    if ( ilnpk->lratio != NULL )
      free_dvector(ilnpk->lratio,1,ilnpk->ipos);
    if ( ilnpk->ts0 != NULL )
      free_dvector(ilnpk->ts0,0,ilnpk->t);    
    if ( ilnpk->work != NULL )
      free_dvector(ilnpk->work,0,ilnpk->t);    
    if ( ilnpk->kpvt != NULL )
      free_ivector(ilnpk->kpvt,0,ilnpk->t);    
 if ( debugging > 2 ) {
        sprintf(gwarn,"In cdlinpakws(), deallocated 1 linpakws node at %p\n",ilnpk);
        MemoryAccount(gwarn);
 }
    free( (char *) ilnpk );
  }
  else if ( cd == 1 ) {
#if defined(MACWARRIOR) || defined(WINWARRIOR)
    olnpk = (linpakws *) malloc((size_t) sizeof(linpakws));
#else
    olnpk = (linpakws *) malloc((unsigned)   sizeof(linpakws));
#endif
 if ( debugging > 2 ) {
        sprintf(gwarn,"In cdlinpakws(), allocated 1 linpakws node at %p\n",olnpk);
        MemoryAccount(gwarn);
 }
    olnpk->xx = dmatrix(1,p,1,n);
    olnpk->xsave = dmatrix(1,p,1,n);
    olnpk->y = dmatrix(0,t,1,n);
    olnpk->rsd = dmatrix(0,t,1,n);
    olnpk->bb = dmatrix(0,t,1,p);
    olnpk->qraux = dvector(1,p);
    olnpk->jpvt = imatrix(0,t,1,p);
    olnpk->s2 = dmatrix(0,t,0,t);
    olnpk->s2i = dmatrix(0,t,0,t);
    olnpk->n = olnpk->ldx = n;
    olnpk->p = olnpk->k = p;
    olnpk->t = t;
    olnpk->work = NULL; /* dvector(0,t);   is this going to be a problem?  looks like it.*/
    olnpk->kpvt = NULL; /* ivector(0,t);*/
    olnpk->pp1 = NULL;
    olnpk->pp2 = NULL;
    olnpk->pv = NULL;
    olnpk->qv = NULL;
    olnpk->estimates = dmatrix(0,t,1,12);
    olnpk->bp = NULL;
    olnpk->wrsd = NULL;
    olnpk->wy = NULL;
    olnpk->samplesize = NULL;
    olnpk->pcnts = NULL;
    olnpk->lratio = NULL;
    olnpk->ipos = 0;
    olnpk->ts0 = NULL;
  }
  return(olnpk);
}

/*
This was written by Chris Basten.  It creates some workspace for MImapqtl.
There is also a lot of extra space that is required for some of the 
subroutines in the QTL Cartographer system.  Windows doesn't like it
when you allocate and deallocate large chunks of memory.  It fragements the
memory and causes general protection faults.  This subroutine will create
workspace if cd==1 and dellocate it if cd==0.  You need to give an
input pointer (which can be NULL if cd==1) and it returns an output
pointer (NULL if cd==0).

The use of this is to allocate the space in the main section of a program, and
to pass the pointer to subroutions that can then access the memory.

ilnpk  input pointer
n      sample size
p      columns of genotypes   
t      traits in multitrait analysis
olnpk  output pointer
*/
linpakws *cd_mimws(linpakws *ilnpk, int n, int p, int t, int cd)
{
  linpakws *olnpk;
  int i;
  
  olnpk = NULL;
  if ( cd == 0 && ilnpk != NULL ) {
    if ( ilnpk->xx != NULL ) {
      for ( i=1; i<=ilnpk->n ; i++ )
        if ( ilnpk->xx[i] != NULL )
          free_dvector( ilnpk->xx[i],1,ilnpk->gtindex[i][0]);
      free_dsvector(ilnpk->xx,1,ilnpk->n);
    }
    if ( ilnpk->xsave != NULL ) {
      for ( i=0; i<=ilnpk->n ; i++ )
        if ( ilnpk->xsave[i] != NULL )
          free_dvector( ilnpk->xsave[i],1,ilnpk->gtindex[i][0]);
      free_dsvector(ilnpk->xsave,0,ilnpk->n);
    }
    if ( ilnpk->y != NULL )
      free_dmatrix(ilnpk->y,0,ilnpk->t,1,ilnpk->n);

    if ( ilnpk->gtindex != NULL ) {
      for ( i=0; i<=ilnpk->n ; i++ )
        if ( ilnpk->gtindex[i] != NULL )
          free_lvector( ilnpk->gtindex[i],0,ilnpk->gtindex[i][0]);
      free_lsvector(ilnpk->gtindex,0,ilnpk->n);
    }
    if ( ilnpk->rsd != NULL )
      free_dmatrix(ilnpk->rsd,1,3,1, (int) MAXQTL);


    if ( ilnpk->s2 != NULL)
      free_dmatrix(ilnpk->s2,0,ilnpk->t,0,ilnpk->t);
    if ( ilnpk->s2i != NULL)
      free_dmatrix(ilnpk->s2i,0,ilnpk->ldx,0,ilnpk->ldx);

    if ( ilnpk->jpvt != NULL)
      free_imatrix(ilnpk->jpvt,0,ilnpk->t,1,(int) MAXQTL );
    if ( ilnpk->samplesize != NULL)
      free_ivector(ilnpk->samplesize,0,ilnpk->t);
      
    if ( ilnpk->pp1 != NULL)
      free_dvector(ilnpk->pp1,0,ilnpk->n);
    if ( ilnpk->pp2 != NULL)
      free_dvector(ilnpk->pp2,0,ilnpk->n);
  
      
 if ( debugging > 2 ) {
        sprintf(gwarn,"In cd_mimws(), deallocated 1 linpakws node at %p\n",ilnpk);
        MemoryAccount(gwarn);
 }
    free( (char *) ilnpk );
  }
  else if ( cd == 1 ) {
#if  defined(MACWARRIOR) || defined(WINWARRIOR)
    olnpk = (linpakws *) malloc((size_t) sizeof(linpakws));
#else
    olnpk = (linpakws *) malloc((unsigned)   sizeof(linpakws));
#endif
 if ( debugging > 2 ) {
        sprintf(gwarn,"In cd_mimws(), allocated 1 linpakws node at %p\n",olnpk);
        MemoryAccount(gwarn);
 }
    olnpk->xx = dsvector(1,n);    /*  dmatrix(1,n,1,p);   Use for posteriors*/
    olnpk->xsave =  dsvector(0,n);    /* dmatrix(1,n,1,p); Use for priors */
    olnpk->xsave[0] = dvector(1,p);
    olnpk->y = dmatrix(0,t,1,n);        /*use for traits*/
    olnpk->rsd = dmatrix(1,3,1, (int) MAXQTL );  /* For MImapqtl, hold the prior probabilities of qq, Qq and QQ genotypes over all loci */
    olnpk->bb = NULL;
    olnpk->qraux = NULL;
    olnpk->gtindex = lsvector(0,n);
    olnpk->gtindex[0] = lvector(0,p);   /* lmatrix(1,n,0,p); Use to index genotypes*/
    olnpk->gtindex[0][0] = p;
    olnpk->jpvt = imatrix(0,t,1,(int) MAXQTL );  /*reuseable vector for a multilocus genotype */
    olnpk->s2 = dmatrix(0,t,0,t);
    olnpk->s2i =NULL;  /*  Will be dmatrix(0, nparams, 0, nparams) */
    olnpk->n = olnpk->ldx = n;
    olnpk->p = olnpk->k = p;
    olnpk->t = t;
    olnpk->work = NULL; /* dvector(0,t);   is this going to be a problem?  looks like it.*/
    olnpk->kpvt = NULL; /* ivector(0,t);*/
    olnpk->pp1 = dvector(0,n);
    olnpk->pp2 = dvector(0,n);
    olnpk->pv = NULL;
    olnpk->qv = NULL;
    olnpk->estimates = NULL;
    olnpk->bp = NULL;
    olnpk->wrsd = NULL;
    olnpk->wy = NULL;
    olnpk->samplesize = ivector(0, t);
    olnpk->pcnts = NULL;
    olnpk->lratio = NULL;
    olnpk->ipos = 0;
    olnpk->ts0 = NULL;
  }
  return(olnpk);
}



/* How many traits in a multitrait analysis? Some of the other subroutines will
use the matrix pointed to by themap->tnames.  If it doesn't exist, create it. */
int how_many_traits(params *theparams,markermap *themap)
{  
  int jj,t,wt;
  wt = theparams->whichtrait;
  t = 0;
  if ( themap->tnames != NULL ) {
    for (jj = 1; jj <= theparams->traits; jj++) 
      if ( (wt == 0 && themap->tnames[jj][0] == '+') || (wt > theparams->traits && themap->tnames[jj][0] != '-') || (wt == jj) ) 
        t = t+1;
  }
  else {
    themap->tnames = cmatrix( 1, themap->traits, 0, MAXNAME);
    for ( jj = 1 ; jj <= themap->traits ; jj++ )
      if ( wt > theparams->traits )
        sprintf(themap->tnames[jj],"+Trait_%d",jj);
      else
        sprintf(themap->tnames[jj],"Trait_%d",jj);
    if ( wt > theparams->traits )
      t = theparams->traits;
    else
      t = 1;  
  }
  return(t);
}

/*  
Here's the rule, where wt is which trait and t is the number of traits:  

wt=0   If we are using no traits but those starting in +,
          then use only those individuals with phenotypic data at 
          ALL traits with names starting with +.
wt>t   If we are using all traits except those starting with -,
          then use only those individuals with phenotypic data at 
          ALL traits with names not starting with a -.
0<wt<t If we are doing single trait analysis, 
          then use only those individuals with phenotypic data.
*/
int actual_trait(params *theparams,markermap *themap,int trait)
{
  int k,jj,thistrait;
  k = 0;

  if ( theparams->whichtrait == 0 ) { /*  use no traits but those with names starting in + */
    for (jj = 1; jj <= theparams->traits; jj++) 
      if ( themap->tnames[jj][0] == '+'  ) {
        k = k+1;
        if ( k == trait )
          thistrait = jj;
      }
  }
  else if (theparams->whichtrait > theparams->traits ) {/*  use all traits except those with names starting in - */
    for (jj = 1; jj <= theparams->traits; jj++) 
      if ( themap->tnames[jj][0] != '-'  ) {
        k = k+1;
        if ( k == trait )
          thistrait = jj;
      }
  }
  else { /*  use  only trait wt ...  put it in row 1*/
    thistrait = theparams->whichtrait;
  }
  return(thistrait);
}

/*
 * Indicate whether 1 or more QTLs reside in each interval.  Use this only after
 * QTLs have been created, and a new genome linked list has been created.
 */

void zplace_qtls(params *theparams, aqtl *qtlptr, genome *first)
{
  int kk, knum;
  genome *gptr;
  if ( qtlptr != NULL ) {
	  knum = 0;
	  for (kk = 1; kk <= qtlptr[1].map->traits; kk++)
	    knum = knum + qtlptr[1].map->knum[kk];
	  for (kk = 1; kk <= knum; kk++) {
	    gptr = first;
	    while (gptr != NULL)
	      if (gptr->chrom == qtlptr[kk].chrm && gptr->markr == qtlptr[kk].mrk ) {
		    if ( qtlptr[kk].trait == theparams->whichtrait )
		      gptr->whichqtl = kk;
		    else
		      gptr->whichqtl = -kk;
		    gptr->pxo = gptr->pos + mapfunc(  qtlptr[kk].c1 , 1 ); 
		    /* rec to Morgans, pxo is position. Needed for pick_markers */
		    gptr->mxo = qtlptr[kk].c2;
		    gptr = NULL;
	      }
	      else
		    gptr = gptr->next;
	  }
  }
}

/*
This was written by Chris Basten.  It creates some workspace for MultiRegress.

 This subroutine will create
workspace if cd==1 and dellocate it if cd==0.  You need to give an
input pointer (which can be NULL if cd==1) and it returns an output
pointer (NULL if cd==0).

The use of this is to allocate the space in the main section of a program, and
to pass the pointer to subroutions that can then access the memory.

This routine defines the space for MultiRegress.

ilnpk  input pointer
n      sample size
p      maximum number of rows in design matrix
t      traits in multitrait analysis
olnpk  output pointer
*/
linpakws *cd_multiregressws(linpakws *ilnpk, int n, int p, int t, int cd)
{
  linpakws *olnpk;

  olnpk = NULL;
  if ( cd == 0 && ilnpk != NULL ) {
    if ( ilnpk->xx != NULL )
      free_dmatrix(ilnpk->xx,1,ilnpk->p,1,ilnpk->n);   
    if ( ilnpk->xsave != NULL )
      free_dmatrix(ilnpk->xsave,1,ilnpk->p,1,ilnpk->n);   
    if ( ilnpk->y != NULL )
      free_dmatrix(ilnpk->y,0,ilnpk->t,1,ilnpk->n);   
    if ( ilnpk->rsd != NULL )
      free_dmatrix(ilnpk->rsd,0,ilnpk->t,1,ilnpk->n);   
    if ( ilnpk->bb != NULL )
      free_dmatrix(ilnpk->bb,0,ilnpk->t,1,ilnpk->p);   
    if ( ilnpk->qraux != NULL )
      free_dvector(ilnpk->qraux,1,ilnpk->p);  
    if ( ilnpk->pp1 != NULL)
      free_dvector(ilnpk->pp1,1,ilnpk->n);   
    if ( ilnpk->pp2 != NULL)
      free_dvector(ilnpk->pp2,1,ilnpk->n);   
    if ( ilnpk->pv != NULL)
      free_dvector(ilnpk->pv,1,ilnpk->n);   
    if ( ilnpk->qv != NULL)
      free_dvector(ilnpk->qv,1,ilnpk->n);  
    if ( ilnpk->wrsd != NULL)
      free_dmatrix(ilnpk->wrsd,0,ilnpk->t,1,ilnpk->n);   
    if ( ilnpk->wy != NULL)
      free_dvector(ilnpk->wy,1,ilnpk->n);   
    if ( ilnpk->samplesize != NULL)
      free_ivector(ilnpk->samplesize,1,ilnpk->ipos);   
    if ( ilnpk->pcnts != NULL )
      free_ivector(ilnpk->pcnts,1,ilnpk->ipos);  
    if ( ilnpk->lratio != NULL )
      free_dvector(ilnpk->lratio,1,ilnpk->ipos);   
    if ( ilnpk->jpvt != NULL )
      free_imatrix(ilnpk->jpvt,0,ilnpk->t,1,ilnpk->p);  
    if ( ilnpk->bp != NULL )
      free_imatrix(ilnpk->bp,1,ilnpk->k,0,ilnpk->n);  
 if ( debugging > 2 ) {
        sprintf(gwarn,"In cdmultiregressws(), deallocated 1 linpakws node at %p\n",ilnpk);
        MemoryAccount(gwarn);
 }
    free( (char *) ilnpk );   
  }
  else if ( cd == 1 ) {
#if defined(MACWARRIOR) || defined(WINWARRIOR)
    olnpk = (linpakws *) malloc((size_t) sizeof(linpakws));
#else
    olnpk = (linpakws *) malloc((unsigned)   sizeof(linpakws));
#endif
 if ( debugging > 2 ) {
        sprintf(gwarn,"In cdlinpakws(), allocated 1 linpakws node at %p\n",olnpk);
        MemoryAccount(gwarn);
 }
    olnpk->xx = dmatrix(1,p,1,n);   /* Design matrix       */
    olnpk->xsave = NULL;
    olnpk->y = NULL;                /* Will be allocated in the program. */
    olnpk->rsd = dmatrix(0,t,1,n);  /* Residuals           */
    olnpk->bb = dmatrix(0,t,1,p);   /* Parameter estimates */
    olnpk->qraux = dvector(1,p);    /* auxiliary matrix for QR decomposition */
    olnpk->jpvt = imatrix(0,t,1,p); /* pivot matrix        */
    olnpk->s2 = NULL;
    olnpk->s2i = NULL;
    olnpk->n = n;    /* sample size   */
    olnpk->ldx = 0;  /* number of rows in design matrix due to categorical traits  */
    olnpk->p = p;    /* maximum number of parameters */
    olnpk->k = 0;    /* number of categorical traits */
    olnpk->t = t;    /* number of traits */
    olnpk->work = NULL;  
    olnpk->kpvt = NULL;  
    olnpk->pp1 = NULL;
    olnpk->pp2 = NULL;
    olnpk->pv = NULL;
    olnpk->qv = NULL;
    olnpk->estimates = NULL;
    olnpk->bp = NULL;  /*  Will have categorical trait data and be imatrix(ilnpk->bp,1,ilnpk->k,0,ilnpk->n)*/
    olnpk->wrsd = NULL;
    olnpk->wy = NULL;
    olnpk->samplesize = NULL; /*  will be ivector(1,lnpk->ipos)*/
    olnpk->pcnts = NULL;      /*  will be ivector(1,lnpk->ipos)*/
    olnpk->lratio = NULL;     /*  will be dvector(1,lnpk->ipos)*/
    olnpk->ipos = 0;          /*  Will be number of sites or analysis positions*/
    olnpk->ts0 = NULL;
  }
  return(olnpk);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MLnpkws.c
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (SRfunc.c) is part of QTL Cartographer
         
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



/*  Stepwise regression subroutines. */


#include "Main.h"


/*
Initialize the matrix used by the  stepwise regression subroutines

	  AA    Aa    aa
  a    1     0    -1
  d  -1/2   1/2  -1/2
  
  If missing, a = pAA - paa
              d = 1/2 - |a|

*/
FPN init_xmyv(FPN **xm, FPN *yv, genome *gptr, genome *lgptr, int n, individual *individs, markermap *themap, params *theparams)
{
  int ii,jj,pp,thismark,cross,chrom,mark,ind,roffset,kk,row;
  genome *tgptr;
  FPN ysum,ysum2,syy;
  
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  chrom = lgptr->chrom;
  mark = lgptr->markr;
  jj = 1;
  ysum=ysum2=(FPN)0.0;
  for ( ii = 1 ; ii <= theparams->nn ; ii++ ) 
    if ( individs[ii].y[theparams->whichtrait] > (FPN) MISS_VAL ) {
      yv[jj] = individs[ii].y[theparams->whichtrait];
      ysum = ysum + yv[jj];
      ysum2 = ysum2 + yv[jj] * yv[jj];
      jj = jj+1;
    }
  syy = ysum2 - ysum*ysum/((FPN) jj-1) ;
  for ( ii = 1 ; ii <= n ; ii++ )
     xm[1][ii] = (FPN)1.0;
  pp = 2;
  for ( tgptr = gptr ; tgptr != NULL ; tgptr = tgptr->next ) {
    if ( tgptr->markr != 0 && tgptr->whichqtl > 0  ) {
      jj = 1;
      for ( ii = 1 ; ii <= theparams->nn ; ii++ ) 
        if ( individs[ii].y[theparams->whichtrait] > (FPN) MISS_VAL ) {
          thismark =  individs[ii].markers[tgptr->chrom][tgptr->markr] ;
	      if (thismark < -1 || thismark > 1)
	        xm[pp][jj] = expect_mark(theparams,individs, ii, tgptr);
          else
	        xm[pp][jj] = (FPN) thismark ;
		  if (  cross == 3 || cross == 4) 
		    xm[pp+1][jj] = (FPN)0.5 - (FPN)fabs(xm[pp][jj]);
		  
		  jj = jj+1;
	   }
      pp = pp + 1;
      if (  cross == 3 || cross == 4 ) 
         pp = pp+1;
    }
  }
  roffset = 0;
  for ( ii = 1 ; ii <= themap->otraits ; ii++ )
    if ( themap->onames[ii][0] == '+' )
      roffset = roffset + themap->otypes[ii] - 1;
  if ( roffset > 0 ) {
    ind = 0;    
    for ( ii = 1 ; ii <= theparams->nn ; ii++ ) {
      row = pp;
      if (  individs[ii].y[theparams->whichtrait] > (FPN) MISS_VAL ) {  
        ind = ind+1; /* The rest of the otraits in rows pp onward Otraits and Otraits*marker */
        for ( kk = 1 ; kk <= themap->otraits ; kk++ ) {
          if ( themap->onames != NULL && themap->onames[kk][0] == '+' ) {
	        for ( jj=row; jj < row + themap->otypes[kk]-1 ; jj++ ) 
	          if ( individs[ii].oyt[3] == jj-row+1 ) {
	            xm[jj][ind] = (FPN)1.0 ;/*  What follows is the Marker by otrait interaction term. */ /*
	            if (*(*((individs + ii)->markers + chrom) + mark) < -1 || *(*((individs + ii)->markers + chrom) + mark) > 1)   
	              *(*(xm+jj+roffset)+ind) = expect_mark(theparams,individs, ii, tgptr);
	            else
	              *(*(xm+jj+roffset)+ind) = (FPN) *(*( (individs+ii)->markers+chrom)+mark);*/	            
	          }
	          else 
	            xm[jj][ind] = (FPN)0.0; /* *(*(xm+jj+roffset)+ind) = 0.0; */
	        row = jj;
	      }	      
        }
    }
    }
  }  
  return(syy);
}

/*
  Code to do the combination forward-backward stepwise regression.

*/
genome *for_back_swr(individual *individs,params *theparams,markermap *themap,genome *first,linpakws *lnpk,int nbp,int explan,int steps)
{

  genome *gptr,*lgptr,*tgptr;
  int n,pp,i,k,fact,error, saverank,*jv,thestep;
  FPN max,**xm,*qv,*bv,*yv,*rv,tol,Syy,Syyp,SSep,pval,tvar;
  i=nbp;

/* 
   Forward stepwise regression.  This will rank the markers.  It will
   rank the first nbp markers, provided nbp is less than n-1 (or n/2 - 1 if an sfx or rfx).
*/

  fact = calc_fact(theparams);
  tol = (FPN)0.0;
  n = calc_samplesize(theparams,individs);
  yv = lnpk->y[1];
  rv = lnpk->rsd[1];
/*  steps = nbp;   Start out by assuming we could use all markers... 
  if ( steps >= (n-1-explan)/fact )  
    steps = (n-1-explan)/fact - 1;     ...but decrease this if that is too many. 
  if ( theparams->verbosity == 1 )
    printf("\n\n\tSetting maximum number of markers to %d\n\n",steps);*/
 

  xm = lnpk->xx;
  qv = lnpk->qraux;
  jv = lnpk->jpvt[1];
  bv = lnpk->bb[1];
  gptr = first;
  for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next )
	 lgptr->whichqtl = 0; /*  No Markers in model... */

  pp = 1 + explan ;	 

  Syyp = init_xmyv(xm,yv,gptr,gptr,n,individs,themap,theparams);
  error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
  SSep = sdot(n,rv,1,rv,1);

  for ( i = 1 ; i <= steps ; i++ ) {
	if ( theparams->verbosity == 1 ) {
	  fprintf(stdout,"\nStep %d     ",i);
	  thestep = Rotator(0);
	}
	for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next )
	  lgptr->pxo = (FPN)0.0;
	pp = 1 + fact * i +  explan ;
	for ( lgptr = gptr; lgptr != NULL ; lgptr = lgptr->next )   /* do each marker in turn */
	  if ( lgptr->whichqtl == 0 && lgptr->markr != 0 ) {
/* This allows you to iconify the program when it runs in MS Windows. */
#if defined(BORLAND)
	    WYield();
#endif
        lgptr->whichqtl = 1; /*Put it in...*/
	  	if ( theparams->verbosity == 1 )
	  	  thestep = Rotator(thestep);
		Syy = init_xmyv(xm,yv,gptr,lgptr,n,individs,themap,theparams);
		error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
		lgptr->pxo = sdot(n,rv,1,rv,1);
		lgptr->mxo = (SSep/lgptr->pxo - (FPN)1.0) * (FPN) (n-pp);
        lgptr->whichqtl = 0; /*Take it out...*/
      }
      max = (FPN)0.0;
	  for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next ) /* First, find max...*/
		if ( lgptr->markr != 0 && lgptr->whichqtl == 0 && lgptr->mxo > max && lgptr->mxo > 0.0) {
		  max = lgptr->mxo;
		  tgptr = lgptr;
		}
     /*  Is the max > srf1?  */
     
      pval = betai( (FPN)0.5*(n-pp), (FPN)0.5* (FPN) fact, (FPN) (n-pp) / (n - pp + tgptr->mxo) );
      if ( pval <= theparams->srf1 ) {
        SSep = tgptr->pxo;
	    tgptr->whichqtl = i;
	  }
	  else {    /* Can't add another variable, time to recheck all... */
	  	if ( theparams->verbosity == 1 ) {
		  fprintf(stdout,"\nRechecking all variables previously entered   ");
		  thestep = Rotator(0);
		} 
	    pp = 1 + fact * (i-1) +  explan ;
	    Syyp = init_xmyv(xm,yv,gptr,gptr,n,individs,themap,theparams);
	    error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
	    SSep = sdot(n,rv,1,rv,1);
        pp = pp-1;
	    for ( lgptr = gptr; lgptr != NULL ; lgptr = lgptr->next )   /* do each marker in turn */
		  if ( lgptr->whichqtl > 0 && lgptr->markr != 0 ) {
#if defined(BORLAND)    /* This allows you to iconify the program when it runs in MS Windows. */
	        WYield();
#endif
            saverank = lgptr->whichqtl;
            lgptr->whichqtl = 0; /*Take it out ...*/
	  	    if ( theparams->verbosity == 1 )
	  	      thestep = Rotator(thestep);
		    Syy = init_xmyv(xm,yv,gptr,lgptr,n,individs,themap,theparams);
		    error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
		    lgptr->pxo = sdot(n,rv,1,rv,1);
		    lgptr->pxo = ((Syyp-SSep)-(Syy-lgptr->pxo)) * ((FPN) (n-pp-1))/SSep;
		    tvar = (FPN) (n-pp-1);
		    tvar = tvar / (tvar + lgptr->pxo); 
		    if (tvar > 0.0 && tvar < 1.0 )
		      lgptr->pxo =  betai( (FPN)0.5*(n-pp-1), (FPN)0.5* (FPN) fact, tvar );
		    else 
		      lgptr->pxo = (FPN)1.0;
            lgptr->whichqtl = saverank; /*Put it in...*/
		  }
	  for ( lgptr = gptr; lgptr != NULL ; lgptr = lgptr->next )   /* Negate those with large pvals */
	    if ( lgptr->whichqtl > 0 && lgptr->markr != 0 && lgptr->pxo > theparams->srb1 )  
          lgptr->whichqtl = -lgptr->whichqtl;
	  i = steps+1;/* this gets us out of the loop.*/
	}
  }
  return(first);
}



/*
  Code to do the forward stepwise regression.

*/

genome *forward_swr(individual *individs,params *theparams,markermap *themap,genome *first,linpakws *lnpk,int nbp,int explan,int steps)
{
genome *gptr,*lgptr,*tgptr;
int n,pp,i,k,fact,error,*jv,thestep;
FPN max,**xm,*qv,*yv,*bv,*rv,tol,Syy,Syyp,SSep;
i=nbp;
/* 
   Forward stepwise regression.  This will rank the markers.  It will
   rank the first nbp markers, provided nbp is less than n-1 (or n/2 - 1 if an sfx or rfx).
*/

  fact = calc_fact(theparams);
  tol = (FPN)0.0;
  n = calc_samplesize(theparams,individs);
  yv = lnpk->y[1];
  rv = lnpk->rsd[1];
/*
  steps = nbp;
  if ( steps >= (n-1-explan)/fact ) {
    steps = (n-1-explan)/fact - 1;
    if ( theparams->verbosity == 1 )
      printf("\n\n\tToo many parameters in the model...Setting maximum number of markers to %d\n\n",steps);
  }
*/
  xm = lnpk->xx;
  qv = lnpk->qraux;
  jv = lnpk->jpvt[1];
  bv = lnpk->bb[1];
  gptr = first;
  for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next )
	 lgptr->whichqtl = 0; /*  No Markers in model... */

  pp = 1 + explan ;	 

  Syyp = init_xmyv(xm,yv,gptr,gptr,n,individs,themap,theparams);
  error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
  SSep = sdot(n,rv,1,rv,1);

  for ( i = 1 ; i <= steps ; i++ ) {
	 if ( theparams->verbosity == 1 ) {
		fprintf(stdout,"\nStep %d    ",i);
		thestep = Rotator(0);
     }
	 for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next )
		lgptr->pxo = (FPN)0.0;
	 pp = 1 + fact * i +   explan ;
	 for ( lgptr = gptr; lgptr != NULL ; lgptr = lgptr->next )   /* do each marker in turn */
		if ( lgptr->whichqtl == 0 && lgptr->markr != 0 ) {
/* This allows you to iconify the program when it runs in MS Windows. */
#if defined(BORLAND)
	 WYield();
#endif
          lgptr->whichqtl = 1; /*Put it in...*/
		  if ( theparams->verbosity == 1 )
		    thestep = Rotator(thestep); /* {
			 fprintf(stdout,".");
			 fflush(stdout);
		  }*/
		  Syy = init_xmyv(xm,yv,gptr,lgptr,n,individs,themap,theparams);
		  error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
		  lgptr->pxo = sdot(n,rv,1,rv,1);
		  lgptr->mxo = (SSep/lgptr->pxo - (FPN)1.0) * (FPN) (n-pp);
          lgptr->whichqtl = 0; /*Take it out...*/
		}

	 max = (FPN)0.0;
	 for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next ) /* First, find max...*/
		if ( lgptr->markr != 0 && lgptr->whichqtl == 0 && lgptr->mxo > max ) {
		  max = lgptr->mxo;
		  tgptr = lgptr;
		}
		
     SSep = tgptr->pxo;
	 tgptr->whichqtl = i;
  }
/*  pp = 1 + fact*lnbp;*/
  return(first);
}


/*
  Code to do the backward stepwise regression.
*/
genome *backward_swr(individual *individs,params *theparams,markermap *themap,genome *first,linpakws *lnpk,int nbp,int explan,int steps)
{
genome *gptr,*lgptr,*tgptr;
int n,pp,i,k,fact,error,lnbp,*jv,thestep;
FPN max,min,**xm,*qv,*bv,*yv,*rv,tol,Syy,Syyp,SSep,SSrp;
i=steps; i=nbp;
/* 
   Backward stepwise regression.  This will rank the markers.  It will
   rank all markers.   
*/

  fact = calc_fact(theparams);
  tol = (FPN)0.0;
  n = calc_samplesize(theparams,individs);
  yv = lnpk->y[1];
  rv = lnpk->rsd[1];
  lnbp = themap->ml;
  pp = 1 + fact*themap->ml + explan;
  if ( pp >= n ) {
    if ( theparams->verbosity == 1 )
      printf("\n\n\tToo many parameters in the model...Now exiting\n\n");
    exit(1);
  }
  xm = lnpk->xx;
  qv = lnpk->qraux;
  jv = lnpk->jpvt[1];
  bv = lnpk->bb[1];
  gptr = first;
  for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next ) /*Put all markers in...*/
	 lgptr->whichqtl = 1;


  for ( i = 1 ; i <= lnbp ; i++ ) {
	 if ( theparams->verbosity == 1 ) {
		fprintf(stdout,"\nStep %d   ",i);
		thestep = Rotator(0);
     }
	 for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next )
		lgptr->pxo =(FPN) 0.0;
     pp = 1 + fact * (lnbp-i+1) + explan;
     Syyp = init_xmyv(xm,yv,gptr,gptr,n,individs,themap,theparams);
     error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
     SSep = sdot(n,rv,1,rv,1);
     SSrp = Syyp - SSep;
     pp = pp-1;
	 for ( lgptr = gptr; lgptr != NULL ; lgptr = lgptr->next )   /* do each marker in turn */
		if ( lgptr->whichqtl == 1 && lgptr->markr != 0 ) {
/* This allows you to iconify the program when it runs in MS Windows. */
#if defined(BORLAND)
	 WYield();
#endif
          lgptr->whichqtl = 0; /*Take it out...*/
		  if ( theparams->verbosity == 1 ) 
		    thestep = Rotator(thestep);
		    /*{
			 fprintf(stdout,".");
			 fflush(stdout);
		  }*/
		  Syy = init_xmyv(xm,yv,gptr,lgptr,n,individs,themap,theparams);
		  error = sqrst(xm,n,n,pp,yv,tol,bv,rv,&k,jv,qv);
		  lgptr->pxo = sdot(n,rv,1,rv,1);
		  lgptr->mxo = (SSrp - (Syy-lgptr->pxo) ) * (FPN) (n-pp-1) / SSep;
		  lgptr->whichqtl = 1; /*Put it in...*/
          
		}
/*  Find the QTL with the minimal effect, that is with the minimum F statistic. */
	 max = (FPN)0.0;
	 for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next )
		if ( lgptr->markr != 0 && lgptr->whichqtl == 1 && lgptr->mxo > max ) {
		  max = lgptr->mxo;
		  tgptr = lgptr;
		}
	 min = max;
	 for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next ) /* Then, find min...*/
		if ( lgptr->markr != 0 && lgptr->whichqtl == 1 && lgptr->mxo < min ) {
		  min = lgptr->mxo;
		  tgptr = lgptr;
		}
/*     SSep = tgptr->pxo;*/
 	 tgptr->whichqtl =  i - themap->ml - 1 ;
  }

  for ( lgptr = gptr ; lgptr != NULL ; lgptr = lgptr->next )
		lgptr->whichqtl = -lgptr->whichqtl;
  return(first);
}

void write_sr_results(genome *agptr,params *theparams,char *srfile,int n,int explan)
{
  FILE *errorf;
  int jj,fact;
  genome *tptr;
  
  fact = calc_fact(theparams);
  
    /*  Write the markers to the output file....*/
    errorf = fileopen(srfile, "a");
    if ( theparams->srm == 3 )
	  fprintf(errorf,"\n#\n-Zmapqtl  analysis converted to ranks for -trait %d",theparams->whichtrait);
    else if ( theparams->srm == 2 )
	  fprintf(errorf,"\n#\n-FB Stepwise regression analysis for -trait %d",theparams->whichtrait);
    else if ( theparams->srm == 1 )
      fprintf(errorf,"\n#\n-Backward Stepwise regression analysis for -trait %d",theparams->whichtrait);
    else  
	  fprintf(errorf,"\n#\n-Forward Stepwise regression analysis for -trait %d",theparams->whichtrait);
    fprintf(errorf,"\n# %d DOF in the numerator: The given DOF is for the denominator...\n",fact);
    for ( jj = 1 ; jj<= 55 ; jj++ )
      fprintf(errorf,"-");
    fprintf(errorf, "\n");
    for ( jj = 1 ; jj<= 55 ; jj++ )
      fprintf(errorf,"-");
    if ( theparams->srm == 3 )
      fprintf(errorf,"\nChromosome    Marker      WhichQTL     LRStat      Model\n");
    else
      fprintf(errorf,"\nChromosome    Marker      Rank         F-Stat      DOF\n");
    for ( jj = 1 ; jj<= 55 ; jj++ )
      fprintf(errorf,"-");
    fprintf(errorf,"              -start");
    if ( theparams->srm == 3 )
      fprintf(errorf," ranks");
    for ( tptr = agptr ; tptr != NULL ; tptr = tptr->next )
	  if ( tptr->whichqtl > 0 ) {
	    if ( theparams->srm == 3 )
	      fprintf(errorf,"\n%10d%10d%10d%15.5f%10d",tptr->chrom, tptr->markr, tptr->whichqtl,tptr->mxo, n );
	    else
	      fprintf(errorf,"\n%10d%10d%10d%15.5f%10d",tptr->chrom, tptr->markr, tptr->whichqtl,tptr->mxo, n-(1+tptr->whichqtl*fact+explan) );
      }
    if ( theparams->srm == 3 )
      fprintf(errorf,"              -end  ranks\n");
    else
      fprintf(errorf,"              -end\n");
    for ( jj = 1 ; jj<= 55 ; jj++ )
      fprintf(errorf,"-");
    fprintf(errorf,"\n");
    for ( jj = 1 ; jj<= 55 ; jj++ )
      fprintf(errorf,"-");
    fprintf(errorf, "\n");
    fileclose(srfile, errorf);
}

/*

  Fact will tell us how many variables per marker.  If there are two 
  genotypic classes, then fact is 1.  If there are three, then fact is
  2.  Backcrosses and recombinant inbreds return a 1, while intercrosses
  return a 2.
  
*/
int calc_fact(params *theparams)
{
  int fact;

  if ( theparams->cross == 1 || theparams->cross == 2 || theparams->cross == 5)
    fact = 1;
  else if ( theparams->cross == 3 || theparams->cross == 4 )
    fact = 2;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    fact = 1;

  return(fact);
}  


int calc_samplesize(params *theparams,individual *individs)
{
  int i,n;
  
  n = 0;
  for ( i = 1 ; i <= theparams->nn ; i++ )
    if ( individs[i].y[theparams->whichtrait] > (FPN) MISS_VAL )
		n = n+1;
  return(n);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file SRfunc.c
------------------------------------------------------------------ */


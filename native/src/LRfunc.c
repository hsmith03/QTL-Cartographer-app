/* ------------------------------------------------------ XCutXCodeXSkip
     This file (LRfunc.c) is part of QTL Cartographer
         
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
/*  Simple linear regression routines.*/


void do_mlr_model(individual *iptr, params *theparams, genome *tgptr, int trait, markermap *map, FPN **beta0, FPN **beta1, FPN **fstat, FPN **lrvalues, linpakws *lnpk)
{
  int ii,jj,kk,ind,error,k,mark,chrom,row,roffset;
  FPN   SSr,SSe,rr,ff,MSe;
  FPN SSeM,SSrM,tol,ysum,ydot;
  mark = tgptr->markr;
  chrom = tgptr->chrom;
  roffset = (lnpk->p-2)/2;
  tol = (FPN) 0.0;
  ind = 0;
  for ( ii = 1 ; ii <= theparams->nn ; ii++ ) 
    if (  iptr[ii].y[trait] > (FPN) MISS_VAL ) { /* put trait in yy */
        ind = ind+1;
        lnpk->y[1][ind] =  iptr[ii].y[trait];
    
	    lnpk->xx[1][ind] = (FPN) 1.0;  /* mean in row 1 */
/* The rest of the otraits in rows 2 through pp-1 */
        row = 2;
        for ( kk = 1 ; kk <= map->otraits ; kk++ ) {
          if ( map->onames != NULL && map->onames[kk][0] == '+' ) {
	        for ( jj = row ; jj < row + map->otypes[kk]-1 ; jj++ ) 
	          if (  iptr[ii].oyt[3] == jj-row+1 ) {
	            lnpk->xx[jj][ind] = (FPN) 1.0 ;
	            if ( iptr[ii].markers[chrom][mark] < -1 ||  iptr[ii].markers[chrom][mark] > 1)   /* This is the Marker by otrait interaction term. */
	              lnpk->xx[jj+roffset][ind] = expect_mark(theparams,iptr, ii, tgptr);
	            else
	              lnpk->xx[jj+roffset][ind] = (FPN)   iptr[ii].markers[chrom][mark];	            
	          }
	          else 
	            lnpk->xx[jj][ind] = lnpk->xx[jj+roffset][ind] = (FPN) 0.0; 
	        row = jj;
	      }	      
        }
	    if ( iptr[ii].markers[chrom][mark] < -1 ||  iptr[ii].markers[chrom][mark] > 1)  /* marker in row pp */
	      lnpk->xx[lnpk->p][ind] = expect_mark(theparams,iptr, ii, tgptr);
	    else
	      lnpk->xx[lnpk->p][ind] = (FPN)   iptr[ii].markers[chrom][mark];
   }
     
/*  This part does the Marker and Marker by otrait interactions... 
    Y = mu + Marker + Sum(Factors + Marker*Factors) + e */      
  for ( ii = 1 ; ii <= lnpk->p ; ii++ )
    scopy(lnpk->n, lnpk->xx[ii], 1, lnpk->xsave[ii], 1);
  ysum = ssum(ind,lnpk->y[1],1);
  ydot = sdot(ind,lnpk->y[1],1,lnpk->y[1],1);
  error = sqrst(lnpk->xx,ind,ind,lnpk->p,lnpk->y[1],tol,lnpk->bb[1],lnpk->rsd[1],&k,lnpk->jpvt[1], lnpk->qraux);
  SSe = sdot(ind,lnpk->rsd[1],1,lnpk->rsd[1],1);
  SSr =  ydot - ysum*ysum/((FPN) ind) - SSe  ;
  MSe = SSe / (FPN) (ind-lnpk->p);
  ff = (SSr/ (FPN) (lnpk->p-1)) / MSe ;  /* F stat that at least one beta not 0 */ 

/*  This part does the Null hypothesis of no marker effect... 
    Y = mu + Sum(Factors + Marker*Factors) + e */
  for ( ii = 1 ; ii < lnpk->p ; ii++ )
    scopy(lnpk->n, lnpk->xsave[ii], 1, lnpk->xx[ii], 1);
  error = sqrst(lnpk->xx,ind,ind,lnpk->p-1,lnpk->y[1],tol,lnpk->bb[1],lnpk->rsd[1],&k,lnpk->jpvt[1],lnpk->qraux);
  SSeM = sdot(ind,lnpk->rsd[1],1,lnpk->rsd[1],1);
  SSrM =  ydot - ysum*ysum/((FPN) ind) - SSeM  ;
  rr = (FPN) 1.0;
  ff = ( (SSr-SSrM)/ rr ) / MSe;
  fstat[chrom][mark] = (FPN) ind * (FPN) log( 1.0 + ff/(FPN)(ind-lnpk->p) );
  if ( ff > (FPN) 0.0 )
    beta0[chrom][mark] =  betai( (FPN) (ind-lnpk->p)/(FPN)2.0, rr/(FPN)2.0, (FPN) (ind-lnpk->p)/((FPN)(ind-lnpk->p) + rr*ff)  );
  else 
    beta0[chrom][mark] = (FPN) 1.0;

/*  This part does the Null hypothesis of no marker by factors effect... 
    Y = mu + Marker + Sum(Factors) + e */
  for ( ii = 1 ; ii <= roffset+1 ; ii++ )
    scopy(lnpk->n, lnpk->xsave[ii], 1, lnpk->xx[ii], 1);
  scopy(lnpk->n, lnpk->xsave[roffset+2], 1, lnpk->xx[lnpk->p], 1);
  error = sqrst(lnpk->xx,ind,ind,roffset+2,lnpk->y[1],tol,lnpk->bb[1],lnpk->rsd[1],&k,lnpk->jpvt[1],lnpk->qraux);
  SSeM = sdot(ind,lnpk->rsd[1],1,lnpk->rsd[1],1);
  SSrM =  ydot - ysum*ysum/((FPN) ind) - SSeM  ;
  rr = (FPN) roffset;
  ff = ( (SSr-SSrM)/ rr ) / MSe;
  if ( ff > (FPN) 0.0 )
    beta1[chrom][mark] =  betai( (FPN) (ind-lnpk->p)/(FPN)2.0, rr/(FPN)2.0, (FPN) (ind-lnpk->p)/((FPN)(ind-lnpk->p) + rr*ff)  );
  else 
    beta1[chrom][mark] = (FPN) 1.0;

/*  This part does the Null hypothesis of no marker or marker by factors effect... 
    Y = mu + Sum(Factors) + e */
  for ( ii = 1 ; ii <= roffset+1 ; ii++ )
    scopy(lnpk->n, lnpk->xsave[ii], 1, lnpk->xx[ii], 1);
  error = sqrst(lnpk->xx,ind,ind,roffset+1,lnpk->y[1],tol,lnpk->bb[1],lnpk->rsd[1],&k,lnpk->jpvt[1],lnpk->qraux);
  SSeM = sdot(ind,lnpk->rsd[1],1,lnpk->rsd[1],1);
  SSrM =  ydot - ysum*ysum/((FPN) ind) - SSeM  ;
  rr = (FPN) (roffset+1);
  ff = ( (SSr-SSrM)/ rr ) / MSe;
  if ( ff > (FPN) 0.0 )
    lrvalues[chrom][mark] = betai( (FPN) (ind-lnpk->p)/(FPN)2.0, rr/(FPN)2.0, (FPN) (ind-lnpk->p)/((FPN)(ind-lnpk->p) + rr*ff)  );
  else 
    lrvalues[chrom][mark] = (FPN) 1.0;




/*  This part does the  hypothesis of no marker by factors effect... 
    Y = mu + Marker + Sum(Factors) + e */
  for ( ii = 1 ; ii <= roffset+1 ; ii++ )
    scopy(lnpk->n, lnpk->xsave[ii], 1, lnpk->xx[ii], 1);
  scopy(lnpk->n, lnpk->xsave[lnpk->p], 1, lnpk->xx[roffset+2], 1);
  error = sqrst(lnpk->xx,ind,ind,roffset+2,lnpk->y[1],tol,lnpk->bb[1],lnpk->rsd[1],&k,lnpk->jpvt[1],lnpk->qraux);
  SSe = sdot(ind,lnpk->rsd[1],1,lnpk->rsd[1],1);
  SSr =  ydot - ysum*ysum/((FPN) ind) - SSe  ;
  MSe = SSe / (FPN) (ind-roffset-2);

/*  This part does the Null hypothesis of no marker or marker by factors effect... 
    Y = mu + Sum(Factors) + e */
  for ( ii = 1 ; ii <= roffset+1 ; ii++ )
    scopy(lnpk->n, lnpk->xsave[ii], 1, lnpk->xx[ii], 1);
  error = sqrst(lnpk->xx,ind,ind,roffset+1,lnpk->y[1],tol,lnpk->bb[1],lnpk->rsd[1],&k,lnpk->jpvt[1],lnpk->qraux);
  SSeM = sdot(ind,lnpk->rsd[1],1,lnpk->rsd[1],1);
  SSrM =  ydot - ysum*ysum/((FPN) ind) - SSeM  ;
  rr = (FPN) 1.0;
  ff = ( (SSr-SSrM)/ rr ) / MSe;

/*
We compare the two models:

    Y = mu + Marker + Sum(Factors) + e  
    Y = mu + Sum(Factors) + e 

This makes sense if there is no interaction between the Marker and the Factors.
    
*/

  fstat[chrom][mark] = (FPN) ind * (FPN) log( 1.0 + ff/(FPN)(ind-roffset-2) );
  if ( ff > (FPN) 0.0 )
    beta0[chrom][mark] =  betai( (FPN) (ind-roffset-2)/(FPN)2.0, rr/(FPN)2.0, (FPN) (ind-roffset-2)/((FPN)(ind-roffset-2) + rr*ff)  );
  else 
    beta0[chrom][mark] = (FPN) 1.0;


}


/*
  Calculate the linear regression statistics for a single marker, but 
  add in any other traits that are indicated.
*/
int calc_mlrstats(params *theparams, markermap *themap, individual *individs, int trait, FPN **beta0, FPN **beta1, FPN **fstat, FPN **lrvalues, genome *gptr, linpakws *lnpk)
{
  int  error;
  genome *tgptr;
  error = 0;
  if (theparams->nn <= 0)
    return (-1);
 
  for ( tgptr = gptr;  tgptr != NULL ; tgptr = tgptr->next ) {
    if ( tgptr->markr > 0 )
      do_mlr_model(individs,theparams,tgptr,trait,themap,beta0,beta1,fstat,lrvalues,lnpk);  
  }

   return (error);
}

int calc_lrstats(params *theparams, markermap *themap, individual *individs, int trait, FPN **beta0, FPN **beta1, FPN **fstat, genome *gptr)
{
  int ii, jj, kk, error,  missingdata;
  FPN xbar, xsqr, xiyi, ysqr, ysum, xij, Sxx, Sxy, Syy, SSe, SSr, xsum;
  genome *tgptr;
  ii=themap->traits;
  error = 0;
  if (theparams->nn <= 0)
    return (-1);

  for ( tgptr = gptr;  tgptr != NULL ; tgptr = tgptr->next ) {
    ii = tgptr->chrom;
    jj = tgptr->markr;
    if ( tgptr->markr > 0 )   {
	  missingdata = 0;
	  ysum = ysqr = (FPN) 0.0;
	  xsum = (FPN) 0.0;
	  xsqr = (FPN) 0.0;
	  xiyi = (FPN) 0.0;
	  ysqr = (FPN) 0.0;
	  for (kk = 1; kk <= theparams->nn; kk++)
	    if ( individs[kk].y[trait] > (FPN) MISS_VAL) {
	      if (individs[kk].markers[ii][jj] < -1 || individs[kk].markers[ii][jj] > 1) 
	        xij = expect_mark(theparams,individs, kk, tgptr);
	      else
	        xij = (FPN) individs[kk].markers[ii][jj];
/*
	      if ( theparams->cross >= 2)
	        xij = xij + 1.0;
*/
	      ysum = ysum + individs[kk].y[trait];
	      ysqr = ysqr + individs[kk].y[trait] * individs[kk].y[trait];

	      xsum = xsum + xij;
	      xsqr = xsqr + xij * xij;
	      xiyi = xiyi + xij * individs[kk].y[trait];
	    }
	    else
	      missingdata = missingdata + 1;

	  Syy = ysqr - ysum * ysum / (FPN) (theparams->nn - missingdata);

	  xbar = xsum / (FPN) theparams->nn;
	  Sxx = xsqr - xsum * xsum / (FPN) (theparams->nn - missingdata);
	  Sxy = xiyi - xsum * ysum / (FPN) (theparams->nn - missingdata);

	  if (Sxx != (FPN) 0.0)
	    beta1[ii][jj] = Sxy / Sxx;
	  else {
	    beta1[ii][jj] = (FPN) 0.0;
	    error = error + 1;
	  }
	  beta0[ii][jj] = ysum / (FPN) (theparams->nn - missingdata) - beta1[ii][jj] * xbar;
	  SSe = Syy - beta1[ii][jj] * Sxy;
	  SSr = beta1[ii][jj] * Sxy;
	  if (SSe != (FPN) 0.0)
	    fstat[ii][jj] = (FPN) (theparams->nn - missingdata - 2) * SSr / SSe;
	  else {
	    fstat[ii][jj] = (FPN) 0.0;
	    error = error + 1;
	  }
    }	
  }
  return (error);
}



void print_lrstats(int nn, int trait, markermap *themap, char *minfile, char *iinfile, char *outfile, FPN **beta0, FPN **beta1, FPN **fstat, int cross)
{
  FILE *outf;
  int ii, jj;
  FPN lr, halfnn, prf;
  ii=cross;
  lr = (FPN) 0.0;
  halfnn = (FPN) (nn - 2) / (FPN)2.0;
  if (*outfile == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  fprintf(outf, "\n\n\tThis output is based on the map in (%s)\n", minfile);
  fprintf(outf, "\tAnd the data in (%s)\n\n", iinfile);
  fprintf(outf, "\n\t\tSample Size............%10d\n\n", nn);
  fprintf(outf, "\nThis analysis fits the data to the simple linear regression model");
  fprintf(outf, "\n\t\t y = b0 + b1 x + e");
  fprintf(outf, "\n The results below give the estimates for b0, b1 and the F statistic");
  fprintf(outf, "\n for each marker.  ");
  fprintf(outf, "\n\n We are interested in whether the marker is linked to a QTL.   We test");
  fprintf(outf, "\n this idea by determining if b1 is significantly different from zero.   The F ");
  fprintf(outf, "\n statistic compares the hypothesis H0: b1 = 0 to an alternative H1: b1 not 0.   ");
  fprintf(outf, "\n The pr(F) is a measure of how much support there is for H0.   A smaller pr(F)  ");
  fprintf(outf, "\n indicates less support for H0 and thus more support for H1.   Significance at");
  fprintf(outf, "\n the 5%%, 1%%, 0.1%% and 0.01%% levels are indicated by *, **, *** and");
  fprintf(outf, "\n ****, respectively.   \n");
  fprintf(outf, "\n\nNote that our Likelihood ratio test statistic compares two nested hypotheses");
  fprintf(outf, "\nand is two times the negative natural log of the ratio of the likelihoods.  For example,");
  fprintf(outf, "\nassume that  hypothesis H0 is nested within H1 and that they have likelihoods L0 and L1 respectively.");
  fprintf(outf, "\nThen, the \"Likelihood Ratio Test Statistic\" is -2ln(L0/L1). \n#");
  if ( themap->tnames != NULL )
    fprintf(outf, "\tThis trait is: %s, and ",themap->tnames[trait]);
  fprintf(outf, "\n-t %4d is the number of trait being analyzed.", trait);
  putline(outf,'-',78);
  putline(outf,'-',78);
  fprintf(outf, "\n\n   Chrom.  Marker     b0         b1      -2ln(L0/L1)  F(1,n-2)      pr(F)");
  putline(outf,'-',78);
  for (ii = 1; ii <= themap->m; ii++) {
    for (jj = 1; jj <= themap->mpc[ii]; jj++) {
      lr = (FPN) (nn - 2) / ((FPN) (nn - 2) + fstat[ii][jj]);
      if ( lr < (FPN) 0.0 || lr > (FPN) 1.0 )
        prf = -(FPN) 1.0;
      else
        prf = betai(halfnn, (FPN)0.5, lr);
      lr = nn * (FPN) log(1.0 + fstat[ii][jj] / (FPN) (nn - 2));
      fprintf(outf, "\n%6d  %6d  %9.3f  %9.3f  %11.3f %11.3f  %11.3f", ii, jj, beta0[ii][jj], beta1[ii][jj], lr, fstat[ii][jj], prf);
      fstat[ii][jj] = lr;
      if (prf < (FPN) 0.05)
	    fprintf(outf, " *");
      if (prf < (FPN) 0.01)
	    fprintf(outf, "*");
      if (prf < (FPN) 0.001)
	    fprintf(outf, "*");
      if (prf < (FPN) 0.0001)
	    fprintf(outf, "*");
    }

  }
  putline(outf,'-',78);
  putline(outf,'-',78);
  fprintf(outf,"\n");
  fileclose(outfile, outf);
}

void print_mlrstats(params *theparams, int trait, markermap *themap, FPN **beta0, FPN **beta1, FPN **fstat, FPN **lrvalue)
{
  FILE *outf;
  int   ii, jj ;
  FPN   prf;

 
  outf = fileopen(theparams->lrfile, "a");
  if (outf == NULL)
    outf = fileopen(theparams->lrfile, "w");
  
  fprintf(outf, "\n\n\tThis output is based on the map in (%s)\n", theparams->map);
  fprintf(outf, "\tAnd the data in (%s)\n\n", theparams->ifile );
  fprintf(outf, "\n\t\tSample Size............%10d\n\n", theparams->nn);
  fprintf(outf, "\nThis analysis fits the data to the simple linear regression model");
  fprintf(outf, "\n\t\t Trait = Mean + Marker ");
  for ( ii = 1 ; ii <= themap->otraits ; ii++ ) 
    if ( themap->onames != NULL && themap->onames[ii][0] == '+' ) {
      fprintf(outf, "%s + ",themap->onames[ii]);
      themap->onames[ii][0] = '*';
      fprintf(outf,"Marker%s ",themap->onames[ii]);
      themap->onames[ii][0] = '+';      
    }  
  fprintf(outf, "+ Error");
  fprintf(outf, "\nVersus the null model ");
  fprintf(outf, "\n\t\t Trait = Mean ");
  for ( ii = 1 ; ii <= themap->otraits ; ii++ ) 
    if ( themap->onames != NULL && themap->onames[ii][0] == '+' ) 
      fprintf(outf, "%s ",themap->onames[ii]);
  fprintf(outf, "+ Error");
  fprintf(outf, "\nThe results in the last two columns below give the Likelihood ratio (LR)");
  fprintf(outf, "\ntest statsitics and the p values for marker assuming no interactions. ");
  fprintf(outf, "\nTo clarify, let T be a trait, M a marker and F the explanatory factors.  ");
  fprintf(outf, "\n M*F are the marker-factor interactions and pr() is a p value comparing two models.   ");
  fprintf(outf, "\n\nNote that our Likelihood ratio test statistic compares two nested hypotheses");
  fprintf(outf, "\nand is two times the negative natural log of the ratio of the likelihoods.  For example,");
  fprintf(outf, "\nassume that  hypothesis H0 is nested within H1 and that they have likelihoods L0 and L1 respectively.");
  fprintf(outf, "\nThen, the \"Likelihood Ratio Test Statistic\" is -2ln(L0/L1). \n#");
  putline(outf,'-',65);
  putline(outf,'-',65);
  fprintf(outf, "\n p-value       Model 1             Model 0       Test of   ");
  putline(outf,'-',65);
  fprintf(outf, "\n pr(M+M*F)   T = M + F + M*F     T = F          M and M*F");
  fprintf(outf, "\n pr(M*F)     T = M + F + M*F     T = M + F      M*F  ");
  fprintf(outf, "\n pr(M)       T = M + F           T = F          M assuming M*F = 0");
  putline(outf,'-',65);
  fprintf(outf, "\nFor for the final p-value,   significance at the 5%%, 1%%, 0.1%% and 0.01%% levels ");
  fprintf(outf, "\nare indicated by *, **, *** and ****, respectively. \n");
  if ( themap->tnames != NULL )
    fprintf(outf, "\tThis trait is: %s, and ",themap->tnames[trait]);
  fprintf(outf, "\n-t %4d is the number of trait being analyzed.", trait);
  putline(outf,'-',78);
  putline(outf,'-',78);
  fprintf(outf, "\n\n   Chrom.  Marker  pr(M+M*F)   pr(M*F)  -2log(L0/L1)    pr(M)");
  putline(outf,'-',78);
  for (ii = 1; ii <= themap->m; ii++) {
    for (jj = 1; jj <= themap->mpc[ii]; jj++) {
      fprintf(outf, "\n%6d  %6d    %9.6f  %9.6f  %11.3f  %9.6f  ", ii, jj, lrvalue[ii][jj], beta1[ii][jj], fstat[ii][jj],beta0[ii][jj]);
      prf = beta0[ii][jj];
      if (prf < (FPN) 0.05)
	    fprintf(outf, " *");
      if (prf < (FPN) 0.01)
	    fprintf(outf, "*");
      if (prf < (FPN) 0.001)
	    fprintf(outf, "*");
      if (prf < (FPN) 0.0001)
	    fprintf(outf, "*");
    }

  }
  putline(outf,'-',78);
  putline(outf,'-',78);
  fprintf(outf, "\n\n");
  fileclose(theparams->lrfile, outf);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file LRfunc.c
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (BTL_LRfunc.c) is part of QTL Cartographer
         
    		Copyright (C) 1999-2005
	Lauren McIntyre and Jun Wu

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




/*  
    Additional copyright belongs to Lauren McIntyre and Jun Wu 
    of Biological Sciences and Statistics at Purdue University


          BTL_LRfunc.c,  basic functions for Binary LRmapqtl.          

       Created by Jun Wu from QTL Cartographer source code.
*/

#include "Main.h"

   

/*
   Translate for BTL backcross, B1 or B2. not for F cross.
   Translate the data read in from the Rcross.out formatted data
   file into a format expected by the rest of the programs.  
   
   GT     Rcross.out Read in as  Translated to
   -------------------------------------------
   AA        2          3             0
   Aa        1          1             1
   aa        0          0             0
   -------------------------------------------   
   
   Note:  the old style had Aa encoded by 1 or 2.  If the data file is in
   the old style, 2's read in by get_the_datapoints are not changed to 3.
   New data files are read in as the above.
   
   Originally, a 1 and a 2 stood for aA and Aa heterozygotes, respectively.
   In the simulation programs, they still do.  Some people thought that
   genotypes would be better encoded via the last column. That's why we have this
   byzantine recoding of data points.
*/
/* translate for B1 and B2 back cross only */
void BTL_trans_data(int **dataptr,markermap *themap)
{
  int ii, jj, uc;
  uc = themap->l;
  for (ii = 1; ii <= themap->m; ii++) {
    if (themap->sigl > (FPN) 0.0)
      uc = themap->mpc[ii];
    for (jj = 1; jj <= uc; jj++)
      switch (dataptr[ii][jj]) {
       case 0:   dataptr[ii][jj] =  0; break;
       case 1:   dataptr[ii][jj] =  1;  break;
       case 2:   dataptr[ii][jj] =  0;  break;
       case 3:   dataptr[ii][jj] =  0;  break;
       default:  dataptr[ii][jj] = -1;  break;
    }
  }
}




int BTL_calc_lrstats(params *theparams,   individual *individs, int trait, FPN **beta0, FPN **beta1, FPN **fstat, genome *gptr)
{
  int ii, jj, kk, error, missingdata;
  FPN xbar, xsqr, xiyi, ysqr, ysum, xij, Sxx, Sxy, Syy, SSe, SSr, xsum;
  genome *tgptr;
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
	  for (kk = 1; kk <= theparams->nn; kk++){
		  if ( individs[kk].y[trait]==1 || individs[kk].y[trait]==0){ 
	      if (individs[kk].markers[ii][jj] ==0 || individs[kk].markers[ii][jj] ==1) 
	     
	        xij = (FPN) individs[kk].markers[ii][jj];
	        ysum = ysum + individs[kk].y[trait];
	        ysqr = ysqr + individs[kk].y[trait] * individs[kk].y[trait];

	        xsum = xsum + xij;
	        xsqr = xsqr + xij * xij;
	        xiyi = xiyi + xij * individs[kk].y[trait];
	    }
	    else
	      missingdata = missingdata + 1;
	  }
	  Syy = ysqr - ysum * ysum / (FPN) (theparams->nn - missingdata);

	  xbar = xsum / (FPN) theparams->nn;
	  Sxx = xsqr - xsum * xsum / (FPN) (theparams->nn - missingdata);
	  Sxy = xiyi - xsum * ysum / (FPN) (theparams->nn - missingdata);
	  printf("Sxx is %6.3f\n",Sxx);
	  printf("Sxy is %6.3f\n",Sxy);

	  if (Sxx != (FPN) 0.0)
	    beta1[ii][jj] = Sxy / Sxx;
	  else {
	    beta1[ii][jj] = (FPN) 0.0;
	    error = error + 1;
		printf("error in , beta1 is 0.....\n");
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

/*Gaussian elimination with partial pivoting for Ax=b*/

FPN Abs(FPN value)
{
	if(value<0)
		value=-value;
	return value;
}
int Gauss_elim(FPN **Matrix,  FPN *Ycolumn,int size,FPN *Betta)
{
	int i,j,k;
	int imax=0;
	FPN max=0;
	FPN buffer=(FPN)1.0*(FPN)exp(-10);
    FPN temp;
    int n=size;
	FPN m;
	FPN r,sum;
	/*initialize Betta[] to zero*/
	for(i=1;i<=n;i++)
	{
		Betta[i]=0;
	}

	for(k=1;k<n;k++)
	{
	  /*Find Maximum element in column*/
      for(i=k;i<=n;i++)
	  {
		  if(Abs(Matrix[i][k])>max)
		  {
			  imax=i;
	          max=Abs(Matrix[i][k]);
		  }

	  }

	  /*Test for zero pivot (singularity). buffer is a small machine dependent number*/
	  if(Abs(Matrix[imax][k])<buffer)
	  {	
		  printf("the matrix is singularity, no unique solution\n");
		  return 0;
	  }

	  /*interchange rows k and imax*/
	  for(j=1;j<=n;j++)
	  {
		  temp=Matrix[k][j];
		  Matrix[k][j]=Matrix[imax][j];
		  Matrix[imax][j]=temp;
	  }
	  /*interchange rows k and imax of ycolumn*/
	  
		  temp=Ycolumn[k];
	      Ycolumn[k]=Ycolumn[imax];
		  Ycolumn[imax]=temp;
	  
	 
	  /* loop over rows*/
	  for(i=k+1;i<=n;i++)
	  {
          /*compute multiplier, save in k-th column*/
          Matrix[i][k]=Matrix[i][k]/Matrix[k][k];
		  m=Matrix[i][k];
		  /*loop over element of row i*/
		  for(j=k+1;j<=n;j++)
		  {
			  Matrix[i][j]=Matrix[i][j]-m*Matrix[k][j];
			  /*Do right side of Row I*/
			  		  
		  }
		  Ycolumn[i]=Ycolumn[i]-m*Ycolumn[k];

	  }
     
	
	}/* end of k-loop*/
 
	/*Back Substitution*/
	/*Loop over Rows (Reverse Order)*/
	for(i=n;i>=1;i--)
	{
		sum=0;
       for(j=i+1;j<=n;j++)
	   {
		   sum+=Matrix[i][j]*Betta[j];
	   }
	   r=Ycolumn[i]-sum;
	   /*Test for Zero Pivot*/
	   if(Abs(Matrix[i][i])>buffer)
		   Betta[i]=r/Matrix[i][i];
	   else if(Abs(r)<=buffer)
		   Betta[i]=0;
	   else
	   {
		   printf("inconsistent system, Matrix is singular\n");
		   return 0;
	   }
	}
	
	
	
	return 1;
}

			   
/*single marker, F2 cross.*/

int BTL_calc_F2(params *theparams,  individual *individs, int trait, FPN **beta0, FPN **beta1, FPN **beta2,FPN **fstat,FPN **fstatx1,FPN **fstatx2, genome *gptr)
{
  FPN **MatrixF2;
  FPN *Ycol;
  FPN *Betta;
  int ii, jj, kk, error, missingdata;
 

  FPN x1bar, x1sqr, x1iyi, ysqr, ysum, x1ij,SSe, SSmodel,x1sum;
  FPN x2bar, x2sqr, x2iyi, x2sum,x2ij,x1x2,ybar;
  FPN n=1;
  FPN RSSE1,RSSE2;
  int nn=theparams->nn;
  genome *tgptr;
  error = 0;
  if (theparams->nn <= 0)
    return (-1);

  /*initialize*/
  MatrixF2 = dmatrix(1, 3 ,1, 3);
  Ycol = dvector(1,3);
  Betta = dvector(1,3);

  for(ii=1;ii<=3;ii++)
      {
		  
		  for(jj=1;jj<=3;jj++)
		  {
			 MatrixF2[ii][jj]=(FPN) 0.0;
		  }

		  Ycol[ii]=(FPN) 0.0;
	  }

	
  for ( tgptr = gptr;  tgptr != NULL ; tgptr = tgptr->next )
  {
    ii = tgptr->chrom;
    jj = tgptr->markr;
    if ( tgptr->markr > 0 )  
	{
	  missingdata = 0;
	  ysum = ysqr = (FPN) 0.0;
	  x1sum = (FPN) 0.0;
	  x2sum = (FPN) 0.0;
	  x1sqr = (FPN) 0.0;
	  x2sqr = (FPN) 0.0;
	  x1iyi = (FPN) 0.0;
	  x2iyi = (FPN) 0.0;
	  x1x2  = (FPN) 0.0;
	  ysqr  = (FPN) 0.0;

	  /* M1/M2(Aa) x1=1, other =0;  M2/M2(AA) x2=1 other =0;  define*/
	  /*original data AA=2, Aa=1,aa=0 . */
	  for (kk = 1; kk <= nn; kk++)
	  {		  
	      if (individs[kk].markers[ii][jj] ==0 || individs[kk].markers[ii][jj] ==1 || individs[kk].markers[ii][jj] ==2 || individs[kk].markers[ii][jj]==3) 
		  {

			  
            if(individs[kk].markers[ii][jj]==2 || individs[kk].markers[ii][jj]==3)
			{
			
				x1ij= (FPN) 0.0;
				x2ij= (FPN) 1.0;
			}
			else if(individs[kk].markers[ii][jj]==1)
			{
	            x1ij = (FPN) 1.0;
                x2ij = (FPN) 0.0;
			}
			else
			{
				x1ij= (FPN) 0.0;
				x2ij= (FPN) 0.0;
			}


	      ysum = ysum + individs[kk].y[trait];
	      ysqr = ysqr + individs[kk].y[trait] * individs[kk].y[trait];

	      x1sum = x1sum + x1ij;
		  x2sum = x2sum + x2ij;

	      x1sqr = x1sqr + x1ij * x1ij;
		  x2sqr = x2sqr + x2ij * x2ij;

	      x1iyi = x1iyi + x1ij * individs[kk].y[trait];
	      x2iyi = x2iyi + x2ij * individs[kk].y[trait];
         
		  x1x2=   x1x2+ x1ij*x2ij;

		  }
	    else
	      missingdata = missingdata + 1;
	  }

	  x1bar = x1sum / (FPN) nn;
	  x2bar = x2sum / (FPN) nn;
      ybar  = ysum  / (FPN) nn;

      MatrixF2[1][1]=(FPN) 1.0;
	  MatrixF2[1][2]=x1bar;
	  MatrixF2[1][3]=x2bar;
	  MatrixF2[2][1]=x1sum;
	  MatrixF2[2][2]=x1sqr;
	  MatrixF2[2][3]=x1x2;
	  MatrixF2[3][1]=x2sum;
	  MatrixF2[3][2]=x1x2;
	  MatrixF2[3][3]=x2sqr;
 
	  Ycol[1]= ybar;
	  Ycol[2]= x1iyi;
      Ycol[3]= x2iyi;

	 

	  /*calculate the Betta coefficient*/
	  
	  Gauss_elim(MatrixF2, Ycol,3,Betta);
      
      beta0[ii][jj] = Betta[1];
	  beta1[ii][jj] = Betta[2];
	  beta2[ii][jj] = Betta[3];

      
      
      SSmodel = Betta[1]*ysum+Betta[2]*x1iyi+Betta[3]*x2iyi-nn*ybar*ybar;
	  SSe     = ysqr-Betta[1]*ysum-Betta[2]*x1iyi-Betta[3]*x2iyi;

      /*calculate the SSE reduce model*/
      MatrixF2[1][1]=(FPN) 1.0;
	  MatrixF2[1][2]=x1bar;
	  MatrixF2[2][1]=x1sum;
	  MatrixF2[2][2]=x1sqr;

	  Ycol[1]=ybar;
	  Ycol[2]=x1iyi;
      
	  Gauss_elim(MatrixF2, Ycol, 2, Betta);
	  RSSE1 = ysqr - Betta[1]*ysum-Betta[2]*x1iyi;

	  MatrixF2[1][1]=(FPN) 1.0;
	  MatrixF2[1][2]=x2bar;
	  MatrixF2[2][1]=x2sum;
	  MatrixF2[2][2]=x2sqr;

	  Ycol[1]=ybar;
	  Ycol[2]=x2iyi;
      
	  Gauss_elim(MatrixF2, Ycol, 2, Betta);
	  RSSE2 = ysqr - Betta[1]*ysum-Betta[2]*x2iyi;



	  
	  /*calculate F value*/
	 
	  if (SSe != (FPN) 0.0)
	  {
	    fstat[ii][jj] = (FPN) ((SSmodel/(3-1))/(SSe/(nn -3)));
		fstatx1[ii][jj] = (FPN) ((RSSE2-SSe)/(SSe/(nn -3)));
		fstatx2[ii][jj] = (FPN) ((RSSE1-SSe)/(SSe/(nn -3)));
	 }
	  else 
	  {
	    fstat[ii][jj] = (FPN) 0.0;
	    error = error + 1;
	  }
    }	
  }
  free_dmatrix(MatrixF2, 1, 3, 1, 3);
  free_dvector(Ycol,1,3);
  free_dvector(Betta,1,3);
  return (error);
}


void BTL_print_lrstats(int nn, int trait, markermap *themap, char *minfile, char *iinfile, char *outfile)
{
  FILE *outf;
  
  FPN lr, halfnn;
  lr = (FPN) 0.0;
  halfnn = (FPN) (nn - 2) / (FPN)2.0;
  if (*outfile == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  fprintf(outf, "\n\n\tThis output is result of binary trait locus computation\n");
  fprintf(outf, "\n\n\tThis output is based on the map in (%s)\n", minfile);
  fprintf(outf, "\tAnd the data in (%s)\n\n", iinfile);
  fprintf(outf, "\n\t\tSample Size............%10d\n\n", nn);
  fprintf(outf, "\nThis analysis fits the data to the simple linear regression model");
  fprintf(outf, "\n\t\t y = b0 + b1 x + e");
  fprintf(outf, "\n The results below give the estimates for b0, b1 and the F statistic");
  fprintf(outf, "\n for each marker. The F statistic is for the hypothesis that the marker");
  fprintf(outf, "\n is unlinked to the quantitative trait. The column headed by PR is the");
  fprintf(outf, "\n probability that the trait is unlinked to the marker.  Significance at");
  fprintf(outf, "\n the 5%%, 1%%, 0.1%% and 0.01%% levels are indicated by *, **, *** and");
  fprintf(outf, "\n ****, respectively.  LR is -2log(L0/L1).\n");
  if ( themap->tnames != NULL )
    fprintf(outf, "\tThis trait is: %s, and ",themap->tnames[trait]);
  fprintf(outf, "\n-t %4d is the number of trait being analyzed.", trait);
 
  fprintf(outf,"\n");
  if (*outfile != '-')
    fileclose(outfile, outf);
}

/*print data calculated by specific p1 and p2 value in B1 or B2 cross*/
void BTL_print_S(int nn,markermap *themap, FILE *outf, FPN **beta0, FPN **beta1,FPN **beta2,FPN **Rmq, FPN p1,FPN p2,FPN **p3, FPN **fstat, FPN **fstatx1, FPN **fstatx2, int cross)
{
  int ii, jj, lastmark;
  FPN lr,lr1,lrb1,lrb2, halfnn, halfnn1,prf,prf1,prfb1,prfb2;
  lr = (FPN) 0.0;
  halfnn = (FPN) (nn - 3) / (FPN)2.0;
  halfnn1 = (FPN) (nn -2) /(FPN)2.0;

  if(cross==1)
  {
	  putline(outf,'-',84);
      putline(outf,'-',84);
	  fprintf(outf, "\n\n   Chrom.  Marker   b0      b1       p1       p3      Rmq       LR   F(1,n-2)  pr(F)");
      putline(outf,'-',84);
  }
  else if(cross==2)
  {
	  putline(outf,'-',84);
      putline(outf,'-',84);
	  fprintf(outf, "\n\n   Chrom.  Marker   b0      b1       p2       p3      Rmq       LR   F(1,n-2)  pr(F)");
      putline(outf,'-',84); 
  }
  else if(cross==3)
  {
	  putline(outf,'-',135);
      putline(outf,'-',135);
	  fprintf(outf, "\n\n   Chrom.  Marker   b0      b1       b2       p1       p2       p3      Rmq    F(2,n-3)  pr(F)  Fb1(1,n-3)  pr(Fb1)  Fb2(1,n-3)  pr(Fb2)");
      putline(outf,'-',135);
  }
  for (ii = 1; ii <= themap->m; ii++) {
    if (themap->sigl > (FPN) 0.0)
      lastmark = themap->mpc[ii];
    else
      lastmark = themap->l;
    for (jj = 1; jj <= lastmark; jj++) {
       lr1=(FPN) (nn - 2) / ((FPN) (nn - 2) + fstat[ii][jj]);
      lr = (FPN) (nn - 3) / ((FPN) (nn - 3) + 2*fstat[ii][jj]);
      lrb1 = (FPN) (nn - 3) / ((FPN) (nn - 3) + fstatx1[ii][jj]);
      lrb2 = (FPN) (nn - 3) / ((FPN) (nn - 3) + fstatx2[ii][jj]);
      
	  if ( lr < (FPN) 0.0 || lr > (FPN) 1.0 )
        prf = -(FPN) 1.0;
      else
	  {
        prf = betai(halfnn, 1, lr);
		prf1 = betai(halfnn1,(FPN)0.5,lr1);
		prfb1= betai(halfnn, (FPN)0.5,lrb1);
		prfb2= betai(halfnn, (FPN)0.5,lrb2);
      }
     if(cross==1)
	  {
		  fprintf(outf, "\n%6d  %6d  %6.4f  %6.4f   %6.3f   %6.3f   %6.3f   %6.3f   %6.4f  %6.3f", ii, jj, beta0[ii][jj], beta1[ii][jj], p1,p3[ii][jj],Rmq[ii][jj], lr1, fstat[ii][jj], prf1);
      }
	  else if(cross==2)
	  {
		  fprintf(outf, "\n%6d  %6d  %6.4f  %6.4f   %6.3f   %6.3f   %6.3f   %6.3f   %6.4f  %6.3f", ii, jj, beta0[ii][jj], beta1[ii][jj], p2,p3[ii][jj],Rmq[ii][jj], lr1, fstat[ii][jj], prf1);
	  }
      else if(cross==3)
	  {
		  fprintf(outf, "\n%6d  %6d  %6.4f  %6.4f   %6.4f   %6.3f   %6.3f   %6.3f   %6.3f     %6.4f  %6.4f   %6.4f    %6.4f    %6.4f    %6.4f", ii, jj, beta0[ii][jj], beta1[ii][jj], beta2[ii][jj],p1,p2,p3[ii][jj],Rmq[ii][jj], fstat[ii][jj], prf, fstatx1[ii][jj],prfb1, fstatx2[ii][jj],prfb2);
    
	  }
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
  if(cross==1)
  {
      putline(outf,'-',84);
      putline(outf,'-',84);
  }
  else if(cross==2)
  {
      putline(outf,'-',84);
      putline(outf,'-',84);
  }
  else if(cross==3)
  {
      putline(outf,'-',135);
      putline(outf,'-',135);
  }
  fprintf(outf,"\n");

}




int calc_Btl(int nn, char *outfile,markermap *themap,int cross,FPN **beta0,FPN **beta1,FPN **beta2,FPN **Rmq,FPN **p3,params *theparams, FPN **fstat, FPN **fstatx1, FPN **fstatx2)
{
	int ii,jj,lastmark;
	FPN kk,ll;
	int count=0;
	int numb=0;
	

	FILE *outf;
 
    FPN lr,lr1,lrb1,lrb2, halfnn,halfnn1, prf,prf1,prfb1,prfb2;
    lr = (FPN) 0.0;
    halfnn = (FPN) (nn - 3) / (FPN)2.0;
	halfnn1= (FPN) (nn - 2) / (FPN)2.0;
    if (*outfile == '-')
		 outf = stdout;
    else {
         outf = fileopen(outfile, "a");
         if (outf == NULL)
            outf = fileopen(outfile, "w");
	}
	if(theparams->Btl_mode==1)
		{
		/*print out the name of parameters	*/
		    if(cross==1 || cross==2)
			{
	    		putline(outf,'-',56);
				putline(outf,'-',56);
				fprintf(outf, "\n\n   Chrom.  Marker   b0      b1      LR   F(1,n-2)  pr(F)");
				putline(outf,'-',56);
			}
			else if(cross==3)
			{
				putline(outf,'-',100);
				putline(outf,'-',100);
				fprintf(outf, "\n\n   Chrom.  Marker   b0      b1      b2     F(2,n-3)  pr(F)  Fb1(1,n-3)  pr(Fb1)  Fb2(1,n-3)  pr(Fb2)");
				putline(outf,'-',100);
			}
		}
	
	for (ii = 1; ii <= themap->m; ii++) 
	{
		if (themap->sigl > (FPN) 0.0)
		    lastmark = themap->mpc[ii];
		else
			lastmark = themap->l;
		


        for (jj = 1; jj <= lastmark; jj++) 
		{
	
			if(theparams->Btl_mode==0 )
			{
				 if(cross==1)
				 {
					p3[ii][jj]=beta1[ii][jj]+2*beta0[ii][jj]-theparams->p1;	
				    Rmq[ii][jj]=(theparams->p1-beta0[ii][jj])/(2*theparams->p1-2*beta0[ii][jj]-beta1[ii][jj]);
				   
				 }
				 else if(cross==2)
				 {
					p3[ii][jj]=beta1[ii][jj]+2*beta0[ii][jj]-theparams->p2;
    				Rmq[ii][jj]=(theparams->p2-beta0[ii][jj])/(2*theparams->p2-2*beta0[ii][jj]-beta1[ii][jj]);

				 }
				 else if(cross==3)
				 {
                    p3[ii][jj]  = (FPN)2*beta0[ii][jj]+beta1[ii][jj]+(FPN)0.5*beta2[ii][jj]-(FPN)0.5*(theparams->p1+theparams->p2);
                    Rmq[ii][jj] = (FPN)0.5*(theparams->p1-theparams->p2+beta2[ii][jj])/(theparams->p1-theparams->p2);
				 }
			}
			else if(theparams->Btl_mode==1)
			{
			
			    lr1 =(FPN) (nn -2) /((FPN)(nn-2) +fstat[ii][jj]);
				lr = (FPN) (nn - 3) / ((FPN) (nn - 3) + 2*fstat[ii][jj]);
                lrb1 = (FPN) (nn - 3) / ((FPN) (nn - 3) + fstatx1[ii][jj]);
                lrb2 = (FPN) (nn - 3) / ((FPN) (nn - 3) + fstatx2[ii][jj]);
      
				if ( lr < (FPN) 0.0 || lr > (FPN) 1.0 )
						prf = -(FPN) 1.0;
				   else
				   {	
					   prf1 = betai(halfnn1, (FPN)0.5, lr1);
				       prf = betai(halfnn, 1, lr);
					   prfb1= betai(halfnn, (FPN)0.5,lrb1);
					   prfb2= betai(halfnn, (FPN)0.5,lrb2); 
				   }
	               if(cross==1 || cross ==2)
				   {
						putline(outf,'-',56);
						fprintf(outf, "\n%6d  %6d  %6.3f  %6.3f   %6.3f   %6.3f  %6.3f\n", ii, jj, beta0[ii][jj], beta1[ii][jj],  lr1, fstat[ii][jj], prf1);
						putline(outf,'-',56);
				   }
				   else if(cross ==3)
				   {
						putline(outf,'-',100);
						fprintf(outf, "\n%6d  %6d  %6.4f  %6.4f   %6.4f   %6.4f    %6.4f  %6.4f     %6.4f     %6.4f    %6.4f\n", ii, jj, beta0[ii][jj], beta1[ii][jj],beta2[ii][jj], fstat[ii][jj], prf,fstatx1[ii][jj],prfb1,fstatx2[ii][jj],prfb2);
						putline(outf,'-',100);
				   }
				if(cross==1)
				{
	     		
					/*print out the p1,p3,rmq table for each beta0 and beta1*/
                   fprintf(outf, "\np1     |     p3     Rmq\n");
				   putline(outf,'-',23);
                   for(kk=theparams->p1_start;kk<=theparams->p1_end;kk+=theparams->p1_step)
				   {
					   p3[ii][jj]=beta1[ii][jj]+2*beta0[ii][jj]-kk;	
				       Rmq[ii][jj]=(kk-beta0[ii][jj])/(2*kk-2*beta0[ii][jj]-beta1[ii][jj]);
					   fprintf(outf,"\n%6.3f | %6.3f  %6.3f",kk,p3[ii][jj], Rmq[ii][jj]);
				   
				   }
				}
				else if(cross==2)
				 {
				  
				   /*print out the p1,p3,rmq table for each beta0 and beta1*/
                   fprintf(outf, "\np2     |     p3     Rmq\n");
				   putline(outf,'-',23);
				   for(kk=theparams->p2_start;kk<=theparams->p2_end;kk+=theparams->p2_step)
				   {
					   p3[ii][jj]=beta1[ii][jj]+2*beta0[ii][jj]-kk;	
				       Rmq[ii][jj]=(kk-beta0[ii][jj])/(2*kk-2*beta0[ii][jj]-beta1[ii][jj]);
					   fprintf(outf,"\n%6.3f | %6.3f  %6.3f",kk,p3[ii][jj], Rmq[ii][jj]);
				   
				   }


				 
				}
				else if(cross==3)
				{
                   fprintf(outf, "\n   |   p1");
				   putline(outf,'_',12);
				   fprintf(outf, "\np2 |(p3,rmq)\n\n");  
				   fprintf(outf,"      |");
				   count=0;
				   for(ll=theparams->p1_start;ll<=theparams->p1_end;ll+=theparams->p1_step)
				   {
					   fprintf(outf,"%14.3f ",ll);
					   count++;
				   }
				   numb=count*15+7;
				   
				   putline(outf,'-',numb);
				 
				   for(kk=theparams->p2_start;kk<=theparams->p2_end;kk+=theparams->p2_step)
				   { 
					   fprintf(outf,"\n%6.3f|",kk);
						for(ll=theparams->p1_start;ll<=theparams->p1_end;ll+=theparams->p1_step)
						{
						    p3[ii][jj]  = (FPN)2*beta0[ii][jj]+beta1[ii][jj]+(FPN)0.5*beta2[ii][jj]-(FPN)0.5*(ll+kk);
                            Rmq[ii][jj] = (FPN)0.5*(ll-kk+beta2[ii][jj])/(ll-kk);
			 	            fprintf(outf,"(%6.3f,%6.3f)",p3[ii][jj], Rmq[ii][jj]);
				   
						}
				   }

				}
			}
			else
			{
					printf("Invalid Cross value! must be 1,2,or 3.\n");
					exit(0);
			}
		}
		
	}
    if(theparams->Btl_mode==0 )
	{
		BTL_print_S(nn,themap,outf, beta0,beta1,beta2,Rmq, theparams->p1,theparams->p2,p3, fstat,fstatx1,fstatx2,cross);
	}
	
  if (*outfile != '-')
    fileclose(outfile, outf);
	
return 0;
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file BTL_LRfunc.c
------------------------------------------------------------------ */


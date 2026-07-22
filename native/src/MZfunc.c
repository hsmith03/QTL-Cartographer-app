/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MZfunc.c) is part of QTL Cartographer
         
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


/*  Subroutines for multiple trait mapping.*/

#include "Main.h"

void WriteMarkerName(FILE *outf,markermap *themap,int chrom, int marker );

/*  
put phenotypes in lnpk->y.  Use only those phenotypes that are > (FPN) MISS_VAL 

Here's the rule, where wt is which trait:  

wt=0   If we are using no traits but those starting in +,
          then use only those individuals with phenotypic data at 
          ALL traits with names starting with +.
wt>n   If we are using all traits except those starting with -,
          then use only those individuals with phenotypic data at 
          ALL traits with names not starting with a -.
0<wt<n If we are doing single trait analysis, 
          then use only those individuals with phenotypic data.
          
This means that the sample size (in phenotypic measurements) will be the same for all traits.

(individs+i)->print_flag now tells whether that individual will be used or not.  y for yes, n for no.
*/
void copy_phenotypes(params *theparams,markermap *themap,individual *individs,linpakws *lnpk)
{
  int k,kk,jj,ii,wt;
  wt = theparams->whichtrait;
  kk = 0;

  if ( wt == 0 ) { /*  use no traits but those with names starting in + */
    for ( ii = 1 ; ii <= theparams->nn; ii++) {
      k = 0;
      for (jj = 1; jj <= theparams->traits; jj++) 
        if ( themap->tnames[jj][0] == '+' && individs[ii].y[jj]  <= (FPN) MISS_VAL )
          k = k+1;
      if ( k == 0 ) { /*Then this individual has all traits */
        individs[ii].print_flag = 'y';
        kk = kk+1;
        k = 0;
        for (jj = 1; jj <= theparams->traits; jj++) 
          if ( themap->tnames[jj][0] == '+' )  {
            k = k+1;
            lnpk->y[k][kk]  = individs[ii].y[jj] ;   
          }   
      }
      else
        individs[ii].print_flag = 'n';
    }
  }
  else if (   wt > theparams->traits ) {/*  use all traits except those with names starting in - */
    for ( ii = 1 ; ii <= theparams->nn; ii++) {
      k = 0;
      for (jj = 1; jj <= theparams->traits; jj++) /* Count the number of missing values. */
        if ( themap->tnames[jj][0] != '-' && individs[ii].y[jj]  <= (FPN) MISS_VAL )
          k = k+1;
      if ( k == 0 ) {
        individs[ii].print_flag = 'y';
        kk = kk+1;
        k = 0;
        for (jj = 1; jj <= theparams->traits; jj++) 
          if ( themap->tnames[jj][0] != '-' )  {
            k = k+1;
            lnpk->y[k][kk] =  individs[ii].y[jj]  ;   
          }   
      }
      else
        individs[ii].print_flag = 'n';
    }
  }
  else { /*  use  only trait wt ...  put it in row wt */
    for ( ii = 1 ; ii <= theparams->nn; ii++) {
      if ( individs[ii].y[wt] <= (FPN) MISS_VAL )
        individs[ii].print_flag = 'n';
      else {
        individs[ii].print_flag = 'y'; 
        kk = kk+1;
        lnpk->y[1][kk] = individs[ii].y[wt];   
      }     
    }
  }

  for ( ii = 1 ; ii <= lnpk->t ; ii++ )
    lnpk->samplesize[ii] = kk;
}
 
/**/
int set_cofactors(params *theparams,markermap *themap,genome *first,int **srranks ) {
  int marks,ii,bestrank,jj;
  genome *gptr;
  
/*  set the traits */  
    for (jj=1; jj<=theparams->traits; jj++ ) 
      if ( theparams->whichtrait < 1 && themap->tnames[jj][0] == '+' )
        srranks[0][jj] = 1;
      else if ( theparams->whichtrait > theparams->traits && themap->tnames[jj][0] != '-' )
        srranks[0][jj] = 1;
      else if ( theparams->whichtrait == jj ) 
        srranks[0][jj] = 1;
  
/*  For each marker, determine if it has at least one rank smaller than nbp
    and count the number of such markers.  */  
    marks = 0;
    for (ii=1; ii<=themap->ml; ii++ )  {
      bestrank = themap->ml;
      for ( jj=1; jj<=theparams->traits; jj++ )  
        if ( srranks[0][jj] > 0  && srranks[ii][jj] > 0 && srranks[ii][jj] < bestrank)
          bestrank =   srranks[ii][jj];
      if ( bestrank <= theparams->nbp ) {
        srranks[ii][0] +=1;    
        marks +=1;
      }
    }

/*  set for which markers to use as cofactors */
    for ( gptr=first; gptr != NULL ; gptr = gptr->next )   /*  binary code:  0 means don't use as cofactor, 1 means use*/
        if (gptr->markr > 0 &&  srranks[themap->ttable[gptr->chrom][gptr->markr]][0]  > 0 )
          gptr->whichqtl = 1;
        else
          gptr->whichqtl = 0;      



  return(marks);
}


/**/
void show_cofactors(FILE *outfile,params *theparams,markermap *themap,genome *first,int **srranks) {
  int ii;
  genome *gptr;
  if (srranks != NULL ) {
    fprintf(outfile,"\n# Cofactor usage           Cofactor ranks (at least one <= %d) for traits",theparams->nbp);
    fprintf(outfile,"\n\tChrom  Mark       Name  ");
    for ( ii=1; ii<=theparams->traits; ii++ )
      if ( srranks[0][ii] > 0 )
        fprintf(outfile," T(%d)",ii);
    for ( gptr=first; gptr != NULL ; gptr = gptr->next )   /*  binary code:  0 means don't use as cofactor, 1 means use*/
        if ( gptr->whichqtl   > 0 ) {
          fprintf(outfile,"\n\t  %3d  %4d %10s",gptr->chrom,gptr->markr,themap->names[themap->ttable[gptr->chrom][gptr->markr]]);
          for ( ii=1; ii<=themap->traits; ii++ )
            if ( srranks[0][ii] > 0 )
              fprintf(outfile," %4d",srranks[themap->ttable[gptr->chrom][gptr->markr]][ii]);
        
        }


    fprintf(outfile,"\n#      End of cofactor table");
  }
}




/*
Go into the SRmapqtl results file and get the ranks of the 
markers.
*/
int **jzget_srresults(params *theparams,markermap *themap)
{
  FILE *infile;
  int  ch,  ii,jj, trait,model,start,mark,chrom,**srranks;
  model = trait = -1;
  ch = 10;
  start = 0;
  infile = fileopen(theparams->srfile, "r");
  if ( infile == NULL )
    return(NULL);
  else
    srranks = imatrix(0,themap->ml,0,theparams->traits);
  for (ii=0; ii<=themap->ml; ii++ )
    for ( jj=0; jj<=theparams->traits; jj++ )
      srranks[ii][jj] = 0;
  
    ch = get_next_token(gname, MAXNAME, infile);
    if ( ch == EOF ) {
      fileclose(theparams->srfile, infile);  
      free_imatrix(srranks,0,themap->ml,0,theparams->traits);
      return(NULL);
    }
  
  start = -3;
  while ( (ch = get_next_token(gname, MAXNAME, infile)) != EOF ) {
      if ( start > 0 ) {
        if ( !strncmp(gname,"-end",4)  )
          start = -3;
        else {
          if ( start == 1 )
            chrom = atoi(gname);
          else if ( start == 2 ) 
            mark = atoi(gname);
          else if ( start == 3 )
            srranks[themap->ttable[chrom][mark]][trait] = atoi(gname);
          else if ( start == 5 )   /* end of a line. */
            start = 0;
          start +=1;
        }
      }
      else if ( start == -1 ) {
        trait = atoi(gname);
        start = 0;
      }
      else if ( gname[0] == '-' ) {
         if ( !strncmp(gname,"-start",6)  && start == 0 )
           start = 1; 
         if ( !strncmp(gname,"-trait",6)  && start == -2 )
           start = -1;
         if ( !strncmp(gname,"-FB",3)  && theparams->srm == 2 )
           start = -2;
         if ( !strncmp(gname,"-Forward",8)  && theparams->srm == 0 )
           start = -2;
         if ( !strncmp(gname,"-Backward",9)  && theparams->srm == 1 )
           start = -2;
         if ( !strncmp(gname,"-Zmapqtl",8)  && theparams->srm == 3 )
           start = -2;
         if ( !strncmp(gname,"-JZmapqtl",9)  && theparams->srm == 4 )
           start = -2;
      
      } 
  
  }  
  
  fileclose(theparams->srfile, infile);  
  return(srranks);
}

/*  Determine the starting and ending points for the analysis. */
int determine_endpoints(params *theparams,markermap *themap,genome **start,genome **end,genome *first,int *do_analysis)
{
  int go_on;
  genome *startptr,*endptr;
  *do_analysis = 1;

/* Do analysis for whole genome or theparams->wchrom */
  endptr = startptr = first;
  while (endptr->next != NULL) /*  This will do the entire genome... */
    endptr = endptr->next;
  if (theparams->wchrom > 0 && theparams->wchrom <= themap->m) {  /*  just do this chrom.*/
    while (startptr->chrom != theparams->wchrom)
      startptr = startptr->next;
    endptr = startptr;
    while (endptr->next != NULL && endptr->next->chrom == theparams->wchrom)
      endptr = endptr->next;
  }
  go_on = 0;
  *start = startptr;
  *end = endptr;
  return(go_on);
}



/*  
lnpk->xsave will contain the static rows of the design matrix.  
Row 1 will have a 1.0 in each column, and stands for the mean in the regression model.
Rows 2 - lnpk->k+1 will contain the otrait information.  
*/
void init_xsave(params *theparams,markermap *themap,individual *individs,linpakws *lnpk)
{
  int k,kk,jj,ii,row;
  k  = 0;

  for ( ii = 1 ; ii <= theparams->nn ; ii++ ) 
    if ( individs[ii].print_flag == 'y' ) {
      k  = k +1;  /* This is the column...that is the kth individual with data */
      row = 1;
      lnpk->xsave[row][k] = (FPN) 1.0;
      for ( jj = 1 ; jj <= themap->otraits ; jj++ )
        if ( themap->onames != NULL && themap->onames[jj][0] == '+' ) 
          for ( kk = 1 ; kk < themap->otypes[jj] ; kk++ ) {
            row = row+1;
            if ( individs[ii].oyt[jj] == kk ) 
              lnpk->xsave[row][k] = (FPN) 1.0;
            else
              lnpk->xsave[row][k] = (FPN) 0.0;              
          }       
    }
}

/*  
	Add the marker to the design matrix lnpk->xx. 
	Assume that row1 is the first row to be initialized.
	If ad = 1, then the additive effect will be added, 
	  whereas, if it is 0, it will be the dominance effect.
	  

Jiang and Zeng use

	  AA    Aa    aa
  b    2     1     0
  d    0     1     0
  
  If missing, b = 2pAA + pAa
              d =  1 - |b|
	  
This is old                     ....
	  AA    Aa    aa
  b    1     0    -1
  d  -1/2   1/2  -1/2
  
  If missing, b = pAA - paa
              d = 1/2 - |b|    ......
	  
	  
	If oti = 1, then the interaction term between the other 
	  traits and the marker will be used, otherwise not.
	Return the integer value of the next empty row.
	If the marker is a virtual marker, then take the 
	  expected values based on the flanking markers.
*/
int add_marker(params *theparams,markermap *themap,individual *individs,linpakws *lnpk,genome *gptr,int row1,int ad,int oti)
{
  FPN exp_val;
  int k,kk,jj,ii,row;
  k  = 0;
/* */
  for ( ii = 1 ; ii <= theparams->nn ; ii++ ) 
    if ( individs[ii].print_flag == 'y' ) {
      k  = k +1;  /* This is the column...that is the kth individual with data */
      row = row1;
      if ( gptr->markr == -10 ) { /* get the expected value for a virtual marker. */
        exp_val = (FPN) 1.0 + expect_mark(theparams,individs,ii,gptr) ;
        if (ad == 1)
          lnpk->xx[row][k] = exp_val;
        else
          lnpk->xx[row][k] = (FPN) 1.0 -  (FPN) 0.5*(FPN) fabs(exp_val);
      }
	  else {
	    if ( ad == 1 ) {
	      switch ( (individs + ii)->markers[gptr->chrom][gptr->markr]  ) {
	        case 1:  lnpk->xx[row][k]  =  (FPN) 2.0 ; break;
  	        case 0:  lnpk->xx[row][k]  =  (FPN) 1.0 ; break;
	        case -1: lnpk->xx[row][k]  =  (FPN)  0.0 ; break;
	        default: lnpk->xx[row][k] = (FPN) 1.0 + expect_mark(theparams,individs,ii,gptr); break;
	      }
	    }
	    else {
	      switch ( (individs + ii)->markers[gptr->chrom][gptr->markr]  ) {
	        case 1:  case -1:  lnpk->xx[row][k]  =  (FPN) 0.0 ; break;
  	        case 0:  lnpk->xx[row][k]  =  (FPN) 1.0 ; break;
	        default: lnpk->xx[row][k] = (FPN) 1.0 - (FPN) 0.5 * ((FPN) 1.0+(FPN) fabs(expect_mark(theparams,individs,ii,gptr))); break;
	      }
	    }
	  }
      if ( oti == 1 ) 
        for ( jj = 1 ; jj <= themap->otraits ; jj++ )
          if ( themap->onames != NULL && themap->onames[jj][0] == '+' ) 
            for ( kk = 1 ; kk < themap->otypes[jj] ; kk++ ) {
              row = row+1;
              if ( (individs+ii)->oyt[jj] == kk ) /* This puts in the marker by otrait interaction term. */
                lnpk->xx[row][k] = (FPN) 1.0 * lnpk->xx[row1][k];
              else
                lnpk->xx[row][k] = (FPN) 0.0;              
            }       
   }
   return(row+1);
}

/* 
   Determine the number of rows in the design matrix X 
   
*/
int how_many_rows(params *theparams,markermap *themap,individual *individs,int *otr,int marks)
{
  int rows,ii;
  switch (theparams->Model ) {
        default:  case 3:
          rows = 1; break;
        case 6:
          rows = 2*marks + 1; break;
        case 7:
          rows = 2*themap->knum[theparams->whichtrait] + 1;
  }
  *otr = 0;
  if ( themap->otraits > 0 ) {
     process_otraits(individs,theparams->nn,themap);
     for ( ii = 1 ; ii <= themap->otraits ; ii++ )
       if ( themap->onames[ii][0] == '+' ) {
         rows = rows +  2*themap->otypes[ii] - 2;
         *otr  = *otr + themap->otypes[ii] - 1;
       }
  }
  return(rows);
}

/*
c******calculates the determinant of a symetric covariance matrix.
c      a---input covariance matrix, unchanged by the routine.
c      b---contains information about the factorization of a, factorized by
c          spofa; if job=1 or 11, contains inverse of a.
c      det---the determinant of a, calculated by spodi.
c      lda---leading dimension of a and b matrix.
c      nt---number of traits.
c      job= 11 both determinant and inverse are computed
c            1 only inverse is computed
c           10 only determinant is computed
c  a(lda,lda), b(lda,lda), det
  a,b = dmatrix(1,lda,1,lda)
   shouldn't lda and nt be the same?  The variance-covariance matrix is square...
   I'm assuming that when this is called, nt = lda.
   Should these be part of lnpkws?  
*/
FPN invdet(FPN **a,FPN **b,int lda,int nt,FPN *z,int *kpvt,int job) 
{
  FPN det[3],rdet,work[3],work1[21]; /*  What is inert?  */
  int   inert[4],i,j,info,newjob;
  FPN   rcond;
  if ( lda == 1 ) {
    rdet = a[1][1];
    if ( a[1][1] != (FPN) 0.0 )
      b[1][1] = (FPN) 1.0 / a[1][1];
    else 
      b[1][1] = a[1][1];
    return(rdet);
  }
  for ( i=1 ; i<=lda ; i++ ) 
    for ( j=1 ; j<=lda ; j++ )
      b[i][j]=a[i][j];
 
  info = spofa(b,lda,nt);
 
  if ( info != 0 ) {
    for ( i=1 ; i<=lda ; i++ ) 
      for ( j=1 ; j<=lda ; j++ )
        b[i][j]=a[i][j];
    info = ssico(b,lda,nt,kpvt,&rcond,z);       
    if  (  rcond == (FPN) 0.0 ) {
      printf("Residual covariance is sigular");
      mypause();
    }
    newjob = job + 100; 
    info = ssidi(b,lda,nt,kpvt,det,inert,work1,newjob);
    rdet = det[1] * (FPN) pow(10.0,det[2]); 
  }
  else {
    info = spodi(b,lda,nt,work,job);
    rdet = work[1] * (FPN) pow(10.0,work[2]);
  }
/* What is this for?  Does it merely put the upper values in the lower part?
   That is what I'm assuming, and changing it so that the upper values are
   put in the lower triangle, and the diagonal is unchanged.  */
  for ( i=1 ; i<=lda ; i++ ) 
    for ( j=i+1 ; j<=lda ; j++ )
        b[i][j]=b[j][i];
  return(rdet);
}



/*
  Pick the background markers.  Reduce to 3, 6 and 7.
           
model =    
          3  => use no markers, just the mean (LB)
           
          6  =>   Specify two parameters:  ns and ws.  ns is the number of
                 background markers.  Pick the top ns of them as determined by a stepwise 
                 regression.  ws is a window size that blocks out markers on either side.
                 it is 10 cM by default.
          7  => Use the results of a previous scan with Zmapqtl to pick the background
                markers.  Eqtl should have been run, and the results summarized in the
                ->eqtl file will be used.  They should be in Rqtl.out format.
*/
void pick_cofactors(linpakws *lnpk,markermap *themap,aqtl *theqtls,params *theparams,genome *tgptr,genome *agptr )
{
  int ii,col;
  int model,chrom,lfm,ns;
  FPN lpos, rpos,wind;
  genome *gptr;
  
  for ( gptr = agptr ; gptr != NULL ; gptr = gptr->next )
    if ( gptr->whichqtl < 0 ) /* Set all whichqtls to > 0...meaning use none but maintain ranks */
      gptr->whichqtl = -1 * gptr->whichqtl;
      
  model = theparams->Model;
  chrom = tgptr->chrom;
  lfm = tgptr->markr;
  ns = theparams->nbp;
  for ( ii = 1;ii <= themap->ml;ii++ )
    lnpk->bp[1][ii]   = lnpk->bp[2][ii]   = 0;
  
  if ( theparams->window > (FPN) 0.0 )
    wind = theparams->window * (FPN) 0.01;
  else
    wind = (FPN) WIN_SIZE * (FPN) 0.01;


  col = 1;
  if ( model == 6 ) {
/* Use up to ns background markers, where a forward stepwise regression was done to
   determine them.  */
	 rpos = (FPN) MAXCHROMLEN;
	 lpos = (FPN) 0.0;
	 for ( gptr = agptr ; gptr != NULL ; gptr = gptr->next ) {
		if ( gptr->chrom == chrom && gptr->markr == lfm )
		  lpos = gptr->pos - wind;
		else if ( gptr->chrom == chrom && gptr->markr == lfm+1 )
		  rpos = gptr->pos + wind;
		gptr->mxo = (FPN) 0.0;
	 }
	 ii = 1;
	 gptr = agptr;
	 do {
	   if ( gptr->markr != 0 && gptr->whichqtl <= ns && gptr->whichqtl > 0 )
		 if (gptr->chrom != chrom || (gptr->chrom == chrom && gptr->pos < lpos) || (gptr->chrom == chrom && gptr->pos > rpos) )
           ii = ii+1;

       gptr = gptr->next;
     } while ( gptr != NULL );
	 lnpk->bp[1][1]   = ii;
	 lnpk->bp[2][1]   = 6;
    ii = 1;
    gptr = agptr;
    do {
		if ( gptr->markr != 0 && gptr->whichqtl <= ns && gptr->whichqtl > 0)
		  if (gptr->chrom != chrom || (gptr->chrom == chrom && gptr->pos < lpos) || (gptr->chrom == chrom && gptr->pos > rpos) ) {
          ii = ii+1;
          lnpk->bp[1][ii]   = gptr->chrom;
          lnpk->bp[2][ii]   = gptr->markr;
          gptr->whichqtl = -1 * gptr->whichqtl;
        }
      gptr = gptr->next;
    } while ( gptr != NULL );    
  }
  else if ( model == 7 ) {
/* Use the results summarized in the Eqtl.out file.
   Use all that are for this trait. 
   1. get eqtl.out results and put into the standard structure ... done in Zmain.c
   2. use those results for the background...we assume that the
      data have been run through Zmapqtl for the trait we are now
      analyzing.     
      
      Use virtual markers for the background.
      
      
  Now what?  do scan for all traits singly.  Run Eqtl.  Get results.  Use all
  of the qtls? 
  
  
  
  Open Question: What if two QTL estimates are in the same interval?  
*/
     ns = 0;
     for ( ii = 1 ; ii <= themap->traits ; ii++ )
       ns = ns + themap->knum[ii];
     for ( ii = 1 ; ii <= ns ; ii++ )
       theqtls[ii].d = -(FPN) 1.0;

	 rpos = (FPN) MAXCHROMLEN;
	 lpos = (FPN) 0.0;
	 for ( gptr = agptr ; gptr != NULL ; gptr = gptr->next ) {
		if ( gptr->chrom == chrom && gptr->markr == lfm )
		  lpos = gptr->pos - wind;
		else if ( gptr->chrom == chrom && gptr->markr == lfm+1 )
		  rpos = gptr->pos + wind;
	 }
	 ii = 1;
	 gptr = agptr;
	 do {
	  if ( gptr->whichqtl <= ns && gptr->whichqtl > 0 ) /* the pxo takes into account the position of the virtual marker. */
		if (gptr->chrom != chrom || ( gptr->pxo < lpos  ||  gptr->pxo > rpos) ) {
          ii = ii+1;
          theqtls[gptr->whichqtl].d = (FPN) 1.0;
          gptr->whichqtl = -1 * gptr->whichqtl;
        }
      gptr = gptr->next;
    } while ( gptr != NULL );
    
    lnpk->bp[1][1]   = 0;
    for ( ii = 1 ; ii <= ns ; ii++ )
       if ( theqtls[ii].d > (FPN) 0.0) {
         lnpk->bp[1][1]   = lnpk->bp[1][1] + 1;
         if ( theparams->verbosity == 2 )
           printf("\n Chrom %d, Marker %d, Trait %d",theqtls[ii].chrm,theqtls[ii].mrk,theqtls[ii].trait);
       }
     lnpk->bp[2][1]   = 7;
    
  }
  else {	/* The default is interval mapping a la  Lander and Botstein (1989). */
	 lnpk->bp[1][1]  = 1;
	 lnpk->bp[2][1] = 3;
  }
  for ( gptr = agptr ; gptr != NULL ; gptr = gptr->next ) /* change sign of all whichqtls */
    gptr->whichqtl = -1 * gptr->whichqtl;
}


void WriteMarkerName(FILE *outf,markermap *themap,int chrom, int marker ) {

      if ( themap->names != NULL ) {
        if ( marker > 0)
          fprintf(outf," %-12s",themap->names[ themap->ttable[chrom][marker] ] );
        else
          fprintf(outf," Telomere.%-2d ",chrom);
      }
      else
        fprintf(outf," Marker%d.%d",chrom,marker);


}


/*
The convention for lnpk->estimate[trait][1-9]:

Col       1    2    3   4   5   6   7     8     9
hypo
1        H1             a               H1:H0
30                  H3      a3      d3  H3:H0
31       H1  H1:H0  H3  a1  a3      d3  H3:H0  H3:H1
32     H2:H0   H2   H3      a3  d2  d3  H3:H0  H3:H2

14       H1   Hge       a               H1:H0  H1:Hge
34            Hge   H3      a3      d3  H3:H0  H3:Hge

Also, the joint parameters are in the 0th row.
*/
void write_position_results(linpakws *lnpk,params *theparams,char *outfile,genome *gptr,FPN abspos,markermap *themap)
{

  FILE *outf;
  int trait,i,realtrait;
/*Print out the joint LR*/
  trait = 0;
    for ( i = 0 ; i < MAXNAME ; i++ )
      gname[i] = '\0';
    sprintf(gname,"%s%d",outfile,trait);
    outf = fileopen(gname, "a");
    if (outf != NULL) {      
	  fprintf(outf,"\n     %-5d      %-5d",gptr->chrom,gptr->markr );
	  WriteMarkerName(outf,themap,gptr->chrom,gptr->markr);
        
	  fprintf(outf," %10.7f %10.4f",abspos,lnpk->estimates[trait][8] );
	  if ( theparams->ihypo == 31 )
	    fprintf(outf," %10.4f %10.4f", lnpk->estimates[trait][9] ,lnpk->estimates[trait][2]);
	  else if ( theparams->ihypo == 32 )
	    fprintf(outf," %10.4f %10.4f", lnpk->estimates[trait][9], lnpk->estimates[trait][1]);
	  else if ( theparams->ihypo == 34 )
	    fprintf(outf," %10.4f %10.7f %10.7f", lnpk->estimates[trait][9], lnpk->estimates[0][5], lnpk->estimates[0][7] );
 	  else if ( theparams->ihypo == 14 )
	    fprintf(outf," %10.4f %10.7f ", lnpk->estimates[trait][9], lnpk->estimates[0][5]  );
   } 
    fileclose(gname, outf);
 /* t files, one for each trait. */
  for ( trait = 1 ; trait <= lnpk->t ; trait++ ) {
    realtrait = actual_trait(theparams,themap,trait);
    for ( i = 0 ; i < MAXNAME ; i++ )
      gname[i] = '\0';
    sprintf(gname,"%s%d",outfile,realtrait);
    outf = fileopen(gname, "a");
    if (outf != NULL) {      
	  fprintf(outf,"\n     %-5d      %-5d",gptr->chrom,gptr->markr );
	  WriteMarkerName(outf,themap,gptr->chrom,gptr->markr);
        
	  fprintf(outf," %10.7f %10.4f",abspos,lnpk->estimates[trait][8] );
	  if ( theparams->ihypo == 1 || theparams->ihypo == 10 )
	    fprintf(outf," %10.7f", lnpk->estimates[trait][4] ); 
	  if ( theparams->ihypo == 14   )
	    fprintf(outf," %10.4f %10.7f",  lnpk->estimates[trait][9] ,lnpk->estimates[trait][5] ); 
	  if ( theparams->ihypo == 31 )
	    fprintf(outf," %10.4f %10.4f", lnpk->estimates[trait][9] ,lnpk->estimates[trait][2]);
	  else if ( theparams->ihypo == 32 )
	    fprintf(outf," %10.4f %10.4f", lnpk->estimates[trait][9], lnpk->estimates[trait][1]);
	  else if ( theparams->ihypo == 34 )
	    fprintf(outf," %10.4f ", lnpk->estimates[trait][9] );
	  if (  theparams->ihypo > 14 )
	    fprintf(outf," %10.7f %10.7f", lnpk->estimates[trait][5] ,lnpk->estimates[trait][7]);
	  if ( theparams->ihypo == 31 )
	    fprintf(outf," %10.7f", lnpk->estimates[trait][4]);
	  else if ( theparams->ihypo == 32 )
	    fprintf(outf," %10.7f", lnpk->estimates[trait][6]);
    } 
    fileclose(gname, outf);
  }
}

void write_jzheader(char *outfile,params *theparams,char *onamae, char *chptr,markermap *themap,int oc,genome *first,int **srranks)
{
  int cross,trait,t,i,wt,ii;
  FILE *outf;
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  t = how_many_traits(theparams,themap); 
  wt = theparams->whichtrait; 
  for ( ii = 0 ; ii < MAXNAME ; ii++ )
    gname[ii] = '\0';
  sprintf(gname,"%s0",outfile);
  if (oc == 1 )
    print_head(onamae,gname,chptr,1,65,theparams);
  outf = fileopen(gname, "a");
  if (outf != NULL) { 
    if ( oc == 1 ) {
	    fprintf(outf, "\n#  The position is from the left telomere on the chromosome");
	    fprintf(outf, "\n-window     %6.2f      Window size for model 6",theparams->window);
	    fprintf(outf, "\n-background %6d      Max. Rank for Background parameters in model 6",theparams->nbp);
	    if ( theparams->Model == 6 ) 
	      show_cofactors(outf,theparams,themap,first,srranks) ;
	    fprintf(outf, "\n-Model      %6d      Model number",theparams->Model);
	    fprintf(outf, "\n-traits     %6d      Number of Traits in Joint Analysis", t);
	    fprintf(outf, "\n-Names ");
        for ( i = 1 ; i <= theparams->traits ; i++ ) 
          if (  (wt == 0 && themap->tnames[i][0] == '+') || (wt > theparams->traits && themap->tnames[i][0] != '-') || (wt == i)  ) 
	        fprintf(outf, "\n\t %d   %s",i,themap->tnames[i]);
    	fprintf(outf, "\n-cross      %6s      Cross", theparams->thecross);
	    fprintf(outf, "\n-ihypo      %6d      Hypothesis test",theparams->ihypo);
    fprintf(outf, "\n#\n#  Note that our Likelihood ratio test statistic compares two nested hypotheses");
    fprintf(outf, "\n#  and is two times the negative natural log of the ratio of the likelihoods.  For example,");
    fprintf(outf, "\n#  assume that  hypothesis H0 is nested within H1 and that they have likelihoods L0 and L1 respectively.");
    fprintf(outf, "\n#  Then, the \"Likelihood Ratio Test Statistic\" is -2ln(L0/L1). \n#");
		fprintf(outf, "\n# Chromosome   Marker   MarkerName  Position   ");
		if ( theparams->ihypo == 1 || theparams->ihypo == 10)
		  fprintf(outf," LR(H1:H0) " );
		else if ( theparams->ihypo == 14 )
		  fprintf(outf," LR(H1:H0)  LR( GxE )       a      " );
		else  if ( theparams->ihypo == 30 )
		  fprintf(outf," LR(H3:H0) " );
		else if ( theparams->ihypo == 31 )
		  fprintf(outf," LR(H3:H0)  LR(H3:H1)  LR(H1:H0) " );
		else if ( theparams->ihypo == 32 )
		  fprintf(outf," LR(H3:H0)  LR(H3:H2)  LR(H2:H0) " );
		else if ( theparams->ihypo == 34 )
		  fprintf(outf," LR(H3:H0)  LR( GxE )      a         d      " );
	    fprintf(outf, "\n-s");
	}
    else
      fprintf(outf,"\n-e\n");
    fileclose(gname, outf);
  }

  
  trait = 0;
  for ( i = 1 ; i <= theparams->traits ; i++ ) {
    if ( ( (wt == 0 && themap->tnames[i][0] == '+') || (wt > theparams->traits && themap->tnames[i][0] != '-') || (wt == i) || (i==0)) ) {
      for ( ii = 0 ; ii < MAXNAME ; ii++ )
        gname[ii] = '\0';
      sprintf(gname,"%s%d",outfile,i);
      if ( i > 0 )
        trait = trait+1;
      if (oc == 1 )
        print_head(onamae,gname,chptr,1,65,theparams);
	  outf = fileopen(gname, "a");
	  if (outf != NULL) {
	    if ( oc == 1 ) {
		    fprintf(outf, "\n#  The position is from the left telomere on the chromosome");
		    fprintf(outf, "\n-window     %6.2f      Window size for model 6",theparams->window);
		    fprintf(outf, "\n-background %6d      Background parameters in model 6",theparams->nbp);
		    fprintf(outf, "\n-Model      %6d      Model number",theparams->Model);
		    fprintf(outf, "\n-trait      %6d      Trait ", i);
		    if ( themap->tnames != NULL && i > 0)
		      fprintf(outf,": %s",themap->tnames[i]);
		    fprintf(outf, "\n-cross      %6s      Cross", theparams->thecross);
		    fprintf(outf, "\n-ihypo      %6d      Hypothesis test",theparams->ihypo);
    fprintf(outf, "\n#\n#  Note that our Likelihood ratio test statistic compares two nested hypotheses");
    fprintf(outf, "\n#  and is two times the negative natural log of the ratio of the likelihoods.  For example,");
    fprintf(outf, "\n#  assume that  hypothesis H0 is nested within H1 and that they have likelihoods L0 and L1 respectively.");
    fprintf(outf, "\n#  Then, the \"Likelihood Ratio Test Statistic\" is -2ln(L0/L1). \n#");
		    fprintf(outf, "\n# Chromosome   Marker   MarkerName  Position   ");
			if ( theparams->ihypo == 1 || theparams->ihypo == 10 )
			  fprintf(outf," LR(H1:H0)       a    " );
			else if ( theparams->ihypo == 14 )
			  fprintf(outf," LR(H1:H0)  LR( GxE )       a " );
			else  if ( theparams->ihypo == 30 )
			  fprintf(outf," LR(H3:H0)       a3         d3     " );
			else if ( theparams->ihypo == 31 )
			  fprintf(outf," LR(H3:H0)  LR(H3:H1)  LR(H1:H0)      a3        d3         a1" );
			else if ( theparams->ihypo == 32 )
			  fprintf(outf," LR(H3:H0)  LR(H3:H2)  LR(H2:H0)      a3        d3         d2" );
			else if ( theparams->ihypo == 34 )
			  fprintf(outf," LR(H3:H0)  LR( GxE )       a3        d3 " );
	        fprintf(outf, "\n-s");
	    }
	    else
	        fprintf(outf, "\n-e\n");
	    fileclose(gname, outf);
	  }
    }
  }

}

/*
Write a header for the new data set, where expected values are placed at each test site.
*/
void write_alttraits(char *outfile,params *theparams, markermap *themap,individual *individs, linpakws *lnpk)
{
  int t,i,j,wt,cntr;
  FILE *outf;
  t = how_many_traits(theparams,themap); 
  wt = theparams->whichtrait; 
  outf = fileopen(outfile, "a");
  if (outf != NULL) { 
	    fprintf(outf, "\n-n          %6d      Sample Size",lnpk->samplesize[1]);
    for ( i = 1 ; i <= theparams->traits ; i++ ) 
      if (  (wt < 1 && themap->tnames[i][0] == '+') || (wt > theparams->traits && themap->tnames[i][0] != '-') || (wt == i)  ) {
	    fprintf(outf, "\n-Trait %d   %s\n",i,themap->tnames[i]);
	    cntr = 1;
	    for ( j=1 ; j<=lnpk->samplesize[1]; j++ ) {
	      fprintf(outf," %12.4g ",lnpk->y[i][j]);
	      cntr +=1;
	      if ( cntr > 5 && j < lnpk->samplesize[1]) {
	        fprintf(outf,"\n");
	        cntr = 1;
	      }
	    }
	  }
	for ( i = 1 ; i <= themap->otraits ; i++ )
      if ( themap->onames != NULL && themap->onames[i][0] == '+' ) {
        
	    fprintf(outf, "\n-Otrait %d   %s\n",i,themap->onames[i]);
	    cntr = 1;
	    for ( j=1 ; j<=lnpk->samplesize[1]; j++ ) {
	      fprintf(outf," %4d ",individs[j].oyt[i]);
	      cntr +=1;
	      if ( cntr > 15 && j < lnpk->samplesize[1]) {
	        fprintf(outf,"\n");
	        cntr = 1;
	      }
	    }        
      }	  
    fileclose(outfile, outf);
  }
}

/*
Write a header for the new data set, where expected values are placed at each test site.
*/
long write_altheader(char *outfile,params *theparams,char *onamae, char *chptr,markermap *themap,int oc)
{
  int cross,t,i,wt,ot;
  long fileposition;
  FILE *outf;
  ot = 0;
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  t = how_many_traits(theparams,themap); 
  wt = theparams->whichtrait; 
  if (oc == 1 )
    print_head(onamae,outfile,chptr,1,66,theparams);
	for ( i = 1 ; i <= themap->otraits ; i++ )
      if ( themap->onames != NULL && themap->onames[i][0] == '+' ) 
        ot +=1;
  outf = fileopen(outfile, "a");
  if (outf != NULL) { 
    if ( oc == 1 ) {
	    fprintf(outf, "\n#  This output of JZmapqtl is meant to be used with MultiRegress ");
	    fprintf(outf, "\n-walk       %6.2f      Interval distance in cM",100.0*theparams->walk);
    	fprintf(outf, "\n-cross      %6s      Cross", theparams->thecross);
	    fprintf(outf, "\n-otraits    %6d      Number of explanatory variables",ot);
	    fprintf(outf, "\n-traits     %6d      Number of Traits ", t);
	    fprintf(outf, "\n-positions ");
	    fileposition = ftell(outf);
	    fprintf(outf, "                             Number of positions");	    
	    
	}
    fileclose(outfile, outf);
  }
  return(fileposition);
}

/*Insert the number of positions*/
void insert_positions(char *outfile,int positions,long fp) {
  FILE *outf;
  outf = fileopen(outfile, "r+");

      fseek(outf,fp,SEEK_SET);
      
      fprintf(outf,"%d",positions);
    fileclose(outfile, outf);

}

int    do_jzexpecteds(params *theparams,genome *startptr,genome *endptr,markermap *themap, individual *individs,char *outfile, linpakws *lnpk)
{
  int cross,t,wt,ipos ;
  genome *gptr;
  FILE *outf;
  FPN excess;
  ipos = 0;
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  t = how_many_traits(theparams,themap); 
  wt = theparams->whichtrait; 
  outf = fileopen(outfile, "a");
  excess = (FPN) MIN_DIST ;
  for ( gptr=startptr ; gptr != endptr->next ; gptr = gptr->next ) 
    if ( gptr->dist > excess )
      excess = expected_int(lnpk,theparams,themap,individs,gptr,outf,&ipos ,excess);
    else
      excess = excess - gptr->dist;
  fileclose(outfile, outf);
  return(ipos);
}

/*
 Calculate the expected value of the marker at the test site for each individual.
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
*/
FPN expected_int(linpakws *lnpk,params *theparams,markermap *themap,individual *individs,genome *gptr,FILE *outf,int *ipos, FPN excess)
{
  FPN abspos,relpos,thetaL,thetaR,evalue;
  FPN aQQ,aQq,aqq,dQQ,dQq,dqq; 
  int i,cntr;
/*  We will go along the entire interval in gptr...*/
  aQQ = (FPN) 1.0;
  aQq = (FPN) 0.0;
  aqq = -(FPN) 1.0;
  dqq = dQQ = -(FPN) 0.5;
  dQq = (FPN) 0.5;
  if ( theparams->crosstype == 1 ) {
    aQQ = (FPN) 0.5;
    aQq = -(FPN) 0.5;
    aqq = (FPN) 0.0;  
  }
  else if ( theparams->crosstype == 2 ) {
    aQq = (FPN) 0.5;
    aqq = -(FPN) 0.5;
    aQQ = (FPN) 0.0;    
  }
   
  relpos = excess;
  while ( relpos < gptr->dist ) {
    *ipos +=1;
    abspos = gptr->pos + relpos;
	if (gptr->markr == 0)
	  thetaL = -(FPN) 2.0;
	else
	  thetaL = relpos;
	if ( gptr->markr == themap->mpc[gptr->chrom] )
	  thetaR = -(FPN) 2.0;
	else
	  thetaR = gptr->dist - relpos;
    calc_priors(theparams,lnpk->pp1,lnpk->pp2,theparams->nn,individs,gptr,thetaL,thetaR);
/*
    lnpk->pp1 = pQQ unless this is a BC2 population, then it is pqq
    lnpk->pp2 = pQq 
    1-lnpk->pp1-lnpk->pp2 = pqq unless BC2, then 1-lnpk->pp1-lnpk->pp2 = pQQ    
*/
    if ( gptr->markr > 0 ) 
      fprintf(outf,"\n-Site %d -parameter additive -chromosome %d  -marker %d -name %s -position %f  -values\n",*ipos,gptr->chrom, gptr->markr,themap->names[ themap->ttable[gptr->chrom][gptr->markr] ], abspos);
    else
      fprintf(outf,"\n-Site %d -parameter additive -chromosome %d  -marker %d -name Telomere%dL -position %f  -values\n",*ipos,gptr->chrom, gptr->markr,gptr->chrom, abspos);

    cntr = 1;
    for ( i = 1 ; i <= lnpk->samplesize[1]; i++ ) {
      if ( theparams->crosstype == 1 ) 
        evalue = aQQ * lnpk->pp1[i] + aQq * lnpk->pp2[i] ; 
      else if ( theparams->crosstype == 2 ) 
        evalue = aqq * lnpk->pp1[i] + aQq * lnpk->pp2[i] ; 
      else  
        evalue = aQQ * lnpk->pp1[i] + aQq * lnpk->pp2[i] + aqq * ((FPN) 1.0 -  lnpk->pp1[i] -  lnpk->pp2[i]); 


	      fprintf(outf," %12.4g ",evalue);
	      cntr +=1;
	      if ( cntr > 5 && i < lnpk->samplesize[1]) {
	        fprintf(outf,"\n");
	        cntr = 1;
	      }
    }
    if ( theparams->crosstype == 3 || theparams->crosstype == 4 ) {
	    fprintf(outf,"\n-Site %d -parameter dominance -chromosome %d  -marker %d -name ",*ipos,gptr->chrom,gptr->markr);
	    WriteMarkerName(outf,themap,gptr->chrom,gptr->markr);
	    fprintf(outf,"-position %f -values\n",  abspos);
	    cntr = 1;
	    for ( i = 1 ; i <= lnpk->samplesize[1]; i++ ) {
	        evalue = dQQ * lnpk->pp1[i] + dQq * lnpk->pp2[i] + dqq * ((FPN) 1.0 -  lnpk->pp1[i] -  lnpk->pp2[i]); 


		      fprintf(outf," %12.4g ",evalue);
		      cntr +=1;
		      if ( cntr > 5 && i < lnpk->samplesize[1]) {
		        fprintf(outf,"\n");
		        cntr = 1;
		      }
	    }
    }
    relpos = relpos + theparams->walk;
  }
  relpos = relpos - gptr->dist;
  if ( gptr->next != NULL && gptr->chrom != gptr->next->chrom )
    relpos = (FPN) MIN_DIST;
  return(relpos);
}



/*
  This is a driver for actually doing the mapping.  It needs to know the start
  and ending points for the analysis, and a file to print the results.  
*/
FPN do_jzanalysis(params *theparams,genome *startptr,genome *endptr,markermap *themap,aqtl *theqtls,individual *individs,char *outfile,genome *agptr,linpakws *lnpk )
{
  FPN  maxe,tmax,ts0;
  int  ipos,nrows,ii;
  genome *gptr;
  ipos = 0;
  maxe = (FPN) 0.0;
  if ( theparams->verbosity == 1 ) {
        fprintf(stdout,"\nThe Current Mapping Interval is \nchromosome after marker named " );      
	    fflush(stdout);
  }
  for ( gptr=startptr ; gptr != endptr->next ; gptr = gptr->next ) 
    if ( gptr->dist > (FPN) 0.0 ) {
  	  if ( theparams->verbosity == 1 ) {
        fprintf(stdout,"\n %5d          %5d",gptr->chrom,gptr->markr);  
        WriteMarkerName(stdout,themap,gptr->chrom,gptr->markr);
	    fflush(stdout);
	  }
	  pick_cofactors(lnpk,themap,theqtls,theparams,gptr,agptr);
      ts0 = fit_null(lnpk,theparams,individs,themap,theqtls,gptr, &nrows);
      lnpk->ts0[0] = -(FPN) 0.5 * (FPN) lnpk->samplesize[1] * ((FPN) log( fabs(ts0))+ (FPN) lnpk->t);
      for ( ii = 1 ; ii <= lnpk->t ; ii++ )
        lnpk->ts0[ii] = -(FPN) 0.5 * (FPN) lnpk->samplesize[1] * ((FPN) log(lnpk->s2[ii][ii])+ (FPN) 1.0 );
      tmax = zmap_int(lnpk,theparams,themap,individs,gptr,outfile,nrows);
      if (tmax > maxe)
        maxe = tmax; 
    }
  if ( theparams->verbosity == 1 ) 
    printf("\n");

  return(maxe);
}

/*
  Do (composite) interval mapping on the interval pointed to by gptr.
*/
FPN zmap_int(linpakws *lnpk,params *theparams,markermap *themap,individual *individs,genome *gptr,char *outfile,int nrows)
{
  FPN tmax,abspos,relpos,thetaL,thetaR;
  int hypo,bail_status1,bail_status2,bail_status3,trait;
  tmax = (FPN) 0.0;
  hypo = 1;
/*  We will go along the entire interval in gptr...*/

  relpos = (FPN) MIN_DIST;
  while ( relpos < gptr->dist ) {
    abspos = gptr->pos + relpos;
	if (gptr->markr == 0)
	  thetaL = -(FPN) 2.0;
	else
	  thetaL = relpos;
	if ( gptr->markr == themap->mpc[gptr->chrom] )
	  thetaR = -(FPN) 2.0;
	else
	  thetaR = gptr->dist - relpos;
    calc_priors(theparams,lnpk->pp1,lnpk->pp2,theparams->nn,individs,gptr,thetaL,thetaR);

    if ( theparams->ihypo == 1 || theparams->ihypo == 10 || theparams->ihypo == 14) {
      bail_status1 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,1,theparams);
      for ( trait = 0 ; trait <= lnpk->t ; trait++ ) 
        if ( bail_status1 == 0 ) 
          lnpk->estimates[trait][8] = (FPN) 2.0 * (lnpk->estimates[trait][1] - lnpk->ts0[trait]);
        else 
          lnpk->estimates[trait][8] =   -(FPN) SIGNOFDEVIL;
      if (theparams->ihypo == 14) {
        bail_status2 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,4,theparams);
        for ( trait = 0 ; trait <= lnpk->t ; trait++ ) 
          if ( bail_status2 == 0 ) 
            lnpk->estimates[trait][9] = (FPN) 2.0 * (lnpk->estimates[trait][1] - lnpk->estimates[trait][2]);
          else 
            lnpk->estimates[trait][9] =   -(FPN) SIGNOFDEVIL;      
      }
    }
    else  {
      bail_status3 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,3,theparams);
      if ( theparams->ihypo == 31 )  
        bail_status1 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,1,theparams);
      if ( theparams->ihypo == 32 ) 
        bail_status2 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,2,theparams);
      if ( theparams->ihypo == 34 )
        bail_status2 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,4,theparams);
      
      for ( trait = 0 ; trait <= lnpk->t ; trait++ ) {
        if ( bail_status3 == 0 ) 
         lnpk->estimates[trait][8] = (FPN) 2.0 * ( lnpk->estimates[trait][3] - lnpk->ts0[trait]);
        else 
         lnpk->estimates[trait][8] = -(FPN) SIGNOFDEVIL;
        if ( theparams->ihypo == 31 ) {
          if ( bail_status1 == 0 )  
            lnpk->estimates[trait][2] = (FPN) 2.0 * ( lnpk->estimates[trait][1] - lnpk->ts0[trait]);
          else 
            lnpk->estimates[trait][2] =  -(FPN) SIGNOFDEVIL;
          if ( bail_status1 == 0 && bail_status3 == 0 )  
             lnpk->estimates[trait][9] = (FPN) 2.0 * ( lnpk->estimates[trait][3] - lnpk->estimates[trait][1]);
          else 
             lnpk->estimates[trait][9] = -(FPN) SIGNOFDEVIL;
        }
        if ( theparams->ihypo == 32 ) {
          if ( bail_status2 == 0 )  
            lnpk->estimates[trait][1] = (FPN) 2.0 * ( lnpk->estimates[trait][2] - lnpk->ts0[trait]);
          else 
            lnpk->estimates[trait][1] =  -(FPN) SIGNOFDEVIL;
          if ( bail_status2 == 0 && bail_status3 == 0 )  
             lnpk->estimates[trait][9] = (FPN) 2.0 * ( lnpk->estimates[trait][3] - lnpk->estimates[trait][2]);
          else 
             lnpk->estimates[trait][9] = -(FPN) SIGNOFDEVIL;
        }
        if ( theparams->ihypo == 34 ) {
          if ( bail_status2 == 0 )  
             lnpk->estimates[trait][9] = (FPN) 2.0 * ( lnpk->estimates[trait][3] - lnpk->estimates[trait][2]);
          else 
             lnpk->estimates[trait][9] = -(FPN) SIGNOFDEVIL;
        }
        
      }
    }
    if ( lnpk->estimates[lnpk->t][1] > tmax )
      tmax = lnpk->estimates[lnpk->t][1];
    write_position_results(lnpk,theparams,outfile,gptr,abspos,themap);
    relpos = relpos + theparams->walk;
  }
  return(tmax);
}

/* Fit the null hypothesis for  the background markers specified in bp. 
 gptr->chrom, gptr->markr,   *(*(lnpk->bp + 1) + 1)*/
FPN fit_null(linpakws *lnpk,params *theparams,individual *individs,markermap *themap,aqtl *theqtls ,genome *first,int *nrows)
{

  int ii,kk,next_row,n,nn,info,*itmp,cross;
  genome *tgptr,lg,rg;
  FPN *tmp,ts0,thetaL,thetaR;
  cross = theparams->cross;
  if ( theparams->tcross ==1 || theparams->tcross == 2 ) 
    cross = theparams->tcross;   /* If this is a test cross to P1 or P2, treat it like a backcross.*/
  if ( theparams->tcross == 12 ) /* otherwise, treat it like an SF2 */
    cross = 3;
  if ( cross == 4 )              /* RFx is like SFx as far as the estimates go. */
    cross = 3;
  
  itmp = NULL;
  tmp = NULL;
  nn=n = lnpk->samplesize[1];
  /* Initialize the design matrix xx. */
  for ( ii = 1; ii <= lnpk->k+1 ; ii++ ) /* copy first k+1 rows into xx from xsave. */
    scopy(n, lnpk->xsave[ii], 1,  lnpk->xx[ii], 1);
  next_row = lnpk->k+2;
  if ( theparams->Model > 3 ) /* now, copy in all markers that have already been selected. */
    for ( tgptr = first ; tgptr != NULL ; tgptr = tgptr->next )
      if ( tgptr->whichqtl > 0 ) { 
        if ( theparams->Model == 7 ) {           
	      thetaL = mapfunc((theqtls+tgptr->whichqtl)->c1,1);
	      thetaR = mapfunc((theqtls+tgptr->whichqtl)->c2,1);
          add_virtual_node(tgptr,&lg,&rg,thetaL,thetaR);
          next_row = add_marker(theparams,themap,individs,lnpk,&rg,next_row,1,1);         
          if ( cross == 3 )
            next_row = add_marker(theparams,themap,individs,lnpk,tgptr,next_row,0,1);
          del_virtual_node(tgptr);
        }
        else {
          next_row = add_marker(theparams,themap,individs,lnpk,tgptr,next_row,1,1);
          if ( cross == 3 )
            next_row = add_marker(theparams,themap,individs,lnpk,tgptr,next_row,0,1);
        }
      }
   *nrows = next_row - 1;
/*
  Zero out the residual matrix and qraux for safety.   
*/
  for ( kk = 1 ; kk <= lnpk->t ; kk++ )
    for (ii = 1; ii <= nn; ii++)  
       lnpk->wrsd[kk][ii] = (FPN) 0.0;
  for (ii = 1; ii <= *nrows ; ii++)
    *(lnpk->qraux + ii) = (FPN) 0.0;
/*  Do an sqrdc on the design matrix.  */
  info = sqrdc(lnpk->xx,n,n,*nrows,lnpk->qraux,itmp,0);
/*  
    For each trait in turn, do an sqrsl on it.  Put the result in the
    row for lnpk->rsd.  Then, the Variance coveriance matrix is the set of 
    dot products.
*/
  for ( kk = 1 ; kk <= lnpk->t ; kk++ ) {
    scopy(n,lnpk->y[kk],1,lnpk->wy,1);
    info = sqrsl(lnpk->xx,n,n,*nrows,*nrows,lnpk->qraux,lnpk->wy,tmp,lnpk->rsd[kk],tmp,lnpk->rsd[kk],tmp,10L);
  }
  for ( kk = 1 ; kk <= lnpk->t ; kk++ ) 
    for (ii = 1; ii <= kk; ii++)  
      lnpk->s2[ii][kk] = lnpk->s2[kk][ii] = sdot(n,lnpk->rsd[kk],1,lnpk->rsd[ii],1)/ (FPN) n;
  ts0 = invdet(lnpk->s2,lnpk->s2i,lnpk->t,lnpk->t,lnpk->work,lnpk->kpvt,11); 
  return(ts0);
}



/*
  Equation numbers are those in Z.-B. Zeng (1994) Genetics 136:1457-1468.

lnpk->estimates  ests is a pointer that will hold the estimates.
lnpk->xw    xw is the output of sqrdc
            y  is the vector of phenotypes
            ldx = nn is the sample size
            pp = kk is the number of rows in xw
lnpk->rsd   rsd is the residual vector from the null hypothesis, the result of sqrdc and sqrsl
lnpk->qraux qraux is aux. info from sqrdc
            gval is the conditional probabilities of the QTL given flanking markers
            marks is the matrix of flanking markers
            table is a way to determine which value of gval should be used, based on marks
            hypo is the hypothesis to be tested
            cross is the type of cross

lnpk is workspace.  

The convention for lnpk->estimates[trait][1-9]:

        ln(Li).........    Estimates
Col       1    2      3   4   5   6   7     
          H1   H2    H3  a1  a3  d2  d3   
Also, the joint parameters are in the 0th row.


The convention for lnpk->estimate[trait][1-9]:

Col       1    2    3   4   5   6   7     8     9
hypo
1        H1             a               H1:H0
30                  H3      a3      d3  H3:H0
31       H1  H1:H0  H3  a1  a3      d3  H3:H0  H3:H1
32     H2:H0   H2   H3      a3  d2  d3  H3:H0  H3:H2

14       H1   Hge       a               H1:H0  H1:Hge
34            Hge   H3      a3      d3  H3:H0  H3:Hge

*/
int ecm_solve(linpakws *lnpk,int nn, int pp, int ihypo,params *theparams)
{
  FPN *tmp;
  FPN test_stat,test_stat0, s2,be,de;
  int ktime,  go_on, info,cross,bailout,aa,dd,trait,col;
  bailout = 0;  /* Just an indicator if one should bail out...0 means no, 1 means yes*/
/*
  With all the new experimental designs, we have to modify this section
  to estimate the effects correctly.  Most notably, the Design III is probably
  not correct.     Otherwise, T(XX)SFx seems to be ok (it's like a backcross).
*/
  aa = 5;
  dd = 7;  
  if ( ihypo == 1 ) {
    col=1;
    aa = 4;
  }
  else if ( ihypo == 2 ) {
    dd = 6;
    col=2;
  }
  else if ( ihypo == 3 ) 
    col=3;
  else if ( ihypo == 4 )
    col = 2;
  else if ( ihypo == -101 || ihypo == -110 ) 
    col=1;
  else if ( ihypo == -301 || ihypo == -310 ) 
    col=3;
  
  
  
  
  
  
  cross = theparams->cross;
  if ( theparams->tcross ==1 || theparams->tcross == 2 ) 
    cross = theparams->tcross;   /* If this is a test cross to P1 or P2, treat it like a backcross.*/
  if ( theparams->tcross == 12 ) /* otherwise, treat it like an SF2 */
    cross = 3;
  if ( cross == 4 )              /* RFx is like SFx as far as the estimates go. */
    cross = 3;
  tmp = NULL;
  scopy(nn, lnpk->pp1, 1, lnpk->pv, 1);  /* copy a priori probabilities into ...*/
  scopy(nn, lnpk->pp2, 1, lnpk->qv, 1);  /* ... a posteriori probs. initially */
  for ( trait = 1 ; trait <= lnpk->t ; trait++ )
    scopy(nn, lnpk->rsd[trait], 1, lnpk->wrsd[trait], 1); /* copy residuals into a residual workspace */
  ktime = 0;
  go_on = 1;
  de = be = (FPN) 0.0;
  while (go_on == 1) {
/*  Do the M step first...*/
	bailout = AssignGeneticParameters(lnpk,nn,cross,ihypo,theparams,aa,dd);
/*  Is this where we reset the genetic parameters to their weighted averages if we want to do GxE for design I? 
	What is s2i if this is the first iteration?   */	
    if ( ihypo == 4 ) 
      GxEParameters(lnpk, cross, theparams,aa,dd);

    for ( trait = 1 ; trait <= lnpk->t ; trait++ ) {  /*Individual trait mapping comes first.*/ 
      DoTheTrick(lnpk,nn,cross,ihypo,theparams,aa,dd,trait);
      info = sqrsl(lnpk->xx,nn,nn,pp,pp,lnpk->qraux,lnpk->wy,tmp,lnpk->wrsd[trait],tmp,lnpk->wrsd[trait],tmp,10L);
      UnDoTheTrick(lnpk,nn,cross,ihypo,theparams,aa,dd,trait);
      bailout = CalculateTraitlnL(lnpk,nn,cross,ihypo,theparams,aa,dd,trait,&go_on);
    }   
    CreateVarCovar(lnpk,nn,cross,ihypo,theparams,aa,dd);             /*Calculate the Variance-Covariance Matrix*/
    s2 = invdet(lnpk->s2,lnpk->s2i,lnpk->t,lnpk->t,lnpk->work,lnpk->kpvt,11);     /*  Calculate V-1 and det(V) */ 
    test_stat0 = lnpk->estimates[0][col];               /*  Calculate the log of L1 for the Joint hypothesis */ 
    if ( s2 == (FPN) 0.0 )
      bailout = 1;  /* if s2 is 0, something's awry and we should bail */
    else           
      lnpk->estimates[0][col] = - (FPN) 0.5 *  (FPN) nn *  (FPN) log(fabs(s2)) ;
    if ( bailout == 0 )  /* Calculate the joint mapping ln Likelihood*/
      CalculateJointlnL(lnpk,nn,cross,ihypo,theparams,aa,dd);
    else if ( bailout == 1 ) {
      if ( theparams->verbosity == 1 )
        printf("\nBailing out of ECM algorithm...");
      go_on = 0;
    }
    if ( ktime >= (int) M_TIME ) /*  Do we go on or not? */
      go_on = 0;
    else if ( ktime >= 2 && bailout == 0 ) {  /* test_stat is the joint ln(L1) */
      if ( test_stat0 > (FPN) 1.0 &&  (FPN) fabs((lnpk->estimates[0][col]-test_stat0) / test_stat0) < (FPN) STOP_EM  ) 
	    go_on = 0;
      else if ( (FPN) fabs(lnpk->estimates[0][col] - test_stat0) < (FPN) STOP_EM )  
	    go_on = 0;
    }
    if (go_on == 1) /*  Here is where the E step is performed. */
      bailout = DoEstep(lnpk,nn,cross,ihypo,theparams,aa,dd,&ktime,s2);
  } /* End of ECM loop */
  if ( bailout == 1 ) {  /*If the ECM algorithm failed, set all estimates to SIGNOFDEVIL */
    for ( trait = 0 ; trait <= lnpk->t ; trait++ ) 
	  test_stat = lnpk->estimates[trait][aa] = lnpk->estimates[trait][dd] = (FPN) SIGNOFDEVIL;
  }
  else if ( cross == 2 ) { /* we've estimated  -b for B2, thus multiply by -1 */
    for ( trait = 0 ; trait <= lnpk->t ; trait++ ) 
	   lnpk->estimates[trait][aa]   = - lnpk->estimates[trait][aa];
  }
  else if ( cross == 5 ) { /* we've estimated  2b for RI, thus multiply by 0.5 */
    for ( trait = 0 ; trait <= lnpk->t ; trait++ ) 
	   lnpk->estimates[trait][aa]   = (FPN) 0.5 * lnpk->estimates[trait][aa];
  }
  return(bailout);
}

/*
  These are equations 6 and 7 from Jiang and Zeng, Genetics 140:1111-1127 
*/
int AssignGeneticParameters(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd) {
  FPN cp,cq,ctemp1,ctemp2 ;
  int bailout,trait,i;
  i=theparams->traits;
    bailout = 0;
    cp = ssum(nn, lnpk->pv, 1);	            /* cp = q2.1  is the Sum of the elements of pv */
    if (cross == 3) 
      cq = ssum(nn, lnpk->qv, 1);          	/* cq = q1.1  is the Sum of the elements of qv */
    for ( trait = 1 ; trait <= lnpk->t ; trait++ ) {   
      ctemp1 = sdot(nn, lnpk->wrsd[trait], 1, lnpk->pv, 1);	/*ctemp1 = q2.(Y-XB) Numerator of Eqn 6 in Jiang and Zeng 95 */

      if (cross == 3) {
        ctemp2 = sdot(nn, lnpk->wrsd[trait], 1, lnpk->qv, 1);	/* ctemp2 = q1.(Y-XB)  */

        if (ihypo == 2 )
	      lnpk->estimates[trait][aa]  = (FPN) 0.0;
        else if ( cp != (FPN) 0.0 )  /* b = q2.(Y-XB) / (2 q2.1) */
	      lnpk->estimates[trait][aa] = ctemp1 / ((FPN) 2.0 * cp);
	    else /* if the denominator is 0, something's awry. */
	      bailout = 1; 


        if (ihypo == 1  )
	      lnpk->estimates[trait][dd] = (FPN) 0.0;
        else if ( cq != (FPN) 0.0 )  /* d = q1.(Y-XB) / (q1.1) - b*/
	      lnpk->estimates[trait][dd] = ctemp2 / cq - lnpk->estimates[trait][aa];
	    else /* if the denominator is 0, something's awry. */
	      bailout = 1; 
      }
       /* else if ( cross == 5 && cp != 0.0 )       Eqn 6    b = q2.(Y-XB)/(2 q2.1) 
        be = lnpk->estimates[trait][aa] = ctemp1 / ((FPN) 2.0 * cp);	 */
      else if (   cp != (FPN) 0.0 )  	/* Eqn 6    b = q2.(Y-XB)/(q2.1)      is this ok for RI lines?  */
        lnpk->estimates[trait][aa] = ctemp1 / cp;
	  else /* if the denominator is 0, something's awry. */
	    bailout = 1;


      if ( ihypo < 0 ) {
        if ( trait == 1 && ( ihypo == -101 || ihypo == -301 ) )
          lnpk->estimates[1][dd]  = lnpk->estimates[1][aa]  = 0.0;
        if ( trait == 2 && ( ihypo == -110 || ihypo == -310 ) )
          lnpk->estimates[2][dd]  = lnpk->estimates[2][aa]  = 0.0;
      
      
      }



	} 
	return(bailout);
}

/*
  These are equations 34 and 35 from Jiang and Zeng, Genetics 140:1111-1127 
*/
void GxEParameters(linpakws *lnpk,  int cross,  params *theparams,int aa,int dd) {
  FPN cp,cq,ctemp1;
  int trait,i;
  i=theparams->traits;
    lnpk->estimates[0][aa] = lnpk->estimates[0][dd]  = (FPN) 0.0;
    cp = cq = (FPN) 0.0;
/* Do we need to recreate s2i? */
    for ( trait=1; trait<=lnpk->t; trait++ )
      lnpk->work[trait] = ssum(lnpk->t, lnpk->s2i[trait], 1);
    ctemp1 = ssum(lnpk->t,lnpk->work,1); 

    for ( trait=1; trait<=lnpk->t; trait++ )
      cp = cp + lnpk->estimates[trait][aa] * lnpk->work[trait];
    lnpk->estimates[0][aa] = cp / ctemp1;  /* weighted average of the bi */
    if ( cross == 3 ) {
      for ( trait=1; trait<=lnpk->t; trait++ )
        cq = cq + lnpk->estimates[trait][dd] * lnpk->work[trait];
       lnpk->estimates[0][dd]  = cq/ctemp1;  /* weighted average of the di */
    }    
}

/*
This essentially does equation 14 (middle) from Jiang and Zeng, Genetics 140:1111-1127 for various
experimental designs.
*/
void DoTheTrick(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int trait) {
  FPN be,de;
  int row;
  row = theparams->traits;
  if ( ihypo == 4 )
    row = 0;
  else
    row = trait;
      be = lnpk->estimates[row][aa]; 
      de = lnpk->estimates[row][dd];
      if (cross == 3 && ihypo == 2 )
        be = (FPN) 0.0;
      if (cross == 3 && ihypo == 1 )
        de = (FPN) 0.0;
      scopy(nn,lnpk->y[trait],1,lnpk->wy,1);
      if ( cross == 1 || cross == 2 || cross == 5) /*   Use (Y - XB -  b.q2) ... */
        saxpy(nn,-be,lnpk->pv ,1,lnpk->wy ,1);
      else if ( cross == 3 ) { /*   (Y - (2 pv + qv) b  -  qv d - XB) = (Y - XB - 2 b.q2 - (b+d).q1)  */
        saxpy(nn,-(FPN) 2.0 * be,lnpk->pv ,1,lnpk->wy ,1);
        saxpy(nn,-(be+de),lnpk->qv ,1,lnpk->wy ,1);
      }
}


/*
This essentially does equation 15 (right hand side) from Jiang and Zeng, Genetics 140:1111-1127 for various
experimental designs.
*/
void UnDoTheTrick(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int trait) {
  FPN be,de;
  int row;
  row = theparams->traits;
  if ( ihypo == 4 )
    row = 0;
  else
    row = trait;
      be = lnpk->estimates[row][aa]; 
      de = lnpk->estimates[row][dd];
      if ( cross == 1 || cross == 2 || cross == 5) /*... then add back  b.q2 */ 
        saxpy(nn,be,lnpk->pv,1,lnpk->wrsd[trait],1);
      else if ( cross == 3 ) { /*... then add back            2 b.q2 + (b+d).q1     */
        saxpy(nn,(FPN) 2.0 * be,lnpk->pv ,1,lnpk->wrsd[trait],1);
        saxpy(nn,(be+de),lnpk->qv ,1,lnpk->wrsd[trait],1);
      }
}

/*
  This calculates the individual trait lnL
*/
int CalculateTraitlnL(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int trait,int *go_on) {
  FPN rssp,s2,be,de,temp,cp,cq;
  int bailout,ii,row,col;
  bailout = 0;
  if ( ihypo == 4 ) {
    row = 0;
    col = 2;
  }
  
  else {
    row = trait;
    if ( ihypo < -300 )
      col = 3;
    else if ( ihypo < -100 )
      col = 1;
    else
      col = ihypo;
  }
    cp = ssum(nn, lnpk->pv, 1);	            /* cp = q2.1  is the Sum of the elements of pv */
    if (cross == 3) 
      cq = ssum(nn, lnpk->qv, 1);          	/* cq = q1.1  is the Sum of the elements of qv */

      be = lnpk->estimates[row][aa]; 
      de = lnpk->estimates[row][dd];
	  /* for each trait, calculate the log of L1 */
      rssp = sdot(nn, lnpk->wrsd[trait], 1, lnpk->wrsd[trait], 1); /* (Y-XB).(Y-XB)    */ 
      if (cross == 3) /*  [(Y-XB).(Y-XB) - 4 (q2.1) (b'.b) - (q1.1)(b+d)'.(b+d)]/n  */ 
        s2 = (rssp - (FPN)4.0 * cp * be * be - cq * (be + de) * (be + de)) / (FPN) nn;
      else            /*  [(Y-XB).(Y-XB) - (q2.1) (b'.b)]/n  */ 
        s2 = (rssp - cp * be * be) / (FPN) nn;	/* This is equation 9 */

      if (s2 > (FPN) 0.0)	 
        lnpk->estimates[trait][col] = -((FPN) nn *  (FPN) log(s2)  + rssp / s2) / (FPN) 2.0;
      else if ( s2 < (FPN) 0.0 ) 
        lnpk->estimates[trait][col] =  -((FPN) nn *   (FPN) log(-s2)  + rssp / s2) / (FPN) 2.0;
      else /* if s2 is 0, something's awry. */
        bailout = 1;
      if ( bailout == 0 ) {
        for (ii = 1; ii <= nn; ii++) {
          if (cross == 3) {
	        temp =  lnpk->pp1[ii] * (FPN)exp(2.0 * be * (lnpk->wrsd[trait][ii] - be) / s2);
	        temp = temp + lnpk->pp2[ii] * (FPN)exp((be + de) * (lnpk->wrsd[trait][ii] - (be + de) / (FPN) 2.0) / s2);
	        temp = temp + (FPN) 1.0 - lnpk->pp1[ii] - lnpk->pp2[ii];
          }
          else if (cross == 5) /*   */
	        temp = lnpk->pp1[ii] * (FPN)exp(be * (2.0 * lnpk->wrsd[trait][ii] - be) / (2.0 * s2)) + (FPN) 1.0-lnpk->pp1[ii];
          else
	        temp = lnpk->pp1[ii] * (FPN)exp(be * (2.0 * lnpk->wrsd[trait][ii] - be) / (2.0 * s2)) + lnpk->pp2[ii];

          if (temp > (FPN) 0.0)	/*********************Something to check on************************************/
	        lnpk->estimates[trait][col]  = lnpk->estimates[trait][col]  + (FPN) log(temp);
        }
      }
      else if ( bailout == 1 ) {
        if ( theparams->verbosity == 1 )
          printf("\nBailing out of ECM algorithm doing single trait Likelihood...");
        *go_on = 0;
      }
  return(bailout);
}

/*
  Create the variance/covariance matrix
*/
void  CreateVarCovar(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd) {
  FPN be1,be2,de1,de2,cp,cq;
  int trait,trait2,row;
  row = theparams->traits;
    cp = ssum(nn, lnpk->pv, 1);	            /* cp = q2.1  is the Sum of the elements of pv */
    if (cross == 3) 
      cq = ssum(nn, lnpk->qv, 1);          	/* cq = q1.1  is the Sum of the elements of qv */
    for ( trait = 1 ; trait <= lnpk->t ; trait++ )  /*  Create V (variance covariance matrix)*/
      for ( trait2 = trait ; trait2 <= lnpk->t ; trait2++ ) { /* This is equation 9 */
        if ( ihypo == 4 ) 
          row = 0;
        else
          row = trait;
        if (ihypo == 2)
	      be1 = (FPN) 0.0;
        else 
	      be1 = lnpk->estimates[row][aa];
        if (ihypo == 1)
	      de1 = (FPN) 0.0;
        else  
	      de1 = lnpk->estimates[row][dd];
        if ( ihypo == 4 ) 
          row = 0;
        else
          row = trait2;
        if (ihypo == 2)
	      be2 = (FPN) 0.0;
        else 
	      be2 = lnpk->estimates[row][aa];
        if (ihypo == 1)
	      de2 = (FPN) 0.0;
        else  
	      de2 = lnpk->estimates[row][dd];
        if (cross == 3)  
          lnpk->s2[trait2][trait] = lnpk->s2[trait][trait2] = (sdot(nn, lnpk->wrsd[trait], 1, lnpk->wrsd[trait2], 1) - (FPN)4.0 * cp * be1 *  be2 - cq * (be1 + de1) * (be2 + de2)) / (FPN) nn;
        else 
          lnpk->s2[trait2][trait] = lnpk->s2[trait][trait2] = (sdot(nn, lnpk->wrsd[trait], 1, lnpk->wrsd[trait2], 1) - cp * be1 * be2) / (FPN) nn;	  
      }
}


/*  

This section should be the last sum in Jiang and Zeng (1995), eqn 10.  
----------------------------------------------------------------------------------         
----------------------------------------------------------------------------------         
First line of Equation (10) from Jiang and Zeng:         
----------------------------------------------------------------------------------         
SFx,RFx    ln(L1) = k - (n/2)ln(|V|) + sum_j ln{  
T(SFy)SFx           p2j exp[-0.5 (yj-2b*-xj B).V-1.(yj-2b*-xj B)']  +
                    p1j exp[-0.5 (yj-b*-d*-xj B).V-1.(yj-b*-d*-xj B)']
                    p0j exp[-0.5 (yj-xj B).V-1.(yj-xj B)']   }
----------------------------------------------------------------------------------         
Modifications for other crosses:
----------------------------------------------------------------------------------         
B1, B1x    ln(L1) = k - (n/2)ln(|V|) + sum_j ln{  
T(B1){S||R}Fx       p2j exp[-0.5 (yj-b*-xj B).V-1.(yj-b*-xj B)']
                    p1j exp[-0.5 (yj-xj B).V-1.(yj-xj B)']   }
B2,B2x     ln(L1) = k - (n/2)ln(|V|) + sum_j ln{  
T(B2){S||R}Fx       p1j exp[-0.5 (yj-b*-xj B).V-1.(yj-b*-xj B)']
                    p0j exp[-0.5 (yj-xj B).V-1.(yj-xj B)']   }
RIx        ln(L1) = k - (n/2)ln(|V|) + sum_j ln{  
                    p2j exp[-0.5 (yj-2b*-xj B).V-1.(yj-2b*-xj B)']  +      !!!! are the 2b*'s right?  
                    p0j exp[-0.5 (yj-xj B).V-1.(yj-xj B)']   }        
----------------------------------------------------------------------------------         
----------------------------------------------------------------------------------         
Doing it:           Cross
           B1     B2      SFx,RFx     RIx          D3?
  --------------------------------------------------------      
  pAA     p(AA)  0.0       p(AA)      p(AA)
  pAa     p(Aa)  p(Aa)     p(Aa)      0.0
  paa      0.0   p(aa)     p(aa)      p(aa)
  --------------------------------------------------------     
       T(B1)SFx T(B2)SFx  T(SFy)SFx
        B1x       B2x
        
        
        
        ihypo = 1, 2, 3, 4, -101, -110, -301, -311
*/
void CalculateJointlnL(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd) {
  FPN tempsum,pAA,pAa,paa,be1,de1,atemp,temp;
  int trait,ii,col,row;
    row = theparams->traits;

  if ( ihypo == 4 ) 
    col = 2;
  else if ( ihypo < -300 )
    col = 3;
  else if ( ihypo < -100 )
    col = 1;
  else 
    col = ihypo;
      tempsum = (FPN) 0.0;
      for (ii = 1; ii <= nn; ii++) {  
        pAA = lnpk->pp1[ii];
        pAa = lnpk->pp2[ii];
        paa = (FPN) 1.0 - pAA - pAa;
        if ( cross == 1 || cross == 2 || cross == 5 ) {
          pAa = lnpk->pp1[ii];
          paa = (FPN) 1.0 - pAa;
          pAA = (FPN) 0.0;
        }
/*         if ( cross == 2 ) { 
          pAa = lnpk->pp1[ii];
          paa = (FPN) 1.0 - pAa;
          pAA = 0.0;
        }
        else if ( cross == 1   ) {
          paa = pAa;
          pAa = pAA;
          pAA = 0.0;
        }  
        else if ( cross == 5 ) {   Need to use first and third terms? 
 
          pAa = pAA;     
          pAA = 0.0;    Need to use second and third??
          paa = 1-pAA; pAa=0.0; 
        } */
        for ( trait = 1 ; trait <= lnpk->t ; trait++ ) {/* = p2j exp(-0.5 (yj-2b-xjB).V-1.(yj-2b-xjB)' */
          if (ihypo == 2)
	        be1 = (FPN) 0.0;
	      else if ( ihypo == 4 )
	        be1 = lnpk->estimates[0][aa];
          else 
	        be1 = lnpk->estimates[trait][aa];
          lnpk->s2i[0][trait]  = lnpk->wrsd[trait][ii] - 2*be1;    /*  (yj-2b-xjB) in s2i[0] */
        }        
        for ( trait  = 1 ;trait  <= lnpk->t ; trait++ ) 
          lnpk->s2[0][trait] = sdot(lnpk->t,lnpk->s2i[0],1,lnpk->s2i[trait],1);            
        atemp  =  sdot(lnpk->t,lnpk->s2[0],1,lnpk->s2i[0],1);  
        temp =  pAA * (FPN)exp(-0.5 * atemp);
        
        for ( trait = 1 ; trait <= lnpk->t ; trait++ ) {/* = p1j exp(-0.5 (yj-b-d-xjB).V-1.(yj-b-d-xjB)' */
          if (ihypo == 2)
	        be1 = (FPN) 0.0;
	      else if ( ihypo == 4 )
	        be1 = lnpk->estimates[0][aa];	      
          else 
	        be1 = lnpk->estimates[trait][aa];
          if (ihypo == 1)
	        de1 = (FPN) 0.0;
	      else if ( ihypo == 4 )
	        de1 = lnpk->estimates[0][dd];
          else  
	        de1 = lnpk->estimates[trait][dd];
          lnpk->s2i[0][trait]  = lnpk->wrsd[trait][ii] - be1 - de1;   /*  (yj-b-d-xjB) in s2i[0] */
        }        
        for ( trait  = 1 ;trait  <= lnpk->t ; trait++ ) 
          lnpk->s2[0][trait] = sdot(lnpk->t,lnpk->s2i[0],1,lnpk->s2i[trait],1);            
        atemp  =  sdot(lnpk->t,lnpk->s2[0],1,lnpk->s2i[0],1);  
        temp = temp + pAa * (FPN)exp(-0.5 * atemp);
        
        for ( trait = 1 ; trait <= lnpk->t ; trait++ ) /* = p0j exp(-0.5 (yj-xjB).V-1.(yj-xjB)' */
          lnpk->s2i[0][trait]  = lnpk->wrsd[trait][ii] ;   /* row zero of s2i has y-xB*/
        for ( trait  = 1 ;trait  <= lnpk->t ; trait++ )  
          lnpk->s2[0][trait] = sdot(lnpk->t,lnpk->s2i[0],1,lnpk->s2i[trait],1); 
             
        atemp  =  sdot(lnpk->t,lnpk->s2[0],1,lnpk->s2i[0],1);  
        temp = temp + paa * (FPN)exp(-0.5 * atemp);
        if (temp > (FPN) 0.0)	/*********************Something to check on********/
	      tempsum = tempsum + (FPN) log(temp);
      }
	  lnpk->estimates[0][col] = lnpk->estimates[0][col] + tempsum;
}

/*
  This is for equation 5 in Jiang and Zeng, Genetics 140:1111-1127
*/
int DoEstep(linpakws *lnpk,int nn, int cross, int ihypo, params *theparams,int aa,int dd,int *ktime,FPN s2) {
	/* Do the E step */
	FPN const2,be1,de1,tc1,tc2,tc3,ctemp1,ctemp2,ctemp3;
	int ii,trait,bailout,row;
	  row = theparams->traits;

	bailout = 0;
      const2 = (FPN) 1.0 / ((FPN) 2.0 * s2);
      *ktime +=1;
      for (ii = 1; ii <= nn; ii++) {
        for ( trait = 1 ; trait <= lnpk->t ; trait++ )  {/*do tc1 first   (r-b).V-1.(r-b) or (r-2b).V-1.(r-2b)  */
            if (ihypo == 2)
	          be1 = (FPN) 0.0;
            else if ( ihypo == 4 )
	          be1 = lnpk->estimates[0][aa];
            else 
	          be1 = lnpk->estimates[trait][aa];
            if ( cross == 3 )
              lnpk->s2i[0][trait] = lnpk->wrsd[trait][ii]-(FPN) 2.0 * be1 ;
            else
              lnpk->s2i[0][trait] = lnpk->wrsd[trait][ii] - be1 ;
        }
        for ( trait = 1 ; trait <= lnpk->t ; trait++ )
          lnpk->s2[0][trait] = sdot(lnpk->t,lnpk->s2i[trait],1,lnpk->s2i[0],1);
        tc1 = sdot(lnpk->t,lnpk->s2i[0],1,lnpk->s2[0],1);        

        for ( trait = 1 ; trait <= lnpk->t ; trait++ ) /*do tc2 second  r.V-1.r , which goes into tc3 if cross == 3 */
          lnpk->s2i[0][trait] =  lnpk->wrsd[trait][ii];
        for ( trait = 1 ; trait <= lnpk->t ; trait++ )
          lnpk->s2[0][trait] = sdot(lnpk->t,lnpk->s2i[trait],1,lnpk->s2i[0],1);
        tc2 = sdot(lnpk->t,lnpk->s2i[0],1,lnpk->s2[0],1);

	    if (cross == 3) {   /*do tc2 last  r.V-1.r , which goes into tc3 if cross == 3 */
	      tc3 = tc2;
          for ( trait = 1 ; trait <= lnpk->t ; trait++ ) {
            if (ihypo == 2)
	          be1 = (FPN) 0.0;
            else if ( ihypo == 4 )
	          be1 = lnpk->estimates[0][aa];
            else 
	          be1 = lnpk->estimates[trait][aa];
            if (ihypo == 1)
	          de1 = (FPN) 0.0;
            else if ( ihypo == 4 )
	          de1 = lnpk->estimates[0][dd];
            else  
	          de1 = lnpk->estimates[trait][dd];
            lnpk->s2i[0][trait] =  lnpk->wrsd[trait][ii] - be1 - de1;
          }
          for ( trait = 1 ; trait <= lnpk->t ; trait++ )
            lnpk->s2[0][trait] = sdot(lnpk->t,lnpk->s2i[trait],1,lnpk->s2i[0],1);
          tc2 = sdot(lnpk->t,lnpk->s2i[0],1,lnpk->s2[0],1);        
        }
        else
          tc3 = (FPN) 0.0;
	    
	    ctemp1 = lnpk->pp1[ii] * (FPN)exp(-tc1/2.0);
	    if ( cross == 5 )
	      ctemp2 = ((FPN) 1.0-lnpk->pp1[ii]) * (FPN)exp(-tc2/2.0);
	    else
	      ctemp2 = lnpk->pp2[ii] * (FPN)exp(-tc2/2.0);
	    if (cross == 3)
	      ctemp3 = ((FPN) 1.0 -lnpk->pp1[ii] - lnpk->pp2[ii]) * (FPN)exp(-tc3/2.0);
	    else
          ctemp3 = (FPN) 0.0;
        ctemp3 = ctemp1 + ctemp2 + ctemp3;
        if ( ctemp3 == (FPN) 0.0 )
          bailout = 1;
        else {
	      lnpk->pv[ii] = ctemp1 /  ctemp3;
	      if (cross == 3)
	       lnpk->qv[ii] = ctemp2 /  ctemp3;
	    }
      }
   return(bailout);
}

  
/*
  This does the scan in joint mapping
*/
void do_scan(params *theparams,markermap *themap, aqtl *theqtls,individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr,int **srranks) {   /*  Scan the genome */
  genome   *startptr,*endptr ;
  int go_on,do_analysis;
  FPN dumm;
    write_jzheader(theparams->zfile, theparams, progname, chptr,themap,1,first,srranks);
  /*  Here's where we do the analysis.  */
    go_on = determine_endpoints(theparams,themap,&startptr,&endptr,first,&do_analysis);
    if (do_analysis == 1) {
	  dumm = do_jzanalysis(theparams,startptr,endptr,themap,theqtls,individs,theparams->zfile,first,lnpk );
	  write_jzheader(theparams->zfile, theparams, progname, chptr,themap,0,first,srranks);
    }
} 



/*
  this will test for Genotype by Environment interactions where the same set of genotypes are raised in different 
  environments.  
  
  
  There are t traits and m sites.   We will do m * c(t,2) tests.  For each site,
  we will do tests for all pairs of traits.
*/
void do_gbye_same(params *theparams,markermap *themap, aqtl *theqtls,individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr,int **srranks) {   /*  Scan the genome */
  genome   *startptr,*endptr ;
  int go_on, do_analysis;
  FPN dumm;
    write_jzheader(theparams->zfile, theparams, progname, chptr,themap,1,first,srranks);
  /*  Here's where we do the analysis.  */
    go_on = determine_endpoints(theparams,themap,&startptr,&endptr,first,&do_analysis);
    if (do_analysis == 1) {
	  dumm = do_jzanalysis(theparams,startptr,endptr,themap,theqtls,individs,theparams->zfile,first,lnpk );
	  write_jzheader(theparams->zfile, theparams, progname, chptr,themap,0,first,srranks);
    }
} 

/*
  this will test for Genotype by Environment interactions where random genotypes are raised in different 
  environments.  
  
  There are t traits, m sites and a categorical trait with 2 classes.  There will thus
  be m*t tests.      
*/
void do_gbye_diff(params *theparams,markermap *themap, aqtl *theqtls,individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr,int **srranks) {   /*  Scan the genome */
  genome   *startptr,*endptr ;
  int go_on, do_analysis;
  FPN dumm;
    if ( themap->otraits > 0 )
      process_otraits(individs,theparams->nn,themap);
    else
      return;
      
    write_jzheader(theparams->zfile, theparams, progname, chptr,themap,1,first,srranks);
  /*  Here's where we do the analysis.  */
    go_on = determine_endpoints(theparams,themap,&startptr,&endptr,first,&do_analysis);
    if (do_analysis == 1) {
	  dumm = do_jzanalysis(theparams,startptr,endptr,themap,theqtls,individs,theparams->zfile,first,lnpk );
	  write_jzheader(theparams->zfile, theparams, progname, chptr,themap,0,first,srranks);
    }
} 


/*
  This is a driver for actually the scanning of pleiotropy.   
  It will do two traits at a time.  For t traits, that means t(t-1)/2 sets of analyses.
  
  The output will be in stem.zp.   
   
*/
void do_pleiotropy(params *theparams,linpakws *lnpk,char *progname,char *chptr,int **srranks)
{
  FILE *outf;
  FPN  maxe,tmax,ts0;
  int  ipos,nrows,ii,jj, t1, t2, *traitstates,wt,doit,savet,elrotate;
  genome *gptr,*iptr,*jptr,*kptr,*tptr;
  char *outfile;
  
  outfile=cvector(0,MAXLINE);
  sprintf(outfile,"%sp",theparams->zfile);

  print_head(progname,outfile,chptr,1,68,theparams);
  wt = theparams->whichtrait;
  savet = lnpk->t;
  /*  Each trait name will start with an X, so that we can insert/delete a + sign.*/
  traitstates = ivector(1,theparams->traits);
  for ( ii=1; ii<=theparams->traits; ii++ ) {  /* store the original values. */
      if (  theparams->themap->tnames[ii][0] == '+' ) 
        traitstates[ii] = 1;
      else if  (  theparams->themap->tnames[ii][0] == '-' ) 
        traitstates[ii] = -1;
      else {
        traitstates[ii] = 0;
        for ( jj= (int) MAXNAME ; jj > 0 ; jj-- ) /* shift the name for the X */
          theparams->themap->tnames[ii][jj] = theparams->themap->tnames[ii][jj-1];
      }
        theparams->themap->tnames[ii][0] = 'X';
  }
  /* Do a scan for each pair of traits. */
  theparams->whichtrait = 0; 
  for (t1=1; t1<=theparams->traits; t1++) 
    for ( t2=t1+1; t2<=theparams->traits; t2++ )  {
      doit = 0;  /*  Decide whether this pair of traits qualifies for analysis. */
      if (  wt < 1 && traitstates[t1] == 1 && traitstates[t2] == 1 ) 
        doit = 1;
      if ( wt > theparams->traits && traitstates[t1] != -1 && traitstates[t2] != -1)
        doit = 1;
      if ( doit == 1 ) {  /* Yes it does, so change the X's to + signs. */
        theparams->themap->tnames[t1][0] = '+'; 
        theparams->themap->tnames[t2][0] = '+';
        copy_phenotypes(theparams, theparams->themap, theparams->thedata, lnpk);
        lnpk->t = 2;  /* trick it into thinking that there are only two traits.*/
        write_pleiotropyheader(outfile,theparams,theparams->themap,theparams->thegenome,srranks) ;   
        ipos = 0;
        maxe = (FPN) 0.0;
        if ( theparams->verbosity == 1 ) 
          fprintf(stdout,"\nPleiotropic scan of trait %s by %s    ",theparams->themap->tnames[t1],theparams->themap->tnames[t2] );      
        elrotate=1;
        for ( gptr=theparams->thegenome ; gptr != NULL ; gptr = gptr->next ) 
          if ( gptr->dist > (FPN) 0.0 ) {
  	        if ( theparams->verbosity == 1 ) 
  	          elrotate=Rotator(elrotate);
	        pick_cofactors(lnpk,theparams->themap,theparams->theqtls,theparams,gptr,theparams->thegenome);
            ts0 = fit_null(lnpk,theparams,theparams->thedata,theparams->themap,theparams->theqtls,gptr, &nrows);
            lnpk->ts0[0] = -(FPN) 0.5 * (FPN) lnpk->samplesize[1] * ((FPN) log( fabs(ts0))+ (FPN) lnpk->t);
            for ( ii = 1 ; ii <= lnpk->t ; ii++ )
              lnpk->ts0[ii] = -(FPN) 0.5 * (FPN) lnpk->samplesize[1] * ((FPN) log(lnpk->s2[ii][ii])+ (FPN) 1.0 );
            tmax = pleiotropy_int(lnpk,theparams,theparams->themap,theparams->thedata,gptr,outfile,nrows);
            if (tmax > maxe)
              maxe = tmax; 
          }
          if ( theparams->verbosity == 1 ) {
            printf("\n");
            for ( gptr=theparams->thegenome ; gptr != NULL ; gptr = gptr->next ) 
              if (gptr->n == 1 )
                printf("(%d,%d)", gptr->chrom,gptr->markr);
            printf("\n");
          }
          lnpk->t = savet;/* restore the value*/
          outf = fileopen(outfile, "a");
          if (outf != NULL)       
	        fprintf(outf,"\n-e\n#\n#\n" );
     
          fileclose(outfile, outf);
 
        }/* end of loop over genome */
        
/*  This is where pleiotropy v close linkage should go.   */        
        for ( gptr=theparams->thegenome ; gptr != NULL ; gptr = gptr->next ) 
          if ( gptr->n == 1 ) {
            iptr=gptr;
            for (tptr=iptr; tptr->chrom==iptr->chrom && tptr->n == 1 && tptr->next !=NULL; tptr=tptr->next);
            gptr=tptr->next;  /*  we will do an analysis on iptr to tptr.  */
            for ( jptr=iptr; jptr!=tptr->next; jptr=jptr->next )
              for ( kptr=iptr; kptr!=tptr->next; kptr=kptr->next ) 
                tmax = PleiotropyCloseLinkage(lnpk,theparams,theparams->themap,theparams->thedata,jptr,kptr, outfile, nrows);
        
          }

        
        
        theparams->themap->tnames[t1][0] = 'X'; 
        theparams->themap->tnames[t2][0] = 'X';
    }  /*  end of dual loop over t1, t2. */
  for ( ii=1; ii<=theparams->traits; ii++ ) {  /* restore the original values. */
      if (  traitstates[ii] == 1 ) 
        theparams->themap->tnames[ii][0] = '+';
      else if  (  traitstates[ii] == -1 ) 
        theparams->themap->tnames[ii][0] =  '-';
      else  {
        for ( jj=0; jj< (int) MAXNAME ;   jj++ )
          theparams->themap->tnames[ii][jj] = theparams->themap->tnames[ii][jj+1];
        theparams->themap->tnames[ii][(int) MAXNAME] = '\0';
      }
  }
  
  free_ivector(traitstates,1,theparams->traits);
  free_cvector(outfile,0,MAXLINE);
}



/*
  This will calculate expected genotypes at each position and write an outputfile
  formatted for MultiRegress.   
*/
void prep_multiregress(params *theparams,markermap *themap ,individual *individs,genome *first,linpakws *lnpk,char *progname,char *chptr) {   /*  Scan the genome */
  genome   *startptr,*endptr ;
  int go_on,do_analysis,ipos;
  long fileposition;
  /*  Convert data to form needed by MultiRegress*/
    if ( themap->otraits > 0 )
      process_otraits(individs,theparams->nn,themap);

    fileposition = write_altheader(theparams->mrinput, theparams, progname, chptr,themap,1);
    write_alttraits(theparams->mrinput, theparams, themap, individs, lnpk);  
    go_on = determine_endpoints(theparams,themap,&startptr,&endptr,first,&do_analysis);
    ipos = do_jzexpecteds(theparams,startptr,endptr,themap, individs,theparams->mrinput, lnpk);
    insert_positions(theparams->mrinput,ipos,fileposition) ;
}

/*
  Do (composite) interval mapping on the interval pointed to by gptr.


The convention for lnpk->estimate[0,1,2][1-9]:

Col       1       2       3      4   5     6     7     8     9
row
 0                       lnL
 1                                   a1          d1
 2                                   a2          d2


Also, the joint parameters are in the 0th row.
*/
FPN pleiotropy_int(linpakws *lnpk,params *theparams,markermap *themap,individual *individs, genome *gptr, char *outfile,int nrows)
{
  FILE *outf;
  FPN tmax,abspos,relpos,thetaL,thetaR, l11, l10, l01, l0,a1a,a2a,d1a,d2a,a1b,d1b,a2c,d2c,lr1,lr2,lr3;
  int hypo,bail_status1,bail_status2,bail_status3;
  tmax = (FPN) 0.0;
  hypo = 1;
/*  We will go along the entire interval in gptr...*/
    l0 = lnpk->ts0[0];
  gptr->n = 0;
  relpos = (FPN) MIN_DIST;
  while ( relpos < gptr->dist ) {
    l11=l10=l01=0.0;
    abspos = gptr->pos + relpos;
	if (gptr->markr == 0)
	  thetaL = -(FPN) 2.0;
	else
	  thetaL = relpos;
	if ( gptr->markr == themap->mpc[gptr->chrom] )
	  thetaR = -(FPN) 2.0;
	else
	  thetaR = gptr->dist - relpos;
    calc_priors(theparams,lnpk->pp1,lnpk->pp2,theparams->nn,individs,gptr,thetaL,thetaR);
    if ( theparams->ngt == 2 ) {
      bail_status1 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,1,theparams);     /*  Once for joint mapping */
      l11 = lnpk->estimates[0][1];
      lr1 = 2.0*(l11-l0);
      a1a = lnpk->estimates[1][4]; d1a=d2a=0.0;
      a2a = lnpk->estimates[2][4];  /* column 4 has the estimate */
      bail_status2 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-110,theparams);    /*  Hypothesis 10     */
      l10 =  lnpk->estimates[0][1];
      a1b = lnpk->estimates[1][5];  /* column 5 has the estimate */
      d1b = 0.0;
      lr2 = 2.0*(l11-l10);
      bail_status3 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-101,theparams);    /* Hypothesis 01       */
      l01 = lnpk->estimates[0][1];
      a2c = lnpk->estimates[2][5];  /* column 4 has the estimate */
      d2c = 0.0;
      lr3 = 2.0*(l11-l01);
    }
    else {
      bail_status1 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,3,theparams);     /*  Once for joint mapping */
      l11 = lnpk->estimates[0][3];
      lr1 = 2.0*(l11-l0);
      a1a = lnpk->estimates[1][5];  /* columns 5 and 7 have the estimates */
      d1a = lnpk->estimates[1][7];
      a2a = lnpk->estimates[2][5];
      d2a = lnpk->estimates[2][7];
      bail_status2 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-310,theparams);    /*  Hypothesis 10     */
      l10 =  lnpk->estimates[0][3];
      a1b = lnpk->estimates[1][5];
      d1b = lnpk->estimates[1][7];
      lr2 = 2.0*(l11-l10);
      bail_status3 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-301,theparams);    /* Hypothesis 01       */
      l01 = lnpk->estimates[0][3];
      a2c = lnpk->estimates[2][5];
      d2c = lnpk->estimates[2][7];
      lr3 = 2.0*(l11-l01);
    }

    if ( lr1 > theparams->siglevel && lr2 > theparams->siglevel && lr3 > theparams->siglevel )
        gptr->n = 1; /*  if pleiotropy is indicated, then we might test close linkage later. */ 
    
    outf = fileopen(outfile, "a");
    if (outf != NULL) {      
	  fprintf(outf,"\n     %-5d      %-5d",gptr->chrom,gptr->markr );
	  WriteMarkerName(outf,themap,gptr->chrom,gptr->markr);
        
	  fprintf(outf," %10.7f",abspos );
      fprintf(outf," %10.4f %10.4f %10.4f        %10.4f",lr1,lr2,lr3,l11);
      fprintf(outf," %10.7f %10.7f %10.7f %10.7f",a1a,d1a,a2a,d2a);
      fprintf(outf," %10.7f %10.7f %10.7f %10.7f",a1b,d1b,a2c,d2c);
   } 
    fileclose(outfile, outf);
    relpos = relpos + theparams->walk;
  }
  return(tmax);
}


/*

   This interval indicates pleiotropy.  Do a test of pleiotropy vs close linkage for all test sites in this interval
   against all in this interval and following pleiotropic intervals.   For example, if this and the next two intervals
   indicate pleiotropy, and each is 10 cM long, the test sites are marked by x:
   
   |x.x.x.x.x.|x.x.x.x.x.|x.x.x.x.x.|
   
   There are 15 sites.  We will do 14+13+12+11+10 tests

   I am assuming that the calling routine will then advance, and do 9+8+7+6+5 tests
   and then again to do 4+3+2+1, meaning that the calling routine will do 15*14/2 = 105 tests.  

*/
FPN PleiotropyCloseLinkage(linpakws *lnpk,params *theparams,markermap *themap,individual *individs,genome *jptr,genome *gptr,char *outfile,int nrows)
{
  FILE *outf;
  FPN tmax,abspos,relpos,thetaL,thetaR, l11, l10, l01, l0,a1a,a2a,d1a,d2a,a1b,d1b,a2c,d2c,lr1,lr2,lr3;
  int hypo,bail_status1,bail_status2,bail_status3;
  tmax = (FPN) 0.0;
  hypo = 1;
/*  We will go along the entire interval in gptr...*/
  l0 = lnpk->ts0[0];
  relpos = (FPN) MIN_DIST;
  while ( relpos < gptr->dist ) {  /*  Loop in the current interval.  */
    l11=l10=l01=0.0;
    abspos = gptr->pos + relpos;
	if (gptr->markr == 0)
	  thetaL = -(FPN) 2.0;
	else
	  thetaL = relpos;
	if ( gptr->markr == themap->mpc[gptr->chrom] )
	  thetaR = -(FPN) 2.0;
	else
	  thetaR = gptr->dist - relpos;
    calc_priors(theparams,lnpk->pp1,lnpk->pp2,theparams->nn,individs,gptr,thetaL,thetaR);
    if ( theparams->ngt == 2 ) {
      bail_status1 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,1,theparams);     /*  Once for joint mapping */
      l11 = lnpk->estimates[0][1];
      lr1 = 2.0*(l11-l0);
      a1a = lnpk->estimates[1][4]; d1a=d2a=0.0;
      a2a = lnpk->estimates[2][4];  /* column 4 has the estimate */
      bail_status2 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-110,theparams);    /*  Hypothesis 10     */
      l10 =  lnpk->estimates[0][1];
      a1b = lnpk->estimates[1][5];  /* column 5 has the estimate */
      d1b = 0.0;
      lr2 = 2.0*(l11-l10);
      bail_status3 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-101,theparams);    /* Hypothesis 01       */
      l01 = lnpk->estimates[0][1];
      a2c = lnpk->estimates[2][5];  /* column 4 has the estimate */
      d2c = 0.0;
      lr3 = 2.0*(l11-l01);
    }
    else {
      bail_status1 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,3,theparams);     /*  Once for joint mapping */
      l11 = lnpk->estimates[0][3];
      lr1 = 2.0*(l11-l0);
      a1a = lnpk->estimates[1][5];  /* columns 5 and 7 have the estimates */
      d1a = lnpk->estimates[1][7];
      a2a = lnpk->estimates[2][5];
      d2a = lnpk->estimates[2][7];
      bail_status2 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-310,theparams);    /*  Hypothesis 10     */
      l10 =  lnpk->estimates[0][3];
      a1b = lnpk->estimates[1][5];
      d1b = lnpk->estimates[1][7];
      lr2 = 2.0*(l11-l10);
      bail_status3 = ecm_solve(lnpk,lnpk->samplesize[1],nrows,-301,theparams);    /* Hypothesis 01       */
      l01 = lnpk->estimates[0][3];
      a2c = lnpk->estimates[2][5];
      d2c = lnpk->estimates[2][7];
      lr3 = 2.0*(l11-l01);
    }

    if ( lr1 > theparams->siglevel && lr2 > theparams->siglevel && lr3 > theparams->siglevel )
        gptr->n = 1; /*  if pleiotropy is indicated, then we might test close linkage later. */ 
    
    outf = fileopen(outfile, "a");
    if (outf != NULL) {      
	  fprintf(outf,"\n     %-5d      %-5d",gptr->chrom,gptr->markr );
	  WriteMarkerName(outf,themap,gptr->chrom,gptr->markr);
        
	  fprintf(outf," %10.7f",abspos );
      fprintf(outf," %10.4f %10.4f %10.4f        %10.4f",lr1,lr2,lr3,l11);
      fprintf(outf," %10.7f %10.7f %10.7f %10.7f",a1a,d1a,a2a,d2a);
      fprintf(outf," %10.7f %10.7f %10.7f %10.7f",a1b,d1b,a2c,d2c);
   } 
    fileclose(outfile, outf);
    relpos = relpos + theparams->walk;
  }
  return(tmax);
}


void write_pleiotropyheader(char *outfile,params *theparams,markermap *themap,genome *first,int **srranks)
{
  int cross,t,i,wt;
  FILE *outf;
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  t = how_many_traits(theparams,themap); 
  wt = theparams->whichtrait; 
  outf = fileopen(outfile, "a");
  if (outf != NULL) { 
    fprintf(outf, "\n#  The position is from the left telomere on the chromosome");
    fprintf(outf, "\n-window     %6.2f      Window size for model 6",theparams->window);
    fprintf(outf, "\n-background %6d      Max. Rank for Background parameters in model 6",theparams->nbp);
    if ( theparams->Model == 6 ) 
      show_cofactors(outf,theparams,themap,first,srranks) ;
    fprintf(outf, "\n-Model      %6d      Model number",theparams->Model);
    fprintf(outf, "\n-traits     %6d      Number of Traits in Joint Analysis", t);
    fprintf(outf, "\n-Names ");
    for ( i = 1 ; i <= theparams->traits ; i++ ) 
      if (  (wt == 0 && themap->tnames[i][0] == '+') || (wt > theparams->traits && themap->tnames[i][0] != '-') || (wt == i)  ) 
        fprintf(outf, "\n\t %d   %s",i,themap->tnames[i]);
	fprintf(outf, "\n-cross      %6s      Cross", theparams->thecross);
    fprintf(outf, "\n-ihypo      %6d      Hypothesis test",theparams->ihypo);
    fprintf(outf, "\n#\n#  Note that our Likelihood ratio test statistic compares two nested hypotheses");
    fprintf(outf, "\n#  and is two times the negative natural log of the ratio of the likelihoods.  For example,");
    fprintf(outf, "\n#  assume that  hypothesis H0 is nested within H1 and that they have likelihoods L0 and L1 respectively.");
    fprintf(outf, "\n#  Then, the \"Likelihood Ratio Test Statistic\" is -2ln(L1/L0). \n#");
	fprintf(outf, "\n# ------------------------ -------------------    -------- Test Statistics ------                     -----------  H11 ests ------------       -----H10 ests-----    ----  H01 ests ----");
	fprintf(outf, "\n# Chromosome   Marker   MarkerName  Position    LR(L11:L0) LR(L11:L10) LR(L11:L01)      ln(L11)       a1        d1         a2         d2          a1         d1          a2        d2  " );
    fprintf(outf, "\n-s");
    fileclose(outfile,outf);
  }
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MZfunc.c
------------------------------------------------------------------ */


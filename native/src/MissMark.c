/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MissMark.c) is part of QTL Cartographer
         
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

    Functions to handle inferring missing marker data. 
    See:
       C. Jiang and Z.-B. Zeng 1997 Mapping quantitative
       trait loci with dominant and missing markers in 
       various crosses from two inbred lines.   Genetica
       101:47-58.   


*/
int denom0;

void  showa_Mv(FPN *a, FPN (*M)[3], FPN *v);


/*
  If marker jj on chromosome ii in individual kk is missing, calculate it's
  expected value.  Need the probabilities of AA, Aa and aa genotypes.  
*/
FPN expect_mark(params *theparams,individual *individs,int kk,genome *gptr)
{
  FPN evalue,pAA,pAa,paa,vAA,vAa,vaa;
  vAA =  (FPN) 1.0;
  vAa =  (FPN) 0.0;
  vaa = -(FPN) 1.0;
  calc_cond_prob(theparams,gptr,individs,kk,&pAA,&pAa,&paa);
  evalue = vAA*pAA + vAa*pAa + vaa*paa;
  return (evalue);
}

/*
  Add a virtual marker node into the genome chain
*/
void add_virtual_node(genome *tgptr,genome *lgptr,genome *rgptr, FPN thetaL, FPN thetaR)
{
        lgptr->prev = tgptr->prev;
	    if ( tgptr->prev != NULL )
	      lgptr->prev->next = lgptr;
   	    rgptr->prev = lgptr;
	    lgptr->next = rgptr;

	    rgptr->next = tgptr->next;
	    if ( tgptr->next != NULL )
	      rgptr->next->prev = rgptr;

    	lgptr->dist = thetaL;
	    rgptr->dist = thetaR;
	    lgptr->chrom = tgptr->chrom;
	    lgptr->markr = tgptr->markr;
	    rgptr->chrom = tgptr->chrom;
	    rgptr->markr = -10;      
}

/*
Rejoin the chain
*/ 
void del_virtual_node(genome *tgptr)
{
        if ( tgptr->prev != NULL )
  	      tgptr->prev->next = tgptr;
	    if ( tgptr->next != NULL )
	      tgptr->next->prev = tgptr;
}


/* Calculate the prior probabilities of the QTL given the flanking markers and
   map distance to the flanking markers.  */

void calc_priors(params *theparams,FPN *pp1,FPN *pp2,int nn,individual *individs,genome *gptr,FPN thetaL,FPN thetaR)
{
  int ii,kk,go_on;
  FPN pAA,pAa,paa;
  genome   lg,rg;
  int chr,mark,cross;
  cross = theparams->cross;
  if ( theparams->tcross ==1 || theparams->tcross == 2 ) 
    cross = theparams->tcross;
  if ( theparams->tcross == 12 )
    cross = 3;
  if ( cross == 4 )
    cross = 3;
  chr = gptr->chrom;
  mark = gptr->markr;
  go_on = 1;
  add_virtual_node(gptr,&lg ,&rg,thetaL,thetaR);  
  kk = 0;
  for (ii = 1; ii <= nn; ii++)
    if ( individs[ii].print_flag == 'y' ) {
      kk = kk + 1;
	  calc_cond_prob(theparams,&rg,individs,ii,&pAA,&pAa,&paa);
      
      pp1[kk] = pAA;
      pp2[kk] = pAa;
      if (  cross == 2 )
        pp1[kk] = paa;
    }
  del_virtual_node(gptr);
}


/*
Assume that r and  a  are 3x3 square matrices.
Put a into r
*/
void AinR(FPN (*rr)[3], FPN (*aa)[3])
{
  int ii,jj;
  for ( ii = 0 ; ii<=2 ;ii++ )
    for ( jj = 0 ; jj<=2 ; jj++ )
      rr[ii][jj] = aa[ii][jj];
}  

/*
Assume that r, a, and b are 3x3 square matrices.
Calculate r = a.b
*/
void AdotB(FPN (*rr)[3], FPN (*aa)[3], FPN (*bb)[3])
{
  int ii,jj,kk;
  for ( ii =0  ; ii<=2  ;ii++ )
    for ( jj = 0 ; jj<=2  ; jj++ ) 
      for ( kk = 0  ; kk <= 2 ; kk++ )
        rr[ii][jj] = rr[ii][jj] + aa[ii][kk]*bb[kk][jj];
}  


/*
   mm is the matrix that holds the product 
   cross is the cross and = 1, 2, or 3
   rec is the recombination fraction of the interval
   rl is an indicator of whether we are to the right (1) or left (0) 
   of the missing marker.
   
   
   This comes from Jiang and Zeng's paper on using Markov chains to 
   get the probablities of different genotypes.
   
   
   This calculates the matrix  Hx.Ix, where Hx is the transition probability
   matix in the Markov chain, and Ix is a matrix to 'zero-out' some columns,
   depending on the marker.
*/

void markov_rec(params *theparams, FPN (*mm)[3], FPN rec, int rl)
{
  FPN mrec[3][3],tot[3][3],sftpm[3][3],trec,twrec;
  int i,ii,jj;
  for ( ii = 0 ; ii<3 ;ii++ )
    for ( jj = 0 ; jj<3 ; jj++ )
      tot[ii][jj] = mrec[ii][jj] = (FPN) 0.0;
      
  if ( theparams->cross == 1 ) { /* B1 */
    mrec[1][1] = mrec[0][0] = (FPN) 1.0 - rec;
    mrec[0][1] = mrec[1][0] = rec;    
  }
  else if ( theparams->cross == 2 ) { /* B2 */
    mrec[1][1] = mrec[2][2] = (FPN) 1.0 - rec;
    mrec[2][1] = mrec[1][2] = rec;
  }
  else if ( theparams->cross == 3 ) { /* SFx or SFx crossed to P1 or P2 or Design III*/
    if ( theparams->tcross == 0 || theparams->tcross == 12 || theparams->tcrosst > theparams->crosst )
      selfed_f_tpm(mrec,rec,theparams->crosst);
    else if (theparams->tcross == 1) /* test cross to P1 */
      bselfed_f_tpm(mrec,rec,theparams->crosst,1);
    else if (theparams->tcross == 2) /* test cross to P2 */
      bselfed_f_tpm(mrec,rec,theparams->crosst,2);

/*  If there were more than one generations of backcrossing from the SFx, then
     we want to modify the G matrix for that number.    */
  if ( theparams->tcross < 3 && theparams->tcrosst > 1 ) { 
    if ( theparams->tcross == 1 ) { /* B1 */
      sftpm[1][1] = sftpm[0][0] = (FPN) 1.0 - rec;
      sftpm[0][1] = sftpm[1][0] = rec;    
    }
    else if ( theparams->tcross == 2 ) { /* B2 */
      sftpm[1][1] = sftpm[2][2] = (FPN) 1.0 - rec;
      sftpm[2][1] = sftpm[1][2] = rec;
    }
      AinR(tot,mrec);
      for ( i = 2 ; i<= theparams->crosst ; i++ ) {
        AdotB(mrec,tot,sftpm);
        AinR(tot,mrec);
      }
    for ( ii = 0 ; ii<3 ;ii++ )
      for ( jj = 0 ; jj<3 ; jj++ )
        tot[ii][jj] = sftpm[ii][jj] = (FPN) 0.0;    
        
  }



  }
  else if ( theparams->cross == 4 ) { /* RFx */
    if ( theparams->tcross == 0 ) { /* no test cross */
      if ( theparams->crosst == 2 )
        trec = rec;
      else
        trec = (FPN) 0.5 * ((FPN) 1.0 - (FPN) pow( 1.0-rec, (FPN)(theparams->crosst-2) ) * ((FPN)1.0-(FPN) (FPN) 2.0*rec) );
      twrec = (FPN) 1.0-trec;
      mrec[0][0] = mrec[2][2] = twrec*twrec;
      mrec[2][0] = mrec[0][2] = trec*trec;
      mrec[1][0] = mrec[1][2] = trec*twrec;
      mrec[0][1] = mrec[2][1] = (FPN)  2.0*trec*twrec;
      mrec[1][1] = twrec*twrec + trec*trec;
    }
    else { /*   test cross to P1 or P2  */
      trec = (FPN) 0.5 * ((FPN) 1.0 - (FPN) pow( 1.0-rec, (FPN)(theparams->crosst-1) ) * ((FPN) 1.0-(FPN) 2.0*rec) );
      twrec = (FPN) 1.0-trec;
      mrec[0][0] = mrec[1][1] = mrec[2][2] = (FPN) 1.0-twrec ;
      mrec[2][0] = mrec[0][2] = (FPN) 0.0;
      mrec[1][0] = mrec[1][2] = mrec[0][1] = mrec[2][1] = trec ;
      if ( theparams->tcross == 1 )
        mrec[2][1] = mrec[2][2] = mrec[1][2] = (FPN) 0.0;
      else if ( theparams->tcross == 2 )
        mrec[0][0] = mrec[0][1] = mrec[1][0] = (FPN) 0.0;
    }
  }
  else if ( theparams->cross == 5 ) {   /* RI lines...0=>Trudy, 1=>selfed RI, 2=>sibbed RI */
    mrec[0][0] = mrec[2][2] = (FPN) 1.0 - rec;
    mrec[0][2] = mrec[2][0] = rec;
    if ( theparams->crosst == 1 ) {
      mrec[0][0] = mrec[2][2] = (FPN) 1.0/((FPN) 1.0+(FPN) 2.0*rec);
      mrec[0][2] = mrec[2][0] = (FPN) 2.0*rec/((FPN) 1.0+(FPN) 2.0*rec);
    }
    else if ( theparams->crosst == 2 ) {
      mrec[0][0] = mrec[2][2] = (FPN) 1.0/((FPN) 1.0+(FPN)6.0*rec);
      mrec[0][2] = mrec[2][0] = (FPN) 4.0*rec/((FPN) 1.0+(FPN)6.0*rec);
    }
  }
/*  If there were more than one generations of backcrossing from the F1, then
     we want to power up the G matrix for that number.    
      
             */
  if ( theparams->cross < 3 && theparams->crosst > 1 ) { 
      AinR(tot,mrec);
      AinR(sftpm,mrec);
      for ( i = 2 ; i<= theparams->crosst ; i++ ) {
        AdotB(mrec,tot,sftpm);
        AinR(sftpm,mrec);
      }
    for ( ii = 0 ; ii<3 ;ii++ )
      for ( jj = 0 ; jj<3 ; jj++ )
        tot[ii][jj] = sftpm[ii][jj] = (FPN) 0.0;    
        
  }
  switch (rl) {
    default :
    case -3 : 
      break;
    case -2 : /*a- genotypes indicate that AA are impossible, so 0 first column. */
      for (i=0;i<3;i++) 
        mrec[i][0] = (FPN) 0.0;
      break;
    case -1 : /*aa genotypes indicate that AA, Aa are impossible, so 0 first 2 columns. */
      for (i=0;i<3;i++) 
        mrec[i][0] = mrec[i][1] = (FPN) 0.0;
      break;
    case 0  : /*Aa genotypes indicate that AA, aa are impossible, so 0 first and third columns. */
      for (i=0;i<3;i++) 
        mrec[i][0] = mrec[i][2] = (FPN) 0.0;
      break;
    case 1  : /*AA genotypes indicate that aa, Aa are impossible, so 0 second 2 columns. */
      for (i=0;i<3;i++) 
        mrec[i][2] = mrec[i][1] = (FPN) 0.0;
      break;
    case 2  : /*A- genotypes indicate that aa are impossible, so 0 last column. */
      for (i=0;i<3;i++) 
        mrec[i][2] = (FPN) 0.0;
      break;
  }
  AdotB(tot,mm,mrec);
  AinR(mm,tot);
}

/*
    For some position, calculate the probability distribution of the
    various genotypes conditioned on the flanking markers.
    
    cross = 1, 2, 3, 4, 5 for B1, B2, SFx, RFx, RI
    crosst = 2, 3, ... for SF, RF, 0, 1 for RI
    tcross = 1, 2, 3, 12 for SFx test crossed to P1, P2, SFx+1 or Design III
    tcrosst = 2, 3, ...  for SFX test crossed SFx...
    
*/
void calc_cond_prob(params *theparams,genome *gptr,individual *individs,int kk,FPN *pAA,FPN *pAa,FPN *paa)
{
  FPN pLuk[3],pRuk[3],Muv[3][3],quk[3],mml[3][3],mmr[3][3],rec;
  genome *tgptr;
  int ii,jj,k,lmark,rmark,mark,go_on;
  denom0=0;
  lmark=rmark=-3;
  for ( ii = 0 ; ii < 3 ; ii++ )
    for ( jj = 0 ; jj < 3 ; jj++ )
      Muv[ii][jj]=mml[ii][jj]=mmr[ii][jj] = (FPN) 0.0;
  for ( ii = 0 ; ii < 3 ; ii++ )
    quk[ii]=Muv[ii][ii]=mml[ii][ii]=(FPN) 1.0;

/* For test crosses.  SFx then either Design III or a few more SF generations. */
  if ( theparams->cross == 3 && (theparams->tcross == 3 || theparams->tcross == 12 ) && theparams->tcrosst > theparams->crosst ) {
     Muv[1][1] = (FPN) 1.0/(FPN) pow(2.0, (FPN) (theparams->tcrosst - theparams->crosst));
     Muv[0][1] = Muv[2][1] = (FPN) 0.5*((FPN) 1.0-Muv[1][1]);
     Muv[0][0] = Muv[2][2] = (FPN) 1.0;  /* This is actually the tranpose of the Muv matrix, which is all that is needed.*/
     
    if (  theparams->tcross == 12 ) {
      if (  (individs + kk)->bc == 1 ) {
        mml[1][0] = mml[1][1] = (FPN) 0.5;
        mml[2][2] = (FPN) 0.0;
        mml[2][1] = (FPN) 1.0;
      }
      else if ( (individs + kk)->bc == 2 ) {
        mml[1][2] = mml[1][1] = (FPN) 0.5;
        mml[1][1] = (FPN) 0.0;
        mml[0][1] = (FPN) 1.0;
      }
      for ( ii = 0 ; ii < 3 ; ii++ )
        for ( jj = 0 ; jj < 3 ; jj++ ) 
          for ( k = 0 ;k<3;k++ )
            mmr[ii][jj] = mmr[ii][jj]+Muv[ii][k]*mml[jj][k];
      for ( ii = 0 ; ii < 3 ; ii++ )
        for ( jj = 0 ; jj < 3 ; jj++ ) 
          Muv[ii][jj] = mmr[jj][ii];/* Muv is transposed. Make sure to test this when Design III's come online.*/
    }
    quk[1] = (FPN) 1.0/(FPN) pow(2.0, (FPN) (theparams->crosst - 1));
    quk[0]=quk[2] = (FPN) 0.5*((FPN) 1.0-quk[1]);
  }
  else
    for ( ii = 0 ; ii < 3 ; ii++ )
      quk[ii]=Muv[ii][ii] =(FPN) 1.0;

/*

The Muv matrix is actually the transpose of the Muv matrix.   That is all that is needed.  
For u=v, Muv is the Identity matrix and equal to its tranpose.
*/
  
  for ( ii = 0 ; ii < 3 ; ii++ )
    for ( jj = 0 ; jj < 3 ; jj++ )
       mml[ii][jj]=mmr[ii][jj] = (FPN) 0.0;
  for ( ii = 0 ; ii < 3 ; ii++ )
    mml[ii][ii]=mmr[ii][ii]=(FPN) 1.0;
  *pAA = *pAa = *paa = (FPN) 0.0;
  if ( gptr->markr > 0 )
    mark =  individs[kk].markers[gptr->chrom][gptr->markr];
  else
    mark = -3;
  switch (mark) {
    case -1 :
      *paa = (FPN) 1.0;
      break;
    case 0 :
      *pAa = (FPN) 1.0;
      break;
    case 1 :
      *pAA = (FPN) 1.0;
      break; 
    case -2 :
      mml[0][0] = (FPN) 0.0;
    case 2 :
      if (mark == 2)
        mml[2][2] = (FPN) 0.0; 
    default :
      tgptr = gptr;
      if ( tgptr->prev != NULL ) {
        tgptr = tgptr->prev;
        go_on = 1;
      }
      else 
        go_on = 0;
        
        
      while ( go_on == 1 ) 
        if ( tgptr->chrom == gptr->chrom ) {
          rec = mapfunc( tgptr->dist , -1);
          if (tgptr->markr > 0) /*  CHECK markers matrix ZZZZ   Checks out, value is -2*/
            lmark =  individs[kk].markers[tgptr->chrom][tgptr->markr];
          else if (tgptr->markr < 0)
            lmark = individs[kk].vqtls[tgptr->chrom][-tgptr->markr];
          else  /*  Marker 0 is a virtual marker too.  Should be set to -10 */
            lmark = -10;

          markov_rec( theparams, mml,  rec  , lmark  );
          if ( lmark > -2 && lmark < 2 )
            go_on = 0;
          if ( tgptr->prev != NULL   ) 
            tgptr = tgptr->prev;
          else 
            go_on = 0;
        }
        else
          go_on = 0;
      a_Mv(pLuk,mml,quk);
            
      tgptr = gptr;
      if ( tgptr->next != NULL ) 
        go_on = 1;
       else 
        go_on = 0;
        
        
      while ( go_on == 1 ) 
        if ( tgptr->next->chrom == gptr->chrom ) {
          rec = mapfunc( tgptr->dist , -1);
          if (tgptr->next->markr >= 0)  /*  CHECK markers matrix ZZZZ   Checks out, value is -2 */
            rmark = individs[kk].markers[tgptr->next->chrom][tgptr->next->markr];
          else if (tgptr->markr < 0)
            rmark = individs[kk].vqtls[tgptr->next->chrom][-tgptr->next->markr];

          markov_rec(theparams,mmr,  rec  , rmark  );
          if ( rmark > -2 && rmark < 2 )
            go_on = 0;
          tgptr = tgptr->next;
          if ( tgptr->next == NULL   ) 
            go_on = 0;
        }
        else
          go_on = 0;
      a_Mv(pRuk,mmr,quk); /* pRuk =  Hzk+1.Hzk+2...Hzl.c */
      
      
      assign_q(quk,theparams);
      cond_prob(theparams,Muv,quk,pRuk,pLuk);
      *pAA = quk[0];
      *pAa = quk[1];
      *paa = quk[2];   
      if ( denom0 == 1 ) {
        sprintf(gwarn, "\ndenom0:  -i %3d -g %3d -c %3d -m %3d", kk, mark, gptr->chrom, gptr->markr);
        LogTheError( theparams->error, gwarn);
      }   
      break;
  }
}

/*
  Calc. conditional probability vector
  
  Return the conditional probability vector
  of the AA, Aa and aa genotypes in quk.
*/

void cond_prob(params *theparams, FPN (*Muv)[3], FPN *quk, FPN *pRuk, FPN *pLuk)
{
  FPN pCp[3],denom;

  int i;
  
  for ( i=0;i<3;i++ )
    pCp[i] = pRuk[i]*pLuk[i];
  for ( i=0;i<3;i++ )
    pRuk[i] = quk[i]*pCp[i];
  denom = 0;
  for ( i=0;i<3;i++ )
    denom = denom + quk[i]*pCp[i];
  if (denom == (FPN) 0.0 ) {
    if (theparams->verbosity == 1 )
      printf("\ndenom = 0.0 in cond_prob...");
    denom = (FPN) 1.0;
    denom0 = 1;
  }
  for ( i=0;i<3;i++ )
    pLuk[i] = pRuk[i]/denom;
  a_Mv(quk,Muv,pLuk);

}

/*
   This does a = M.v 
*/
void a_Mv(FPN *a, FPN (*M)[3], FPN *v)
{
int i,j;
  for ( i=0;i<3;i++ ) {
    a[i] = (FPN) 0.0;
    for ( j=0;j<3;j++ ) 
      a[i] = a[i] + M[i][j]*v[j];
  }
}

void showa_Mv(FPN *a, FPN (*M)[3], FPN *v) {
int i;
for ( i=0; i<3; i++ )  
  printf("\n  %6.4f  | %6.4f %6.4f %6.4f | %6.4f", a[i], M[i][0], M[i][1],M[i][2], v[i] );
printf("\n\n");
mypause();
}


/*

This assigns the vector qk' which is written explicitly for an F2
population in equation 3 of Jiang and Zeng's paper on dominanant and
missing markers.

*/
void assign_q(FPN *ql,params *theparams)
{
  int ii;
  for ( ii = 0 ; ii < 3 ; ii++ )
    ql[ii] = (FPN) 0.0;
  if ( theparams->cross == 1 ) {  /*  pr(AA) = 1-pr(Aa);  pr(Aa) = pow(.5,t) for B1t */
    ql[1] = (FPN) pow(0.5, (FPN) theparams->crosst) ;
    ql[0] = (FPN) 1.0 - ql[1];
  }
  else if ( theparams->cross == 2 ) {  /*  pr(aa) = 1-pr(Aa);  pr(Aa) = pow(.5,t) for B2t */
    ql[1] = (FPN) pow(0.5, (FPN) theparams->crosst) ;
    ql[2] = (FPN) 1.0 - ql[1];
  }
  else if ( theparams->cross == 3 ) {  /*  pr(AA) = pr(aa) = .25 for B1 */
    if ( theparams->tcross == 1 )
      ql[0]=ql[1]=(FPN) 0.5;
    else if  ( theparams->tcross == 2 )
      ql[2]=ql[1]=(FPN) 0.5;
    else  {
      ql[1] = (FPN) 1.0/(FPN) pow(2.0, (FPN) (theparams->crosst-1));
      ql[0] = ql[2] = (FPN) 0.5*((FPN) 1.0-ql[1]);
    } 
  }
  else if ( theparams->cross == 4 ) {
    if ( theparams->tcross == 1 )
      ql[0]=ql[1]=(FPN) 0.5;
    else if  ( theparams->tcross == 2 )
      ql[2]=ql[1]=(FPN) 0.5;
    else  {
      ql[0] = ql[2] = (FPN)0.25 ;
      ql[1] = (FPN) 0.5;
    }  
  }
  else if ( theparams->cross == 5 ) 
    ql[0] = ql[2] = (FPN) 0.5 ;
}

/*
  Suppose we are producing Selfed Ft's.  This will calculate the
  transition matrix.
*/
void selfed_f_tpm(FPN (*sftpm)[3], FPN r, int t)
{
  FPN twotm2,twotm1, c1, c2, c3;

  if ( t != 2 )
    twotm2 = (FPN) pow( 2.0, (FPN) t - 2.0 );
  else
    twotm2 = (FPN) 1.0;
  twotm1 = twotm2*(FPN) 2.0;
  c1 = (FPN) 1.0+2*r;
  c2 = (FPN) pow( (0.5-r), (FPN) t ) ;
  c3 = (FPN) pow( 0.5-r*((FPN) 1.0-r) , (FPN) t - (FPN) 1.0 );
  
  sftpm[0][0] = sftpm[2][2] =   (twotm1/c1       - (FPN) 1.0 - twotm1*c2/c1 + twotm2*c3)/(twotm1-(FPN) 1.0);
  sftpm[0][1] = sftpm[2][1] =   (                  (FPN) 1.0 -                twotm1*c3)/(twotm1-(FPN) 1.0);
  sftpm[0][2] = sftpm[2][0] =   ((FPN) 2.0*r*twotm1/c1 - (FPN) 1.0 + twotm1*c2/c1 + twotm2*c3)/(twotm1-(FPN) 1.0);
  sftpm[1][0] = sftpm[1][2] =   (FPN) 0.5 -  twotm2*c3 ;
  sftpm[1][1] =                 twotm1*c3 ;

}

void bselfed_f_tpm(FPN (*sftpm)[3], FPN r, int t, int bc)
{
  FPN  c1, c2;


  c1 = (FPN) pow( (0.5-r), (FPN) t);
  c2 = (FPN) 1.0+(FPN) 2.0*r;
   
  sftpm[0][0] = sftpm[1][1] = sftpm[2][2] =  c1 + ((FPN) 1.0-c1)/c2;
  sftpm[0][1] = sftpm[1][0] = sftpm[0][2] = sftpm[2][1] =  (c2 - (FPN) 1.0 + c1)/c2 - c1;
  if (bc==1) 
    sftpm[0][2] = sftpm[1][2] = sftpm[2][2] = sftpm[2][0] = sftpm[2][1] =  (FPN) 0.0;
  if (bc==2)
    sftpm[0][0] = sftpm[0][1] = sftpm[0][2] = sftpm[1][0] = sftpm[2][0] =  (FPN) 0.0;
}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MissMark.c
------------------------------------------------------------------ */


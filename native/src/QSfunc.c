/* ------------------------------------------------------ XCutXCodeXSkip
     This file (QSfunc.c) is part of QTL Cartographer
         
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
/*  Functions to calculate basic statistics on the phenotypes and to 
    summarize missing markers and calculate segregation distortion.*/


/*
  A vector holds the stats.   Here is what is stored where.
  
tptr = dvector(1,nind)
thestats = dvector(0,kk)
(only kk = 15 useful)

kk
 0         n
 1       Sum(x)
 2       Sum(x^2)
 3       Sum(x^3)
 4       Sum(x^4)
 5       Mean(x) = Sum(x) / n
 6       Var(x)  = (Sum(x^2) - Mean(x)^2)/ (n-1)
 7       StDev(x) = Sqrt(Var(x)) = s
 8       Skew(x) = Sum(x-Mean(x)^3 / (n s^3)
 9       Kurtosis(x) = Sum(x-Mean(x)^4 / (n s^3)   - 3
10       Skw(x) = n^2  ( Sum(x^3)/n - 3 Sum(x^2)/n + 2 Mean(x)^3 ) / ( (n-1)(n-2) )
11       Kur(x) =  n^2 (n+1) ( Sum(x^4)/n - 4 Sum(x^3)/n + 6 Sum(x^2)/n  - 3 Mean(x)^4) / ( (n-1)(n-2)(n-3) )
12       k3 = Skw(x)/ Var(x)^(3/2)
13       k4 = Kur(x)/Var(x)^2   - 3
14       S =  nk3^2/6 + nk4^2/24
15       AveDev(x) = Sum( |x-Mean(x)| ) /n



*/
void calc_qstats(FPN *tptr,int nind,FPN *thestats,int kk)
{
  FPN  z1,z2,z3,z4,n;
  int ii;
  for ( ii=0; ii<=kk; ii++ )
    thestats[ii] = (FPN) 0.0;
  thestats[0] = nind;
  if ( nind <= 3 )
    return;
  for (ii = 1; ii <= nind ; ii++) {  
    thestats[1] = thestats[1] + tptr[ii];
    thestats[2] = thestats[2] + tptr[ii] * tptr[ii];
    thestats[3] = thestats[3] + (FPN)pow(tptr[ii], 3.0 );
    thestats[4] = thestats[4] + (FPN)pow(tptr[ii], 4.0 );  
  }
  thestats[5] = thestats[1] / (FPN) thestats[0];
  thestats[6] = (thestats[2] - thestats[5] * thestats[1]) / (FPN) (thestats[0]-1);
  if ( thestats[6] > (FPN) 0.0 )
    thestats[7] = (FPN) sqrt(thestats[6]);  
  else
    return;
  
  n = thestats[0]; 
  z1 = thestats[1]/n;
  z2 = thestats[2]/n;
  z3 = thestats[3]/n; 
  z4 = thestats[4]/n; 
  
  for (ii = 1; ii <= nind ; ii++) {  /* Average deviation, Skew, Kurtosis */
    thestats[15] = thestats[15] + (FPN) fabs(tptr[ii]-thestats[5]) ;
    thestats[8] = thestats[8] + (FPN) pow( (tptr[ii]-thestats[5]) , 3.0);
    thestats[9] = thestats[9] + (FPN) pow( (tptr[ii]-thestats[5]) , 4.0);    
  }
  thestats[15] = thestats[15] / n; 
  thestats[8] = thestats[8] / (n * (FPN)pow(thestats[7], 3.0)); 
  thestats[9] = ( thestats[9] / ( n * (FPN)pow(thestats[7], 4.0) ) )  - (FPN)3.0; 
/* Lynch and Walsh Skw and Kur*/  
  thestats[10] = (n * n) * ( z3- (FPN)3.0 * z2 * z1 + (FPN)2.0 * (FPN)pow(z1,3.0) ) /  ((n-1) * (n-2));
  thestats[11] = (n * n * (n+(FPN)1.0)) * ( z4  - (FPN)4.0 * z3 * z1 + (FPN)6.0 * z2 * z1 * z1 - (FPN)3.0 * (FPN)pow(z1,4.0) ) / ((n-(FPN)1.0) * (n-(FPN)2.0) * (n-(FPN)3.0));
/*Lynch and Walsh, k3 and k4*/    
  thestats[12] = thestats[10] / (FPN)pow(thestats[6], 1.5) ; 
  thestats[13] = thestats[11] / (FPN)pow(thestats[6], 2.0)  - (FPN)3.0;   
/*Lynch and Walsh, S (Eqn 11.5, page 296)*/
  thestats[14] =  n * ((FPN)4.0 * thestats[12] * thestats[12] + thestats[13] * thestats[13]) / (FPN)24.0;
}

int move_phenotypes( FPN *trptr,individual *individs,int trait, int nn)  
{
  int ii,nind;
  nind = 0;
  for (ii = 1; ii <= nn; ii++)
    if ( individs[ii].y[trait] > (FPN) MISS_VAL ) {
      nind +=1;
	  trptr[nind] = individs[ii].y[trait] ;
	}
  return(nind);
}

void print_qstats(FPN *stats,params *theparams,int trait,markermap *themap)
{
  FILE *outf;
  if (theparams->qstat[0] == '-')
    outf = stdout;
  else {
    outf = fileopen(theparams->qstat, "a");
    if (outf == NULL)
      outf = fileopen(theparams->qstat, "w");
  }
  if (outf != NULL) {
    fprintf(outf, "\n#\n#\tThis output is based on the map in (%s)\n#", theparams->map);
    fprintf(outf, "\tAnd the data in (%s)\n#\n#", theparams->ifile);
    putline(outf, '-',54);
    putline(outf, '-',54);
    if ( themap->tnames != NULL )
      fprintf(outf, "\n\tThis is for -trait %d called %s",trait,themap->tnames[trait]);
    else
      fprintf(outf, "\n\tThis is for -trait number %d",trait);
    putline(outf, '-',54);
    fprintf(outf, "\n\t\tSample Size................%14d", (int) stats[0] );
    fprintf(outf, "\n\t\tM(1).......................%14.4f", stats[1]/stats[0]);
    fprintf(outf, "\n\t\tM(2).......................%14.4f", stats[2]/stats[0]);
    fprintf(outf, "\n\t\tM(3).......................%14.4f", stats[3]/stats[0]);
    fprintf(outf, "\n\t\tM(4).......................%14.4f", stats[4]/stats[0]);
    fprintf(outf, "\n\t\tMean Trait Value...........%14.4f", stats[5]);
    fprintf(outf, "\n\t\tVariance...................%14.4f", stats[6]);
    fprintf(outf, "\n\t\tStandard Deviation.........%14.4f", stats[7]);
    fprintf(outf, "\n\t\tCoefficient of Variation...%14.4f", stats[7]/stats[5]);
/*
    fprintf(outf, "\n\t\tSkew...................%14.4f", stats[8]);
    fprintf(outf, "\n\t\tKurtosis...............%14.4f", stats[9]);
*/
    fprintf(outf, "\n\t\tAverage Deviation..........%14.4f", stats[15]);
    fprintf(outf, "\n\t\tSkw..LW(24)................%14.4f", stats[10]);
    fprintf(outf, "\n\t\t.....Sqrt(6/n).............%14.4f", sqrt(6.0/stats[0]) );
    fprintf(outf, "\n\t\tKur..LW(29)................%14.4f", stats[11]);
    fprintf(outf, "\n\t\t.....Sqrt(24/n)............%14.4f", sqrt(24.0/stats[0]) );
    fprintf(outf, "\n\t\tk3...LW(24)................%14.4f", stats[12]);
    fprintf(outf, "\n\t\tk4...LW(28)................%14.4f", stats[13]);
    fprintf(outf, "\n\t\tS (5%%: 5.99, 1%%: 9.21).....%14.4f", stats[14]);
    putline(outf, '-',54);
    putline(outf, '-',54);
    fprintf(outf, "\n#\n#");
    fileclose(theparams->qstat, outf);
  }
}
/*
  For each marker in turn, test whether it conforms to Hardy-Weinberg
  proportions.
*/
void  marker_segregation(individual *individs,markermap *themap,params *theparams) {
  FILE *fptr;
  FPN freqs[3],test_stat,opAA,opAa,opaa,epAA,epAa,epaa , test_stat2;
  int ii, chrom, mark, marks, nm, nmt,counts[3];

  assign_q(freqs,theparams);



  fptr = fileopen(theparams->qstat, "a");
  if (fptr == NULL)
    fptr = fileopen(theparams->qstat, "w");
  fprintf(fptr, "\n                      Summary of marker segregation ");
  putline(fptr, '-',80);
  putline(fptr, '-',80);
  fprintf(fptr, "\n   Chrom  Mark     Name            type    n(m)           Chi2             LR ");
  putline(fptr, '-',80);
  fprintf(fptr, "                  -begin segregation ");
    
    
  for (chrom = 1; chrom <= themap->m; chrom++) {
    for (mark = 1; mark <= themap->mpc[chrom]; mark++) {
      counts[0]=counts[1]=counts[2]=0;
      test_stat = (FPN) 0.0;
    
      if (themap->ttable != NULL)
	    marks =  themap->ttable[chrom][mark];
      else
	    marks = marks + 1;
      if ( themap->names != NULL )
        fprintf(fptr, "\n %5d %5d   %-20s", chrom, mark, themap->names[marks]);
      else
        fprintf(fptr, "\n %5d %5d                     ", chrom, mark);
      if ( themap->types != NULL )
        switch ( themap->types[chrom][mark] ) {
          case -1 :
            fprintf(fptr, " a-  ");   
            break;  
          case  0 :
            fprintf(fptr, " co  ");   
            break;  
          case  1 :
            fprintf(fptr, " A-  ");   
            break;  
          case  -99 :
            fprintf(fptr, "insuf");   
            break;  
          default :
            fprintf(fptr, "     ");   
            break;  
        }
      else
        fprintf(fptr, "     ");     
      nm = nmt = 0;
      for ( ii = 1 ; ii <= theparams->nn ; ii++ )       
        switch ( individs[ii].markers[chrom][mark] ) {
          case -2 :
          case -1 :
            counts[0] +=1;   
            break;  
          case  0 :
            counts[1] +=1;   
            break;  
          case  1 :
          case  2 :
            counts[2] +=1;   
            break;  
        }
      nm = counts[0] + counts[1] + counts[2];
      opAA = (FPN) counts[2]  ;
      opAa = (FPN) counts[1]  ;
      opaa = (FPN) counts[0]  ;
      epAA = freqs[0] * (FPN) nm;
      epAa = freqs[1] * (FPN) nm;
      epaa = freqs[2] * (FPN) nm;
      test_stat = (FPN) 0.0;
      if ( themap->types != NULL )
        switch ( themap->types[chrom][mark] ) {
          case -1 :  /* a dominant system*/
            if ( epAa+epaa > (FPN) 0.0 )
              test_stat = (opAa+opaa-epAa-epaa)*(opAa+opaa-epAa-epaa)/(epAa+epaa);
            if ( epAA > (FPN) 0.0 )
              test_stat = test_stat + (opAA-epAA)*(opAA-epAA)/epAA ;   
            break;  
          case  0 :   /* codominant system*/
            if ( epAA > (FPN) 0.0 )
              test_stat =  (opAA-epAA)*(opAA-epAA)/epAA ;   
            if ( epAa > (FPN) 0.0 )
              test_stat = test_stat + (opAa-epAa)*(opAa-epAa)/epAa ;   
            if ( epaa > (FPN) 0.0 )
              test_stat = test_stat  + (opaa-epaa)*(opaa-epaa)/epaa ;   
            break;  
          case  1 : /* A dominant  system*/
            if ( epAa+epAA > (FPN) 0.0 )
              test_stat = (opAa+opAA-epAa-epAA)*(opAa+opAA-epAa-epAA)/(epAa+epAA);
            if ( epaa > (FPN) 0.0 )
              test_stat = test_stat  + (opaa-epaa)*(opaa-epaa)/epaa ;   
            break;  
          default :
            test_stat = (FPN)-1.0;   
            break;  
        }
      else
        test_stat = (FPN) 0.0;
      test_stat2 = (FPN) 0.0;
      if ( themap->types != NULL )
        switch ( themap->types[chrom][mark] ) {
          case -1 :  /* a dominant system*/
            if ( epAa+epaa > (FPN) 0.0 && opAa+opaa > 0 )
              test_stat2 = test_stat2  + (opAa+opaa)* (FPN) log((opAa+opaa)/(epAa+epaa))  ;
            if ( epAA > (FPN) 0.0 && opAA > 0)
              test_stat2 = test_stat2 +  opAA*(FPN)log(opAA/epAA)  ;   
            break;  
          case  0 :   /* codominant system*/
            if ( epAA > (FPN) 0.0 && opAA > 0)
              test_stat2 =  test_stat2  +  opAA*(FPN)log(opAA/epAA) ;   
            if ( epAa > (FPN) 0.0 && opAa > 0)
              test_stat2 = test_stat2 +  opAa*(FPN)log(opAa/epAa) ;   
            if ( epaa > (FPN) 0.0 && opaa > 0)
              test_stat2 = test_stat2  +  opaa*(FPN)log(opaa/epaa) ;   
            break;  
          case  1 : /* A dominant  system*/
            if ( epAa+epAA > (FPN) 0.0 && opAa+opAA > 0)
              test_stat2 = test_stat2  +  (opAa+opAA)*(FPN)log((opAa+opAA)/(epAa+epAA));
            if ( epaa > (FPN) 0.0 && opaa > 0)
              test_stat2 = test_stat2  + opaa*(FPN)log(opaa/epaa) ;   
            break;  
          default :
            test_stat2 = (FPN)-0.5;   
            break;  
        }
      else
        test_stat2 = (FPN) 0.0;
      test_stat2 = (FPN)2.0*test_stat2;
      fprintf(fptr, " %5d       %10.4f      %10.4f", nm ,test_stat,test_stat2);


    }
  }
  putline(fptr, '-',80);
  putline(fptr, '-',80);
  fprintf(fptr,"                  -end segregation\n\n");
  fileclose(theparams->qstat, fptr);

}

/*
  For each marker in turn, determine the amount of missing data.
*/
void    miss_mark_summary(individual *individs,markermap *themap,params *theparams,int trait)
{
  FILE *fptr;

  int ii, chrom, mark, marks,ss, nm, nmt;

  fptr = fileopen(theparams->qstat, "a");
  if (fptr == NULL)
    fptr = fileopen(theparams->qstat, "w");
  marks = 0;
  ss = 0;
  for ( ii = 1 ; ii <= theparams->nn ; ii++ )
    if ( individs[ii].y[trait] > (FPN) MISS_VAL )
      ss = ss+1;

  if ( themap->tnames != NULL )
    fprintf(fptr, "\n\nSummary of missing data for trait number %d (%s) \n\t with %d individuals\n\n",trait,themap->tnames[trait],ss);
  else
    fprintf(fptr, "\n\nSummary of missing data for trait %d \n\t where %d individuals have trait data\n",trait,ss);

  fprintf(fptr, "n(m) individuals have marker data, n(m+t) have trait and marker data\n\n");
  putline(fptr, '-',65);
  putline(fptr, '-',65);
  fprintf(fptr, "\n   Chrom  Mark     Name            type    n(m)    n(m+t)  %%(m+t)  ");
  putline(fptr, '-',65);
  fprintf(fptr, "                  -begin missmark  %d",trait);
  for (chrom = 1; chrom <= themap->m; chrom++) {
    for (mark = 1; mark <=  themap->mpc[chrom]; mark++) {
      if (themap->ttable != NULL)
	marks =  themap->ttable[chrom][mark];
      else
	marks = marks + 1;
      if ( themap->names != NULL )
        fprintf(fptr, "\n %5d %5d   %-20s", chrom, mark, themap->names[marks]);
      else
        fprintf(fptr, "\n %5d %5d                     ", chrom, mark);
      if ( themap->types != NULL )
        switch ( themap->types[chrom][mark] ) {
          case -1 :
            fprintf(fptr, " a-  ");   
            break;  
          case  0 :
            fprintf(fptr, " co  ");   
            break;  
          case  1 :
            fprintf(fptr, " A-  ");   
            break;  
          case  -99 :
            fprintf(fptr, "insuf");   
            break;  
          default :
            fprintf(fptr, "     ");   
            break;  
        }
      else
        fprintf(fptr, "     ");     
      nm = nmt = 0;
      for ( ii = 1 ; ii <= theparams->nn ; ii++ ) {
        if ( individs[ii].y[trait] > (FPN) MISS_VAL &&  individs[ii].markers[chrom][mark] >= 0 )
          nmt = nmt+1;
        if (  individs[ii].markers[chrom][mark]  >= 0 )
          nm = nm+1;
      }
      fprintf(fptr, " %5d    %5d   %6.2f", nm, nmt, 100.0 * (FPN) nmt / (FPN) ss );


    }
  }
  putline(fptr, '-',65);
  putline(fptr, '-',65);
  fprintf(fptr, "                  -end missmark   \n\n");
  fileclose(theparams->qstat, fptr);
}

/*
  For each individual in turn, summarize the proportion of missing marker, trait and
  categorical trait data.
*/
void    ind_data_summary(individual *individs,markermap *themap,params *theparams)
{
  FILE *fptr;
  int mark,trait,otrait,imark,chrom,itrait,ii;
  FPN pmark,ptrait,potrait;
#if defined(BORLAND)
	char dummy[15];
#endif
  
  fptr = fileopen(theparams->qstat, "a");
  if (fptr == NULL)
    fptr = fileopen(theparams->qstat, "w");
  fprintf(fptr,"\n\n  Here is a summary of missing data for each individual.\n");
  fprintf(fptr,"\nThere are %d markers, %d traits and %d categorical traits.\n",themap->ml,themap->traits,themap->otraits);
  fprintf(fptr,"The table below lists the raw number and the percentage of data \n\tpoints for each individual.\n");
  putline(fptr, '-',60);
  putline(fptr, '-',60);
  fprintf(fptr, "\nIndividual      Markers         Traits         Cat. Traits");
  if ( individs[1].name == NULL )
    fprintf(fptr, "\n  Number        #      %%        #     %%        #     %% ");
  else
    fprintf(fptr, "\n   #   Name     #      %%        #     %%        #     %% ");
  putline(fptr, '-',60);
  fprintf(fptr, "                  -begin missind");
  for ( ii = 1 ; ii <= theparams->nn ; ii++ ) {
    mark = trait = otrait = 0;
    pmark=ptrait=potrait= (FPN) 0.0;
    for (chrom = 1; chrom <= themap->m; chrom++) {
      for (imark = 1; imark <= themap->mpc[chrom]; imark++) 
        if (  individs[ii].markers[chrom][imark] >= 0 )
          mark = mark+1;
    }
    if ( themap->ml > 0 )
      pmark = (FPN) (100.0 * (FPN) mark / (FPN) themap->ml);      
    for (itrait = 1; itrait <= themap->traits; itrait++) {
		if ( individs[ii].y[itrait] > (FPN) MISS_VAL )
			 trait = trait+1;
	}
	
    if ( themap->traits > 0 )
      ptrait = (FPN) (100.0 *  (FPN) trait / (FPN) themap->traits);
    if ( themap->otraits > 0 ) {
      for (itrait = 1; itrait <= themap->otraits; itrait++) 
        if (  individs[ii].oy[itrait][0] != '\0' )
          otrait = otrait+1; 
      potrait = (FPN) (100.0 *  (FPN) otrait / (FPN) themap->otraits);
    }  
    if ( individs[ii].name == NULL )  
      fprintf(fptr,"\n   %4d       %4d  %6.2f   %4d  %6.2f   %4d  %6.2f",ii,mark,pmark,trait,ptrait,otrait,potrait);
    else
      fprintf(fptr,"\n %4d %7s %4d  %6.2f   %4d  %6.2f   %4d  %6.2f",ii,individs[ii].name,mark,pmark,trait,ptrait,otrait,potrait);
    
  }  
  putline(fptr, '-',60);
  putline(fptr, '-',60);
  fprintf(fptr,"                  -end missind\n\n");
  fileclose(theparams->qstat, fptr);
}


void CalcMarkProbs(params *theparams) {
  genome *gptr;
  FILE *fptr;
  int i,marks;
  FPN pQQ, pQq, pqq,ea,ed;

  
  fptr = fileopen(theparams->qstat, "a");
  if (fptr == NULL)
    fptr = fileopen(theparams->qstat, "w");
  fprintf(fptr,"\n\n  This is a section to infer missing marker data. Each marker in turn is analyzed.");
  fprintf(fptr,"\n  After the individual number follows p(QQ), p(Qq), p(qq), E(a) and E(d).  Zeros and ");
  fprintf(fptr,"\n  Ones indicatenon-missing data.   Distribution vi Jiang and Zeng 1997.  The expected ");
  fprintf(fptr,"\n  marker values are E(a) = p(QQ)-p(qq) and E(d) = [p(Qq)-p(QQ)-p(qq)]/2. \n\n");
 

  fprintf(fptr, "\n-begin markerprobs");
  for (gptr=theparams->thegenome; gptr!= NULL; gptr=gptr->next ) {

      if (theparams->themap->ttable != NULL)
	    marks =  theparams->themap->ttable[gptr->chrom][gptr->markr];

      if ( theparams->themap->names != NULL )
        fprintf(fptr, "\n-markername %-20s  -chromosome %5d -marker %5d", theparams->themap->names[marks], gptr->chrom, gptr->markr);
      else
        fprintf(fptr, "\n %5d %5d",  gptr->chrom, gptr->markr );
    fprintf(fptr,"\n-h   ind.   p(QQ)   p(Qq)    p(qq)   E(a)     E(d)");
    for ( i=1; i<=theparams->nn; i++ ) {
      if ( theparams->thedata[i].markers[gptr->chrom][gptr->markr] == 1 )
       fprintf(fptr, "\n\t%5d  1       0       0        1.0     -0.5", i  );
      else if ( theparams->thedata[i].markers[gptr->chrom][gptr->markr] == 0 )
       fprintf(fptr, "\n\t%5d  0       1       0        0.0      0.5", i  );
      else if ( theparams->thedata[i].markers[gptr->chrom][gptr->markr] == -1 )
       fprintf(fptr, "\n\t%5d  0       0       1       -1.0     -0.5", i  );
      else {    
        calc_cond_prob(theparams,gptr,theparams->thedata,i,&pQQ,&pQq,&pqq);
        ea = pQQ-pqq;
        ed = 0.5*(pQq-pQQ-pqq);
        fprintf(fptr, "\n\t%5d  %6.4f  %6.4f  %6.4f  %7.4f  %7.4f", i, pQQ, pQq, pqq, ea, ed );

      }   
    
    
    }

  
  
  }  
  fprintf(fptr,"\n-end markerprobs\n\n");
  fileclose(theparams->qstat, fptr);

}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file QSfunc.c
------------------------------------------------------------------ */


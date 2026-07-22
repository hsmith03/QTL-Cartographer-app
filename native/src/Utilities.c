/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Utilities.c) is part of QTL Cartographer
         
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


/*  Utility functions used by all the programs. */

#include "Main.h"



char time_buffer[MAXNAME+1];
char gbuffer[MAXLINE+1];      /*Reusable global buffer space*/
char gname[MAXNAME+1];        /*Reusable global name space*/
char gwarn[MAXNAME+1];
FPN ranf(long int inix);
long ix, Ix;
int whosemf,debugging;
FPN mapparam;


/*
   Given a token and a file stream, go to a matching token.
   Return 0 if found, -1 if none.  

*/
int MoveToToken(FILE *fptr, char *token, int usecase) {
  int go_on;
  int ch;
  char *chptr;
  
  if ( usecase == 0 )
     chptr = strlwr(token);
  do {

    ch = get_next_token(gname, MAXNAME, fptr);  /*  first token.  */
    if ( usecase == 0)
      chptr = strlwr(gname);
     
    if (!strcmp(token,gname)   ) 
      go_on = 0;
    else if ( ch == EOF ) 
      go_on = -1;  
    else
      go_on = 1;
  }  while ( go_on == 1); 
  return(go_on);
}


void putline(FILE *fptr, char ch, int len) {
  int i;
  fputc('\n',fptr);
  for ( i=1; i<=len; i++ ) 
    fputc(ch,fptr);
}

void MemoryAccount(char *buffer)
{
  LogTheError("memacct.txt",buffer);
}

void LogTheError(char *filename, char *buffer)
{
  FILE *fptr;

  fptr = fileopen(filename,"a");
  if (fptr == NULL)
    fptr = fileopen(filename,"w");
  fprintf(fptr,"%s",buffer);
  fileclose(filename,fptr);
}

/*
lrtolod converts the likelihood ratio (LR) to the log of the odds (LOD) score.

 LOD = -log10( exp( -LR/2 ) )
 
Because  
  LR  -2ln(L0/L1) 
and 
  LOD = -log10(L0/L1)
*/
FPN lrtolod(FPN lr)
{
  FPN lod;
  lod = (FPN) - log10( exp( -lr/2.0 ) );  
  return(lod);
}

/*
lrtolod converts the likelihood ratio (LR) to the log of the odds (LOD) score.

 LOD = -log10( exp( -LR/2 ) )
 
Because  
  LR  -2ln(L0/L1) 
and 
  LOD = -log10(L0/L1)
*/
FPN lodtolr(FPN lod)
{
  FPN lr;
  lr  = (FPN) 2.0 * lod * (FPN) log(10);
  return(lr);
}
/*
  return 0 if ok, -1 if not an number.
*/
int is_number(char *buffer)
{
  int len,not_int,ii,start;
  not_int = 0;
  len = (int) strlen(buffer);
  if ( buffer[0] == '-' ) 
    start = 1;
  else 
    start = 0;
  for ( ii = start ; ii < len ; ii++ )
    if ( !isdigit(buffer[ii])  )
      not_int = -1;

  return(not_int);
}

int is_pinteger(char *buffer)
{
  int len,not_int,ii;
  not_int = 0;
  len = (int) strlen(buffer);
  for ( ii = 0 ; ii < len ; ii++ )
    if ( !isdigit(buffer[ii]) )
      not_int = -1;
  if ( not_int == 0 )
    not_int = atoi(buffer);
  return(not_int);
}

/*  Assume that buffer contains a string.  
    Assume that names = cmatrix(1,n,0,somelength)   contains a set of strings
    Look for row in which buffer resides and return index.
    If no there, return -1.
    Check at most len chars.  
*/
int find_token(int n, int len, char *buffer,char **names) {
  int l;
  for ( l=1; l<=n; l++ )
    if ( ! strncmp( buffer, names[l], (size_t) len ) )
      return(l);      
  return(-1);

}
int myfgets(char *buffer, int nn, FILE *fptr)
{
  int ii,k;
  char *chptr;
  k=0;
  for (ii = 0; ii < nn; ii++)
    buffer[ii] = '\0';
  chptr = fgets(buffer, nn, fptr);
  for (ii = 0; ii < nn-1 && buffer[ii] != '\n'; ii++) k+=1;
    if (ii < nn)
      buffer[ii] = '\0';
  return(ii);
}


#if defined(BORLAND)
#define BIG_NUMBER DBL_MAX
#define LITTLE_NUMBER DBL_MIN



void WYield()
{
  MSG msg;
  while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
    TranslateMessage(&msg);
    DispatchMessage(&msg);
  }
}
#endif

#define A 16807.0
#define P 2147483647.0
#define a 16807
#define b15 32768
#define b16 65536
#define p 2147483647
#define xnorm 4.656612875E-10


/*
  Return a random integer from 1 to in, using ranf(ix) as the
  source of uniform deviates.
     Fortran source from John Monahan,
     and translated into c by Chris Basten 19 January 1994.
*/
long iran(long int xix, long int in)
{
  FPN rv;
  rv = ranf(xix);
  return ((long) (rv * (FPN) in) + 1L);
}
/*
 * This subroutine will compute an array of uniformly
 * distributed random numbers.  You need to specify the
 * size of the array, and the numbers will be FPNs.
 */

FPN *ran_arry(int size)
{
  FPN *arry;
  int ii;
  long I, bound;
  if (size < 0) {
    bound = -size;
    Ix = (long) ((FPN) P * ranf(bound));
  }
  else
    bound = size;
  arry = dvector(0, bound);
  for (ii = 0; ii <= bound; ii++) {
    I = (long) ((Ix * A) / P);
    Ix = (long)  (Ix * A - I * P);
    arry[ii] = Ix * (FPN) xnorm;
  }
  return arry;
}

FPN ranf(long int inix)
{
  long xhi, xalo, leftlo, fhi, k;
  FPN xx;
  if (inix < 0)
    ix = -1 * inix;
  xhi = ix / b16;
  xalo = (ix - xhi * b16) * a;
  leftlo = xalo / b16;
  fhi = xhi * a + leftlo;
  k = fhi / b15;
  ix = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
  if (ix < 0)
    ix = ix + p;
  xx = (FPN) ix * (FPN) xnorm;
  return xx;
}

#undef A  
#undef P  
#undef a  
#undef b15  
#undef b16  
#undef p  
#undef xnorm  

/*
 * The following come from Numerical Recipes in C, and
 * handle the allocation of memory.  The character
 * allocators were modified by me.
 */


void nrerror(char *error_text)
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

FPN *vector(int nl, int nh)
{
  FPN *v;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  v = (FPN *) malloc((size_t) (nh - nl + 1) * sizeof(FPN));
#else
  v = (FPN *) malloc((unsigned) (nh - nl + 1) * sizeof(FPN));
#endif  
  if (!v)
    nrerror("allocation failure in vector()");
  else if ( debugging > 2 ) {

    sprintf(gwarn,"In vector(), allocated %d FPNs at %p\n",nh-nl+1,v);
    MemoryAccount(gwarn);
  } 
  return v - nl;
}
void free_vector(FPN *v, int nl, int nh)
{
  
  if ( debugging > 2 ) {
    sprintf(gwarn,"In free_vector(), deallocated %d FPNs at %p\n",nh-nl+1,v+nl);
    MemoryAccount(gwarn);
  } 
  free((char *) (v + nl));
}


FPN *dvector(int nl, int nh)
{
  int ii;
  FPN *v;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  v = (FPN *) malloc((size_t) (nh - nl + 1) * sizeof(FPN));
#else
  v = (FPN *) malloc((unsigned) (nh - nl + 1) * sizeof(FPN));
#endif  
  if (!v)
    nrerror("allocation failure in dvector()");
  else if ( debugging > 2 ) {

    sprintf(gwarn,"In dvector(), allocated %d FPNs at %p\n",nh-nl+1,v);
    MemoryAccount(gwarn);
  } 
  for ( ii = 0; ii <= nh-nl ; ii++ )
	 v[ii] = (FPN) 0.0;
  return v - nl;
}
void free_dvector(FPN *v, int nl, int nh)
{
   
  if ( debugging > 2 ) {
    sprintf(gwarn,"In free_dvector(), deallocated %d FPNs at %p\n",nh-nl+1,v+nl);
    MemoryAccount(gwarn);
  } 
  free((char *) (v + nl));
}

int *ivector(int nl, int nh)
{
  int *v,ii;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  v = (int *) malloc((size_t) (nh - nl + 1) * sizeof(int));
#else
  v = (int *) malloc((unsigned) (nh - nl + 1) * sizeof(int));
#endif
  if (!v)
    nrerror("allocation failure in ivector()");
  else if ( debugging > 2 ) {

    sprintf(gwarn,"In ivector(), allocated %d ints at %p\n",nh-nl+1,v);
    MemoryAccount(gwarn);
  } 
  for ( ii = 0; ii <= nh-nl ; ii++ )
	 v[ii] = 0;
  return v - nl;
}
void free_ivector(int *v, int nl, int nh)
{
  
  if ( debugging > 2 ) {
    sprintf(gwarn,"In free_ivector(), deallocated %d ints at %p\n",nh-nl+1,v+nl);
    MemoryAccount(gwarn);
  } 
  free((char *) (v + nl));
}

char *cvector(int nl, int nh)
{
  char *v;
  int ii;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  v = (char *) malloc((size_t) (nh - nl + 1) * sizeof(char));
#else
  v = (char *) malloc((unsigned) (nh - nl + 1) * sizeof(char));
#endif
  if (!v)
    nrerror("allocation failure in cvector()");
  else if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In cvector(), allocated %d chars at %p\n",nh-nl+1,v);
    MemoryAccount(gwarn);
  } 
  for ( ii = 0; ii <= nh-nl ; ii++ )
	 v[ii] = '\0';
  return v - nl;
}
void free_cvector(char *v, int nl, int nh)
{
  
  if ( debugging > 2 ) {
    sprintf(gwarn,"In free_cvector(), deallocated %d chars at %p\n",nh-nl+1,v+nl);
    MemoryAccount(gwarn);
  } 
  free((char *) (v + nl));
}

long *lvector(int nl, int nh)
{
  int ii;
  long *v;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  v = (long *) malloc((size_t) (nh - nl + 1) * sizeof(long));
#else
  v = (long *) malloc((unsigned) (nh - nl + 1) * sizeof(long));
#endif
  if (!v)
    nrerror("allocation failure in lvector()");
  else if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In lvector(), allocated %d longs at %p\n",nh-nl+1,v);
    MemoryAccount(gwarn);
  } 
  for ( ii = 0; ii <= nh-nl ; ii++ )
	 v[ii] = 0L;
  return v - nl;
}
void free_lvector(long int *v, int nl, int nh)
{
  
  if ( debugging > 2 ) {
    sprintf(gwarn,"In free_lvector(), deallocated %d longs at %p\n",nh-nl+1,v+nl);
    MemoryAccount(gwarn);
  } 
  free((char *) (v + nl));
}


long **lsvector(int nrl, int nrh )
{
  int i;
  long **m;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  m = (long **) malloc((size_t) (nrh - nrl + 1) * sizeof(long *));
#else
  m = (long **) malloc((unsigned) (nrh - nrl + 1) * sizeof(long *));
#endif
  
  if (!m)
    nrerror("allocation failure 1 in lsvector()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)  
    m[i] = NULL;
  if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In lsvector(), allocated %d long vector of pointers  at %p\n",nrh-nrl+1, m+nrl);
    MemoryAccount(gwarn);
  } 
  return m;
}
void free_lsvector(long **m, int nrl, int nrh )
{
  
  if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In free_lsvector(), deallocated %d long vector of pointers at %p\n",nrh-nrl+1,m+nrl);
    MemoryAccount(gwarn);
  } 
  free((char *) (m + nrl));
}

FPN **dsvector(int nrl, int nrh )
{
  int i;
  FPN **m;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  m = (FPN **) malloc((size_t) (nrh - nrl + 1) * sizeof(FPN *));
#else
  m = (FPN **) malloc((unsigned) (nrh - nrl + 1) * sizeof(FPN *));
#endif

  if (!m)
    nrerror("allocation failure 1 in dsvector()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)  
    m[i] = NULL;
  if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In dsvector(), allocated %d FPN vector of pointers at %p\n",nrh-nrl+1, m+nrl);
    MemoryAccount(gwarn);
  } 
  return m;
}
void free_dsvector(FPN **m, int nrl, int nrh )
{
   
  if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In free_dsvector(), deallocated %d FPN vector of pointers %p\n",nrh-nrl+1,m+nrl);
    MemoryAccount(gwarn);
  } 
  free((char *) (m + nrl));
}

FPN **svector(int nrl, int nrh )
{
  int i;
  FPN **m;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  m = (FPN **) malloc((size_t) (nrh - nrl + 1) * sizeof(FPN *));
#else
  m = (FPN **) malloc((unsigned) (nrh - nrl + 1) * sizeof(FPN *));
#endif
  if (!m)
    nrerror("allocation failure 1 in svector()");
  m -= nrl;
  for (i = nrl; i <= nrh; i++)  
    m[i] = NULL;
  if ( debugging > 2 ) {
    sprintf(gwarn,"In svector(), allocated %d FPN vector of pointers at %p\n",nrh-nrl+1, m+nrl);
    MemoryAccount(gwarn);
  } 
  return m;
}
void free_svector(FPN **m, int nrl, int nrh )
{
  if ( nrh > nrl ) {
	  if ( debugging > 2 ) {
	    sprintf(gwarn,"In free_svector(), deallocated FPN vector of pointers %p\n",m+nrl);
	    MemoryAccount(gwarn);
	  } 
	  free((char *) (m + nrl));
  }
}

char **csvector(int nrl, int nrh )
{
  int i;
  char **m;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  m = (char **) malloc((size_t) (nrh - nrl + 1) * sizeof(char *));
#else
  m = (char **) malloc((unsigned) (nrh - nrl + 1) * sizeof(char *));
#endif
  if (!m)
    nrerror("allocation failure 1 in csvector()");
  m -= nrl;
  for (i = nrl; i <= nrh; i++)  
    m[i] = NULL;
  if ( debugging > 2 ) {
    sprintf(gwarn,"In csvector(), allocated %d char vector of pointers at %p\n",nrh-nrl+1, m+nrl);
    MemoryAccount(gwarn);
  } 
  return m;
}
void free_csvector(char **m, int nrl, int nrh )
{
  if ( debugging > 2 ) {
    sprintf(gwarn,"In free_csvector(), deallocated %d char vector of pointers %p\n",nrh-nrl+1,m+nrl);
    MemoryAccount(gwarn);
  } 
  free((char *) (m + nrl));
}

int **isvector(int nrl, int nrh )
{
  int i;
  int **m;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  m = (int **) malloc((size_t) (nrh - nrl + 1) * sizeof(int *));
#else
  m = (int **) malloc((unsigned) (nrh - nrl + 1) * sizeof(int *));
#endif
  if (!m)
    nrerror("allocation failure 1 in isvector()");
  m -= nrl;
  for (i = nrl; i <= nrh; i++)  
    m[i] = NULL;
  if ( debugging > 2 ) {
    sprintf(gwarn,"In isvector(), allocated %d int vector of pointers at %p\n",nrh-nrl+1, m+nrl);
    MemoryAccount(gwarn);
  } 
  return m;
}
void free_isvector(int **m, int nrl, int nrh )
{
  if ( debugging > 2 ) {
    sprintf(gwarn,"In free_isvector(), deallocated %d int vector of pointers %p\n",nrh-nrl+1,m+nrl);
    MemoryAccount(gwarn);
  } 
  free((char *) (m + nrl));
}

char **cmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  char **m;
  m = csvector(nrl,nrh);
  for (i = nrl; i <= nrh; i++)
    m[i] = cvector(ncl,nch);
  return(m);
}
int **imatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  int **m;
  m = isvector(nrl,nrh);
  for (i = nrl; i <= nrh; i++)
    m[i] = ivector(ncl,nch);
  return(m);
}
long **lmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  long **m;
  m = lsvector(nrl,nrh);
  for (i = nrl; i <= nrh; i++)
    m[i] = lvector(ncl,nch);
  return(m);
}
FPN **matrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  FPN **m;
  m = svector(nrl,nrh);
  for (i = nrl; i <= nrh; i++)
    m[i] = vector(ncl,nch);
  return(m);
}
FPN **dmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  FPN **m;
  m = dsvector(nrl,nrh);
  for (i = nrl; i <= nrh; i++)
    m[i] = dvector(ncl,nch);
  return(m);
}

void free_dmatrix(FPN **m, int nrl, int nrh, int ncl, int nch)
{
  int i;  
  for (i = nrh; i >= nrl; i--)
    free_dvector(m[i],ncl,nch);
  free_dsvector(m,nrl,nrh);
}
void free_matrix(FPN **m, int nrl, int nrh, int ncl, int nch)
{
  int i;  
  for (i = nrh; i >= nrl; i--)
    free_vector(m[i],ncl,nch);
  free_svector(m,nrl,nrh);
}
void free_lmatrix(long **m, int nrl, int nrh, int ncl, int nch)
{
  int i;  
  for (i = nrh; i >= nrl; i--)
    free_lvector(m[i],ncl,nch);
  free_lsvector(m,nrl,nrh);
}
void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch)
{
  int i;  
  for (i = nrh; i >= nrl; i--)
    free_cvector(m[i],ncl,nch);
  free_csvector(m,nrl,nrh);
}
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
  int i;  
  for (i = nrh; i >= nrl; i--)
    free_ivector(m[i],ncl,nch);
  free_isvector(m,nrl,nrh);
}

/*****************************************************************************
 *        This is a general Gamma variate generator with the shape
 *     parameter beta > 0.25. Coded from the Algorithm GBH of Cheng and Feast
 *     (1980 Communications of the ACM 23:389-394) by Zhao-Bang Zeng on Jan. 30,
 *     1992.
 *     This was translated into c by Chris Basten, 19 Jan. 1994.
 *****************************************************************************/
FPN gamgbh(FPN beta, long int xix)
{
  FPN aa, bb, c, d, t, h1, h2, u, u1, u2, w, tmp;

  if (beta <= 0.25)
    nrerror("Error in gamgbh:  beta <= 0.25...");
  aa = beta - (FPN) 0.25;
  bb = beta / aa;
  c = (FPN) 2.0 / aa;
  d = c + (FPN) 2.0;
  t = (FPN) 1.0 / (FPN) sqrt(beta);
  h1 = ((FPN) 0.4417 + (FPN) 0.0245 * t / beta) * t;
  h2 = ((FPN) 0.222 - (FPN) 0.043 * t) * t;
  do {
    u = ranf(xix);
    u1 = ranf(xix);
    u2 = u1 + h1 * u - h2;
    if (u2 <= (FPN) 0.0 || u2 >= (FPN) 1.0)
      tmp = (FPN) 1.0;
    else {
      w = bb * (FPN) pow((u1 / u2), (FPN) 4.0);
      if (c * u2 - d + w + (FPN) 1.0 / w <= (FPN) 0.0)
	return (aa * w);
      else
	tmp = c * (FPN) log(u2) - (FPN) log(w) + w - (FPN) 1.0;
    }
  } while (tmp >= (FPN) 0.0);
  return (aa * w);
}

/****************************************************************************
      This gamma generator generates gamma variates by composition-
   rejection from the Weibull density for the shape parameter beta
   smaller than 1 (0 < beta < 1). Coded from I. Vaduva (1977 Math.
   Operationsforsch. Statist., Ser. Statistics, Vol. 8:545-576) by
   Zhao-Bang Zeng on Feb. 5 1992.   (Best one of this kind)

   Translated into c by Chris Basten on 19 January 1994.
 *****************************************************************************/
FPN gamnl1(FPN beta, FPN aa, FPN bb, FPN pp, long int xix)
{
  FPN gamnl1r, u1, u2, u0;

  if (beta >= 1.0 || beta <= 0.0)
    nrerror("beta in gamnl1 is not in (0,1)");

/*  Change > to < in the following line suggested by Ian Painter. 21 Oct. 1996*/
  if (ranf(xix) < pp) {
    u0 = ranf(xix);
    gamnl1r = u0 = (FPN) pow(u0, aa);
    u1 = ranf(xix);
    while (u0 >= u1) {
      u2 = ranf(xix);
      if (u1 < u2) {
	u0 = ranf(xix);
	gamnl1r = u0 = (FPN) pow(u0, aa);
      }
      else
	u0 = u2;
      u1 = ranf(xix);
    }
  }
  else
    do {
      u0 = ranf(xix);
      u0 = (FPN) pow(u0, bb);
      gamnl1r = (FPN) 1.0 - (FPN) log(ranf(xix));
    } while (gamnl1r >= u0);
  return (gamnl1r);
}


/*
For:
  MIN_INT < lb < ub < MAX_INT  and ub - lb < MAX_INT
Do:
  Shuffle  v = ivector(lb,ub)
Where:
  ranf() gives a uniform random number in [0.0, 1.0)

foreach i in [lb,ub) in order, pick a random j from (i,ub]
and switch v(i), v(j)

*/
void shuffle_ivector(int *v, int lb, int ub)
{
  int i, j, k;
  for (i = lb; i < ub; i++) {
    j = i + (int) ((FPN) (ub - i) * ranf(i)) + 1;
    k = v[i];
    v[i] = v[j]  ;
    v[j]   = k;
    if (j > ub || j <= i)
      printf("\nlb = %d, i = %d, j = %d, ub = %d\n", lb, i, j, ub);
  }
}



/* return 0 if no file 1 if file exists */
int isfile(char *filename)
{
  int isf;
  FILE *fptr;
  isf = 0;
  if (*(filename) != '\0') {
    fptr = fopen(filename, "r");
    if (fptr != NULL) {
      isf = 1;
      fclose(fptr);
    }
  }
  return (isf);
}


int get_int(void)
{
  char buffer[15], ch;
  int ii, ans,k;
  k=0;
  for (ii = 0; ii < 15; ii++)
    buffer[ii] = '\0';
  while (isspace(ch = (char) getchar())) k+=1;
  buffer[0] = ch;
  for (ii = 1; ii < 15 && buffer[ii-1]  != '\n'; ii++)
    buffer[ii] = (char) getchar();
  if (ii < 15) {
    buffer[ii-1]   = '\0';
    ans = atoi(buffer);
  }
  else
    ans = -1;
  return (ans);
}

/*
This is a function to convert between recombination
frequencies and distance in Morgans or centiMorgans.

value is what will be converted.
flag indicates how it will be converted:
flag =
   -2  => value cM to rvalue Rec. Freq.
   -1  => value  M to rvalue Rec. Freq.
    0  => value is returned unchanged.
    1  => value Rec. Freq. to rvalue M
    2  => value Rec. Freq. to rvalue cM


The global variable whosemf determines which mapfunction to
use.

  whosemf =
             1 => Haldane (1919)  
             2 => Kosambi (1944)  
             3 => Morgan (1928) (Fixed)  
             4 => Carter and Falconer (1951)
             5 => Rao et al (1977)
             6 => Sturt (1976)
             7 => Felsenstein (1979)
             8 => Karlin (1984)
See Ben Hui Liu (1998) "Statistical Genomics: Linkage, Mapping and QTL Analysis" p319
CRC Press

Send Morgans or r to the following functions:


Kosambi, iKosambi
Haldane, iHaldane
Morgan, iMorgan
CarterFalconer, iCarterFalconer
Rao, iRao
Sturt, iSturt
Felsenstein, iFelsenstein
Karlin, iKarlin


Check that 
  0.0 < m 
  0.0 < r < 0.5
*/
FPN mapfunc(FPN value, int flag)
{
  FPN rvalue;
  FPN mval,rval;
  /* 0.0 < value
     0.0 < value < 0.5 if value is a recombination frequency */
  if ( value < (FPN) 0.0 )
    return((FPN) -1.0);  /* -1 for a negative distance/probability */
  else if ( value == (FPN) 0.0 ) 
    return((FPN) 0.0);    
  if ( flag > 0 && value >= (FPN) 0.5 ) /* -2 for a rec. prob. >= 1/2 */
    return((FPN) -2.0);
  
    
  mval = rval = value;
  if (flag == -2 )
    mval = mval* (FPN) 0.01;
  else if ( flag == 0)
    return (value);
   
  if ( flag < 0 ) {
    switch(whosemf) {
      default: case 1: rval = iHaldane(mval);  break;
      case 2: rval = iKosambi(mval);  break;
      case 3: rval = iMorgan(mval);  break;
      case 4: rval = iCarterFalconer(mval);  break;
      case 5: rval = iRao(mval);  break;
      case 6: rval = iSturt(mval);  break;
      case 7: rval = iFelsenstein(mval);  break;
      case 8: rval = iKarlin(mval);  break;
    }
  
  }
  else if ( flag > 0 ) {
    switch(whosemf) {
      default: case 1: mval = Haldane(rval);  break;
      case 2: mval = Kosambi(rval);  break;
      case 3: mval = Morgan(rval);  break;
      case 4: mval = CarterFalconer(rval);  break;
      case 5: mval = Rao(rval);  break;
      case 6: mval = Sturt(rval);  break;
      case 7: mval = Felsenstein(rval);  break;
      case 8: mval = Karlin(rval);  break;
    }  
  }
  
    
  if (flag == 2 )  /*change Morgans to Centimorgans*/
    rvalue = mval* (FPN) 100.0;
  else if ( flag == 1 )
    rvalue = mval;
  else if ( flag == -1 || flag == -2)
    rvalue = rval;
  else 
    rvalue = value;
  
  return (rvalue);
}

/*Inverse of Kosambi mapping function*/
FPN iKosambi(FPN mm)
{
  FPN rr;   
  rr = (FPN) 0.5 * ((FPN) 1.0 - (FPN) exp( -4.0 * mm)) / ((FPN) 1.0 + (FPN) exp( -4.0 * mm));
  return (rr);
}
/*Kosambi mapping function*/
FPN Kosambi(FPN rr)
{
  FPN mm;   
  mm = (FPN) 0.25 * (FPN) log((1.0 + 2.0 * rr) / (1.0 - 2.0 * rr));
  return (mm);
}

/*Inverse of Morgan mapping function*/
FPN iMorgan(FPN mm)
{
  FPN rr;   
  rr = mm;
  return (rr);
}
/*Morgan mapping function*/
FPN Morgan(FPN rr)
{
  FPN mm;   
  mm = rr;
  return (mm);
}

/*Inverse of Hadane mapping function*/
FPN iHaldane(FPN mm)
{
  FPN rr;   
  rr = (FPN) 0.5 * ((FPN) 1.0 - (FPN) exp(-2.0 * mm));;
  return (rr);
}
/*Hadane mapping function*/
FPN Haldane(FPN rr)
{
  FPN mm;   
  mm = (FPN) -0.5 * (FPN) log(1.0 - 2.0 * rr);
  return (mm);
}
/*Inverse of CarterFalconer mapping function*/
FPN iCarterFalconer(FPN mm)
{
  FPN rr,ru,rl,delta,mt;   
      delta = (FPN) MAPDELTA;
      rl = delta;
      ru = (FPN) 0.5-delta;
      do {
        rr = (FPN) 0.5*(rl+ru);
        mt = CarterFalconer(rr);
        if ( mt > mm )
          ru = rr;
        else
          rl = rr;
      } while ( fabs(mt-mm) > delta );
  return (rr);
}
/*CarterFalconer mapping function*/
FPN CarterFalconer(FPN rr)
{
  FPN mm;    
  mm = (FPN) 0.5*( (FPN) atan(2.0*rr) + (FPN) 0.5* (FPN) log( (1.0+2.0*rr)/(1.0-2.0*rr) ) ) ;
  return (mm);
}
/*Inverse of Felsenstein mapping function*/
FPN iFelsenstein(FPN mm)
{
  FPN rr,kk;   
  kk = mapparam;   
  rr = ((FPN) 1.0- (FPN) exp(2.0*(kk-2.0)*mm))/((FPN) 2.0*((FPN) 1.0-(kk- (FPN) 1.0)* (FPN) exp(2.0*(kk-2.0)*mm)));
  return (rr);
}
/*Felsenstein mapping function*/
FPN Felsenstein(FPN rr)
{
  FPN mm,kk;
  kk = mapparam;   /* kk can't be 2 */
  mm = (FPN) log( (1.0-2.0*rr)/(1.0-2.0*(kk-1.0)*rr) )/((FPN) 2.0*(kk-(FPN) 2.0)) ;
  return (mm);
}
/*Inverse of Karlin mapping function*/
FPN iKarlin(FPN mm)
{
  FPN rr,nn;
  nn = mapparam;   
  rr = (FPN) 0.5*((FPN) 1.0 - (FPN) pow((1.0-2.0*mm/nn),nn) );
  return (rr);
}
/*Karlin mapping function*/
FPN Karlin(FPN rr)
{
  FPN mm,nn;
  nn = mapparam;      
  mm = (FPN) 0.5 * nn *((FPN) 1.0- (FPN) pow((1.0-2.0*rr),1.0/nn));
  return (mm);
}
/*Inverse of Rao mapping function*/
FPN iRao(FPN mm)
{
  FPN rr,ru,rl,delta,mt;   
      delta = (FPN) MAPDELTA;
      rl = delta;
      ru = (FPN) 0.5-delta;
      do {
        rr = (FPN) 0.5*(rl+ru);
        mt = Rao(rr);
        if ( mt > mm )
          ru = rr;
        else
          rl = rr;
      } while ( fabs(mt-mm) > delta );
  return (rr);
}
/*Rao mapping function*/
FPN Rao(FPN rr)
{
  FPN mm,pp;
  pp = mapparam;     
  mm = pp*((FPN)2.0*pp-(FPN)1.0)*((FPN)1.0-(FPN)4.0*pp)* (FPN) log(1.0-2.0*rr)/(FPN)6.0 + ((FPN)8.0*pp*(pp-(FPN)1.0)*((FPN)2.0*pp-(FPN)1.0)*(FPN) atan(2.0*rr) +  pp*((FPN)1.0-pp)*((FPN)4.0*pp+(FPN)1.0) * (FPN) log((1.0+2.0*rr)/(1.0-2.0*rr)))/(FPN)3.0 + ((FPN)1.0-pp)*((FPN)1.0-(FPN)2.0*pp)*((FPN)1.0-(FPN)4.0*pp)*rr;
  return (mm);
}
/*Inverse of Sturt mapping function*/
FPN iSturt(FPN mm)
{
  FPN rr,ll;
  ll = (FPN) mapparam;   
  if ( mm < ll )
    rr = (FPN) 0.5*((FPN) 1.0-((FPN) 1.0-mm/ll)*(FPN) exp(mm*(1.0-2.0*ll)/ll) );
  else
  	rr = (FPN) 0.5;
  return (rr);
}
/*Sturt mapping function*/
FPN Sturt(FPN rr)
{
  FPN mm,ml,mu,rt,delta;   
  ml = delta = (FPN) MAPDELTA;
  do {
    mu = ml+ (FPN) 1.0;
    rt = iSturt(mu);
    if (rt < rr )
      ml = ml + (FPN) 1.0;
  
  } while ( rt < rr );
  do {
        mm = (FPN) 0.5*(ml+mu);
        rt = iSturt(mm);
        if ( rt > rr )
          mu = mm;
        else
          ml = mm;
  } while ( fabs(rt-rr) > delta );
  return (mm);
}



int get_next_line(char *buffer, int nn, FILE *fileptr)
{
  int ii, ch;
  for (ii = 0; ii < nn; ii++)
    buffer[ii] = '\0';

  for (ii = 0; ((ch = fgetc(fileptr)) != EOF) && (ch != '\n') && (ii < nn); ii++)
    buffer[ii] = (char) ch;
  if (ch != EOF)
    ch = buffer[0];
  return (ch);
}





int get_next_token(char *xtemp, int nn, FILE *fileptr)
{
  int ii, ch;
  for (ii = 0; ii < nn; ii++)
    xtemp[ii] = '\0';
  ch = fgetc(fileptr);
  while (isspace(ch) && ch != EOF )
	 ch = fgetc(fileptr);
  if ( ch != EOF )
  for (ii = 0; ii < nn && !isspace(ch) && ch != EOF; ii++) {
    xtemp[ii] = (char) ch;
    ch = fgetc(fileptr);
  }
  return (ch);
}

void gnuplot_values(FPN **xx, FPN **yy, int lr, int ur, int lc, int uc, char *thefile, char *title, char *xaxis, char *yaxis)
{
  char *stem,   *outputf, *gnufile;
  FILE *output;
  int ii, jj;
  stem = cvector(0,MAXNAME);
  outputf = cvector(0,MAXNAME);
  gnufile = cvector(0,MAXNAME);
  for (ii = 0; (thefile[ii] != '.' && thefile[ii] != '\0') && ii < MAXNAME - 4; ii++)
    stem[ii] = thefile[ii];
  stem[ii] = '.';
  strcpy(gnufile, stem);
  strcat(gnufile, "gnuplot");
  output = fileopen(gnufile, "w");
  fprintf(output, "# For use in gnuplot...\n");
  fprintf(output, "\nset autoscale");
  fprintf(output, "\nset title \"%s\"", title);
  fprintf(output, "\nset xlabel \"%s\"", xaxis);
  fprintf(output, "\nset ylabel \"%s\"", yaxis);
  fileclose(gnufile, output);

  for (ii = lr; ii <= ur; ii++) {
    for (jj = 0; jj < MAXNAME; jj++)
      outputf[jj] = '\0';
    sprintf(outputf,"%s%d",stem,ii);
    output = fileopen(outputf, "w");
    for (jj = lc; xx[ii][jj] >= 0.0 && jj <= uc; jj++)
      fprintf(output, "%9.4f %9.4f\n", xx[ii][jj], yy[ii][jj]);
    fileclose(outputf, output);
    output = fileopen(gnufile, "a");
    fprintf(output, "\nplot \'%s\' with lines", outputf);
    fprintf(output, "\npause -1 \"Hit return to continue\" \n");
    fileclose(gnufile, output);
  }
  free_cvector( stem,0,MAXNAME);
  free_cvector( outputf,0,MAXNAME);
  free_cvector( gnufile,0,MAXNAME);

}



/*
  This function reads the first line of the file, and if there is
  a '#' as the very first character, it gets the long number immediately
  following, and returns it.  This should be a unique number characterizing
  the file.
*/
long get_identifier(char *thefile)
{
  FILE *infile;
  int ii, ch;
  long identifier;
  identifier = 0;
  infile = fileopen(thefile, "r");
  if (infile == NULL)
    return (-1);
  for (ii = 0; ii < MAXLINE; ii++)
    gbuffer[ii] = '\0';
  for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
    gbuffer[ii] = (char) ch;
  if (*gbuffer == '#') {
    get_field(2, gwarn, gbuffer);
    identifier = atol(gwarn);
  }
  fileclose(thefile, infile);
  return (identifier);
}

/*

y = dvector(1,n) is a vector of FPNs.  
n is the number of FPNs
intervals is the number of intervals in the histogram


print out a histogram of the y vector.     Determine min and max values and
divide the  range into intervals intervals.   Scale it so that the ordinate goes
from 0 to half the number of intervals.  

This had a bug in it that was fixed on 27 April 2001.   

*/
void print_histogram(FPN *y, int n, int intervals, char *outfile)
{
  FPN ymax, ymin, ydelta,*yvect,q1,q2,q3;
  int ii, jj,  histmax, height, *hist,element ;
  FILE *outf;
  height = intervals / 2;
  hist = ivector(1, intervals);
  yvect = dvector(1,n);
  for ( ii=1; ii<=n; ii++)  
    yvect[ii] = y[ii];
  sort(n, yvect);

  if (outfile[0] == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  if (outf != NULL) {
    fprintf(outf, "\n\nHere is a histogram for the quantitative trait...\n");

  for (ii = 1; ii <= intervals; ii++)
        hist[ii] = 0;

    ymin =  yvect[1];
    ymax =  yvect[n];
 
    ydelta = (ymax - ymin) / (FPN) intervals;
    for (ii = 1 ; ii <= n; ii++)  {
		  element =  ((int) ((y[ii] - ymin) / ydelta) + 1);
		  if ( element == intervals+1 )
			 element = intervals;
		  if ( element > 0 && element <= intervals )
			hist[element] = hist[element]+1;
		}



    histmax = hist[1];
    for (ii = 2; ii <= intervals; ii++)
      if (histmax < hist[ii])
	    histmax = hist[ii];

    for (ii = 1; ii <= intervals; ii++)
      hist[ii] = (int) ((FPN) height * (FPN) hist[ii] / (FPN) histmax) + 1;

    for (jj = height; jj > 0; jj--) {
      fprintf(outf, "\n  |");
      for (ii = 1; ii <= intervals; ii++)
	if (hist[ii] > jj)
	  fprintf(outf, "*");
	else
	  fprintf(outf, " ");
    }
    fprintf(outf, "\n  -+");
    for (ii = 2; ii <= intervals; ii++)
      if (ii == intervals / 2 || ii == intervals)
	fprintf(outf, "+");
      else
	fprintf(outf, "-");
    fprintf(outf, "\n %5.2f", ymin);
    for (ii = 1; ii <= intervals / 2 - 6; ii++)
      fprintf(outf, " ");
    fprintf(outf, "%5.2f", (ymax + ymin) / 2.0);
    for (ii = 1; ii <= intervals / 2 - 5; ii++)
      fprintf(outf, " ");
    fprintf(outf, "%5.2f\n", ymax);

    fprintf(outf, "  min ");
    for (ii = 1; ii <= intervals / 2 - 6; ii++)
      fprintf(outf, " ");
    fprintf(outf, "  Y  ");
    for (ii = 1; ii <= intervals / 2 - 5; ii++)
      fprintf(outf, " ");
    fprintf(outf, " max\n\n");
    
    q1 = (yvect[ (int) floor(0.25* (FPN) n )   ] + yvect[(int) ceil(0.25* (FPN) n ) ])/(FPN) 2.0;
    q2 = (yvect[ (int) floor(0.5* (FPN) n )   ] + yvect[(int) ceil(0.5* (FPN) n ) ])/(FPN) 2.0;
    q3 = (yvect[ (int) floor(0.75* (FPN) n )   ] + yvect[(int) ceil(0.75* (FPN) n ) ])/(FPN) 2.0;
    
    
    
    fprintf(outf, "\n--------------------");
    fprintf(outf, "\n      n = %9d",n);
    fprintf(outf, "\n Min(Y) = %9.6g",ymin);
    fprintf(outf, "\n     Q1 = %9.6g", q1);
    fprintf(outf, "\n Median = %9.6g",q2);
    fprintf(outf, "\n     Q3 = %9.6g",q3);
    fprintf(outf, "\n Max(Y) = %9.6g", ymax);
    fprintf(outf, "\n--------------------\n\n");


    if (*outfile != '-')
      fileclose(outfile, outf);
  }
	 free_ivector(hist, 1, intervals);
	 free_dvector(yvect, 1, n);
}






long get_a_seed(void)
{
  time_t tptr;
  time(&tptr);
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  tptr = tptr / 2L;
#endif
  return ((long) tptr);
}

char *
     asctime2(void)
{
/*  static char time_buffer[MAXNAME];*/
  int i;
  time_t tptr;
  struct tm *tms;
  size_t len;
  for ( i=0;i<MAXNAME;i++ )
    time_buffer[i] = '\0';
  time(&tptr);
  tms = localtime(&tptr);
  len = strftime(time_buffer, MAXNAME, "%H:%M:%S on %A, %d %B %Y\n", tms);
#if defined(MACWARRIOR) || defined(WINWARRIOR) 
    return time_buffer;
#else
  if (len == 0)
    return NULL;
  else
    return time_buffer;
#endif
}


#if defined(DSIGN)
FPN dsign(FPN val1, FPN val2)
{
  FPN xx;
  xx = (FPN) fabs(val1);
  if (val2 < 0.0)
    xx = -xx;
  return (xx);
}
#endif



/*
 * dtranspose does an arbitrary transpose of one matrix
 * onto another. You must have allocated the memory that
 * mm1 and mm2 have pointed to.  lr, and lc are the initial
 * row and column while ur and uc are the final row and
 * column for the transpositon.  Note that mm1 =
 * dmatrix(lr,ur,lc,uc); mm2 = dmatrix(lc,uc,lr,ur); at
 * least.  These can be bigger, so that one can transpose
 * an arbitrary part of one matrix onto another.
 *
 * By Chris Basten, January 1994
 */
int dtranspose(FPN **mm1, FPN **mm2, int lr, int lc, int ur, int uc)
{
  int ii, jj;
  if (lr >= ur || lc >= uc)
    return (1);
  for (ii = lr; ii <= ur; ii++)
    for (jj = lc; jj <= uc; jj++)
      mm2[jj][ii] = mm1[ii][jj];
  return (0);
}



#if defined(ITOA)
/*
 * I think that every compiler has these now.  
 */
char *strupr(char *s)
{
  int len, ii;
  if ( s == NULL )
    return(s);
  len = (int) strlen(s);
  for (ii = 0; ii < len; ii++)
    if ( islower(s[ii]) )
      s[ii] = toupper(s[ii]);
  return s;
}

char *strlwr(char *s)
{
  int len, ii;
  if ( s == NULL )
    return(s);
  len = (int) strlen(s);
  for (ii = 0; ii < len; ii++) 
    if ( isupper(s[ii]) )
      s[ii] = tolower(s[ii]);
  return(s);
}

#endif


#if defined(DIVT)

div_t
div(int numer, int denom);
{
  div_t x;
  x.quot = 0;
  x.rem = 0;
  if (denom != 0) {
    x.quot = numer / denom;
    x.rem = numer - x.quot * denom;
  }
  else
    fprintf(stderr, "\n\nSorry, divide by zero in _div...\n");
  return x;
}

ldiv_t
ldiv(long numer,long denom);
{
  ldiv_t x;
  x.quot = 0;
  x.rem = 0;
  if (denom != 0) {
    x.quot = numer / denom;
    x.rem = numer - x.quot * denom;
  }
  else
    fprintf(stderr, "\n\nSorry, divide by zero in _ldiv...\n");
  return x;
}
#endif

FILE *
     fileopen(char *name, char *mode)
{
 /* Open files, writing an error message if it fails. */
  FILE *fp = fopen(name, mode);
  if (fp == NULL) {
    printf("Can't open %s for ", name);
    switch (mode[0]) {
     case 'r':
      printf("reading\n");
      break;
     case 'w':
      printf("writing\n");
      break;
     case 'a':
      printf("appending\n");
      break;
     default:
      printf("some strange mode\n");
      break;
    }
  }
  return fp;
}

void fileclose(char *name, FILE *fp)
{
 /* Close files, writing an error message if it fails. */
  if (fp != NULL && fclose(fp) == EOF)
    printf("Error closing %s\n", name);
}

void mypause(void)
{
  int k;
  char answer;
  k = 0;
  printf("\nType a <CR> to go on...\n");
  while ((answer = (char) getchar()) != '\n') k+=1;
}


void get_field(int xfield, char *xtemp, char *xbuffer)
{
  int i, j,k, field, look;
  k=0;
  for (i = 0; i < MAXNAME; i++)
    xtemp[i] = '\0';
  for (i = 0; xbuffer[i] != '\0'; i++) k+=1;
  xbuffer[i] = ' ';
  field = 0;
  look = 1;
  for (i = 0; field < xfield; i++) {
    if (xbuffer[i] != ' ' && look == 1) {
      field = field + 1;
      look = 0;
    }
    if (xbuffer[i] == ' ' && look == 0)
      look = 1;
  }
  for (j = i - 1; !isspace(xbuffer[j]); j++)
    xtemp[j-i+1] = xbuffer[j];
}

int Rotator(int thestep) {
  
      if ( thestep == 0 )
        printf(" ");
      else if ( thestep == 1 )
        printf("\b|");
      else if ( thestep == 2 )   
        printf("\b/");
      else if ( thestep == 3 )   
        printf("\b-");
      else if ( thestep == 4 )   
        printf("\b\\");
      else if ( thestep == 5 )   
        printf("\b|");
      else if ( thestep == 6 )   
        printf("\b/");
      else if ( thestep == 7 )   
        printf("\b-");
      else if ( thestep == 8 )  { 
        printf("\b\\");
        thestep = 0;
      } 
        
      fflush(stdout);
      return(thestep+1);
}





void ShowDVector(FPN *invec, int lb, int ub, int cols) {

  int cntr,i;
  
  cntr = 0;
  printf("\n\n");
  for ( i=lb; i<=ub; i++ ) {
    cntr +=1;
    printf(" %8.4g",invec[i]);
    if ( cntr == cols ) {
      printf("\n");
      cntr = 0;
    }
  
  }
}

void ShowIVector(int *invec, int lb, int ub, int cols) {

  int cntr,i;
  
  cntr = 0;
  printf("\n\n");
  for ( i=lb; i<=ub; i++ ) {
    cntr +=1;
    printf(" %8d",invec[i]);
    if ( cntr == cols ) {
      printf("\n");
      cntr = 0;
    }
  
  }
}

void ShowDMatrix(FPN **mm, int lr, int ur, int lc, int uc ) {

  int i,j;
  
  printf("\n\n");
  for ( i=lr; i<=ur; i++ ) {
    printf("\n");
    for ( j=lc; j<=uc; j++ )
      printf(" %6.4g",mm[i][j]);
  
  }
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Utilities.c
------------------------------------------------------------------ */


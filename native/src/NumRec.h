/* ------------------------------------------------------ XCutXCodeXSkip
     This file (NumRec.h) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */

/*
  NumRec.h header for NumRec.c
  
  All the subroutines in this file come from


    Numerical Recipes in C, The Art of Scientific Computing
  by W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling.
  1990, Cambridge University Press, Cambridge, England
  ISBN 0-521-35465-X  (the book)
  ISBN 0-521-35466-8  (MS-Dos diskette with C sources)


  They have been modified ony slightly by Chris J. Basten, on
  19 January 1994.  The modifications should not effect the behavior
  of the programs in any way.  The source has been reformatted to
  the liking of C.J. Basten, and the function definitions have been
  moved to "NumRec.h". All 'FPN' declarations were changed to 'FPN'.
  Finally, the uniform random number generator
  ranf() was substituted for ran1().  ranf() is an implementation of
  of a FORTRAN subroutine written by John Monohan at North Carolina
  State University.
*/




/*  The following are for xnumrec.c...*/
#define ITMAX 100
#define EPS 3.0e-7

void   indexx(int n, FPN arrin[],int indx[]);
void   moment(FPN data[],int n,FPN  *ave,FPN  *adev,FPN  *sdev,FPN  *svar,FPN  *skew,FPN  *curt);
void   sort(int n,FPN ra[]);
FPN gammln(FPN xx);
FPN gammp(FPN la,FPN  x);
FPN gasdev(int *idum);
void   gcf( FPN *gammcf,FPN la, FPN x, FPN *gln);
void   gser( FPN *gamser, FPN la, FPN x, FPN *gln);
FPN poidev(FPN xm,int *idum);
FPN betai(FPN la,FPN  b,FPN  x);
FPN beta(FPN z,FPN  w);
FPN betacf(FPN la,FPN  b,FPN  x);
FPN gammq( FPN la,FPN x);
FPN  chiprob( int     dof, FPN  value);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file NumRec.h
------------------------------------------------------------------ */


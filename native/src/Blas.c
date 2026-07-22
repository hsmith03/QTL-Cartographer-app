/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Blas.c) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */




#include "Blas.h"


/*


  The subroutines below were originally written in
   Fortran. They appeared in Dongarra, J.J, Moler, C.B.,
   Bunch, J.R. and Stewart G.W. 1979 LINPACK Users' Guide,
   SIAM, Philadelphia.
  
   The C translations were done by Christopher J. Basten at
   North Carolina State University in December of 1993.
   These are simple ports in that they do not unroll loops
   when the increments are both one.
  
 */




int srotg(FPN *sa, FPN *sb, FPN *cc, FPN *ss)

/*
 * construct givens plane rotation.  use pointers here
 * because the values change.  use function as: ind =
 * srotg(&sa,&sb,&cc,&ss);
 *                      TEST THIS ONE SOME MORE.  IT SEEMS TO WORK OK.
 * Since there are no arrays, this can be used as is.
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  FPN roe, scale, rr, zz;
  roe = *sb;
  if (fabs(*sa) > fabs(*sb))
    roe = *sa;
  scale = (FPN) (fabs(*sa) + fabs(*sb));
  if (scale != 0.0) {
    rr = scale * (FPN)  sqrt((*sa / scale) * (*sa / scale) + (*sb / scale) * (*sb / scale));
    if (roe < 0.0)
      rr = - (FPN)  1.0 * rr;
    *cc = *sa / rr;
    *ss = *sb / rr;
  }
  else {
    *cc = (FPN) 1.0;
    *ss = (FPN) 0.0;
    rr = (FPN) 0.0;
  }
  zz = (FPN) 1.0;
  if (fabs(*sa) > fabs(*sb))
    zz = *ss;
  if (fabs(*sb) >= fabs(*sa) && *cc != 0.0)
    zz = (FPN) 1.0 / *cc;
  *sa = rr;
  *sb = zz;
  return (0);
}



int srot(int nn, FPN *sx,int  isx, FPN *sy,int  isy, FPN cc, FPN ss)
/*
 * apply a plane rotation to sx and sy return(1) => error
 * return(0) => aok
 *                      TEST THIS ONE SOME MORE.  IT SEEMS TO WORK OK.
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  FPN stemp;
  int ii, jx, iy;
  if (nn == 0)
    return (1);
  if (isx < 0)
    jx = -(nn + 1) * isx + 1;
  else
    jx = 1;
  if (isy < 0)
    iy = -(nn + 1) * isy + 1;
  else
    iy = 1;
  for (ii = 1; ii <= nn; ii++) {
    stemp = cc * sx[jx]   + ss *  sy[iy] ;
    sy[iy] = cc * sy[iy] - ss * sx[jx];
   sx[jx] = stemp;
    jx = jx + isx;
    iy = iy + isy;
  }
  return (0);
}

int isamax(int nn, FPN *sx,int isx)
/* TESTED 1.11.94, SEEMS AOK.
 * Find index of element of sx with largest magnitude
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  int ii, jx, imax = 0;
  FPN smax;
  if (nn < 1)
    return (imax);
  imax = 1;
  if (nn == 1)
    return (imax);
  jx = 1;
  smax = (FPN) fabs(sx[jx]);
  jx = jx + isx;
  for (ii = 2; ii <= nn; ii++) {
    if (fabs(sx[jx]) > smax) {
      imax = jx;
      smax = (FPN) fabs(sx[jx]);
    }
    jx = jx + isx;
  }
  return (imax);
}

FPN sasum(int nn, FPN *sx, int isx)

/*  TESTED 1.11.94, SEEMS AOK.
 * sum the absolute values of vector sx. return(-1.0) =>
 * error return(x>=0.0) => aok
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  FPN stemp = (FPN) 0.0;
  int ii, nisx;
  if (nn <= 0)
    return ( (FPN)  -1.0 );
  nisx = nn * isx;
  for (ii = 1; ii <= nisx; ii = ii + isx)
    stemp = stemp + (FPN) fabs(sx[ii]);
  return (stemp);
}

FPN ssum(int nn, FPN *sx,int isx)
/*  TESTED 1.11.94, SEEMS AOK.
 * sum the   values of vector sx. return(-1.0) =>
 * error return(x>=0.0) => aok
 *
 * Based on a subroutine sasum written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  FPN stemp = (FPN)  0.0;
  int ii, nisx;
  if (nn <= 0)
    return ( (FPN) -1.0);
  nisx = nn * isx;
  for (ii = 1; ii <= nisx; ii = ii + isx)
    stemp = stemp + sx[ii];
  return (stemp);
}

int saxpy(int nn, FPN sa, FPN *sx,int isx, FPN *sy,int isy)
/* TESTED 1.11.94, SEEMS AOK.
 * do sy = a sx + sy (constant times a vector plus a
 * vector) putting the result where sy points to. return(1)
 * => error return(0) => aok
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  int ii, jx, iy;
  if (nn == 0)
    return (1);
  if (isx < 0)
    jx = -(nn + 1) * isx + 1;
  else
    jx = 1;
  if (isy < 0)
    iy = -(nn + 1) * isy + 1;
  else
    iy = 1;
  for (ii = 1; ii <= nn; ii++) {
    sy[iy] = sy[iy] + sa * sx[jx];
    jx = jx + isx;
    iy = iy + isy;
  }
  return (0);
}

FPN ssaxpy(int nn, FPN sa, FPN *sx,int isx, FPN *sy,int isy)
/* TESTED 1.11.94, SEEMS AOK.
 * do sum ( a sx + sy ) (sum of a constant times a vector plus a
 * vector) returning the result
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  int ii, jx, iy;
  FPN sum;
  sum = (FPN) 0.0;
  if (nn == 0)
    return (1);
  if (isx < 0)
    jx = -(nn + 1) * isx + 1;
  else
    jx = 1;
  if (isy < 0)
    iy = -(nn + 1) * isy + 1;
  else
    iy = 1;
  for (ii = 1; ii <= nn; ii++) {
    sum = sum + sy[iy] + sa * sx[jx];
    jx = jx + isx;
    iy = iy + isy;
  }
  return (sum);
}

int scopy(int nn, FPN *sx,int isx, FPN *sy,int isy)
/* TESTED 1.11.94, SEEMS AOK.
 * copy the vector sx into sy return(1) => error return(0)
 * => aok
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  int ii, jx, iy;
  if (nn == 0)
    return (1);
  if (isx < 0)
    jx = -(nn + 1) * isx + 1;
  else
    jx = 1;
  if (isy < 0)
    iy = -(nn + 1) * isy + 1;
  else
    iy = 1;
  for (ii = 1; ii <= nn; ii++) {
    sy[iy] = sx[jx];
    jx = jx + isx;
    iy = iy + isy;
  }
  return (0);
}

FPN sdot(int nn, FPN *sx,int isx, FPN *sy,int isy)
/* TESTED 1.11.94, SEEMS AOK.
 * calculate the dot product of two vectors pointed to by
 * sx and sy return(0.0) => error return(else) => aok
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  FPN stemp = (FPN) 0.0;
  int ii, jx, iy;
  if (nn == 0)
    return ( (FPN)  0.0);
  if (isx < 0)
    jx = -(nn + 1) * isx + 1;
  else
    jx = 1;
  if (isy < 0)
    iy = -(nn + 1) * isy + 1;
  else
    iy = 1;
  for (ii = 1; ii <= nn; ii++) {
    stemp = stemp + sx[jx] * sy[iy];
    jx = jx + isx;
    iy = iy + isy;
  }
  return (stemp);
}


FPN shadamer(int nn, FPN *sx, int isx, FPN *sy, int isy, FPN *sr,int isr)
/*  
 * calculate the Hadamer product of two vectors pointed to by
 * sx and sy.  Put the product into the vector pointed to by
   sr.   

   return(0.0)  => error ( nn = 0)
   return(else) => aok

This one was written by Chris Basten 30 April 1996.
 */
{
  int ii, jx, iy,ir;
  if (nn == 0)
    return ((FPN)  0.0);
  if (isx < 0)
    jx = -(nn + 1) * isx + 1;
  else
    jx = 1;
  if (isy < 0)
    iy = -(nn + 1) * isy + 1;
  else
    iy = 1;
  if (isr < 0)
    ir = -(nn + 1) * isr + 1;
  else
    ir = 1;
  for (ii = 1; ii <= nn; ii++) {
    sr[ir] =   sx[jx] * sy[iy];
    jx = jx + isx;
    iy = iy + isy;
    ir = ir + isr;
  }
  return ((FPN)  1.0);
}


FPN snrm2(int nn, FPN *sx, int isx)
/* TESTED 1.11.94, SEEMS AOK.
 * Calculate the Euclidean norm of sx.  Scaling in cases
 * where the sum would be too small or large is accounted
 * for. return(0.0) => nn <= 0 return(-1.0) => isx < 1
 * return(else) => aok
 *
 * Based on a subroutine written by C.L. Lawson in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  int ii, ub, max_val;
  FPN absmax, sfctr = (FPN) 1.0, total = (FPN) 0.0;
  if (nn <= 0)
    return ((FPN) 0.0);
  else if (isx < 1)
    return (- (FPN)  1.0);
  max_val = isamax(nn, sx, isx);
  absmax = (FPN) fabs(sx[max_val]);
  if (absmax < sqrt(DBL_MIN / DBL_EPSILON))
    sfctr = absmax;
  else if (absmax > sqrt(DBL_MAX) / (FPN) nn)
    sfctr = (FPN) 1.0 / absmax;
  ub = isx * nn;
  for (ii = 1; ii <= ub; ii = ii + isx)
    total = total + sfctr * (FPN)  pow(sx[ii], 2.0);
  if (total > 0.0)
    total = (FPN)  sqrt(total) / sfctr;
  return (total);
}

int sscal(int nn, FPN sa, FPN *sx, int isx)
/* TESTED 1.11.94, SEEMS AOK.
 * scale the FPN vector pointed to by sx by the factor
 * sa return(0) is no errors. return(1) indicates someting
 * went wrong.
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  int ii, ub;
  if (nn <= 0 || sa == 0.0)
    return (1);
  ub = isx * nn;
  for (ii = 1; ii <= ub; ii = ii + isx)
    sx[ii] = sa * sx[ii];
  return (0);
}


int sswap(int nn, FPN *sx,int isx, FPN *sy,int isy)
/* TESTED 1.11.94, SEEMS AOK.
 * a simple swap of the elements of the FPN vectors
 * pointed to by sx and sy.  the increments can be any
 * integers. return(0) is no errors. return(1) indicates
 * someting went wrong.
 *
 * Based on a subroutine written by Jack Dongarra in Fortran.
 * The original is in Dongarra, J.J, Moler, C.B., Bunch,
 * J.R. and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 */
{
  FPN stemp;
  int ii, jx, iy;
  if (nn == 0)
    return (1);
  if (isx < 0)
    jx = -(nn + 1) * isx + 1;
  else
    jx = 1;
  if (isy < 0)
    iy = -(nn + 1) * isy + 1;
  else
    iy = 1;
  for (ii = 1; ii <= nn; ii++) {
    stemp = sx[jx];
    sx[jx] = sy[iy];
    sy[iy] = stemp;
    jx = jx + isx;
    iy = iy + isy;
  }
  return (0);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Blas.c
------------------------------------------------------------------ */


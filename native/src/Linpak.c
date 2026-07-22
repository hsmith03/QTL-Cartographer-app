/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Linpak.c) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */


#include "Main.h"


/*   
 *
 * Comments in lower case are from the Fortran source.
 * Those in upper case are questions to myself.  Upper case
 * comments should be removed at some point.
 *
 * The subroutines below were originally written in Fortran.
 * They appeared in Dongarra, J.J, Moler, C.B., Bunch, J.R.
 * and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.
 *
 * The C translations were done by Christopher J. Basten at
 * North Carolina State University in January of 1994.
 *
 * These C translations work with the transpose of the matrices
 * relative to the original Fortran subroutines.  The comments are
 * a bit mixed up as a result.  For  solutions to the system
     xx * bb = yy  (1)
 * use the transpose of xx and all should work out fine.
 */


/*
 * sqrdc uses Householder transformations to compte the QR
 * factorization of an pp by nn matrix yy. Row pivoting
 * based on the 2-norms of the reduced rows may be
 * performed at the users option.
 *
 * On entry
 *
 * yy = dmatrix(1,pp,1,ldx) requires ldx >= nn.  It contains
 *      the matrix whose decompostion is to be computed.
 *
 * ldx  is the leading dimension of yy
 *
 * nn is the number of rows of yy
 *
 * pp is the number of rows of yy
 *
 * jpvt = ivector(1,pp) contains intergers that control the
 *        selection of the pivot rows.  The kk-th row of yy
 *        is placed in one of three classes according to the value
 *        of *(jpvt+kk).
 *
 *    *(jpvt+kk) > 0      kk-th row is an initial row
 * if *(jpvt+kk) = 0 then kk-th row is a free row
 *    *(jpvt+kk) < 0      kk-th row is a final row
 *
 * Before the decomposition is computed, initial rows are
 * moved to the beginning of the matrix xx and final
 * rows to the end.  Both initial and final rows are
 * frozen in place during the compution and only free
 * rows are moved.  At the kk-th stage of the reduction,
 * if the kk-th row is a free one it is interchanged
 * with the free row of largest reduced norm.  jpvt is
 * not referenced if job = 0
 *
 *
 * job is an integer that initiates row pivoting.
 * if job == 0 no pivoting
 * if job != 0 pivot
 *
 * On return
 *
 * yy  contains in its upper triangle the upper triangular
 *     matrix R of the QR factorization. Below its diagonal xx
 *     contains information from which the orthogaonal part to
 *     the decomposition can be recovered.  Note that if
 *     pivoting has been requested, the decomposition is not
 *     that of the origial matrix yy but that of yy with its
 *     rows permuted as described by jpvt.
 *
 * qraux = dvector(1,pp) contains further information required
 *         to recover the orthogonal part of the decomposition.
 *
 * jpvt   the kk-th element contains the index of the row
 *        of the original matrix that has been interchanged into
 *        the kk-th row, if pivoting was requested.
 *
 * LINPACK.  This version dated 08/14/78. G.W. Stewart,
 * University of Maryland, Argonne National Lab. Translated
 * into C by Chris Basten 1/3/94 at North Carolina State
 * University, Departmen of Statistics.
 *
 * The following are used: from blas.c:  saxpy, sdot, sscal,
 * sswap, snrm2
 */

int sqrdc(FPN **yy,int ldx,int nn,int pp, FPN *qraux,int *jpvt,int job)
{
/* work = dvector(1,pp) is a work vector.  work is not
 * referenced if job = 0. */
  int j, jj, jp, ll, lp1, lup, maxj, pl, pu, error;
  FPN maxnrm, t, tt, nrmxl, *work;
  int negj, swapj;	/* these will be logicals, i.e. = 0 (false) or = 1 (true) */
  int doofus;
      doofus=ldx;
  pl = 1;
  pu = 0;

  work = dvector(1, pp);

  if (job != 0) {	/* pivot rows according to jpvt as requested */
    for (jj = 1; jj <= pp; jj++) {
      if ( jpvt[jj] > 0)
	    swapj = 1;
      else
	    swapj = 0;
      if (jpvt[jj] < 0)
	    negj = 1;
      else
	    negj = 0;
      jpvt[jj] = jj;
      if (negj == 1)
	    jpvt[jj] = -jj;
      if (swapj == 1) {
	    if (jj != pl)
	      error = sswap(nn, yy[pl], 1,  yy[jj], 1);
	    jpvt[jj] = jpvt[pl];
	    jpvt[ pl]  = jj;
	    pl = pl + 1;
      }
    }
    pu = pp;
    for (jj = 1; jj <= pp; jj++) {
      j = pp - jj + 1;
      if (jpvt[j] < 0) {
	    jpvt[j] = -jpvt[j];
	    if (j != pu) {
	       error = sswap(nn,  yy[pu], 1, yy[j], 1);
	       jp = jpvt[pu ] ;
	       jpvt[pu]  = jpvt[j] ;
	       jpvt[j]  = jp;
	    }
	    pu = pu - 1;
      }
    }
  }
  if (pu >= pl)	/* compute norms of the free rows */
    for (jj = pl; jj <= pu; jj++)
      work[jj] = qraux[jj] = snrm2(nn, yy[jj], 1);
 /* do the Householder transformation */
  if (nn < pp)
    lup = nn;
  else
    lup = pp;
  for (ll = 1; ll <= lup; ll++) {
    if (ll >= pl && ll < pu) {
      maxnrm = (FPN) 0.0;
      maxj = ll;
      for (j = ll; j <= pu; j++)
	    if (qraux[j] > maxnrm) {
	      maxnrm = qraux[j];
	      maxj = j;
	    }
      if (maxj != ll) {
	    error = sswap(nn, yy[ll], 1, yy[maxj], 1);
	    qraux[maxj] = qraux[ll];
	    work[maxj] = work[ll];
	    jp = jpvt[maxj];
	    jpvt[maxj] = jpvt[ll];
	    jpvt[ll] = jp;
      }
    }
    qraux[ll] = (FPN)  0.0;
    if (ll != nn) {	/* Compute the Householder transformation for row ll. */
      nrmxl = snrm2(nn - ll + 1, yy[ll] + ll - 1, 1);
      if (nrmxl != (FPN) 0.0) {
	    if (yy[ll][ll] != (FPN) 0.0)
	      nrmxl = dsign(nrmxl, yy[ll][ll]);
	    error = sscal(nn - ll + 1, (FPN) 1.0 / nrmxl, yy[ll] + ll - 1, 1);

	    yy[ll][ll] = (FPN) 1.0 + yy[ll][ll];
      /* Apply the transformation to the remaining rows,
         updating the norms */
	lp1 = ll + 1;
	if (pp >= lp1)
	  for (j = lp1; j <= pp; j++) {
	    tt = -sdot(nn - ll + 1, yy[ll] + ll - 1, 1, yy[j] + ll - 1, 1) / yy[ll][ll];
	    error = saxpy(nn - ll + 1, tt, yy[ll] + ll - 1, 1, yy[j] + ll - 1, 1);
	    if (j >= pl && j <= pu)
	      if (qraux[j] != (FPN) 0.0) {
		tt = yy[j][ll] / qraux[j];
		tt = 1 - tt * tt;
		if (tt < (FPN) 0.0)
		  tt = (FPN) 0.0;
		t = tt;
		tt = qraux[j] / work[j];
		tt = (FPN) 1.0 + (FPN) 0.05 * t * tt * tt;
		if (tt == (FPN) 1.0) {
		  qraux[j] = snrm2(nn - ll, yy[j] + ll, 1);
		  work[j] = qraux[j];
		}
		else
		  qraux[j] = qraux[j] * (FPN)sqrt(t);
	      }
	  }
	qraux[ll] =  yy[ll][ll];
	yy[ll][ll] = -nrmxl;
      }
    }
  }
  free_dvector(work, 1, pp);
  return (0);
}




/*
 * sqrsl applies the output of sqrdc to compute coordinate
 * transformations, projections, and least squares
 * solutions.  for kk <= min(pp,nn), let xk be the matrix
 *
 * xk = [ *(yy+jpvt(1)), *(yy+jpvt(2)), ..., *(yy+jpvt(kk)) ]
 *
 * formed from the rows jpvt of the original pp by nn
 * matrix yy that was input to sqrdc (if no pivoting was
 * done, xk consists of the first kk rows of yy in their
 * original order.  sqrdc produces a factored orthogonal
 * matrix Q and an upper triangular matrix R such that
 *
 *
 * transpose(xk) = Q * (R)
                       (0)
 *
 *
 * this information is contained in coded form in the
 * arrays yy and quaux. xk = dmatrix(1,kk,1,nn);
 *
 *
 * On entry
 *
 * yy    = dmatrix(1,pp,1,ldx) contains the output of sqrdc.
 *
 * pp      is the number of rows in matrix yy.
 *
 * ldx     is the leading dimension of the array yy.
 *
 * nn      is the number of columns of the matrix xk.  It must
 *         have the same value as nn in sqrdc.
 *
 * kk      is the number of rows of the matrix xk.  kk must
 *         not be greater than min(nn,pp), where pp is the same as
 *         in the calling sequence ot sqrdc.
 *
 * qraux = dvector(1,pp) contains the auxiliary output from sqrdc.
 *
 * y     = dvector(1,nn) contains an nn-vector that is to be
 *         manipulated by sqrsl.
 *
 * job     specifies what is to be computed. job has the
 *         decimal expansion ABCDE, with the following meaning:
 *
 *
    If           then      Compute
  A != 0                    qy
B,C,D || E != 0             qty
  C != 0                    bb
  D != 0                    rsd
  E != 0                    xb
 *
 * Note that a request to compute bb, rsd, or xb
 * automatically triggers the computation of qty, for which
 * an array must be provided in the calling sequence.
 *
 *
 * On return
 *
 * qy   = dvector(1,nn) contans qq*yy, if its computation has
 *        been requested.
 *
 * qty  = dvector(1,nn) contains trans(qq)*yy if its
 *        computation has been requested. Here trans(qq) is the
 *        transpose of the matrix qq.
 *
 * bb   = dvector(1,kk) contains the soution of the least
 *        squares problem Minimize Norm2(yy - xk*bb) If its
 *        computation has been requested.  (Note that if pivoting
 *        was requested in sqrdc, the j-th component of bb will be
 *        associated with row *(jpvt+j) of the original matrix
 *        xx that was input into sqrdc.)
 *
 * rsd  = dvector(1,nn) contains the least squares residual yy
 *        - xk*bb, if its computation has been requested.  rsd is
 *        also the orthogonal projection of yy onto the orthogonal
 *        complement of the row space of xk.
 *
 * xb   = dvector(1,nn)  contains the least squares
 *        aproximation xk*bb, if its computation has been
 *        requested.  xb is also the orthogonal projection of yy
 *        onto the row space of xx.
 *
 * info   is zero unless the computation of bb has been
 *        requested and rr is exactly singular. In this case, info
 *        is the index of the first zero diagonal element of rr
 *        and bb is left unaltered.
 *
 * The parameteres qy, qty, bb rsd and xb are not referencd if
 * their computation is not requested and in this case can
 * be replaced by dummy variables in the calling program.
 * To save storage, the user may in some cases use the smae
 * array for different parameters in the calling sequence.
 * A frequently occuring example is when one wishes to
 * compute any of bb, rsd or xb and does not need yy or
 * qty.  In this case one may identify yy, qty and one of
 * bb, rsd or xb while provideing separate arrays for
 * anything else that is to be computed. Thus the calling
 * sequence
 *
 * info = sqrsl(xx,ldx,nn,kk,qrraux,y,dum,y,bb,y,dum,110);
 *
 * will result in the computation of bb and rsd, with rsd
 * overwriting yy.  More generally, each item in the
 * following list contains groups of permissible
 * identifications of a single calling sequence.
 *
 *
   1. (y, qty, bb) (rsd) (xb) (qy)
   2. (y, qty, rsd) (bb) (xb) (qy)
   3. (y, qty, xb) (bb) (rsd) (qy)
   4. (y, qy) (qty, bb) (rsd) (xb)
   5. (y, qy) (qty, rsd) (bb) (xb)
   6. (y, qy) (qty, xb) (bb) (rsd)
 *
 *
 * In any group the value returned in the array allocated
 * to the group corresponds to the last member of the
 * group.
 *
 *
 * Translated by Chris Basten into c in January 1994 form
 * Linpak, version dated 08/14/78 G.W. Stewart, University
 * of Maryland, Argonne National Lab.
 *
 * sqrsl uses the following functions
 *
 * blas.c  saxpy, scopy, sdot
 */

int sqrsl(FPN **yy,int ldx,int nn,int pp,int kk, FPN *qraux, FPN *y, FPN *qy, FPN *qty, FPN *bb, FPN *rsd, FPN *xb,long job)
{
  ldiv_t test;
  long tmp;
  int i, j, jj, ju, kp1, info, error;
  FPN t, temp;
  int cb, cqy, cqty, cr, cxb;	/* these will be logicals, i.e. = 0 (false) or = 1 (true) */
  int doofus;
      doofus=pp;
      doofus=ldx;
  info = 0;	/* set info flag */
 /* Determine what is to be computed */
  cb = cqy = cqty = cr = cxb = 0;
  tmp = job;
  test = ldiv(tmp, 10);
  if (test.rem != 0)
    cxb = cqty = 1;
  tmp = tmp / 10;
  test = ldiv(tmp, 10);
  if (test.rem != 0)
    cr = cqty = 1;
  tmp = tmp / 10;
  test = ldiv(tmp, 10);
  if (test.rem != 0)
    cb = cqty = 1;
  tmp = tmp / 10;
  test = ldiv(tmp, 10);
  if (test.rem != 0)
    cqty = 1;
  tmp = tmp / 10;
  test = ldiv(tmp, 10);
  if (test.rem != 0)
    cqy = 1;

/*
  printf("\nNow in sqrsl...\n job = %3d ",job);
  if ( cxb == 1 ) printf("cxb ");
  if ( cqty == 1 ) printf("cqty ");
  if ( cr == 1 ) printf("cr ");
  if ( cb == 1 ) printf("cb ");
  if ( cqy == 1 ) printf("cqy ");
  printf("\n\n");
*/
  ju = nn - 1;
  if (ju > kk)
    ju = kk;
 /* nn = 1 is a special case. */
  if (ju == 0) {
    if (cqy == 1)
      qy[1] = y[1];
    if (cqty == 1)
      qty[1] = y[1];
    if (cxb == 1)
      xb[1] = y[1];
    if (cb == 1) {
      if (yy[1][1] == (FPN) 0.0)
	    info = 1;
      else
	    bb[1] = y[1] / yy[1][1];
    }
    if (cr == 1)
      rsd[1] = (FPN) 0.0;
    return (info);
  }
  if (cqy == 1)
    error = scopy(nn, y, 1, qy, 1);
  if (cqty == 1)
    error = scopy(nn, y, 1, qty, 1);
  if (cqy == 1)	/* Compute qy */
    for (jj = 1; jj <= ju; jj++) {
      j = ju - jj + 1;
      if (qraux[j] != (FPN) 0.0) {
	    temp = yy[j][j];
	    yy[j][j] = qraux[j];
	    t = -sdot(nn - j + 1, *(yy + j) + j - 1, 1, (qy + j - 1), 1) / yy[j][j];
	    error = saxpy(nn - j + 1, t, *(yy + j) + j - 1, 1, (qy + j - 1), 1);
	    yy[j][j] = temp;
      }
    }
  if (cqty == 1)	/* Compute trans(q)*y */
    for (j = 1; j <= ju; j++)
      if (qraux[j] != (FPN) 0.0) {
	temp = yy[j][j];
	yy[j][j] = qraux[j];
	t = -sdot(nn - j + 1, *(yy + j) + j - 1, 1, (qty + j - 1), 1) / yy[j][j];
	error = saxpy(nn - j + 1, t, *(yy + j) + j - 1, 1, (qty + j - 1), 1);
	yy[j][j] = temp;
      }
 /* Set up to compute bb, rsd, or xb */
  if (cb == 1)
    error = scopy(kk, qty, 1, bb, 1);
  kp1 = kk + 1;
  if (cxb == 1)
    error = scopy(kk, qty, 1, xb, 1);
  if (cr == 1 && kk < nn)
    error = scopy(nn - kk, (qty + kk), 1, (rsd + kk), 1);
  if (cxb == 1 && kp1 <= nn)
    for (i = kp1; i <= nn; i++)
      xb[i] = (FPN) 0.0;
  if (cr == 1)
    for (i = 1; i <= kk; i++)
      rsd[i] = (FPN) 0.0;
  if (cb == 1) {	/* compute bb */
    jj = 1;
    while (jj <= kk) {
      j = kk - jj + 1;
      if (yy[j][j] != (FPN) 0.0) {
	bb[j] = bb[j] / yy[j][j];
	if (j != 1) {
	  t = - bb[j];
	  error = saxpy(j - 1, t, *(yy + j), 1, bb, 1);
	}
      }
      else {
	info = j;
	jj = kk;
      }
      jj = jj + 1;
    }
  }
  if (cr == 1 || cxb == 1)	/* Compute rsd or xb as required */
    for (jj = 1; jj <= ju; jj++) {
      j = ju - jj + 1;
      if (qraux[j] != (FPN) 0.0) {
	temp = yy[j][j];
	yy[j][j] = qraux[j];
	if (cr == 1) {
	  t = -sdot(nn - j + 1, *(yy + j) + j - 1, 1, (rsd + j - 1), 1) / yy[j][j];
	  error = saxpy(nn - j + 1, t, *(yy + j) + j - 1, 1, (rsd + j - 1), 1);
	}
	if (cxb == 1) {
	  t = -sdot(nn - j + 1, *(yy + j) + j - 1, 1, (xb + j - 1), 1) / yy[j][j];
	  error = saxpy(nn - j + 1, t, *(yy + j) + j - 1, 1, (xb + j - 1), 1);
	}
	yy[j][j] = temp;
      }
    }
  return (info);
}








/*********************************************************************************

         The following was translated into c by Chris Basten in January, 1994.

C
C     SPODI computes the determinant and inverse of a certain real
C     symmetric positive definite matrix (see below) using the factors
C     computed by SPOCO, SPOFA or SQRDC.
C
C     On entry:
C
C        a       = dmatrix(1,lda,1,n)
C                The output A from SPOCO or SPOFA or the output X from
C                SQRDC.
C
C        lda     INTEGER.
C                The leading dimension of the array a.
C
C        n       INTEGER.
C                The order of the matrix A.
C
C        job     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On return:
C
C        a       If SPOCO or SPOFA was used to factor A then SPODI
C                produces the upper half of inverse(A). If SQRDC was
C                used to decompose X then SPODI produces the upper half
C                of inverse(trans(X)*X) where trans(X) is the transpose.
C                Elements of A below the diagonal are unchanged. If the
C                units digit of JOB is zero, A is unchanged.
C
C        det     dvector(1,2)
C                Determinant of A or of trans(X)*X if requested.
C                Otherwise not referenced. Determinant =
C                DET(1)*10.0**DET(2) with 1.0 <= DET(1) < 10.0
C                or DET(1) = 0.0.
C
C     Error condition
C
C        A division by zero will occur if the unput factor contains a
C        zero on the diagonal and the inverse is requested. It will not
C        occur if the subroutines are called correctly and if SPOCO or
C        SPOFA has set INFO.EQ.0.
C
C     LINKPACK. This version dated 08/14/78.
C     CLEVE MOLER, University of New Mexico, Argonne National Lab.
C
C     Subroutines and functions.
C
C     BLAS SAXPY, SSCAL
C     FORTRAN MOD
C
*/

int spodi(FPN **yy, int lda, int n, FPN *det, int job)
{
  div_t test;
  int error;
  FPN t, s;
  int i, j, jm1, k, kp1,doofus;
      doofus=lda;

  test = div(job, 10);
  if (test.quot != 0) {	/* Compute determinant */
    det[1] = (FPN) 1.0;
    det[2] = (FPN) 0.0;
    s = (FPN)10.0;
    for (i = 1; i <= n; i++) {
      det[1] = (FPN)pow(yy[i][i], 2.0) * det[1];
      if (det[1] == (FPN) 0.0)
	    break;	/* ...exit */
      while (det[1] < (FPN) 1.0) {
	    det[1] = s * det[1];
	    det[2] = det[2] - (FPN) 1.0;
      }
      while (det[1] >= s) {
	    det[1] = det[1] / s;
	    det[2] = det[2] + (FPN) 1.0;
      }
    }
  }
  if (test.rem != 0) {	/* Compute inverse(R) */
    for (k = 1; k <= n; k++) {
      yy[k][k] = (FPN) 1.0 / yy[k][k];
      t = -yy[k][k];
      error = sscal(k - 1, t, *(yy + k), 1);
      kp1 = k + 1;
      if (n >= kp1) {
	for (j = kp1; j <= n; j++) {
	  t = yy[j][k];
	  yy[j][k] = (FPN) 0.0;
	  error = saxpy(k, t, *(yy + k), 1, *(yy + j), 1);
	}

      }
    }
    for (j = 1; j <= n; j++) {	/* Form inverse(R)*trans(inverse(R)) */
      jm1 = j - 1;
      if (jm1 >= 1)
	for (k = 1; k <= jm1; k++) {
	  t = *(*(yy + j) + k);
	  error = saxpy(k, t, *(yy + j), 1, *(yy + k), 1);
	}
      t = *(*(yy + j) + j);
      error = sscal(j, t, *(yy + j), 1);
    }
  }
  return (error);
}



/****************************************************************************


C
C     STRSL solves systems of the form
C            T * X = B
C     or     trans(T) * X = B
C     where T is a triangular matrix of order N. Here trans(T) denotes
C     the transpose of the matrix T.
C
C     On entry:
C
C      **t       = dmatrix(ldt,n)
C                t contains the matrix of the system. The zero elements
C                of matrix are not referenced, and the corresponding
C                elements of the array can be used to store other
C                information.
C
C        ldt     INTEGER.
C                ldt is the leading dimension of the array t.
C
C        n       INTEGER.
C                n is the order of the system.
C
C       *b       dvector(1,n)
C                b points to the right hand side of the system.
C
C        job     INTEGER
C                job specifies what kind of system is to be solved. If
C                job is
C                  00   solve T*X=B, T lower triangular,
C                  01   solve T*X=B, T upper triangular,
C                  10   solve trans(T)*X=B, T lower triangular,
C                  11   solve trans(T)*X=B, T upper triangular.
C
C     On return:
C
C       *b       b points to the solution, if INFO == 0. Otherwise B is
C                unaltered.
C
C        info    INTEGER
C                info contains zero if the system is nonsingular.
C                Otherwise INFO contains the index of the first zero
C                diagonal element of T.
C
C     LINKPACK. This version dated 08/14/78.
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C
C     Subroutines and functions.
C
C     BLAS SAXPY, SDOT
C     FORTRAN MOD
C
*/
int strsl(FPN **yy, int ldt, int n, FPN *b, int job)
{
  FPN temp;
  int j, jj, info, error,doofus;
      doofus=ldt;
 /* Begin block permitting. Check for zero diagonal
    elements, and if any, Exit */
  for (info = 1; info <= n; info++)
    if (*(*(yy + info) + info) == (FPN) 0.0)
      return (info);
  info = 0;

  if (job == 0) {	/* Solve T*X=B for T lower triangular */
    *(b + 1) = *(b + 1) / *(*(yy + 1) + 1);
    if (n >= 2)
      for (j = 2; j <= n; j++) {
	temp = -*(b + j - 1);
	error = saxpy(n - j + 1, temp, *(yy + j - 1) + j - 1, 1, b + j - 1, 1);
	*(b + j) = *(b + j) / *(*(yy + j) + j);
      }
  }
  else if (job == 1) {	/* Solve T*X=B for T upper triangular. */
    *(b + n) = *(b + n) / *(*(yy + n) + n);
    if (n >= 2)
      for (jj = 2; jj <= n; jj++) {
	j = n - jj + 1;
	temp = -*(b + j + 1);
	error = saxpy(j, temp, *(yy + j + 1), 1, b, 1);
	*(b + j) = *(b + j) / *(*(yy + j) + j);
      }
  }
  else if (job == 10) {	/* Solve trans(T)*X=B for T lower triangular. */
    *(b + n) = *(b + n) / *(*(yy + n) + n);
    if (n >= 2)
      for (jj = 2; jj <= n; jj++) {
	j = n - jj + 1;
	*(b + j) = *(b + j) - sdot(jj - 1, *(yy + j) + j, 1, (b + j), 1);
	*(b + j) = *(b + j) / *(*(yy + j) + j);
      }
  }
  else if (job == 11) {	/* Solve trans(T)*X=B for T upper triangular. */
    *(b + 1) = *(b + 1) / *(*(yy + 1) + 1);
    if (n >= 2)
      for (j = 2; j <= n; j++) {
	*(b + j) = *(b + j) - sdot(j - 1, *(yy + j), 1, b, 1);
	*(b + j) = *(b + j) / *(*(yy + j) + j);
      }
  }
  return (info);
}


/*
 * sqrst is a function to compute the least squares
 * solutions to the system
     xx * bb = yy  (1)
 * which may be
 * either under-determined or over-determined.  The user
 * may supply a tolerance to limit the rows of xx used
 * in computing the solution.  In effect, a set of rows
 * with a condition number approximately bounded by 1/tol
 * is tused, the other components of bb being set to zero.
 *
 * On entry
 *
 * xx  = dmatrix(1,pp,1,ldx) is a pointer to the
 *       FPN matrix xx(pp,ldx), where ldx >= nn. xx contains
 *       the pp by nn coefficient matrix of the system (1), and
 *       xx is destroyed by sqrst.
 *
 * ldx is the leading dimension of xx.
 *
 * nn  is the number of columns in xx.
 *
 * pp  is the number of rows in xx.
 *
 * yy  = dvector(1,nn), is the RHS of (1)
 *
 * tol is the nonnegative tolerance used to determine the
 *     subset of rows of xx included in the solution.  if
 *     tol is zero, a full complement of min(nn,pp) rows is
 *     used.
 *
 * jpvt = ivector(1,pp)  is a vector used by sqrdc.
 *
 * qraux = dvector(1,pp) is a vector used by sqrdc and sqrsl.
 *
 *
 * On return
 *
 * xx points the output array from sqrdc.
 *
 * bb = dvector(1,pp) contains the solution vector.
 *      components corresponding to rows not used are set to
 *      zero.
 *
 * rsd = dvector(1,nn) contains the residual vector yy-xx*bb.
 *
 * kk  contains the number of rows used in the solution.
 *     use &kk so that this will be changed.
 *
 * jpvt contains the pivot information from sqrdc
 *
 * qraux contains the vector output by sqrdc.
 *
 * On return the arrays pointed to by xx, jpvt and qraux
 * contain the usual output from sqrdc, so that the QR
 * decompositon of xx with pivoting is fully available to
 * the user.  In particular, rows (ROWS) *(jpvt+1),
 * *(jpvt+2), ...,*(jpvt+kk) were used in the solution, and
 * the condition number associated with those rows is
 * estimated by fabs(*(*(xx+1)+1)/    *(*(xx+kk)+kk)).
 *
 * sqrst uses the linpack functions sqrdc and sqrsl.  In
 * essence, it is a driver for these subroutines.
 */
int sqrst(FPN **xx, int ldx, int nn, int pp, FPN *yy, FPN tol, FPN *bb, FPN *rsd, int *kk, int *jpvt, FPN *qraux)
{
  int info, ii, jj, error, mm, k;
  FPN tt;
  for (jj = 1; jj <= pp; jj++)
    *(jpvt + jj) = 0;	/* clearout jpvt */
  error = sqrdc(xx, ldx, nn, pp, qraux, jpvt, 1);	/* reduce xx */
  if (error == 1)
    return (error);
  ii = 1;
  *kk = 0;
  if (nn < pp)
    mm = nn;
  else
    mm = pp;
  for (ii = 1; ii <= mm; ii++) {
    if (fabs(*(*(xx + ii) + ii)) <= tol * fabs(*(*(xx + 1) + 1)))
      break;
    *kk = ii;
  }
  k = *kk;
  if (k != 0)
    info = sqrsl(xx, ldx, nn, pp, k, qraux, yy, rsd, rsd, bb, rsd, rsd, 110L);
  if (info > 0)
    return (info);
  for (jj = 1; jj <= pp; jj++) {
    *(jpvt + jj) = *(jpvt + jj) * -1;
    if (jj > *kk)
      *(bb + jj) = (FPN) 0.0;
  }	/* set unused components of bb to zero and initialize jpvt for unscrambling */
  for (jj = 1; jj <= pp; jj++)
    if (*(jpvt + jj) <= 0) {
      *(jpvt + jj) = *kk = -*(jpvt + jj);
      while (*kk != jj) {
	tt = *(bb + jj);
	*(bb + jj) = *(bb + *kk);
	*(bb + *kk) = tt;
	*(jpvt + *kk) = -*(jpvt + *kk);
	*kk = *(jpvt + *kk);
      }
    }
  return (0);
}


/*c********************************

Chapter 3, page 3.

c      spofa factors a real symmetric positive definite matrix.
c
c      spofa is usually called by spoco, but it can be called derectly with 
c      a saving in time if recond is not needed.
c      (time for spoco)=(1+18/n)*(time for spofa).
c
c      on entry:
c
c      a     real(lda,n)
c            The symmetric matrix to be factored. Only the diagonal and
c            upper triangle are used.
c
c      lda   integer.
c            The leading dimension of the array a.
c
c      n     integer.
c            The order of the matrix a.
c
c      on return:
c
c      a     An lower triangular matrix r so that a=trans(r)*r where 
c            trans(r) is the transpose.
c            The strict upper triangle is unaltered. If info.ne.0, the
c            factorization is not complete.
c
c      info  integer
c            = 0 for normal return.
c            = k signal an error condition. The leading minor of order k
c              is not positive definite.
c
c      linpack. This version dated 08/14/78.
c      Cleve Moler, University of New Mexico, Argonne National Lab.
c
c      subroutines and functions
c
c      blas, sdot
c      fortran sqrt
c
c      internal variables
c
*/
int spofa(FPN **a, int lda, int n)
{
  FPN t,s;
  int  j,jm1,k,info,doofus;
      doofus=lda;
/*          */
  for ( j = 1 ; j <= n ; j++ ) {
    info = j;
    s = (FPN) 0.0;
    jm1 = j-1;
    if ( jm1 >= 1) {
      for ( k = 1 ; k <= jm1 ; k++) {
        t = a[k][j] - sdot(k-1,a[k],1,a[j],1);
        t = t/a[k][k];
        a[j][k] = t;
        s = s+t*t;
      }
    }
    s = a[j][j] - s;
    if ( s <= (FPN) 0.0 ) 
      j = n+2;
    else  
      a[j][j] = (FPN)sqrt(s);
  }
  if ( j == n+1 )
     info = 0;
  return(info);
}



/*
This explanation comes from page 3 of Chapter 5 of the LINPAK  Users Guide

INPUT RETURN

a = dmatrix(1,n,1,lda)  contains in its upper triangle (including diagonal) the information
                        necessary to construct a matrix U and block diagonal matrix D
                        so that a = UDU'
kpvt = ivector(1,n)
work = dvector(1,n)     work array.  If a is close to singular, then ||a z|| = rcond ||a|| ||z||
recond is an estimate of the reciprocal condition.  

If SSICO set rcond==0.0, or SSIF set info to nonzero, then no inverse should be requested.       

*/
int ssico(FPN **a, int lda, int n, int *kpvt, FPN *rcond, FPN *z)
{
  FPN ak,akm1,bk,bkm1,denom,ek,t;
  FPN anorm,s,ynorm;
  int i,info,j,jm1,k,kp,kps,ks;
/*
     find norm of a, using only lower half
*/
  for ( j = 1 ; j <= n ; j++ ) {
    z[j] = sasum( j,a[j],1);
    jm1 = j-1;
    if ( jm1 >= 1 ) 
      for ( i=1 ; i <= jm1 ; i++ )
            z[i] = z[i] + (FPN) fabs(a[j][i]);
  }
  anorm = (FPN) 0.0;
  for ( j=1 ; j<=n ; j++ )
    anorm = amax1(anorm,z[j]);     
    info = ssifa(a,lda,n,kpvt);
    ek = (FPN) 1.0;
    for ( j=1 ; j<= n ; j++ )
         z[j]=(FPN) 0.0;
    k = n;
    
    if ( k != 0 ) {
      ks = 1;
      if ( kpvt[k] < 0) 
        ks = 2;
      kp = abs(kpvt[k]);
      kps = k+1-ks;
      if ( kp != kps ) {  
        t = z[kps];
        z[kps] = z[kp];
        z[kp] = t;
      }
      if ( z[k] != (FPN) 0.0) 
        ek = dsign(ek,z[k]);   
      z[k] = z[k] + ek;
      info = saxpy(k-ks,z[k],a[k],1,z ,1);
      if ( ks != 1) {
        if ( z[k-1] != (FPN) 0.0) 
          ek = dsign(ek,z[k-1]);
        z[k-1] = z[k-1] + ek;
        info = saxpy(k-ks,z[k-1],a[k-1],1,z,1);
      }
      if ( ks == 2 ) {
        ak = a[k][k]/a[k][k-1];
        akm1 = a[k-1][k-1]/a[k][k-1];
        bk = z[k]/a[k][k-1];
        bkm1 = z[k-1]/a[k][k-1];
        denom = ak*akm1-(FPN) 1.0;
        z[k] = (akm1*bk-bkm1)/denom;
        z[k-1] = (ak*bkm1-bk)/denom;
      }
      else {
        if ( fabs(z[k]) > fabs(a[k][k]) ) {
          s=(FPN)fabs(a[k][k])/(FPN)fabs(z[k]);
          info = sscal(n,s,z,1);
          ek=s*ek;
        }
        if (a[k][k] != (FPN) 0.0) 
          z[k]=z[k]/a[k][k];
        else
          z[k]=(FPN) 1.0;
      }
      k = k-ks;
    }
    s = (FPN) 1.0 /sasum(n,z,1);
    info =  sscal(n,s,z,1);
/*
     solve trans(u)*y=w
*/
  k=1;
  while ( k < n) {
    ks = 1;
    if ( kpvt[k] < 0 ) 
       ks = 2;
    if ( k !=1 ) {
      z[k] = z[k] + sdot(k-1,a[k],1,z,1);
      if (ks == 2 )
        z[k+1] = z[k+1] + sdot(k-1,a[k+1],1,z,1);
      kp = abs(kpvt[k]);  
      if (kp != k) {
        t = z[k];
        z[k] = z[kp];
        z[kp] = t;
      }
    }
    k=k+ks;
  }
  s=(FPN) 1.0 /sasum(n,z,1);
  info =  sscal(n,s,z,1);
  ynorm=(FPN) 1.0;
/*
     solve u*d*v=y
*/
  k=n;
  while ( k != 0) {  
    ks=1;
    if ( kpvt[k] < 0 ) 
      ks=2;
    if ( k != ks) {
      kp = abs(kpvt[k]); 
      kps = k+1-ks;
      if ( kp != kps ) {
        t = z[kps];
        z[kps] = z[kp];
        z[kp] = t;
      }
      info = saxpy(k-ks,z[k],a[k],1,z,1);
      if ( ks == 2 ) 
        info = saxpy(k-ks,z[k-1],a[k-1],1,z,1);
    }
    if ( ks == 2) {
       ak = a[k][k]/a[k][k-1];
       akm1 = a[k-1][k-1]/a[k][k-1];
       bk = z[k]/a[k][k-1];
       bkm1 = z[k-1]/a[k][k-1];
       denom = ak*akm1-(FPN) 1.0;
       z[k] = (akm1*bk-bkm1)/denom;
       z[k-1] = (ak*bkm1-bk)/denom;
    }
    else {
      while (fabs(z[k]) > fabs(a[k][k]) ) {
        s=(FPN)fabs(a[k][k])/(FPN)fabs(z[k]);
        info = sscal(n,s,z,1);
        ynorm=s*ynorm;
      }
      if ( a[k][k] != (FPN) 0.0 ) 
        z[k]=z[k]/a[k][k];
      else
        z[k]=(FPN) 1.0;
    }
    k=k-ks;
  }
  s = (FPN) 1.0 / sasum(n,z,1);
  info = sscal(n,s,z,1);
  ynorm = s * ynorm;
/*
     solve trans(u)*z=v
*/
  k = 1;
  while (k <= n) {
    ks = 1;
    if ( kpvt[k] < 0) 
      ks = 2;
      if ( k != 1) {
        z[k] = z[k] + sdot(k-1,a[k],1,z,1);
        if ( ks == 2 )
          z[k+1] = z[k+1] + sdot(k-1,a[k+1],1,z,1);
        kp = abs(kpvt[k]);
        if ( kp != k) {
          t = z[k];
          z[k] = z[kp];
          z[kp] = t;
        }
      }
    k=k+ks;
  }
/*c     make znorm=1.0*/
  s = (FPN) 1.0/sasum(n,z,1);
  info = sscal(n,s,z,1);
  ynorm = s*ynorm;
  if ( anorm != (FPN) 0.0 ) 
    *rcond=ynorm/anorm;
  else
    *rcond = (FPN) 0.0;      
  return(info);
}

/*
  Chapter 5 page 4
  
  ENTRY
   a = dmatrix(1,n,1,lda) is contains lower symmetric matrix whose facotrization is to be 
   computed.  Only the lower triangle is used.
   
  RETURN
   a contains the information necessary to construct a matrix U and block diagonal matrix D 
   so that A' = UDU'
   kpvt contains the pivot information
   info is indicator 0 => all well
*/
int ssifa(FPN **a, int lda, int n, int *kpvt)
{
      FPN ak,akm1,bk,bkm1,denom,mulk,mulkm1,t;
      FPN absakk,alpha,colmax,rowmax;
      int  info,imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,swap,doofus;
      doofus=lda;
/*    initialize */
      alpha=( (FPN) 1.0 + (FPN)sqrt( 17.0 ) )/(FPN)8.0 ;
      info=0;
/*     main loop on k, which goes from n to 1. */
  k = n;
  while ( k != 0 ) {
    if  ( k <= 1) {
      kpvt[1] = 1;
      if ( a[1][1] == (FPN) 0.0 ) 
        info=1;
    }
    else {
      km1 = k-1;
      absakk = (FPN)fabs(a[k][k]);
      imax = isamax(k-1,a[k],1);
      colmax = (FPN)fabs(a[imax][k]);
      if ( absakk >= alpha*colmax)  {
        kstep = 1;
        swap = 0;
      } 
      else {
        rowmax=(FPN) 0.0;
        imaxp1 = imax+1;
        for ( j = imaxp1 ; j <= k ; j++ ) 
          rowmax =  amax1(rowmax,(FPN) fabs(a[imax][j]));  
        if ( imax != 1) {
          jmax = isamax(imax-1,a[imax],1);
          rowmax = amax1(rowmax,(FPN) fabs(a[jmax][imax]));  
        }
        if ( (FPN)fabs(a[imax][imax]) >= alpha*rowmax ) {
          kstep=1;
          swap=1;
        }
        else {
          if ( absakk < alpha*colmax*(colmax/rowmax) ) {
             kstep = 2;
             if ( imax != km1 )
               swap = 1;
             else
               swap = 0;
          }
          else {
            kstep = 1;
            swap = 0;
          }
        }
      }
      if ( amax1(absakk,colmax) == (FPN) 0.0 ) {  
        kpvt[k] = k;
        info = k;
      }
      else {
       if ( kstep != 2 ) {
         if ( swap == 1) {
           info = sswap(imax,a[imax],1,a[k],1);
           for (  jj = imax ; jj <= k  ; jj++ ) {
             j = k + imax - jj;
             t = a[k][j];
             a[k][j] = a[j][imax];
             a[j][imax] = t;
           }
         }
         for ( jj=1 ; jj<= km1 ; jj++ ) {
           j = k-jj;
           mulk = -a[k][j]/a[k][k];
           t = mulk;
           info = saxpy(j,t,a[k],1,a[j],1);
           a[k][j] = mulk;
         } 
         kpvt[k] = k;
         if ( swap == 1 ) 
           kpvt[k] = imax;
        }
        else {
          if (  swap == 1 ) {
            info = sswap(imax,a[imax],1,a[k-1],1);
            for (  jj = imax ; jj <= km1 ; jj++ ) {
              j = km1+imax-jj;
              t = a[k-1][j];
              a[k-1][j] = a[j][imax];
              a[j][imax] = t;
            }
            t = a[k][k-1];
            a[k][k-1] = a[k][imax];
            a[k][imax] = t;
          }
          km2 = k-2;
          if ( km2 != 0 ) {
            ak = a[k][k]/a[k][k-1];
            akm1 = a[k-1][k-1]/a[k][k-1];
            denom = (FPN) 1.0 - ak*akm1;
            for ( jj=1 ; jj <= km2 ; jj++ ) {
              j = km1 - jj;
              bk = a[k][j]/a[k][k-1];
              bkm1 = a[k-1][j]/a[k][k-1];
              mulk = (akm1 * bk - bkm1)/denom;
              mulkm1 = (ak * bkm1 - bk)/denom;
              t = mulk;
              info = saxpy(j,t,a[k],1,a[j],1);
              t = mulkm1;
              info = saxpy(j,t,a[k-1],1,a[j],1);
              a[k][j] = mulk;
              a[k-1][j] = mulkm1;
            }
          }
          kpvt[k] = 1-k;
          if ( swap == 1 ) 
            kpvt[k] = -imax;
          kpvt[k-1] = kpvt[k];
        }
      }
      k = k-kstep;
    }
  } /*200*/
  return(info);
}
 

/*
 * They appeared in Dongarra, J.J, Moler, C.B., Bunch, J.R.
 * and Stewart G.W. 1979 LINPACK Users' Guide, SIAM,
 * Philadelphia.

This explanation comes from page 6 of Chapter 5 of the LINPAK  Users Guide

INPUT

a = dmatrix(1,n,1,lda)
kpvt = ivector(1,n)
work = dvector(1,n)

job indicates what to compute
   = 111  inertia determinant inverse
     110  inertia determinant
     101  inertia             inverse
     100  inertia
      11          determinant inverse
      10          determinant
       1                      inverse

If SSICO set rcond==0.0, or SSIF set info to nonzero, then no inverse should be requested.       

RETURN

a   contains upper triangular inverse of the original matrix in its upper triangle.
det = dvector(1,2) contains the determinant of a in the form ||a|| = det[1] 10^((int) det[2])
     1.0 <= |det[1]| < 10.0
inert = dvector(1,3) contains the inertia of a.  Gives the number of positive, negative and
        zero eigenvalues of a, respectively.
*/
int ssidi(FPN **a, int lda, int n, int *kpvt, FPN *det, int *inert, FPN *work, int job)
{

      FPN akkp1;
      FPN ten,d,t,ak,akp1;
      int  j,jb,k,km1,ks,kstep,info,doofus;
      int noinv,nodet,noert;
      doofus=lda;
      noinv = nodet = noert = 0;
      if ( job == 10 || job == 100 || job == 110 )
        noinv =  1;      
      if ( job == 101 || job == 100 || job == 1 )
        nodet = 1;
      if ( job == 1 || job == 10  || job == 11 )
        noert = 1;

      if ( nodet == 0 ||  noert == 0 ) {
         if ( noert == 0 ) {
            inert[1] = 0;
            inert[2] = 0;
            inert[3] = 0;
         }
         if (nodet == 0)  {
            det[1] = (FPN) 1.0;
            det[2] = (FPN) 0.0;
            ten = (FPN)10.0;
         }
         t = (FPN) 0.0;
         for ( k=1 ; k<=n ; k++ ) {
            d=a[k][k];

            if ( kpvt[k] <= 0) {

            if ( t == (FPN) 0.0 ) {
               t= (FPN)fabs(a[k+1][k]);
               d= (d/t)*a[k+1][k+1]-t;
            }
            else {
               d=t;
               t=(FPN) 0.0;
            }
         }

         if ( noert == 0 ) {
            if ( d > (FPN) 0.0 ) 
              inert[1]=inert[1]+1;
            else if ( d < (FPN) 0.0 ) 
              inert[2]=inert[2]+1;
            else
              inert[3]=inert[3]+1;
         }  
 /*60       continue*/
         if ( nodet == 0 ) {
            det[1]=d*det[1];
            if ( det[1] != (FPN) 0.0) {  
               while ( fabs(det[1]) < (FPN) 1.0 ) {
                  det[1]=ten*det[1];
                  det[2]=det[2]-(FPN) 1.0;
               }
               while ( fabs(det[1]) >= ten) {
                     det[1]=det[1]/ten;
                     det[2]=det[2]+(FPN) 1.0;
               }
            }
         }
      }
/*  130      continue    140   continue  */
      }
  if ( noinv == 0 ) {
    k=1;
    while ( k <= n ) {
        km1=k-1;
        if ( kpvt[k] >= 0) {
          a[k][k]=(FPN) 1.0/a[k][k];
          if ( km1  >= 1) {
            info = scopy(km1,a[k],1,work,1);
            for ( j=1 ; j<= km1 ; j++ ) {
              a[k][j] = sdot(j,a[j],1,work,1);
              info = saxpy(j-1,work[j],a[j],1,a[k],1);
            }
            a[k][k] = a[k][k]+sdot(km1,work,1,a[k],1);
          }
          kstep=1; 
        }
        else {
          t = (FPN)fabs(a[k][k+1]);
          ak = a[k][k]/t;
          akp1 = a[k+1][k+1]/t;
          akkp1 = a[k+1][k]/t;
          d = t*(ak*akp1-(FPN) 1.0);
          a[k][k] = akp1/d;
          a[k+1][k+1] = ak/d;
          a[k][k+1] = -akkp1/d;
          if ( km1 >= 1 ) {
            info = scopy(km1,a[k+1],1,work,1);
            for ( j=1 ; j<= km1 ; j++) {
              a[k+1][j]=sdot(j,a[j],1,work,1);
              info = saxpy(j-1,work[j],a[j],1,a[k+1],1); /* check this */
            }
            a[k+1][k+1]=a[k+1][k+1]+sdot(km1,work,1,a[k+1],1);
            a[k+1][k]=a[k+1][k]+sdot(km1,a[k],1,a[k+1],1);
            info = scopy(km1,a[k],1,work,1);
            for ( j=1 ; j<=km1 ; j++ ) {
              a[k][j]=sdot(j,a[j],1,work,1);
              info = saxpy(j-1,work[j],a[j],1,a[k],1);
            }
            a[k][k] = a[k][k] + sdot(km1,work,1,a[k],1);
          }
          kstep=2;
        }
        ks = abs(kpvt[k]);
        if ( ks != k) {
          info = sswap(ks,a[ks],1,a[k],1);
          for ( jb=ks ; jb<=k ; jb++ ) {
            j=k+ks-jb;
            dcswap(&a[k][j],&a[j][ks]);
          }
          if ( kstep != 1 ) 
            dcswap(&a[k+1][ks],&a[k+1][k]); 
        }
        k=k+kstep;
      }
  }
  return(info);
}

void dcswap(FPN *a, FPN *b)
{
  FPN tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

FPN amax1(FPN x, FPN y)
{
  if ( x > y )
    return(x);
  else
    return(y);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Linpak.c
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (NumRec.c) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */


#include "Main.h"

/*

  All the subroutines in this file are modified versions of those from


  Numerical Recipes in C, The Art of Scientific Computing
  by W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling.
  1990, Cambridge University Press, Cambridge, England
  ISBN 0-521-35465-X  (the book)
  ISBN 0-521-35466-8  (MS-Dos diskette with C sources)


  They have been modified by Chris J. Basten. The modifications should not effect the behavior
  of the programs in any way.  The source has been reformatted to
  the liking of C.J. Basten, and the function definitions have been
  moved to "NumRec.h". All 'FPN' declarations were changed to 'FPN'.
  Finally, the uniform random number generator
  ranf() was substituted for ran1().  ranf() is an implementation of
  of a FORTRAN subroutine written by John Monahan at North Carolina
  State University.
*/


FPN gammln(FPN xx)
{
  FPN x, tmp, ser;
  static FPN cof[6] = {(FPN) 76.18009173,(FPN) -86.50532033,(FPN) 24.01409822, (FPN) -1.231739516,(FPN) 0.120858003e-2,(FPN) -0.536382e-5};
  int j;
  x = xx - (FPN)1.0;
  tmp = x + (FPN)5.5;
  tmp -= (x + (FPN)0.5) * log(tmp);
  ser = (FPN)1.0;
  for (j = 0; j <= 5; j++) {
    x += (FPN)1.0;
    ser += cof[j] / x;
  }
  return -tmp + (FPN) log( 2.50662827465 * ser);
}



FPN gammp(FPN la,FPN  x)
{
  FPN gamser, gammcf, gln;

  if (x < 0.0 || la <= 0.0)
    nrerror("Invalid arguments in routine GAMMP");
  if (x < (la + 1.0)) {
    gser(&gamser, la, x, &gln);
    return gamser;
  }
  else {
    gcf(&gammcf, la, x, &gln);
    return (FPN) 1.0 - gammcf;
  }
}



FPN gasdev(int *idum)
{
  static int iset = 0;
  static FPN gset;
  FPN fac, r, v1, v2;

  if (iset == 0) {
    do {
      v1 = (FPN) 2.0 * ranf(*idum) - (FPN) 1.0;
      v2 = (FPN) 2.0 * ranf(*idum) - (FPN) 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1.0 || r == 0.0);
    fac = (FPN) sqrt(-2.0 * log(r) / r);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  }
  else {
    iset = 0;
    return gset;
  }
}



void gcf(FPN *gammcf,FPN  la,FPN  x,FPN  *gln)
{
  int n;
  FPN gold = (FPN) 0.0, g, fac = (FPN) 1.0, b1 = (FPN) 1.0;
  FPN b0 = (FPN) 0.0, anf, ana, an, a1, a0 = (FPN) 1.0;

  *gln = gammln(la);
  a1 = x;
  for (n = 1; n <= ITMAX; n++) {
    an = (FPN) n;
    ana = an - la;
    a0 = (a1 + a0 * ana) * fac;
    b0 = (b1 + b0 * ana) * fac;
    anf = an * fac;
    a1 = x * a0 + anf * a1;
    b1 = x * b0 + anf * b1;
    if (a1) {
      fac = (FPN) 1.0 / a1;
      g = b1 * fac;
      if (fabs((g - gold) / g) < EPS) {
	*gammcf = (FPN) exp(-x + la * log(x) - (*gln)) * g;
	return;
      }
      gold = g;
    }
  }
  nrerror("a too large, ITMAX too small in routine GCF");
}



void gser(FPN *gamser,FPN la,FPN x,FPN *gln)
{
  int n;
  FPN sum, del, ap;

  *gln = gammln(la);
  if (x <= 0.0) {
    if (x < 0.0)
      nrerror("x less than 0 in routine GSER");
    *gamser = (FPN) 0.0;
    return;
  }
  else {
    ap = la;
    del = sum = (FPN) 1.0 / la;
    for (n = 1; n <= ITMAX; n++) {
      ap += 1.0;
      del *= x / ap;
      sum += del;
      if (fabs(del) < fabs(sum) * EPS) {
	*gamser = sum * (FPN) exp(-x + la * log(x) - (*gln));
	return;
      }
    }
    nrerror("a too large, ITMAX too small in routine GSER");
    return;
  }
}



void indexx(int n,FPN arrin[],int indx[])
{
  int l, j, ir, indxt, i;
  FPN q;

  for (j = 1; j <= n; j++)
    indx[j] = j;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1)
      q = arrin[(indxt = indx[--l])];
    else {
      q = arrin[(indxt = indx[ir])];
      indx[ir] = indx[1];
      if (--ir == 1) {
	indx[1] = indxt;
	return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && arrin[indx[j]] < arrin[indx[j + 1]])
	j++;
      if (q < arrin[indx[j]]) {
	indx[i] = indx[j];
	j += (i = j);
      }
      else
	j = ir + 1;
    }
    indx[i] = indxt;
  }
}

void moment(FPN data[],int n,FPN *ave,FPN *adev,FPN *sdev,FPN *svar,FPN *skew,FPN *curt)
{
  int j;
  FPN s, lp;


  if (n <= 1)
    nrerror("n must be at least 2 in MOMENT");
  s = (FPN) 0.0;
  for (j = 1; j <= n; j++)
    s += data[j];
  *ave = s / n;
  *adev = (*svar) = (*skew) = (*curt) = (FPN) 0.0;
  for (j = 1; j <= n; j++) {
    *adev += fabs(s = data[j] - (*ave));
    *svar += (lp = s * s);
    *skew += (lp *= s);
    *curt += (lp *= s);
  }
  *adev /= n;
  *svar /= (n - 1);
  *sdev = (FPN) sqrt(*svar);
  if (*svar) {
    *skew /= (n * (*svar) * (*sdev));
    *curt = (*curt) / (n * (*svar) * (*svar)) - (FPN) 3.0;
  }
  else
    nrerror("No skew/kurtosis when variance = 0 (in MOMENT)");
}

void sort(int n,FPN ra[])
{
  int l, j, ir, i;
  FPN rra;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1)
      rra = ra[--l];
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
	ra[1] = rra;
	return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j + 1])
	++j;
      if (rra < ra[j]) {
	ra[i] = ra[j];
	j += (i = j);
      }
      else
	j = ir + 1;
    }
    ra[i] = rra;
  }
}

FPN poidev(FPN xm,int *idum)
{
  static FPN sq, alxm, g, oldm; /* = (-1.0);*/
  FPN em, t, y;

  oldm = (FPN) -1.0;
  
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm = (FPN) xm;
      g = (FPN) exp(-   xm);
    }
    em = -1;
    t = (FPN) 1.0;
    do {
      em += (FPN) 1.0;
      t *= ranf(*idum);
    } while (t > g);
  }
  else {
    if (xm != oldm) {
      oldm = (FPN) xm;
      sq = (FPN) sqrt(2.0 *  xm);
      alxm = (FPN) log( xm);
      g =   (xm * alxm - gammln(  xm + (FPN) 1.0) );
    }
    do {
      do {
	y = (FPN) tan(PI * ranf(*idum));
	em = sq * y + xm;
      } while (em < 0.0);
      em = (FPN) floor(em);
      t = (FPN) 0.9 * ((FPN) 1.0 + y * y) * (FPN) exp(em * alxm - gammln(em + (FPN) 1.0) - g);
    } while (ranf(*idum) > t);
  }
  return em;
}

FPN beta(FPN z,FPN  w)
{
  return  (FPN) exp(gammln(z) + gammln(w) - gammln(z + w)) ;
}


FPN betacf(FPN la,FPN  b,FPN  x)
{
  FPN qap, qam, qab, em, tem, d;
  FPN bz, bm = (FPN)1.0, bp, bpp;
  FPN az = (FPN) 1.0, am = (FPN) 1.0, ap, app, aold;
  int m;

  qab = la + b;
  qap = la + (FPN)1.0;
  qam = la - (FPN)1.0;
  bz = (FPN)1.0 - qab * x / qap;
  for (m = 1; m <= ITMAX; m++) {
    em = (FPN) m;
    tem = em + em;
    d = em * (b - em) * x / ((qam + tem) * (la + tem));
    ap = az + d * am;
    bp = bz + d * bm;
    d = -(la + em) * (qab + em) * x / ((qap + tem) * (la + tem));
    app = ap + d * az;
    bpp = bp + d * bz;
    aold = az;
    am = ap / bpp;
    bm = bp / bpp;
    az = app / bpp;
    bz = (FPN) 1.0;
    if (fabs(az - aold) < (EPS * fabs(az)))
      return(az);
  }
  nrerror("a or b too big, or ITMAX too small in BETACF");
  return((FPN) 0.0);
}

FPN betai(FPN la,FPN  b,FPN  x)
{
  FPN bt;

  if (x < (FPN)0.0 || x > (FPN)1.0)
    nrerror("Bad x in routine BETAI"); 
  if (x == (FPN)0.0 || x == (FPN)1.0)
    bt = (FPN)0.0;
  else
    bt = (FPN) exp(gammln(la + b) - gammln(la) - gammln(b) + la * (FPN) log(x) + b * (FPN) log(1.0 - x));
  if (x < (la + (FPN)1.0) / (la + b + (FPN)2.0))
    return bt * betacf(la, b, x) / la;
  else
    return (FPN) 1.0 - bt * betacf(b, la, (FPN) 1.0 - x) / b;
}

/**************************************************************************/
/*       Function  prob = chiprob(dof,value)                              */
/*                                                                        */
/*  Calculates the chi-squared probability for 'dof' degrees of           */
/* freedom for the 'value' of the X-square.              */
/**************************************************************************/

FPN  chiprob( int     dof, FPN  value)
{
  FPN ddof,prob;
  ddof = (FPN) dof;
  if ( value > 0.0 )
    prob =  gammq(ddof,value/(FPN) 2.0);
  else
    prob = (FPN) 1.0;

  return prob;
}

/*
This is a Numerical Recipes in C function  
*/
FPN gammq( FPN la,FPN x)
{
	FPN gamser,gammcf,gln;

	if (x < 0.0 || la <= 0.0) nrerror("Invalid arguments in routine GAMMQ");
	if (x < (la+1.0)) {
		gser(&gamser,la,x,&gln);
		return (FPN) 1.0-gamser;
	} else {
		gcf(&gammcf,la,x,&gln);
		return gammcf;
	}
}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file NumRec.c
------------------------------------------------------------------ */


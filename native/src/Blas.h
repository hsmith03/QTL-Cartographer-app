/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Blas.h) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */

/*  
 *
 * The subroutines below were originally written in
 * Fortran. They appeared in Dongarra, J.J, Moler, C.B.,
 * Bunch, J.R. and Stewart G.W. 1979 LINPACK Users' Guide,
 * SIAM, Philadelphia.
 *
 * The C translations were done by Christopher J. Basten at
 * North Carolina State University in December of 1993.
 * These are simple ports in that they do not unroll loops
 * when the increments are both one.
 *
 */
   /*  You can change this to float if you want floating point numbers to be floats rather than doubles.*/
#ifndef FPN
#include "LocalD.h"
#endif

#include <stdio.h>
#include <math.h>

#if defined(DBLC)
#define DBL_MAX 1.7976931348623157e+308
#define DBL_EPSILON 2.2204460492503131e-16
#define DBL_MIN 2.2250738585072014e-308
#else 
#include <float.h>
#endif



int sswap(int nn, FPN *sx, int isx, FPN *sy, int isy);
int sscal(int nn, FPN sa, FPN *sx, int isx);
FPN snrm2(int nn, FPN *sx, int isx);
FPN sdot(int nn, FPN *sx, int isx, FPN *sy, int isy);
int scopy(int nn, FPN *sx, int isx, FPN *sy, int isy);
int saxpy(int nn, FPN sa, FPN *sx, int isx, FPN *sy, int isy);
FPN sasum(int nn, FPN *sx, int isx);
int isamax(int nn, FPN *sx, int isx);
int srot(int nn, FPN *sx, int isx, FPN *sy, int isy, FPN cc, FPN ss)
;
int srotg(FPN *sa, FPN *sb, FPN *cc, FPN *ss);

/*My functions...(Chris Basten)*/
FPN ssum(int nn, FPN *sx, int isx);
FPN ssaxpy(int nn, FPN sa, FPN *sx, int isx, FPN *sy, int isy);
FPN shadamer(int nn, FPN *sx, int isx, FPN *sy, int isy, FPN *sr, int isr);



/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Blas.h
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Linpak.h) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */

/*
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




int sqrst(FPN **xx, int ldx, int nn, int pp, FPN *yy, FPN tol, FPN *bb, FPN *rsd, int *kk, int *jpvt, FPN *qraux);
int sqrdc(FPN **yy, int ldx, int nn, int pp, FPN *qraux, int *jpvt, int job);
int sqrsl(FPN **yy, int ldx, int nn, int pp, int kk, FPN *qraux, FPN *y, FPN *qy, FPN *qty, FPN *bb, FPN *rsd, FPN *xb, long int job);
int spodi(FPN **yy, int lda, int n, FPN *det, int job);
int strsl(FPN **yy, int ldt, int n, FPN *b, int job);

void dcswap(FPN *a, FPN *b);
int ssidi(FPN **a, int lda, int n, int *kpvt, FPN *det, int *inert, FPN *work, int job);
int ssifa(FPN **a, int lda, int n, int *kpvt);
int ssico(FPN **a, int lda, int n, int *kpvt, FPN *rcond, FPN *z);
int spofa(FPN **a, int lda, int n);

FPN amax1(FPN x, FPN y);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Linpak.h
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (f2c.h) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */

/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef VERSION
#include "LocalD.h"
#endif

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
/* typedef long long longint; */ /* system-dependent */

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long flag;
typedef long ftnlen;
typedef long ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	FPN d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

typedef long Long;	/* No longer used; formerly in Namelist */

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define llabs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (FPN)abs(x)
#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif

#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

#define dmin(a,b) (FPN)min(a,b)
#define dmax(a,b) (FPN)max(a,b)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef FPN (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef FPN (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef FPN E_f;	/* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif

void s_copy(register char *a, register char *b, ftnlen la, ftnlen lb);
double pow_di(FPN *ap, integer *bp);
double d_sign(FPN *a, FPN *b);
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb);
void s_cat(char *lp, char *rpp[], ftnlen rnp[], ftnlen *np, ftnlen ll);
 
/*  converted blas functions*/
/* Subroutine */ int dcopy_(integer *n, FPN *dx, integer *incx, 	FPN *dy, integer *incy);
/* Subroutine */ int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, FPN *alpha, FPN *a, integer *lda, 	FPN *b, integer *ldb, FPN *beta, FPN *c, integer *ldc);
/* Subroutine */ int dgemv_(char *trans, integer *m, integer *n, FPN *alpha, FPN *a, integer *lda, FPN *x, integer *incx, FPN *beta, FPN *y, integer *incy);
/* Subroutine */ int dger_(integer *m, integer *n, FPN *alpha, FPN *x, integer *incx, FPN *y, integer *incy, FPN *a, integer *lda);
/* Subroutine */ FPN dnrm2_(integer *n, FPN *x, integer *incx);
/* Subroutine */ int dscal_(integer *n, FPN *da, FPN *dx, integer *incx);
/* Subroutine */ int dtrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, FPN *alpha, FPN *a, integer *lda, FPN *b, integer *ldb);
/* Subroutine */ int dtrmv_(char *uplo, char *trans, char *diag, integer *n, FPN *a, integer *lda, FPN *x, integer *incx);
/* Subroutine */ int dtrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, FPN *alpha, FPN *a, integer *lda, FPN *b, integer *ldb);
/* Subroutine */ logical lsame_(char *ca, char *cb);
/* Subroutine */ int xerbla_(char *srname, integer *info);

 
 
 /*  converted lapack functions*/
/* Subroutine */ int dgeqrf_(integer *m, integer *n, FPN *a, integer *	lda, FPN *tau, FPN *work, integer *lwork, integer *info);
/* Subroutine */ int dgeqr2_(integer *m, integer *n, FPN *a, integer *lda, FPN *tau, FPN *work, integer *info);
/* Subroutine */ FPN dlamch_(char *cmach);
/* Subroutine */ int dlamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
/* Subroutine */ int dlamc2_(integer *beta, integer *t, logical *rnd, 	FPN *eps, integer *emin, FPN *rmin, integer *emax, 	FPN *rmax);
/* Subroutine */ FPN dlamc3_(FPN *a, FPN *b);
/* Subroutine */ int dlamc4_(integer *emin, FPN *start, integer *base);
/* Subroutine */ int dlamc5_(integer *beta, integer *p, integer *emin, 	logical *ieee, integer *emax, FPN *rmax);
/* Subroutine */ FPN dlapy2_(FPN *x, FPN *y);
/* Subroutine */ int dlarf_(char *side, integer *m, integer *n, FPN *v, integer *incv, FPN *tau, FPN *c, integer *ldc, FPN *work);
/* Subroutine */ int dlarfb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, FPN *v, integer *ldv, FPN *t, integer *ldt, FPN *c, integer *ldc, FPN *work, integer *ldwork);
/* Subroutine */ int dlarfg_(integer *n, FPN *alpha, FPN *x, 	integer *incx, FPN *tau);
/* Subroutine */ int dlarft_(char *direct, char *storev, integer *n, integer *k, FPN *v, integer *ldv, FPN *tau, FPN *t, integer *ldt);
/* Subroutine */ integer ilaenv_(integer *ispec, char *name, char *opts, integer *n1, integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len);
/* Subroutine */ int dtrtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, FPN *a, integer *lda, FPN *b, integer *ldb, integer *info);
/* Subroutine */ int dormqr_(char *side, char *trans, integer *m, integer *n, integer *k, FPN *a, integer *lda, FPN *tau, FPN *c, integer *ldc, FPN *work, integer *lwork, integer *info);
/* Subroutine */ int dorm2r_(char *side, char *trans, integer *m, integer *n, integer *k, FPN *a, integer *lda, FPN *tau, FPN *c, integer *ldc, FPN *work, integer *info);
/* Subroutine */ int dorgqr_(integer *m, integer *n, integer *k, FPN *a, integer *lda, FPN *tau, FPN *work, integer *lwork, integer *info);
/* Subroutine */ int dorg2r_(integer *m, integer *n, integer *k, FPN *a, integer *lda, FPN *tau, FPN *work, integer *info);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file f2c.h
------------------------------------------------------------------ */


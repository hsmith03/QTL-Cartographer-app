/* ------------------------------------------------------ XCutXCodeXSkip
     This file (ranlib.h) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */

#ifndef VERSION
#include "LocalD.h"
#endif
/* Prototypes for all user accessible RANLIB routines */

extern void advnst(long k);
extern FPN genbet(FPN aa,FPN bb);
extern FPN genchi(FPN df);
extern FPN genexp(FPN av);
extern FPN genf(FPN dfn, FPN dfd);
extern FPN gengam(FPN a,FPN r);
extern void genmn(FPN *parm,FPN *x,FPN *work);
extern void genmul(long n,FPN *p,long ncat,long *ix);
extern void genmul2(long n,long ncat,long *ix);
extern FPN gennch(FPN df,FPN xnonc);
extern FPN gennf(FPN dfn, FPN dfd, FPN xnonc);
extern FPN gennor(FPN av,FPN sd);
extern void genprm(long *iarray,int larray);
extern FPN genunf(FPN low,FPN high);
extern void getsd(long *iseed1,long *iseed2);
extern void gscgn(long getset,long *g);
extern long ignbin(long n,FPN pp);
extern long ignnbn(long n,FPN p);
extern long ignlgi(void);
extern long ignpoi(FPN mu);
extern long ignuin(long low,long high);
extern void initgn(long isdtyp);
extern long mltmod(long a,long s,long m);
extern void phrtsd(char* phrase,long* seed1,long* seed2);
extern FPN RANF(void);
extern void setall(long iseed1,long iseed2);
extern void setant(long qvalue);
extern void setgmn(FPN *meanv,FPN *covm,long p,FPN *parm);
extern void setsd(long iseed1,long iseed2);
extern FPN sexpo(void);
extern FPN sgamma(FPN a);
extern FPN snorm(void);
extern long lennob( char *str );
extern void gsrgs(long getset,long *qvalue);
extern void gssst(long getset,long *qset);
extern FPN fsign( FPN num, FPN sign );
extern void inrgcm(void);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file ranlib.h
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (f.c) has been adapted for use with QTL Cartographer
     The copyright is owned by someone else, as indicated below.
------------------------------------------------------ XCutXCodeXUnSkip */

/*=========================s_copy.c============================*/

#include "f2c.h"

/* assign strings:  a = b 

#ifdef KR_headers
VOID s_copy(a, b, la, lb) register char *a, *b; ftnlen la, lb;
#else
void s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
#endif*/
void s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
{
register char *aend, *bend;

aend = a + la;

if(la <= lb)
	while(a < aend)
		*a++ = *b++;

else
	{
	bend = b + lb;
	while(b < bend)
		*a++ = *b++;
	while(a < aend)
		*a++ = ' ';
	}
}



/*=========================s_cat.c============================*/

/*
#ifdef KR_headers
VOID s_cat(lp, rpp, rnp, np, ll) char *lp, *rpp[]; ftnlen rnp[], *np, ll;
#else
VOID s_cat(char *lp, char *rpp[], ftnlen rnp[], ftnlen *np, ftnlen ll)
#endif*/

VOID s_cat(char *lp, char *rpp[], ftnlen rnp[], ftnlen *np, ftnlen ll)
{
ftnlen i, n, nc;
char *f__rp;

n = (int)*np;
for(i = 0 ; i < n ; ++i)
	{
	nc = ll;
	if(rnp[i] < nc)
		nc = rnp[i];
	ll -= nc;
	f__rp = rpp[i];
	while(--nc >= 0)
		*lp++ = *f__rp++;
	}
while(--ll >= 0)
	*lp++ = ' ';
}




/*=========================pow_di.c============================


#ifdef KR_headers
double pow_di(ap, bp) doublereal *ap; integer *bp;
#else
#endif*/
double pow_di(FPN *ap, integer *bp)
{
double pow, x;
integer n;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for( ; ; )
		{
		if(n & 01)
			pow *= x;
		if(n >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}


/*=========================d_sign.c============================


#ifdef KR_headers
double d_sign(a,b) doublereal *a, *b;
#else
double d_sign(doublereal *a, doublereal *b)
#endif*/
double d_sign(FPN *a, FPN *b)
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}


/*=========================s_cmp.c============================


 * compare two strings  

#ifdef KR_headers
integer s_cmp(a0, b0, la, lb) char *a0, *b0; ftnlen la, lb;
#else
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
#endif*/
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
register unsigned char *a, *aend, *b, *bend;
a = (unsigned char *)a0;
b = (unsigned char *)b0;
aend = a + la;
bend = b + lb;

if(la <= lb)
	{
	while(a < aend)
		if(*a != *b)
			return( *a - *b );
		else
			{ ++a; ++b; }

	while(b < bend)
		if(*b != ' ')
			return( ' ' - *b );
		else	++b;
	}

else
	{
	while(b < bend)
		if(*a == *b)
			{ ++a; ++b; }
		else
			return( *a - *b );
	while(a < aend)
		if(*a != ' ')
			return(*a - ' ');
		else	++a;
	}
return(0);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file f.c
------------------------------------------------------------------ */


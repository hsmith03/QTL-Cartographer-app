/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Utilities.h) is part of QTL Cartographer
         
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


/*
  Header file for basic utilities.  
*/

extern int whosemf, debugging;
extern long     ix, Ix;
extern FPN mapparam;
extern FILE *fileopen(char *name, char *mode);  /* open a file with error messages */
extern void fileclose(char *name, FILE *fp);  /* close a file with error messages */

int MoveToToken(FILE *fptr, char *token, int usecase);
void putline(FILE *fptr, char ch, int len);
void MemoryAccount(char *buffer);
void LogTheError(char *filename, char *buffer);
FPN lrtolod(FPN lr);
FPN lodtolr(FPN lod);
int is_number(char *buffer);
int is_pinteger(char *buffer);
int find_token(int n, int len, char *buffer,char **names);
int myfgets(char *buffer, int nn, FILE *fptr);
void shuffle_ivector(int *v, int lb, int ub);
/*  Map functions*/
FPN mapfunc(FPN value, int flag);
FPN iKosambi(FPN mm); 
FPN Kosambi(FPN rr); 
FPN Haldane(FPN rr); 
FPN iHaldane(FPN mm); 
FPN Morgan(FPN rr); 
FPN iMorgan(FPN mm); 
FPN CarterFalconer(FPN rr); 
FPN iCarterFalconer(FPN mm); 
FPN Rao(FPN rr); 
FPN iRao(FPN mm); 
FPN Sturt(FPN rr); 
FPN iSturt(FPN mm); 
FPN Felsenstein(FPN rr); 
FPN iFelsenstein(FPN mm); 
FPN Karlin(FPN rr); 
FPN iKarlin(FPN mm); 

int isfile(char *filename);
int get_int(void); 
int get_next_line(char *buffer, int nn, FILE *fileptr);
int get_next_token(char *xtemp, int nn, FILE *fileptr);
int dtranspose(FPN **mm1, FPN **mm2, int lr, int lc, int ur, int uc);      /*  Transpose a matrix of FPNs 
                           **mm1, **mm2, lr,lc,ur,uc     */
long get_a_seed(void);
char *asctime2(void); 
void print_histogram(FPN *y, int n, int intervals, char *outfile);
long get_identifier(char *thefile);
void gnuplot_values(FPN **xx, FPN **yy, int lr, int ur, int lc, int uc, char *thefile, char *title, char *xaxis, char *yaxis);

#if defined(DSIGN)
FPN dsign(FPN val1, FPN val2);        /* v1 = dsign(v1,v2) transfers sign from v2 to v1 */
#endif

#if defined(DIVT)
typedef struct {
  int     quot;
  int     rem;
}       div_t;

typedef struct {
  long     quot;
  long     rem;
}       ldiv_t;
div_t   div(int numer, int denom);
ldiv_t  ldiv(long numer,long denom);
#endif


#if defined(ITOA) 
char   *strupr(char *s);   /* transforms a string to upper case */
char   *strlwr(char *s);   /* transforms a string to lower case */
#endif


void mypause(void);             /* read characters until a <CR> */
void get_field(int xfield, char *xtemp, char *xbuffer);         /* gets the specified field in a string */


FPN ranf(long int inix);
FPN *ran_arry(int size);

int     *ivector(int nl, int nh);
long    *lvector(int nl, int nh);
FPN   *vector(int nl, int nh);
FPN  *dvector(int nl, int nh);
char    *cvector(int nl, int nh);
int     **imatrix(int nrl, int nrh, int ncl, int nch);
long    **lmatrix(int nrl, int nrh, int ncl, int nch);
char    **cmatrix(int nrl, int nrh, int ncl, int nch);
FPN   **matrix(int nrl, int nrh, int ncl, int nch);
FPN  **dmatrix(int nrl, int nrh, int ncl, int nch);
FPN  **dsvector(int nrl, int nrh );
long    **lsvector(int nrl, int nrh );
char    **csvector(int nrl, int nrh );
FPN   **svector(int nrl, int nrh );
int     **isvector(int nrl, int nrh );
void    free_svector(FPN **m, int nrl, int nrh );
void    free_csvector(char **m, int nrl, int nrh );
void    free_isvector(int **m, int nrl, int nrh );
void    free_lsvector(long **m, int nrl, int nrh );
void    free_dsvector(FPN **m, int nrl, int nrh );
void    free_ivector(int *v, int nl, int nh);
void    free_lvector(long int *v, int nl, int nh);
void    free_vector(FPN *v, int nl, int nh);
void    free_dvector(FPN *v, int nl, int nh);
void    free_cvector(char *v, int nl, int nh);
void    free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void    free_lmatrix(long **m, int nrl, int nrh, int ncl, int nch);
void    free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch);
void    free_dmatrix(FPN **m, int nrl, int nrh, int ncl, int nch);
void    free_matrix(FPN **m, int nrl, int nrh, int ncl, int nch);
void    nrerror(char *error_text);

long iran(long int xix, long int in);
FPN gamnl1(FPN beta, FPN aa, FPN bb, FPN pp, long int xix); 
FPN gamgbh(FPN beta, long int xix);

int Rotator(int thestep);
/*  Two subroutines used only for debugging.*/
void ShowDVector(FPN *invec, int lb, int ub, int cols);
void ShowIVector(int *invec, int lb, int ub, int cols);
void ShowDMatrix(FPN **mm, int lr, int ur, int lc, int uc );

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Utilities.h
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Idatain.h) is part of QTL Cartographer
         
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

void untrans_data(int **dataptr,markermap *themap);
void trans_data(int **dataptr, markermap *themap);
void GetTheData(params *theparams, char *infile,int translate) ;
void get_TransTable(FILE *fptr, char **TransTable);
void get_marki(FILE *fptr, markermap *themap, individual *iptr, int nn, char **TransTable, int iscase);
int get_traiti(FILE *fptr, markermap *themap, individual *iptr, int nn, int traitst, int iscase, char *miss);
int get_otraiti(FILE *fptr, markermap *themap, individual *iptr, int nn, int otraitst, int iscase, char *miss);
void get_imark(FILE *fptr, markermap *themap, individual *iptr, int nn, char **TransTable, int iscase);
int get_itrait(FILE *fptr, markermap *themap, individual *iptr, int nn, int traitst, int iscase, char *miss);
int get_iotrait(FILE *fptr, markermap *themap, individual *iptr, int nn, int otraitst, int iscase, char *miss);
individual *get_std_data(params *theparams,   char *infile);
int  trans_mark(char *xtemp, char **TransTable);
int which_mark(char *xtemp, markermap *themap, int *cc, int *mm, int iscase, int mapmaker);
int which_ind(char *xtemp, individual *iptr, int nn);

int     get_the_nn(int *n, int *t, char *inputf);
int     get_the_datapoints(individual *datapoints, char *inputf, int nn);
int recheck_file_type(char *infile);

individual *get_mm_data(params *theparams, char *infile );
int get_next_mark(FILE *fptr);
int get_mm_markers(FILE *fptr, individual *iptr, int cc, int mm, int nn, char **TransTable, int col, int *trans);
int get_mm_traits(FILE *fptr, FPN **trait_v, int wtrait, int traits, int nn);
int do_trait_equation(FILE *fptr, FPN **trait_v, int wtrait, int traits, int nn);

void determine_markers(individual *individs, int nn);

individual *convert_plabqtl(params *theparams);
int convert_mark(int ch) ;
int MoveMissPhenotypes(int nn,int wt, individual *ind1);
void cp_individ(individual *ind1, individual *ind2);
individual *indvector(int n, markermap *themap, aqtl *qtlptr);
void free_indvector(individual *indptr, int n);

void print_individuals_mcd(params *theparams,individual *indptr,int n, char *outfile);
void print_individuals(params *theparams, individual *indptr, int n,  char *outfile);
void print_individuals_std(params *theparams, individual *indptr, int n,   char *outfile);
void print_individuals_mm(params *theparams, individual *indptr, int n,   char *outfile);
void print_individuals_R(params *theparams,individual *indptr,int n, char *outfile);
void print_individuals_SAS(params *theparams,individual *indptr,int n, char *outfile);
void print_plabqtl(params *theparams,individual *indptr,int n, char *outfile);
void write_procANOVA(char *x, char *y,FILE *outf);
void write_proc(char *x, char *y, FILE *outf);
void print_marker_names(markermap *, FILE *);


/*
markermap *get_plabqtl_map(params *theparams);
int get_plabqtl_names(params *theparams,markermap *themap);
int get_plabqtl_params(params *theparams);
individual *get_plabqtl_matrix(params *theparams,markermap *themap);
*/

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Idatain.h
------------------------------------------------------------------ */


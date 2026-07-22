/* ------------------------------------------------------ XCutXCodeXSkip
     This file (params.h) is part of QTL Cartographer
         
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


void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag);
void update_params(char **opt_v, params *theparams);
void set_workdir(params *theparams);

void write_trailer(params *theparams,char *chptr,int pfrw);
void shift_fn(char *inbuff);
void insert_wd(char *buff,char *wd,char *fn);
void unset_workdir(params *theparams);
void renew_workdir(params *theparams,  char *xtemp,  char **opt, char  **opt_v,  char **opt_e,  int nopts);
void renew_resource(params *theparams,  char *xtemp, char   **opt, char  **opt_v, char  **opt_e,  int nopts);
int check_params(params *theparams, markermap *themap, int prog);


void renew_stem(char *xtemp,params *theparams);
void rewrite_param_file(params *theparams,char *filename);
void get_param_file(params *theparams,char *qtlrc);
params *create_params(params *aparams,int cd,char *qtlrc);

int show_opts(FILE *fptr,char *tptr,char  *prog,char  *purpose,char  **opt,char  **opt_v,char  **opt_e,int nopts,int  oflag, params *theparams);
void create_opts(char **opt, char **opt_v, char **opt_e, int nopts);
void destroy_opts(char **opt, char **opt_v, char **opt_e, int nopts);
int parse_cross(params *theparams,  char *xtemp);
int     get_cross(char *inputf, params *theparams);
int process_arguments(int argc, char **argv, char *thetime, char *purpose,  int nopts,  params *theparams);

int get_file_type(char *filename);
int file_to_int(char *xtemp);
void write_file_type(int ft,  FILE *fptr);
void print_head(char *prog,char *filename,char *chptr,int outmode,int filetype,  params *theparams);
void check_directory(char *chptr);
void quit_banner(char *stri);
void append_seed( params *theparams,char *filename);

extern int writeseed;
extern int whichprogram;

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file params.h
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (QTLcart.c) is part of QTL Cartographer
         
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

#include "Main.h"
/*  
          QTLcart.c,  main function driver for front end (potentially) to the
          QTL Cartographer system.   Right now, it is pretty much useless.          
*/

int gnopts;

int do_int_options();

#if defined(MACWARRIOR)
#include <console.h>
char **environ;
#endif

int main(argc, argv)
int argc;
char *argv[];
{
  FILE *errorf;
  char *chptr, *purpose;
  char **opt, **opt_v, **opt_e;
  params *theparams;
  int nopts;
  int jj, ii, automatic;

#if  defined(MACWARRIOR)
 /* make environ safe, and give the console window a new
    title */
  unsigned char abuffer[25] = "\pQTL Cartographer\0";
  *environ++ = "";
  *environ = (char *) 0;
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif

 /**/
  purpose = cvector(0,MAXNAME);
  strcpy(purpose, "A simple driver for the QTL Cartographer Suite");
  theparams = NULL;
  automatic = 0;
  gnopts = nopts = 16;
  opt = cmatrix(1, nopts, 0, 5);
  opt_e = cmatrix(1, nopts, 0, MAXNAME);
  opt_v = cmatrix(1, nopts, 0, MAXNAME);
  create_opts(opt, opt_v, opt_e, nopts);
  theparams = create_params(theparams, 1, NULL);
  update_opts(opt, opt_v, opt_e, nopts, theparams, 1);
  chptr = asctime2();

  for (ii = 1; ii < argc; ii++)
    if (*(*(argv + ii) + 0) == '-')
      switch (*(*(argv + ii) + 1)) {
       case 'h':
	automatic = show_opts(stdout, chptr, argv[0], purpose, opt, opt_v, opt_e, nopts, 1, theparams);
	exit(0);
       case 's':
	strcpy(theparams->stem, argv[ii + 1]);
	break;
       case 'f':
	strcpy(theparams->resource, argv[ii + 1]);
	break;
       case 'V':
	theparams->verbosity = 0;
	break;
       case 'A':
	automatic = 1;
	break;
       case 'R':
	strcpy(theparams->resource, argv[ii + 1]);
	get_param_file(theparams, theparams->resource);
	break;
       default:
	fprintf(stderr, "\n\tThere is no such argument [%s ]...\n", argv[ii]);
	break;
      }

  jj = 0;
  if (automatic == 0)
    jj = show_opts(stdout, chptr, argv[0], purpose, opt, opt_v, opt_e, nopts, 2, theparams);

  if (jj == 2)
    exit(2);

  errorf = fileopen(theparams->error, "a");
  if (errorf == NULL)
    errorf = fileopen(theparams->error, "w");
  if (errorf != NULL) {
    fprintf(errorf, "\n");
    for (ii = 1; ii <= 80; ii++)
      fprintf(errorf, ".");
    fprintf(errorf, "\n");
    jj = show_opts(errorf, chptr, argv[0], purpose, opt, opt_v, opt_e, nopts, 1, theparams);
    for (ii = 1; ii <= 80; ii++)
      fprintf(errorf, ".");
    fprintf(errorf, "\n");
    fileclose(theparams->error, errorf);
  }
  if (theparams->verbosity == 1)
    jj = show_opts(stdout, chptr, argv[0], purpose, opt, opt_v, opt_e, nopts, 1, theparams);


  do_int_options( opt_v, opt_e, nopts, theparams, argv[0]);

  free_cvector( purpose,0,MAXNAME);


  return(0);
}

void update_params(char **opt_v,params *theparams)
{
  int i,k;
  k=theparams->traits;
  for ( i=1; i<= gnopts; i++ )
    if ( opt_v[i][0] == '\0' )
      k+=1;
}

void update_opts( char **opt,char **opt_v,char **opt_e,int nopts, params *theparams,int flag)
{
  int ii, jj;
  ii=theparams->traits;
  if (flag == 1) {
    strcpy(opt[1], "-m");
    strcpy(opt[2], "-q");
    strcpy(opt[3], "-c");
    strcpy(opt[4], "-Q");
    strcpy(opt[5], "-l");
    strcpy(opt[6], "-z");
    strcpy(opt[7], "-p");
    strcpy(opt[8], "-P");
    strcpy(opt[9], "-e");
    strcpy(opt[10], "-E");
    strcpy(opt[11], "-b");
    strcpy(opt[12], "-S");
    strcpy(opt[13], "-J");
    strcpy(opt[14], "-M");
    strcpy(opt[15], "-G");
    strcpy(opt[16], "-B");


    strcpy(opt_e[1], "Create or Translate a Genetic Linkage Map");
    strcpy(opt_e[2], "Create or Translate a Genetic Model");
    strcpy(opt_e[3], "Create or Translate a Data Set");
    strcpy(opt_e[4], "Do Some Basic Statistics on the Trait");
    strcpy(opt_e[5], "Map using Linear Regression");
    strcpy(opt_e[6], "Map Using (Composite) Interval Mapping");
    strcpy(opt_e[7], "Prune a dataset, or do a Bootstrap");
    strcpy(opt_e[8], "Process Results for Gnuplot");
    strcpy(opt_e[9], "Estimate Postions and Effects of QTLs");
    strcpy(opt_e[10], "Estimate a Linkage Map from Marker Data");
    strcpy(opt_e[11], "Map QTLs for Binary Traits");
    strcpy(opt_e[12], "Map QTLs by Stepwise Regression");
    strcpy(opt_e[13], "Perform Joint Composite Interval Mapping");
    strcpy(opt_e[14], "Perform Multiple-Interval Mapping");
    strcpy(opt_e[15], "Perform Multiple Regression");
    strcpy(opt_e[16], "Perform Bayesian QTL Mapping");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      *(*(opt_v + ii) + jj) = '\0';

  strcpy(opt_v[1], "Rmap");
  strcpy(opt_v[2], "Rqtl");
  strcpy(opt_v[3], "Rcross");
  strcpy(opt_v[4], "Qstats");
  strcpy(opt_v[5], "LRmapqtl");
  strcpy(opt_v[6], "Zmapqtl");
  strcpy(opt_v[7], "Prune");
  strcpy(opt_v[8], "Preplot");
  strcpy(opt_v[9], "Eqtl");
  strcpy(opt_v[10], "Emap");
  strcpy(opt_v[11], "BTmapqtl");
  strcpy(opt_v[12], "SRmapqtl");
  strcpy(opt_v[13], "JZmapqtl");
  strcpy(opt_v[14], "MImapqtl");
  strcpy(opt_v[15], "MultiRegress");
  strcpy(opt_v[16], "Bmapqtl");

}

int do_int_options(  char     **opt_v, char     **opt_e,int nopts, params *theparams, char *program_path)
{
  char  *tpath;
  char *buffer;
  int ii, i, go_on, ans,error;
#if defined(WIN32)
  char *last_sep;
#endif
  buffer = cvector(0,MAXLINE);
  tpath = cvector(0, MAXNAME);
  for (ii = 0; ii <= MAXLINE; ii++)
    buffer[ii] = '\0';
  for (ii = 0; ii <= MAXNAME; ii++)
    tpath[ii] = '\0';
#if defined(WIN32)
  strncpy(tpath, program_path, MAXNAME);
  tpath[MAXNAME] = '\0';
  last_sep = strrchr(tpath, '\\');
  if (last_sep == NULL)
    last_sep = strrchr(tpath, '/');
  if (last_sep != NULL)
    *last_sep = '\0';
  else {
    tpath[0] = '.';
    tpath[1] = '\0';
  }
#endif
  go_on = 1;
  while (go_on == 1) {
    fprintf(stdout, "\nNo.          Function                                    Program:");
    fprintf(stdout, "\n 0. Continue with these parameters");
    for (i = 1; i <= nopts; i++)
      fprintf(stdout, "\n%2d. %-50s  %-15s", i, opt_e[i], opt_v[i]);
    fprintf(stdout, "\n%2d. Turn off Verbosity", i);
    fprintf(stdout, "\n%2d. Set path for executables                       %s", i + 1, tpath);
    fprintf(stdout, "\n%2d. Help", i + 2);
    fprintf(stdout, "\n%2d. Quit", i + 3);
/*    fprintf(stdout, "\n%2d. Send Mail Message", i + 4);*/
    fprintf(stdout, "\n\n\tPlease enter a number to run a program... ");
    ans = get_int();
    if (ans == i)
      theparams->Inter = -1;
    else if (ans == i + 1) {
      fprintf(stdout, "\nChange executable path from [ %s ] to ... ? ", tpath);
      ii = myfgets(tpath, MAXNAME, stdin);
      i = (int) strlen(tpath);
      while (i > 0 && (tpath[i-1] == '\\' || tpath[i-1] == '/'))
        tpath[--i] = '\0';
    }
    else if (ans == i + 2)
      printf("\n Read the manual...\n");
    else if (ans == i + 3) {
      free_cvector(tpath, 0, MAXNAME);
      return(ans);
    }
/*
   else if (ans == i + 4) {
	strcpy(buffer, "/usr/bin/mail basten </dev/null");
#if defined(MACWARRIOR) || defined(WINWARRIOR)
      printf("\n  Sorry, this OS doesn't have a system call.  Can't run \n [ %s ] \n", buffer);
#else
      error = system(buffer);
      if (error == 0)
	printf("\nCompleted task without error.\n");
      else
	printf("\nCouldn't complete task for some reason...errorlevel = %d\n", error);
#endif

    }
*/
    else if (ans > 0 && ans <= nopts) {
      for (ii = 0; ii <= MAXLINE; ii++)
	buffer[ii] = '\0';
#if defined(WIN32)
      if (strlen(tpath) > 0) {
        if (strlen(theparams->resource) > 0)
          _snprintf(buffer, MAXLINE, "\"%s\\%s.exe\" -R \"%s\"%s", tpath,
                    opt_v[ans], theparams->resource,
                    theparams->Inter == -1 ? " -V" : "");
        else
          _snprintf(buffer, MAXLINE, "\"%s\\%s.exe\"%s", tpath, opt_v[ans],
                    theparams->Inter == -1 ? " -V" : "");
      }
      else {
        if (strlen(theparams->resource) > 0)
          _snprintf(buffer, MAXLINE, "%s.exe -R \"%s\"%s", opt_v[ans],
                    theparams->resource, theparams->Inter == -1 ? " -V" : "");
        else
          _snprintf(buffer, MAXLINE, "%s.exe%s", opt_v[ans],
                    theparams->Inter == -1 ? " -V" : "");
      }
      buffer[MAXLINE] = '\0';
#else
      if ((i = (int) strlen(tpath)) > 0) {
	strcpy(buffer, tpath);
	strcat(buffer, opt_v[ans]);
      }
      else
	strcpy(buffer, opt_v[ans]);

      if ((i = (int) strlen(theparams->resource)) > 0) {
	strcat(buffer, " -R ");
	strcat(buffer, theparams->resource);
      }

      if (theparams->Inter == -1)
	strcat(buffer, " -V");
#endif
#if defined(MACWARRIOR) || defined(WINWARRIOR)
      printf("\n  Sorry, this OS doesn't have a system call.  Can't run \n [ %s ] \n", buffer);
      error = 1;
#else
      error = system(buffer);
      if (error == 0)
	printf("\nCompleted task without error.\n");
      else
	printf("\nCouldn't complete task for some reason...errorlevel = %d\n", error);
#endif
      mypause();

    }
    else
      go_on = 0;
  }
  free_cvector(buffer ,0,MAXLINE);

  free_cvector(tpath, 0, MAXNAME);
 return(ans);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file QTLcart.c
------------------------------------------------------------------ */


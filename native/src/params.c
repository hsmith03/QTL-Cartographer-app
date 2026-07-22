/* ------------------------------------------------------ XCutXCodeXSkip
     This file (params.c) is part of QTL Cartographer
         
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
          params.c, general subroutines to manage parameter values
          and the resource file.
*/


#include "Main.h"

int whichprogram=0;
int writeseed=0;


/*
  Write a standard trailer to the log file and the screen.
*/
void write_trailer(theparams,chptr,pfrw)
params *theparams;
char *chptr;
int pfrw;
{
  FILE *errorf;
  
  int i;
/* A note on asctime2.  It writes to a global buffer chptr.  For an accurate gauge
of when the program started and ended, it shouldn't be used between the
first invocation and now. I am saving the presumed start time to gbuffer and
getting the present time with these two calls. */
  strcpy(gbuffer,chptr);  
  chptr = asctime2();
  if (theparams->verbosity == 1) {
    fprintf(stdout,"\n");
    putline(stdout,'.',(int) HLINE);
    fprintf(stdout, "\nProgram started: %s", gbuffer);
    fprintf(stdout, "and finished:    %s", chptr);
  }
  if ( theparams->error != NULL ) {
	 errorf = fileopen(theparams->error, "a");
     fprintf(errorf, "\nProgram started: %s", gbuffer);
	 fprintf(errorf, "and finished:    %s", chptr);
     putline(errorf,'#',(int) HLINE);
	 fileclose(theparams->error, errorf);
  }
  if ( (i=(int) strlen(theparams->workdir)) > 0 )
	 unset_workdir(theparams);
  if ( pfrw == 1 )
    rewrite_param_file(theparams, theparams->resource);
  else {
    if ( theparams->verbosity == 1 ) 
      printf( " \nA Reminder:  This program doesn't rewrite the resource file.");
  }
  if ( theparams->verbosity == 1 ) 
      quit_banner("\n Aside from deallocating memory, this program is finished.\n");
}

/*
  Generalized method to process options to the programs 
    go through each argument.  If it matches an option, update the options value.  If not,
    check if it is a global.  If not, print an error message. 
*/ 
int process_arguments(argc,argv,thetime,purpose,nopts,theparams)
  int argc;
  char **argv;
  char *thetime, *purpose;
  int nopts;
  params *theparams;
{
  FILE *errorf;
  int automatic,burnin;
  int i,ii,j,jj,hit,error; 
  FPN ss,*tmp;
  char **opt, **opt_v, **opt_e;
#if defined(MACWARRIOR)
  SIOUXSettings.asktosaveonclose = FALSE;
  SIOUXSettings.autocloseonquit = FALSE;
  SIOUXSettings.rows = 50;
  SIOUXSettings.columns = 80;
  SIOUXSettings.showstatusline = FALSE;
#endif
  debugging = 0;
  if ( theparams->verbosity > 9 ) {
    automatic = 1;
    theparams->verbosity = theparams->verbosity - 10;
  }
  else
    automatic = 0;
  opt = cmatrix(1, nopts, 0, 5);
  opt_e = cmatrix(1, nopts, 0, MAXNAME);
  opt_v = cmatrix(1, nopts, 0, MAXNAME);
  create_opts(opt, opt_v, opt_e, nopts);
  update_opts(opt, opt_v, opt_e, nopts, theparams, 1);
  burnin = 0; /*  a variable to burn-in the random-number generator. */

  for (ii = 1; ii < argc; ii++)
    if (  argv[ii][0] == '-' ) {
      hit = 0; /*  Determine if you have hit on one of the possible options */
      for ( j = 1 ; j <= nopts ; j++ )
        if ( !strcmp(opt[j],argv[ii]) ) {
          hit = 1;
          for ( jj = 0 ; jj < MAXNAME  ; jj++ )
            opt_v[j][jj] = '\0';
          strcpy(opt_v[j], argv[ii+1]);
        }
      if ( hit == 0 ) {  /* No hit, therefore, must be one of the global options.*/      
	      switch ( argv[ii][1]) {
	       case 'h':
		     automatic = show_opts(stdout,thetime,argv[0],purpose,opt,opt_v,opt_e,nopts,1,theparams);
             if ( theparams->verbosity == 1 ) 
               quit_banner("\n Now exiting program without doing any calculations.");
		     exit(0);
	       case 'V':
		     theparams->verbosity = 0;
		     break;
	       case 'D':
	         if ( argv[ii][2] == '1' ) 
		       debugging = 1;
		     else if ( argv[ii][2] == '2' )
		       debugging = 2;
		     else if ( argv[ii][2] == '3' )
		       debugging = 3;
		     else
		       debugging = 0;
		     break;
	       case 'A':
		     automatic = 1;
		     break;
	       case 'B':
	         burnin = atoi( argv[ii+1] );  
	         break;
	       case 'R':
	         renew_resource(theparams , argv[ii+1],opt,opt_v,opt_e,nopts);
	         break;
	       case 'W':
	         renew_workdir(theparams,argv[ii+1],opt,opt_v,opt_e,nopts);
	         break;
	       case 'X':
             update_params(opt_v, theparams);
	         renew_stem(argv[ii+1],theparams);
	         update_opts(opt,opt_v,opt_e,nopts,theparams,0);
	         break;
	       default:
		     fprintf(stderr, "\n\tThere is no such argument [%s ]...\n", argv[ii]);
		     break;
	      }
    }
  }
  
  if (automatic == 0)
    automatic = show_opts(stdout, thetime, argv[0], purpose, opt, opt_v, opt_e, nopts, 2, theparams);
  if (automatic == 2)
    exit(2); /*if we did help, or chose to quit.*/
  update_params(opt_v,  theparams);
 /* initialize the random number generators */
  ss = ranf(-theparams->seed);
  for ( i=1; i<=burnin; i++) 
    ss = ranf(i);
  tmp = ran_arry(-2);
  free_dvector(tmp, 0, 2);
/*  set up so that files are put in the working subdirectory... */
  if ( (i=(int) strlen(theparams->workdir)) > 0 ) {
    set_workdir(theparams);
    update_opts(opt,opt_v,opt_e,nopts,theparams,0);
  }
  error = check_params(theparams,NULL,2);
  update_opts(opt,opt_v,opt_e,nopts,theparams,0);
/* Print a message to the error file... */

  if ( (ii=isfile(theparams->error)) == 0  )  
     print_head(argv[0],theparams->error,thetime,0,1,theparams);
  if ( theparams->perms != 1  ) {  
	  errorf = fileopen(theparams->error, "a");
	  if (errorf == NULL)
	    errorf = fileopen(theparams->error, "w");
	  if (errorf != NULL) {
	    fprintf(errorf, "\n");
	    automatic = show_opts(errorf, thetime, argv[0], purpose, opt, opt_v, opt_e, nopts, 1, theparams);
	    putline(errorf,'.',(int) HLINE);
	    fprintf(errorf, "\n");
	    fileclose(theparams->error, errorf);
	  }
  }
  /* Show all options   */
  if (theparams->verbosity == 1)
    automatic = show_opts(stdout, thetime, argv[0], purpose, opt, opt_v, opt_e, nopts, 1, theparams);
  destroy_opts(opt, opt_v, opt_e, nopts);
  return(automatic);  
}




/*
  Check all the parameter values to see if they are ok. 
  prog indicates the program.
           1. Rmap
           2. Rqtl
           3. Rcross
           4. Prune
    prog = 5. Qstats
           6. LRmapqtl
           61. SRmapqtl
           7. Zmapqtl
           71. JZmapqtl
           72. MImapqtl
           73. MultiRegress
           8. Eqtl
           9. Preplot  
          10. MImapqtl 
          11. Bmapqtl 
          12. Emap
           
           Any parameters that are out of range are setset to defaults.
*/
int check_params(params *theparams, markermap *themap, int prog) {
  int error,ii;
  error = 0;
  if ( theparams->verbosity == 1 ) {
	    putline(stdout,'+',(int) HLINE);
      printf("\n\tParameter checking..." );
  }

  if  ( theparams->seed < 0 ) {
    theparams->seed = get_a_seed();
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset random number seed to %ld",theparams->seed);
  }
  if ( theparams->mapfunc != whosemf )
    theparams->mapfunc = whosemf;
  if ( theparams->mapfunc < 1 || theparams->mapfunc > 8 ) {
    whosemf = theparams->mapfunc = 1;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset map function to Haldane");
  }
  if ( prog == 1 &&  ( theparams->gout < 1 || theparams->gout > 3 ) ) {
    theparams->gout = 1;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset Output flag to 1" );
  }
  else if ( prog == 3 &&  ( theparams->gout < 0 || theparams->gout > 7 ) ) {
    theparams->gout = 0;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset Output flag to 0" );
  }
  if ( theparams->ihypo < 1  ) {
     theparams->ihypo = 10;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset hypothesis test flag to 10" );
  }

  if ( themap == NULL ) {
	  if ( theparams->chrom < 1 ) {
	    theparams->chrom = (int) DEFm;
        error=error+1;
        if ( theparams->verbosity == 1 )
          printf("\n\t\tReset chromosome number to %d",theparams->chrom );
      }
	  if ( theparams->mark < 1 ) {
	    theparams->mark = (int) DEFl;
        error=error+1;
        if ( theparams->verbosity == 1 )
          printf("\n\t\tReset average marker number to %d",theparams->mark );
      }
	  if ( theparams->vmark < 0.0 ) {
	    theparams->vmark = (FPN) DEFsigl;
        error=error+1;
        if ( theparams->verbosity == 1 )
          printf("\n\t\tReset marker variance to %f",theparams->vmark );
      }
	  if ( theparams->dist <= 0.0 ) {
	    theparams->dist = (FPN) DEFs;
        error=error+1;
        if ( theparams->verbosity == 1 )
          printf("\n\t\tReset intermarker distance to %f",theparams->dist );
      }
	  if ( theparams->vdist < 0.0 ) {
	    theparams->vdist = (FPN) DEFsigs;
        error=error+1;
        if ( theparams->verbosity == 1 )
          printf("\n\t\tReset variance of intermarker distance to %f",theparams->vdist );
      }
	  if ( theparams->tail < 0.0 ) {
	    theparams->tail = (FPN) DEFbrdrs;
        error=error+1;
        if ( theparams->verbosity == 1 )
          printf("\n\t\tReset flanking DNA to %f",theparams->tail );
      }
  }
  else {
    if ( theparams->verbosity == 1 )
      printf("\n\tChecking map parameters found in %s",theparams->map );
    theparams->chrom = themap->m;
    theparams->mark = themap->l;
    theparams->vmark = themap->sigl;
    theparams->dist = themap->s;
    theparams->vdist = themap->sigs;
    theparams->tail = themap->brdrs;
  }  
  
    
  if ( theparams->wchrom < 0 || theparams->wchrom > theparams->chrom ) {
    theparams->wchrom = 0;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset which chromosome to analyze to 0 (meaning all of them)");
  }

  if ( theparams->qnum < 0 )  {
    theparams->qnum = (int) DEFqtl;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset average number of QTL per trait to %d",theparams->qnum);
  }
  if ( theparams->dom < 1 || theparams->dom > 4 ) {
    theparams->dom = 1;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset dominance to 1 (no dominance at QTL)");
  }
  if ( theparams->beta <= 0.0 ) {
    theparams->beta = (FPN) DEFbeta;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset additive parameter beta to %f",theparams->beta);
  }
  if ( theparams->beta1 <= 0.0 ) {
    theparams->beta1 = (FPN) DEFbeta;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset dominance parameter beta1 to %f",theparams->beta1);
  }
  if ( theparams->beta2 <= 0.0 ) {
    theparams->beta2 = (FPN) DEFbeta;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset dominance parameter beta2 to %f",theparams->beta2);
  }
  if ( theparams->traits < 1 ) {
    theparams->traits = 1;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset the number of traits to %d",theparams->traits);
  }
    
/*    
  if ( theparams->whichtrait < 1 || theparams->whichtrait > theparams->traits ) {
    theparams->whichtrait = 1;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset which trait to analyze to 1");
  }
*/
  if ( theparams->reps < 0 ) {
    theparams->reps = (int) REPS;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset repetitions to %d",theparams->reps);
  }
  if ( theparams->boots < 0 ) {
    theparams->boots = (int) REPS;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset bootstraps to %d",theparams->boots);
  }
  if ( theparams->perms < 0 ) {
    theparams->perms = (int) REPS;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset permutations to %d",theparams->perms);
  }
  if (theparams->Herit < 0.0 || theparams->Herit > 1.0 ) {
    theparams->Herit = (FPN) DEFheritability;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset heritability to %f",theparams->Herit);
  }
/*
  theparams->Environ = -1.0;
*/

  if ( theparams->cross < 1 || theparams->cross > 6 ) {
    theparams->cross = 1;
    theparams->crosst = 0;
    theparams->tcross = 0;
    theparams->tcrosst = 0;  
    for (ii = 0; ii <= MAXNAME; ii++) 
      *(theparams->thecross + ii) = '\0';
    strcpy(theparams->thecross, "B1");
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset the type of cross to B1 (Backcross to P1)");
  }
  if ( theparams->nn <= 0 ) {
    theparams->nn = (int) NN;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset sample size for simulations to %d",theparams->nn);
  }
  if ( theparams->Model < 1 || theparams->Model > 10 ) {
    theparams->Model = 3;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset analysis model to 3 (Interval Mapping)" );
  }
  if ( theparams->walk <= 0.0 ) {
    theparams->walk = (FPN) DELTAX;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset chromosome walking jump distance to %f centimorgans",theparams->walk);
  }


  if ( theparams->boot < 0 ) {
    theparams->boot = 0;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset number of bootstraps to 0");
  }
  
  if ( theparams->siglevel < 0.0 ) {
    theparams->siglevel = (FPN) SIG_LEVEL;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset significance level to %f",theparams->siglevel);
  }
  if ( theparams->nbp < 0 ) {
    theparams->nbp = (int) NUM_SIG;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset number of background markers to %d",theparams->nbp);
  }
  if ( theparams->window < 0.0 ) {
    theparams->window = (FPN) WIN_SIZE;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset window size to %f",theparams->window);
  }
  if ( theparams->maxqtl > (int) MAXQTL ) { 
    theparams->maxqtl = (int) MAXQTL ;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset maximal number of QTL to fit in MImapqtl to %d",theparams->maxqtl);
  }
  if ( theparams->maxepistatics > (int) MAXEPISTATICS ) { 
    theparams->maxepistatics = (int) MAXEPISTATICS;
    error=error+1;
    if ( theparams->verbosity == 1 )
      printf("\n\t\tReset maximal number of interaction terms to fit in MImapqtl to %d",theparams->maxepistatics);
  }
  
  if ( theparams->lodflag < 0 || theparams->lodflag > 1 )
    theparams->lodflag = 0;

  if ( theparams->verbosity == 1 ) {
      printf("\n\tFinished checking parameters.  Found %d bad parameter values",error);
	    putline(stdout,'+',(int) HLINE);
      printf("\n");
  }
  return(error);
}




/*
  make sure that any directory name has a trailing file separator character.
*/
void check_directory(chptr)
  char *chptr;
{
 int ii;
 ii = (int) strlen(chptr);
 if ( ii > 0 )
   if ( chptr[ii-1] != (char) FILESEP )
     chptr[ii] = (char) FILESEP;
}


/*  
Renew the working directory
*/
void renew_workdir(theparams,xtemp,opt,opt_v,opt_e,nopts)
  params *theparams;
  char *xtemp, **opt, **opt_v, **opt_e;
  int nopts;
{
  int jj;
  update_params(opt_v,theparams);
/*  if ( theparams->workdir == NULL )
    theparams->workdir = cvector(0,MAXNAME);
  else */
    for ( jj = 0; jj <= MAXNAME ; jj++ )
      theparams->workdir[jj] = '\0';
  strcpy( theparams->workdir, xtemp );
  check_directory(theparams->workdir);
  update_opts(opt,opt_v,opt_e,nopts,theparams,0);
}


/*  
Renew the resource file
*/
void renew_resource(theparams,xtemp,opt,opt_v,opt_e,nopts)
  params *theparams;
  char *xtemp, **opt, **opt_v, **opt_e;
  int nopts;
{
  int jj;
  update_params(opt_v,theparams);
  for ( jj = 0; jj <= MAXNAME ; jj++ )
      theparams->resource[jj] = '\0';
  strcpy( theparams->resource, xtemp );
  get_param_file(theparams,theparams->resource);
  update_opts(opt,opt_v,opt_e,nopts,theparams,0);
}

/*
Get the type of cross from the datafile.
*/
int get_cross(inputf,theparams)
char *inputf;
params *theparams;
{
  FILE *infile;
  int ch,ii, jj, nyes, cross;
  nyes = 0;
  infile = fileopen(inputf, "r");
  if (infile == NULL) 
    return (-1);
  jj = 0;
  do {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if (gbuffer[0] == '-')
      switch (gbuffer[1]) {
       case 'c':
	     get_field(2, gname, gbuffer);
	     cross = parse_cross(theparams,gname);

	     nyes = 1;
	     break;
       case 'q':
	     fileclose(inputf, infile);
	     printf("\nYou did not specify the type of cross...\n");
	     return (-1);
	     break;
       default:
	     break;
      }
    if (nyes == 1) {
      fileclose(inputf, infile);
      return (cross);
    }
  } while (ch != EOF);
  if (ch == EOF)
    printf("\nTime to close the input file after an EOF...\n");
  fileclose(inputf, infile);
  return (-1);
}


/*

  Determine the type of cross.  Fill in the parameter structure to define it.
  
  theparams->cross  is the initial cross:
  
    1  F1 x P1             B1
    2  F1 x P2             B2
    3  F1 x F1 via selfing SFx   (xth generation, x = 2, 3, ...)
    4  F1 x F1 via random mating  RFx   (xth generation, x = 2, 3, ...)
    5  Ri via selfing or sib mating
    6  New design  ----> does not exist  (SC)

  theparams->crosst is the generation
    2 -> x for SF or RF
    0, 1, 2 for FPNd haploid (0) or Ri vi selfed (1)  or sib-mated (2)
    
  
  If an initial cross is then testcrossed, we will need: 
  
  theparams->tcross indicates the test cross type, either 1, 2 or 3 for B1, B2 or SF
  theparams->tcrosst  is the generation
  
  theparams->crosstype is used by MImapqtl to determine the design matrix.  It is the
    final cross and takes on values 1, 2, 3, 4, 5 (same as initial cross values).
    It only codes for what genotypes are present in the data set.  Theparams->ngt indicates
    the number of marker genotypic classes:  2 or 3.   
*/
int parse_cross(theparams,xtemp)
  params *theparams;
  char *xtemp;
{
  int cross,i,j,k;
  char temp[MAXNAME];
  k=0;
  theparams->crosstype = 1;
  theparams->ngt = 2;
  cross = 0;
  for ( i = 0 ; i<MAXNAME; i++ )
    temp[i] = theparams->thecross[i] = '\0';
  strcpy(theparams->thecross, xtemp);

	  switch (xtemp[0]) {
        case '1': case '2': case '3':  /* Originally, 1, 2 , and 3 for B1, B2 and F2 */
	      theparams->crosstype = cross = atoi(xtemp);
	      if (xtemp[0] == '3' ) {
	        theparams->crosst = 2;
	        theparams->ngt = 3;
	      }
	      else 
	        theparams->crosst = 1;
	      theparams->tcross = theparams->tcrosst = 0;
          break;
	    case 'B' :  /*  Backcrosses:   B1 or B2*/
	    case 'b' :
          if ( xtemp[1] == '1' )
            cross = 1;
          else
            theparams->crosstype = cross = 2;
          if (  (i=isdigit(xtemp[2])) == 0 )
            theparams->crosst = 1;
          else /*  Backcrosses:   B1x or B2x  for repeated backcrossing. */
            theparams->crosst = atoi( (xtemp+2) );
	      theparams->tcross = theparams->tcrosst = 0;
	 	    break;
	    case 'S' : /*  Fx line created via selfing the Fx-1 st line. */
	    case 's' :
	      if ( xtemp[1] == 'F' ||  xtemp[1] == 'f' ) {
            theparams->crosstype = cross = 3;  theparams->ngt = 3;
            for ( i = 2 ; i < MAXNAME ; i++ )
              xtemp[i-2] = xtemp[i];
            theparams->crosst = atoi(xtemp);
          }
          else if ( xtemp[1] == 'C' ||  xtemp[1] == 'c' ) {
            theparams->crosstype = cross = 6;  /*  Unless it is Sue Carson's special line. */
            for ( i = 2 ; i < MAXNAME ; i++ )
              xtemp[i-2] = xtemp[i];
            theparams->tcross = theparams->crosst = atoi(xtemp);
          }
	      break;
	    case 'T' :   /* Test crosses.   */
	    case 't' :
	      for ( i=0;xtemp[i]!=')';i++) k+=1; /*find the final paren.*/
	      for ( j=i+1;j<MAXNAME;j++) temp[j-i-1]=xtemp[j];	      
	      if ( temp[0] == 'S' || temp[0] == 'S' ) {  /* Base cross is SFx*/
	        cross = theparams->cross = 3;
            for ( i = 2 ; i < MAXNAME ; i++ )
              temp[i-2] = temp[i];
            theparams->crosst = atoi(temp);	
            theparams->tcrosst = 1;      
            if ( xtemp[3] == '1' )     /* Backcross SFx to P1*/
              theparams->tcross = 1;
            else if ( xtemp[3] == '2' )/* Backcross SFx to P2*/
              theparams->crosstype = theparams->tcross = 2;
            else if ( xtemp[3] == '3' ) {  /* Design III?  */
              theparams->tcross = 12;
              theparams->tcrosst = theparams->crosst+1;
            }
            else {    /*  SFy+x  from SFx   */
              theparams->crosstype = theparams->tcross = 3;
              theparams->ngt = 3;
              for (i=0;i<MAXNAME; i++) temp[i] = '\0';
              for (i=4;xtemp[i]!=')';i++) temp[i-4]=xtemp[i];
              theparams->tcrosst=atoi(temp);
            }
	      }
	      else {  /* Base cross is RFx*/
	        cross = theparams->cross = 4;
	        if ( xtemp[3] == '1' )    /* There are now two options:  Backcross to P1 or to P2*/
	          theparams->tcross = 1;
	        else
	          theparams->crosstype = theparams->tcross = 2;	      
            for ( i = 2 ; i < MAXNAME ; i++ )
              temp[i-2] = temp[i];
            theparams->crosst = atoi(xtemp);
	      }	      
 	    break;
	    case 'R' :
	    case 'r' :
          if ( xtemp[1] == 'F' || xtemp[1] == 'f' ) {  /*  Fx line created via randomly mating the Fx-1 st line. */
            theparams->crosstype = cross = 4;
            theparams->ngt = 3;    
          }
          else  /*  Recombinant inbred line via selfing (1), sib-mating (2) or a FPNd haploid (0) */
            theparams->crosstype = cross = 5;
          for ( i = 2 ; i < MAXNAME ; i++ )
            xtemp[i-2] = xtemp[i];
          theparams->crosst = atoi(xtemp);
	    break;
	    default :  /*Default is backcross to P1. */
	      cross = 1;
	      theparams->crosst = 1;
	      theparams->tcross = theparams->tcrosst = 0;
	      strcpy(theparams->thecross, "B1\0");
	      break;
	  
	  }
	theparams->cross = cross;
	return(cross);
}


/*
  take the working directory off of the filenames.
*/
void unset_workdir(theparams)
params *theparams;
{
    
    shift_fn(theparams->resource);
    shift_fn(theparams->stem);
    shift_fn(theparams->error);
    shift_fn(theparams->map);
    if ( theparams->mapin[0] != '\0' )
      shift_fn(theparams->mapin);
    shift_fn(theparams->qtl);
    if ( theparams->qtlin[0] != '\0' )
      shift_fn(theparams->qtlin);
    shift_fn(theparams->eqtl);
    shift_fn(theparams->ifile);
    if ( theparams->iinfile[0] != '\0' )
      shift_fn(theparams->iinfile);
    shift_fn(theparams->qstat);
    shift_fn(theparams->lrfile);
    shift_fn(theparams->srfile);
    shift_fn(theparams->zfile);
    shift_fn(theparams->mimfile);
    shift_fn(theparams->mqtfile);
    shift_fn(theparams->bayesfile);
    shift_fn(theparams->mrinput);
    shift_fn(theparams->mroutput);
}



void shift_fn(inbuff)
  char *inbuff;
{
  int ii,jj,kk;
  ii = (int) strlen(inbuff);
  for ( jj = ii ; jj >0 && inbuff[jj] != (char) FILESEP ; jj-- ) ;
    if ( jj>0 ) {
      for ( kk = jj+1 ; kk<ii ; kk++ )
	    inbuff[kk-jj-1] = inbuff[kk];
      for ( kk = ii-jj-1 ; kk < ii ; kk++ )
	    inbuff[kk] = '\0';
  }

}

void insert_wd(buff,wd,fn)
  char *buff,*wd,*fn;
{
  unsigned long slen;
  int scmp;
  
  slen =    strlen(wd);
  scmp = (int)  strncmp(wd,fn,  slen); 
/*  if the first slen characters of fn are the same as wd, then
    wd is already in the string and does not need to be put in again. */  
  if (  scmp !=  0 ) {
    strcpy(buff,wd);
    strcat(buff,fn);
    strcpy(fn, buff);
  }
}

void set_workdir(theparams)
params *theparams;
{
  int ii;


  if ( (ii=(int) strlen(theparams->workdir)) > 0   ) {
 
    check_directory(theparams->workdir);

    insert_wd(gname,theparams->workdir, theparams->resource);
    insert_wd(gname,theparams->workdir, theparams->stem);
    insert_wd(gname,theparams->workdir, theparams->error);
    insert_wd(gname,theparams->workdir, theparams->map);
    insert_wd(gname,theparams->workdir, theparams->qtl);
    insert_wd(gname,theparams->workdir, theparams->eqtl);
    insert_wd(gname,theparams->workdir, theparams->ifile);
    insert_wd(gname,theparams->workdir, theparams->qstat);
    insert_wd(gname,theparams->workdir, theparams->zfile);
    insert_wd(gname,theparams->workdir, theparams->mimfile);
    insert_wd(gname,theparams->workdir, theparams->lrfile);
    insert_wd(gname,theparams->workdir, theparams->srfile);
    insert_wd(gname,theparams->workdir, theparams->mqtfile);
    insert_wd(gname,theparams->workdir, theparams->bayesfile);
    insert_wd(gname,theparams->workdir, theparams->mrinput);
    insert_wd(gname,theparams->workdir, theparams->mroutput);

    if ( theparams->mapin[0] != '\0' ) 
      insert_wd(gname,theparams->workdir, theparams->mapin);
    if ( theparams->qtlin[0] != '\0' ) 
      insert_wd(gname,theparams->workdir, theparams->qtlin);
    if ( theparams->iinfile[0] != '\0' ) 
      insert_wd(gname,theparams->workdir, theparams->iinfile);



  }
}


void create_opts(opt, opt_v, opt_e, nopts)
char **opt, **opt_v, **opt_e;
int nopts;
{
  int ii, jj;
  for (ii = 1; ii <= nopts; ii++) {
    for (jj = 0; jj <= 5; jj++)
      opt[ii][jj] = '\0';
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_e[ii][jj] = opt_v[ii][jj] = '\0';
  }
}

void destroy_opts(opt, opt_v, opt_e, nopts)
char **opt, **opt_v, **opt_e;
int nopts;
{
  free_cmatrix(opt, 1, nopts, 0, 5);
  free_cmatrix(opt_e, 1, nopts, 0, MAXNAME);
  free_cmatrix(opt_v, 1, nopts, 0, MAXNAME);
}

/*

whichprogram is set in the main function for each program.  
at this time, only 1-3  and 8 are needed.
           1. Rmap
           12. Emap
           2. Rqtl
           3. Rcross
           4. Prune
    prog = 5. Qstats
           6. LRmapqtl
           61. SRmapqtl
           7. Zmapqtl
           71. JZmapqtl
           72. MImapqtl
           73. MultiRegress
           74. HKmapqtl  (or is this in MultiRegress?)
           8. Eqtl
           9. Preplot  
           
           returns 0 (all done, it is time to run) or 2 (done and should quit)
*/
int show_opts(fptr, tptr, prog, purpose, opt, opt_v, opt_e, nopts, oflag, theparams)
FILE *fptr;
char *tptr, *prog, *purpose, **opt, **opt_v, **opt_e;
int nopts, oflag;
params *theparams;
{
  int ii,jj, i, go_on, ans, last,simmers,infilelength;
  last = nopts;
  go_on = 1;
  strcpy(gname,VERSION);
  putline(fptr,'=',(int) HLINE);
  fprintf(fptr, "\n\t%s",gname);
#if defined( MACWARRIOR )
  fprintf(fptr, "\tfor Macintosh");
#elif defined(WINWARRIOR )
  fprintf(fptr, "  for Win32");
#else
  fprintf(fptr, "  for Unix");
#endif
  fprintf(fptr, "\n\tCopyright (C) 1996-2001 C. J. Basten, B. S. Weir and Z.-B. Zeng.");
  fprintf(fptr, "\n\tQTL Cartographer comes with ABSOLUTELY NO WARRANTY.");
  fprintf(fptr, "\n\tThis is free software, and you are welcome to redistribute it");
  fprintf(fptr, "\n\tunder certain conditions. For details see the file COPYING.\n");
  if (oflag == 1) {
    putline(fptr,'=',(int) HLINE);

    fprintf(fptr, "\nTIME:    %s", tptr);
    fprintf(fptr, "PROGRAM: %s", prog);
    fprintf(fptr, "\nPURPOSE: %s", purpose);
    fprintf(fptr, "\nUSAGE:   %s ", prog);
    for (i = 1; i <= 3; i++)
      fprintf(fptr, " [%-3s x] ", opt[i]);

    fprintf(fptr, "...\nDEFAULTS:");
    for (i = 1; i <= nopts; i++)
      fprintf(fptr, "\n  [%3s  %15s ] %-50s", opt[i], opt_v[i], opt_e[i]);
    fprintf(fptr, "\n\n  Also: [-h] for help, [-A] for automatic,  [-V] for non-Verbose");
    fprintf(fptr, "\n  [-W path] for a working directory, [-R file] to specify a resource");
    fprintf(fptr, "\n  file and [-X stem] to specify a filename stem.\n");	 
  }
  else if (oflag == 2) {
    simmers = 0;
    if ( whichprogram < 4 && whichprogram > 0 )
      simmers = whichprogram;
      
    while (go_on == 1) {
      infilelength = (int) strlen(opt_v[1]);
	  if ( infilelength > 0 && whichprogram == 1 )  
	      last = 7;
	  else if ( infilelength > 0 && whichprogram == 2 )  
	      last = 5;
	  else if ( infilelength > 0 && whichprogram == 3 )  
	      last = 7;
    
    
      putline(fptr,'=',(int) HLINE);
      fprintf(fptr, "\nNo.          Options                                    Values:");
      putline(fptr,'-',(int) HLINE);
      fprintf(fptr, "\n 0. Continue with these parameters");
      for (i = 1; i <= last; i++)
	    fprintf(fptr, "\n%2d. %-50s  %-15s", i, opt_e[i], opt_v[i]);
      putline(fptr,'-',(int) HLINE);
      fprintf(fptr, "\n%2d. Specify Resource File                               %-15s", i, theparams->resource);
      fprintf(fptr, "\n%2d. Change Filename stem                                %-15s", i + 1, theparams->stem);
      fprintf(fptr, "\n%2d. Change Working Directory: %s", i + 2, theparams->workdir);
      fprintf(fptr, "\n%2d. Quit", i + 3);
      fprintf(fptr, "\n%2d. Quit, but update the Resource File", i + 4);
      putline(fptr,'=',(int) HLINE);
      fprintf(fptr, "\n\n\tPlease enter a number... ");
      ans = get_int();
      if (ans == i) {
	    update_params(opt_v,  theparams);
	    fprintf(fptr, "\nChange resource file (%s)  to ... ? ", theparams->resource);
        ii = myfgets(gname, MAXNAME, stdin);
	    if (ii > 0) 
	      renew_resource(theparams,gname,opt,opt_v,opt_e,nopts);
      }
      else if (ans == i + 1) {
	    update_params(opt_v,  theparams);
	    fprintf(fptr, "\nChange filename stem (%s) to ... ? ", theparams->stem);
          ii = myfgets(gname, MAXNAME, stdin);
	    if (ii > 0) {
	      renew_stem(gname, theparams);
	      update_opts(opt, opt_v, opt_e, nopts, theparams, 0);
	    }
      }
      else if (ans == i + 2) {
	    fprintf(fptr, "\nChange working directory (%s) to ... ? ", theparams->workdir);
        ii = myfgets(gname, MAXNAME, stdin);
        if ( ii > 0 ) 
           renew_workdir(theparams,gname,opt,opt_v,opt_e,nopts);
      }
      else if (ans == i + 3) {
        if ( theparams->verbosity == 1 )  
          quit_banner("\n Now exiting program without doing any calculations.");
	    go_on = 2;
	  }
      else if (ans == i + 4) {
        if ( theparams->verbosity == 1 ) 
          quit_banner("\n Will exit after updating resource file.");
	    update_params(opt_v,  theparams);
        rewrite_param_file(theparams, theparams->resource);
	    go_on = 2;
	  }
      else if (ans == 1 && simmers > 0 ) {
	    fprintf(fptr, "\nChange %s from [ %s ] to ... ? ", opt_e[ans], opt_v[ans]);
        ii = myfgets(gname, MAXNAME, stdin);
        if ( ii == 1 && gname[0] == '.' ) {
          last = nopts;
		  for (  jj = 0 ; jj <= MAXNAME ; jj++ )
            opt_v[ans][jj] = '\0';      
        }  
	    else if (ii > 0) {
	      strcpy(opt_v[ans], gname);
	    }
      }
      else if (ans > 0 && ans <= nopts) {
	    fprintf(fptr, "\nChange %s from [ %s ] to ... ? ", opt_e[ans], opt_v[ans]);
        ii = myfgets(gname, MAXNAME, stdin);
        if ( ii == 1 && gname[0] == '.' )
		  for (  jj = 0 ; jj <= MAXNAME ; jj++ )
            opt_v[ans][jj] = '\0';        
	    else if (ii > 0) {
	      strcpy(opt_v[ans], gname);
	      update_params(opt_v,  theparams);
	      update_opts(opt, opt_v, opt_e, nopts, theparams, 0);
	    }
      }
      else
	    go_on = 0;
    }

  }
 

  if ( theparams->resource[0] == '\0' )
    strcpy( theparams->resource, "qtlcart.rc" );
  return (go_on);
}

void quit_banner(char *stri) 
{
   putline(stdout,'=',(int) HLINE);
  printf("%s",stri);
#if defined(WINWARRIOR)
          printf("\n  You may need to click on the upper right X to close this screen");
#endif
#if defined(MACWARRIOR)
          printf("\n  Command+Q will close this screen");  
#endif
  putline(stdout,'=',(int) HLINE);
  printf("\n");
}


/*
  A new stem has been given.  Change the names of the
  appropriate files to stem. + extension.
*/
void renew_stem(xtemp, theparams)
char *xtemp;
params *theparams;
{
  int jj, ch;
  FILE *fptr;
  strcpy(theparams->stem, xtemp);
/*  Clean out other filenames and reassign... */
  for (jj = 0; jj <= MAXNAME; jj++) {
    theparams->map[jj] = '\0';
    theparams->qtl[jj] = '\0';
    theparams->ifile[jj] = '\0';
    theparams->lrfile[jj] = '\0';
    theparams->srfile[jj] = '\0';
    theparams->zfile[jj] = '\0';
    theparams->error[jj] = '\0';
    theparams->qstat[jj] = '\0';
    theparams->mimfile[jj] = '\0';
    theparams->mqtfile[jj] = '\0';
    theparams->bayesfile[jj] = '\0';
    theparams->mrinput[jj] = '\0';
    theparams->mroutput[jj] = '\0';
  }
  sprintf(theparams->error, "%s.log",xtemp);
  sprintf(theparams->map, "%s.map",xtemp);
  sprintf(theparams->qtl, "%s.qtl",xtemp);
  sprintf(theparams->eqtl, "%s.eqt",xtemp);
  sprintf(theparams->ifile, "%s.cro",xtemp);
  sprintf(theparams->lrfile, "%s.lr",xtemp);
  sprintf(theparams->srfile, "%s.sr",xtemp);
  sprintf(theparams->qstat, "%s.qst",xtemp);
  sprintf(theparams->zfile, "%s.z",xtemp);
  sprintf(theparams->bayesfile, "%s.bys",xtemp);
  sprintf(theparams->mrinput, "%s.zr",xtemp);
  sprintf(theparams->mroutput, "%s.mr",xtemp);
  if ( theparams->mimphase == 0 ) {
    sprintf(theparams->mimfile, "%s.mim",xtemp);
    sprintf(theparams->mqtfile, "%si.mqt",xtemp);
    sprintf(theparams->tfile, "%so.mqt",xtemp);
  }
  else {
    sprintf(theparams->mimfile, "%sPhase%d.mim",xtemp,theparams->mimphase-1);
    sprintf(theparams->mqtfile, "%sPhase%d.mqt",xtemp,theparams->mimphase-1);
    sprintf(theparams->tfile, "%sPhase%d.mqt",xtemp,theparams->mimphase);
  }



  if (theparams->map[0] != '\0') {

    if (  (jj = (int) strlen(theparams->workdir)) > 0 ) {
      if ( theparams->workdir[jj-1] != (char) FILESEP )
        theparams->workdir[jj] =  (char) FILESEP;
      insert_wd(gname,theparams->workdir,theparams->map);
    }
    jj = isfile(theparams->map);
    if ( jj == 1 ) {
      fptr = fileopen(theparams->map, "r");
      if (fptr != NULL) {
        do {
	  ch = get_next_token(gname, MAXNAME, fptr);
        } while ((strcmp("-c", gname) != 0) && ch != EOF);
        if (ch != EOF)
	  ch = get_next_token(gname, MAXNAME, fptr);
        if (ch != EOF)
	  theparams->chrom = atoi(gname);
        fileclose(theparams->map, fptr);
      }
    }
    if ( theparams->workdir != NULL )
      shift_fn(theparams->map);

  }

}



/*
 1 => Allocate space for the parameter structure, read parameters from .qtlrc
 0 => Deallocate space for the parameter structure.
*/
params *create_params(aparams, cd, qtlrc)
params *aparams;
int cd;
char *qtlrc;
{
  int  isresource,verbosity;
  params *theparams;
  if (cd == 0) {
    verbosity = aparams->verbosity;
    theparams = aparams;
    if ( theparams->thedata != NULL ) {
      if ( whichprogram == 10 )
        free_indvector(theparams->thedata, theparams->nn+1);
      else
        free_indvector(theparams->thedata, theparams->nn);
    }
    if (theparams->theqtls != NULL)
      free_qtlvector(theparams->theqtls,theparams->themap);
    if (theparams->thegenome != NULL)
      clear_genome(theparams->thegenome);
    if ( theparams->themap != NULL )
        deallocate_markermap(theparams->themap);
	
    free_cvector(theparams->resource, 0, MAXNAME);
    free_cvector(theparams->error, 0, MAXNAME);
    free_cvector(theparams->map, 0, MAXNAME);
    free_cvector(theparams->mapin, 0, MAXNAME);
    free_cvector(theparams->qtl, 0, MAXNAME);
    free_cvector(theparams->eqtl, 0, MAXNAME);
    free_cvector(theparams->qtlin, 0, MAXNAME);
    free_cvector(theparams->ifile, 0, MAXNAME);
    free_cvector(theparams->qstat, 0, MAXNAME);
    free_cvector(theparams->iinfile, 0, MAXNAME);
    free_cvector(theparams->lrfile, 0, MAXNAME);
    free_cvector(theparams->srfile, 0, MAXNAME);
    free_cvector(theparams->zfile, 0, MAXNAME);
    free_cvector(theparams->mimfile, 0, MAXNAME);
    free_cvector(theparams->mqtfile, 0, MAXNAME);
    free_cvector(theparams->bayesfile, 0, MAXNAME);
    free_cvector(theparams->stem, 0, MAXNAME);
    free_cvector(theparams->thecross, 0, MAXNAME);
    free_cvector(theparams->mimwork, 0, 8);
    free_cvector(theparams->tfile, 0,MAXNAME);
    free_cvector(theparams->mrinput, 0,MAXNAME);
    free_cvector(theparams->mroutput, 0,MAXNAME);
    if ( theparams->workdir != NULL )
      free_cvector(theparams->workdir, 0, MAXNAME);
    if ( theparams->term != NULL )
      free_cvector(theparams->term, 0, MAXNAME);
    if ( debugging > 2 ) {
        sprintf(gwarn,"In create_params(), deallocated 1 param node  at %p\n",theparams);
        MemoryAccount(gwarn);
    }

    free((char *) theparams);
    return (NULL);
  }

/* Allocate space for theparams, and assign default values...*/

#if defined(MACWARRIOR) || defined(WINWARRIOR)  
  theparams = (params *) malloc((size_t) sizeof(params));
#else
  theparams = (params *) malloc((unsigned) sizeof(params));
#endif
    if ( debugging > 2 ) {
        sprintf(gwarn,"In create_params(), allocated 1 param node  at %p\n",theparams);
        MemoryAccount(gwarn);
    }
  
  theparams->resource = cvector(0, MAXNAME);
  theparams->error = cvector(0, MAXNAME);
  theparams->map = cvector(0, MAXNAME);
  theparams->mapin = cvector(0, MAXNAME);
  theparams->qtl = cvector(0, MAXNAME);
  theparams->eqtl = cvector(0, MAXNAME);
  theparams->qtlin = cvector(0, MAXNAME);
  theparams->ifile = cvector(0, MAXNAME);
  theparams->iinfile = cvector(0, MAXNAME);
  theparams->qstat = cvector(0, MAXNAME);
  theparams->lrfile = cvector(0, MAXNAME);
  theparams->srfile = cvector(0, MAXNAME);
  theparams->zfile = cvector(0, MAXNAME);
  theparams->mimfile = cvector(0, MAXNAME);
  theparams->stem = cvector(0, MAXNAME);
  theparams->thecross = cvector(0, MAXNAME);
  theparams->term = cvector(0, MAXNAME);
  theparams->mqtfile =  cvector(0, MAXNAME);
  theparams->bayesfile =  cvector(0, MAXNAME);
  theparams->workdir = cvector(0, MAXNAME);
  theparams->tfile = cvector(0, MAXNAME);
  theparams->mrinput = cvector(0, MAXNAME);
  theparams->mroutput = cvector(0, MAXNAME);
  theparams->mimwork = cvector(0,8);

  strcpy(theparams->stem, "qtlcart");
  renew_stem(theparams->stem,theparams);
  strcpy(theparams->thecross, "B1");
  if ( whichprogram == 8 ) 
    strcpy(theparams->mimwork,"ZMjbp\0");
  else if ( whichprogram == 5 )
    strcpy(theparams->mimwork,"no\0");
  else
    strcpy(theparams->mimwork,"smprtSeC");
    
  strcpy(theparams->term, "x11");

#if defined(MACWARRIOR)
  strcpy(theparams->term, "mac");
#endif
#if defined(WINWARRIOR)
  strcpy(theparams->term, "windows");
#endif

#if defined(AQUA)
  strcpy(theparams->term, "aqua");
#endif

  theparams->seed = get_a_seed();
  theparams->chrom = (int) DEFm;
  theparams->wchrom = 0;
  theparams->mark = (int) DEFl;
  theparams->vmark = (FPN) DEFsigl;
  theparams->dist = (FPN) DEFs;
  theparams->vdist = (FPN) DEFsigs;
  theparams->tail = (FPN) DEFbrdrs;
  theparams->qnum = (int) DEFqtl;
  theparams->dom = 1;
  theparams->beta = (FPN) DEFbeta;
  theparams->beta1 = (FPN) DEFbeta;
  theparams->beta2 = (FPN) DEFbeta;
  theparams->traits = 1;
  theparams->whichtrait = 1;
  theparams->reps = (int) REPS;
  theparams->Herit = (FPN) DEFheritability;
  theparams->Environ = (FPN) -1.0;
  theparams->cross = 1;
  theparams->crosst = 0;
  theparams->tcross = 0;
  theparams->tcrosst = 0;
  theparams->nn = (int) NN;
  theparams->Model = 3;
  theparams->walk = (FPN) DELTAX;
  theparams->Inter = 0;
  theparams->mapfunc = whosemf = (int) MAPFUNCTION;
  theparams->verbosity = 1;
  theparams->gout = 1;
  theparams->Rmode = 0;
  theparams->emethod = 10;
  theparams->emapobj = 0;
  theparams->linksize = (FPN) 0.25;
  theparams->segsize = (FPN) 0.01;
  theparams->ihypo = 10;
  theparams->boot = 0;
  theparams->siglevel = (FPN) SIG_LEVEL;
  theparams->size = (FPN) 0.05;
  theparams->nbp = (int) NUM_SIG;
  theparams->window = (FPN) WIN_SIZE;
  theparams->boots = (int) REPS;
  theparams->perms = (int) REPS;
  theparams->lodflag = 0;
  theparams->srm = 2;
  theparams->srupper = 100;
  theparams->srf1 = (FPN) 0.05;
  theparams->srb1 = (FPN) 0.05;
  theparams->maxeffect = theparams->maxlr = (FPN) 0.0;
  theparams->rwd = 0;
  theparams->mapparam = (FPN) 0.0;
  theparams->mimphase = 0;
  theparams->maxqtl = (int) MAXQTL ;
  theparams->maxepistatics = (int) MAXEPISTATICS  / 10 ;
  theparams->whoseic = 1;
  theparams->mimlod = (FPN) 0.0;
  theparams->null_sse = (FPN) 0.0;
/*  BTmapqtl*/
  theparams->p1=(FPN)0.5;
  theparams->p2=(FPN)0.5;
  theparams->p1_start=(FPN)0.1;
  theparams->p1_end=(FPN)0.9;
  theparams->p1_step=(FPN)0.1;
  theparams->p2_start=(FPN)0.1;
  theparams->p2_end=(FPN)0.9;
  theparams->p2_step=(FPN)0.1;

  theparams->themap = NULL;
  theparams->theqtls = NULL;
  theparams->thegenome = NULL;
  theparams->thedata = NULL;
  renew_stem(theparams->stem,theparams);
/*
  Now, open up the resource file and see if there are any overrides...
*/
  if (qtlrc == NULL) {
    isresource = isfile("qtlcart.rc");
    strcpy(theparams->resource, "qtlcart.rc");
    if (isresource == 1) 
      get_param_file(theparams, theparams->resource);
  }
  else {
    strcpy(theparams->resource, qtlrc );
    get_param_file(theparams, qtlrc);
  }

  return (theparams);
}

void get_param_file(theparams, qtlrc)
params *theparams;
char *qtlrc;
{
  int ch, ii;
  FILE *fileptr;


  fileptr = fileopen(qtlrc, "r");
  if (fileptr != NULL) {
    do {
      for (ii = 0; ii < MAXLINE; ii++)	/* clear gbuffer */
	    gbuffer[ii] = '\0';
      for (ii = 0; ((ch = fgetc(fileptr)) != EOF) && (ch != '\n'); ii++)	/* get a line */
	    gbuffer[ii] = (char) ch;
      if (gbuffer[0] == '-') {	/* if first character in line is - process it */
	    get_field(2, gname, gbuffer);
	    
      if ( !strncmp(gbuffer,"-chrom",6)  )      theparams->chrom = atoi(gname)  ;
      if ( !strncmp(gbuffer,"-cross",6)  )      theparams->cross = parse_cross(theparams,gname);
      if ( !strncmp(gbuffer,"-wchrom",7) )      theparams->wchrom = atoi(gname) ;
      if ( !strncmp(gbuffer,"-mark",5)   )      theparams->mark = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-vmark",6)  )      theparams->vmark = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-dist",5)  )       theparams->dist = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-vdist",6)  )      theparams->vdist = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-tail",5)  )       theparams->tail = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-qnum",5)  )       theparams->qnum = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-dom",4)  )        theparams->dom = atoi(gname)   ;
/*  Things that start with -beta*/
      if ( !strncmp(gbuffer,"-beta1",6)  )      theparams->beta1 = (FPN) atof(gname) ;
      else if ( !strncmp(gbuffer,"-beta2",6))   theparams->beta2 = (FPN) atof(gname) ;
      else if ( !strncmp(gbuffer,"-beta",5) )   theparams->beta = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-traits",7)  )     theparams->traits = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-whichtrait",11) ) theparams->whichtrait = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-reps",5)  )       theparams->reps = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-Hypothesis",11))  theparams->ihypo =   atoi(gname) ;
      if ( !strncmp(gbuffer,"-Herit",6)  )      theparams->Herit = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-Environ",8)  )    theparams->Environ = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-cross",6)  )      theparams->cross = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-nn",3)  )         theparams->nn = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-Model",6)  )      theparams->Model = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-walk",5)  )       theparams->walk = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-window",7)  )     theparams->window = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-nbp",4)  )        theparams->nbp = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-gout",5)  )       theparams->gout = atoi(gname)   ;
	  if ( !strncmp(gbuffer,"-Rmode",6)  )      theparams->Rmode = atoi(gname)   ;
	  if ( !strncmp(gbuffer,"-lodflag",8)  )    theparams->lodflag = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-srF",4)  )        theparams->srf1 = (FPN) atof(gname) ;
	  if ( !strncmp(gbuffer,"-srB",4)  )        theparams->srb1 = (FPN) atof(gname) ;
	  if ( !strncmp(gbuffer,"-srM",4)  )        theparams->srm = atoi(gname)   ;
	  if ( !strncmp(gbuffer,"-srupper",8)  )    theparams->srupper = atoi(gname)   ;
	  if ( !strncmp(gbuffer,"-size",5)  )       theparams->size = (FPN) atof(gname) ;
	  if ( !strncmp(gbuffer,"-siglevel",9)  )   theparams->siglevel = (FPN) atof(gname) ;
	  if ( !strncmp(gbuffer,"-LRmax",6)  )      theparams->maxlr = (FPN) atof(gname) ;
      if ( !strncmp(gbuffer,"-Effectmax",10) )  theparams->maxeffect = (FPN) atof(gname);
	  if ( !strncmp(gbuffer,"-maxqtl",7)  )     theparams->maxqtl = atoi(gname)   ;
	  if ( !strncmp(gbuffer,"-maxepis",8)  )    theparams->maxepistatics = atoi(gname)   ;
	  if ( !strncmp(gbuffer,"-lodmim",7)  )     theparams->mimlod = (FPN) atof(gname) ;
	  if ( !strncmp(gbuffer,"-ic",3)  )         theparams->whoseic = atoi(gname)   ;
	  if ( !strncmp(gbuffer,"-seed",5)  )  {     theparams->seed = atol(gname)   ; writeseed = 1; }
	  if ( !strncmp(gbuffer,"-Verbosity",10)  ) theparams->verbosity = atoi(gname);/* Script added 0, 1, 10 or 11 */
/* Verbosity can be 0 (non-verbose, non-automatic), 1 (verbose, non-automatic), 10 (non-verbose, automatic) or 11 (verbose, automatic) */
	  if ( !strncmp(gbuffer,"-rwd",4)  )        theparams->rwd = 1;                /* Script added rwd flag */
	  if ( !strncmp(gbuffer,"-phasemim",9)  )   theparams->mimphase = atoi(gname)   ;
      if ( !strncmp(gbuffer,"-workdir",8)  )     strcpy(theparams->workdir,gname);
      if ( !strncmp(gbuffer,"-stem",5)  )        renew_stem(gname, theparams);
      if ( !strncmp(gbuffer,"-Thecross",9)  )   theparams->cross = parse_cross(theparams,gname);
      if ( !strncmp(gbuffer,"-error",6)  )       strcpy(theparams->error,gname);
/*  Things that start with -map */            
      if ( !strncmp(gbuffer,"-mapfunc",8)  )    theparams->mapfunc = atoi(gname);
      else if ( !strncmp(gbuffer,"-mapin",6) )   strcpy(theparams->mapin,gname);
      else if ( !strncmp(gbuffer,"-map",4)  )    strcpy(theparams->map,gname);
/*  Things that start with -qtl */      
      if ( !strncmp(gbuffer,"-qtlin",6)  )       strcpy(theparams->qtlin,gname);
      else if ( !strncmp(gbuffer,"-qtl",4)  )    strcpy(theparams->qtl,gname);      
      if ( !strncmp(gbuffer,"-ifile",6)  )       strcpy(theparams->ifile,gname);
      if ( !strncmp(gbuffer,"-iinfile",8)  )     strcpy(theparams->iinfile,gname);
      if ( !strncmp(gbuffer,"-lrfile",7)  )      strcpy(theparams->lrfile,gname);
      if ( !strncmp(gbuffer,"-srfile",7)  )      strcpy(theparams->srfile,gname);
      if ( !strncmp(gbuffer,"-qstat",6)  )       strcpy(theparams->qstat,gname);
      if ( !strncmp(gbuffer,"-zfile",6)  )       strcpy(theparams->zfile,gname);
      if ( !strncmp(gbuffer,"-mimfile",8)  )     strcpy(theparams->mimfile,gname);
      if ( !strncmp(gbuffer,"-mqtfile",8)  )     strcpy(theparams->mqtfile,gname);
      if ( !strncmp(gbuffer,"-bayesfile",10)  )  strcpy(theparams->bayesfile,gname);
      if ( !strncmp(gbuffer,"-mrinput",8)  )     strcpy(theparams->mrinput,gname);
      if ( !strncmp(gbuffer,"-mroutput",9)  )    strcpy(theparams->mroutput,gname);
      if ( !strncmp(gbuffer,"-eqtlfile",9)  )    strcpy(theparams->eqtl,gname);
/*  Things that start with -emap */
      if ( !strncmp(gbuffer,"-emaplink",9)  )        theparams->linksize = (FPN) atof(gname);
      else if ( !strncmp(gbuffer,"-emapseg",8)  )     theparams->segsize = (FPN) atof(gname);      
      else if ( !strncmp(gbuffer,"-emapmeth",9)  )     theparams->emethod = atoi(gname);      
      else if ( !strncmp(gbuffer,"-emapobj",8)  )     theparams->emapobj = atoi(gname);      
	}
    } while (ch != EOF);
    fileclose(qtlrc, fileptr);
  }


}
/*
  This is only invoked by do_a_bootstrap and do_a_permutation in Prune.
  append the last random number to the qtlcart.rc file.  Also, create a file
  to save this same number.   Then, the next program will have it, and the next
  invocation of Prune will get it.   
*/
void append_seed( params *theparams,char *filename)
{ 
  FILE *fileptr;
  int i;
  i=0;
  if ( theparams->nn > 0 )
    i +=1;
  if (filename[0] != '\0') {
  
    fileptr = fileopen(filename, "a");
    if (fileptr != NULL) {
      fprintf(fileptr, "\n-seed         %ld  # (Seed for random number generator)", (long) ix );
      fprintf(fileptr, "\n#\n#  fooled you! \n#\n"); 
  
      fileclose(filename, fileptr);
    }
  }
}


/*
  write a new parameter file...
*/
void rewrite_param_file(theparams, filename)
params *theparams;
char *filename;
{
  FILE *fileptr;
  char *chptr;
  int i;
  if (filename[0] != '\0') {
    chptr = asctime2();

    fileptr = fileopen(filename, "w");
    if (fileptr != NULL) {
      fprintf(fileptr, "#\n#  This is the parameter file for QTLcartographer.");
      fprintf(fileptr, "\n#  It was regenerated %s#\n#        Some basic parameters\n#", chptr);

      if (theparams->chrom > 0)
	    fprintf(fileptr, "\n-chrom        %8d  # (Number of chromosomes: Rmap)", theparams->chrom);
      if (theparams->wchrom > 0)
	    fprintf(fileptr, "\n-wchrom        %8d  # (Chromosome to analyze)", theparams->wchrom);
      if (theparams->mark > 0)
	    fprintf(fileptr, "\n-mark         %8d  # (Average number of markers per chromosome: Rmap)", theparams->mark);
      if (theparams->vmark > 0.0)
	    fprintf(fileptr, "\n-vmark        %8.4f  # (Standard Deviation in the number of markers/chromosome: Rmap)", theparams->vmark);
      if (theparams->dist > 0.0)
	    fprintf(fileptr, "\n-dist         %8.4f  # (Average intermarker distance, in cM: Rmap)", theparams->dist);
      if (theparams->vdist > 0.0)
	    fprintf(fileptr, "\n-vdist        %8.4f  # (Standard Deviation of intermarker distance: Rmap)", theparams->vdist);
      if (theparams->tail > 0.0)
	    fprintf(fileptr, "\n-tail         %8.4f  # (Average amount of telomeric DNA in cM: Rmap)", theparams->tail);
      if (theparams->qnum > 0)
	    fprintf(fileptr, "\n-qnum         %8d  # (Average number of QTLs per trait: Rqtl)", theparams->qnum);
      if (theparams->dom > 0)
	    fprintf(fileptr, "\n-dom          %8d  # (type of dominance: Rqtl)", theparams->dom);
      if (theparams->beta > 0.0)
	    fprintf(fileptr, "\n-beta         %8.4f  # (Value of additive beta: Rqtl)", theparams->beta);
      if (theparams->beta1 > 0.0)
	    fprintf(fileptr, "\n-beta1        %8.4f  # (Value of dominance beta1: Rqtl)", theparams->beta1);
      if (theparams->beta2 > 0.0)
	    fprintf(fileptr, "\n-beta2        %8.4f  # (Value of dominance beta2: Rqtl)", theparams->beta2);
      if (theparams->traits > 0)
	    fprintf(fileptr, "\n-traits       %8d  # (Number of quantitative traits: Rqtl)", theparams->traits);
      if (theparams->whichtrait > 0)
	    fprintf(fileptr, "\n-whichtrait   %8d  # (The trait to be analyzed: Zmapqtl)", theparams->whichtrait);
      if (theparams->reps > 0)
	    fprintf(fileptr, "\n-reps         %8d  # (Number of replications in simulation: Rcross)", theparams->reps);
      if (theparams->Herit > 0.0)
	    fprintf(fileptr, "\n-Herit        %8.4f  # (Heritability: Rcross)", theparams->Herit);
      if (theparams->Environ > 0.0)
	    fprintf(fileptr, "\n-Environ      %8.4f  # (Environmental Variance: Rcross)", theparams->Environ);
      if (theparams->cross > 0)
	    fprintf(fileptr, "\n-cross        %8d  # (Type of cross: Rcross)", theparams->cross);
      if (theparams->nn > 0)
	    fprintf(fileptr, "\n-nn           %8d  # (Sample size: Rcross)", theparams->nn);
      if (theparams->Model > 0)
	    fprintf(fileptr, "\n-Model        %8d  # (Model for analysis: Zmapqtl)", theparams->Model);
      if (theparams->walk > 0.0)
	    fprintf(fileptr, "\n-walk         %8.4f  # (Walking speed along chromosomes, in cM: Zmapqtl)", theparams->walk);
      if (theparams->window > 0.0)
	    fprintf(fileptr, "\n-window       %8.4f  # (Window width, in cM: Zmapqtl)", theparams->window);
      if (theparams->nbp > 0)
	    fprintf(fileptr, "\n-nbp         %8d  # (Number of background parameters: Zmapqtl)", theparams->nbp);
      if (theparams->mapfunc > 0)
	    fprintf(fileptr, "\n-mapfunc      %8d  # (Map Function)", theparams->mapfunc);
      if (theparams->gout > 0)
	    fprintf(fileptr, "\n-gout         %8d  # (Output flag: Rmap)", theparams->gout);
	  fprintf(fileptr, "\n-Rmode        %8d  # (Simulation Mode: Rmap)", theparams->Rmode);
      fprintf(fileptr, "\n-emaplink     %8.4f  # (Size of test for linkage: Emap)", theparams->linksize );
      fprintf(fileptr, "\n-emapseg      %8.4f  # (Size of test for segregation: Emap)", theparams->segsize );      
      fprintf(fileptr, "\n-emapmeth     %8d  #  (Emap method flag)",  theparams->emethod );      
      fprintf(fileptr, "\n-emapobj      %8d  #  (Emap objective function, )",  theparams->emapobj );      
	  fprintf(fileptr, "\n-lodflag      %8d  # (1 => LOD scores in Preplot, Eqtl)", theparams->lodflag);
       /* if (theparams->verbosity > 0)
	  fprintf(fileptr, "\n-verbosity    %8d  # (Verbosity flag)", theparams->verbosity); */
	  fprintf(fileptr, "\n-srF          %8.4f  # (p(Fin): SRmapqtl)", theparams->srf1);
	  fprintf(fileptr, "\n-srB          %8.4f  # (p(Fout): SRmapqtl)", theparams->srb1);
	  fprintf(fileptr, "\n-srM          %8d  # (Regression type: SRmapqtl)", theparams->srm);
	  fprintf(fileptr, "\n-srupper      %8d  # (Maximun number of steps in stepwise regression: SRmapqtl)", theparams->srupper);
	  fprintf(fileptr, "\n-size         %8.4f  # (Size, or alpha)", theparams->size);
	  fprintf(fileptr, "\n-siglevel     %8.4f  # (Size, or alpha)", theparams->siglevel);
	  fprintf(fileptr, "\n-LRmax        %8.4f  # Maximum LR or LOD score", theparams->maxlr);
	  fprintf(fileptr, "\n-Effectmax    %12.6f  # Maximum additive effect ", theparams->maxeffect);
	  fprintf(fileptr, "\n-maxqtl       %8d  # (Maximum QTL to fit in MImapqtl)", theparams->maxqtl);
	  fprintf(fileptr, "\n-maxepis      %8d  # (Maximum Epistatic effect to fit in MImapqtl)", theparams->maxepistatics);
	  fprintf(fileptr, "\n-lodmim       %8.4f  # (LOD for adding/dropping parameters in MImapqtl)", theparams->mimlod);
	  fprintf(fileptr, "\n-ic           %8d  # (Code for IC criterion in MImapqtl)", theparams->whoseic);
	  fprintf(fileptr, "\n-phasemim     %8d  # (Phase of analysis MImapqtl)", theparams->mimphase);
	  fprintf(fileptr, "\n-Hypothesis   %8d  # (Hypothesis test to do/process)", theparams->ihypo);
      if ( writeseed == 1)
        fprintf(fileptr, "\n-seed         %ld  # (Seed for random number generator)", (long) ix);

      fprintf(fileptr, "\n#\n# These are the default filenames...\n#");
      if ( (i=(int) strlen(theparams->workdir)) > 0    )
	    fprintf(fileptr, "\n-workdir %50s  # (The working directory)", theparams->workdir);

      if (theparams->stem[0] != '\0')
	    fprintf(fileptr, "\n-stem      %25s  # (Stem for filenames)", theparams->stem);
      if (theparams->thecross[0] != '\0')
	    fprintf(fileptr, "\n-Thecross  %25s  # (Type of cross)", theparams->thecross);
      if (theparams->error[0] != '\0')
	    fprintf(fileptr, "\n-error     %25s  # (Log and error file: All)", theparams->error);
      if (theparams->map[0] != '\0')
	    fprintf(fileptr, "\n-map       %25s  # (Rmap ouput file, linkage map)", theparams->map);
      if (theparams->mapin[0] != '\0')
	    fprintf(fileptr, "\n-mapin     %25s  # (Rmap input file)", theparams->mapin);
      if (theparams->qtl[0] != '\0')
	    fprintf(fileptr, "\n-qtl       %25s  # (Rqtl output file)", theparams->qtl);
      if (theparams->qtlin[0] != '\0')
	    fprintf(fileptr, "\n-qtlin     %25s  # (Rqtl input file)", theparams->qtlin);
      if (theparams->ifile[0] != '\0')
	    fprintf(fileptr, "\n-ifile     %25s  # (Rcross output file, individual data file)", theparams->ifile);
      if (theparams->iinfile[0] != '\0')
	    fprintf(fileptr, "\n-iinfile   %25s  # (Rcross input file, individual data file)", theparams->iinfile);
      if (theparams->lrfile[0] != '\0')
	    fprintf(fileptr, "\n-lrfile    %25s  # (Results of Linear Regression analysis)", theparams->lrfile);
      if (theparams->srfile[0] != '\0')
	    fprintf(fileptr, "\n-srfile    %25s  # (Results of Stepwise Regression analysis)", theparams->srfile);
      if (theparams->qstat[0] != '\0')
	    fprintf(fileptr, "\n-qstat     %25s  # (Results of Qstats)", theparams->qstat);
      if (theparams->zfile[0] != '\0')
	    fprintf(fileptr, "\n-zfile     %25s  # (Results of (Composite) Interval Mapping Analysis)", theparams->zfile);
      if (theparams->eqtl[0] != '\0')
	    fprintf(fileptr, "\n-eqtlfile     %25s  # (Eqtl output file)", theparams->eqtl);
      if (theparams->mimfile[0] != '\0')
	    fprintf(fileptr, "\n-mimfile   %25s  # (Results of Multiple Interval Mapping Analysis)", theparams->mimfile);
      if (theparams->mqtfile[0] != '\0')
	    fprintf(fileptr, "\n-mqtfile   %25s  # (Multiple Interval Mapping Analysis Model)", theparams->mqtfile);
      if (theparams->bayesfile[0] != '\0')
	    fprintf(fileptr, "\n-bayesfile %25s  # (Results of Bayesian Analysis)", theparams->bayesfile);
      if (theparams->mrinput[0] != '\0')
	    fprintf(fileptr, "\n-mrinput %25s  # (Output of JZmapqtl for input to MultiRegress)", theparams->mrinput);
      if (theparams->mroutput[0] != '\0')
	    fprintf(fileptr, "\n-mroutput %25s  # (Results of MultiRegress)", theparams->mroutput);
      fprintf(fileptr, "\n#\n#  end of file\n#\n");
      fileclose(filename, fileptr);
    }
  }
}

/*

Determine filetype from first line, if possible, otherwise go deeper.

Return isf 


isf          Meaning           Stem Code         Notes
-----------------------------------------
-1           no file
 0           unknown type
 1           QTLcart.log       stem.log
 2           qtlcart.rc        stem.rc
 3           qtlcart.hlp       stem.hlp
 4           qtlcart.mcd       stem.mcd
 10          Rmap.out          stem.map
 11           map.inp          stemm.inp             Emap produces one in Stage 1
 12           mapmaker.maps    stem.maps or stem.mps
 13           Chrom.plt        stemmap.plt
 20          Rqtl.out          stem.qtl or stem.mqt
 21           qtls.inp         stemq.inp
 30          Rcross.out        stem.cro
 31           cross.inp        stemc.inp             Emap produces one in Stage 1
 32           mapmaker.raw     stem.raw
 33           RSplus.inp       stem.r
 34           qtlcart.sas      stem.sas
 40          Qstats.out        stem.qst
 50          LRmapqtl.out      stem.lr
 51          SRmapqtl.out      stem.sr
 60          Zmapqtl.out       stem.zi      +
 61          ZipermE.out       stem.zie     |______ i can be 1-7
 62          ZipermC.out       stem.zic     |          In practice, i can only be 3
 63          Ziboot.out        stem.zib     |
 63          Zijack.out        stem.zMi     +
 65          JZmapqtl.out      stem.zi     t = 0,traits; i = 1-7
 66          JZmapqtl.zr       stem.zr     markers replaced by expected values at regular intervals
 67          MultiRegress.out  stem.mr     MultiRegress output
 70          Eqtl.out          stem.eqt       
 71          Eqtl.tex          stem.tex       the tex and html output formats are on the 
 72          Eqtl.html         stem.htm       wish list.
 73          Ziboots.out       stem.zid     summary file for bootstraps
 74          Zijacks.out       stem.zMi     summary file for jackknives
 80          Preplot.plt       stem.plt
 81          c#t#.?#           c#t#.?# 
 90          HKmapqtl.out      stem.hi      +  these don't exist yet, but I'm hoping they will.
 91          HKipermE.out      stem.hie     |______ i can be 1-7
 92          HKipermC.out      stem.hic     |
 93          HKiboot.out       stem.hib     +
100          MImapqtl.out      stem.mim   
110          Bmapqtl.out       stem.bys  
*/

int get_file_type(filename)
  char *filename;
{

  int isf,ch,cntr;

  FILE *fptr;
  cntr = 0;
  isf = -1;
  if (filename[0] != '\0') {
    fptr = fopen(filename, "r");
    if ( fptr != NULL ) {
      do { /* read tokens on first line.  If we come accross an EOF or a \n, then get out*/
        ch=get_next_token(gname, MAXNAME, fptr);
        if ( ch == '\n' ) 
          cntr = cntr+1;
        else if ( gname[0] == '-' ) {/*If we come accross a -filetype,  put next token in gname and then get out*/
          if ( !strcmp("-filetype",gname) ) {
            get_next_token(gname, MAXNAME, fptr);
            isf = 0;
          }
        }
        else if ( !strcmp("#FileID",gname )  )
          isf = file_to_int("qtlcart.mcd");
        
      } while (ch != EOF && isf < 0 && cntr < 100 );
      fclose(fptr);
      if ( ch !=EOF && cntr < 100 )
        isf = file_to_int(gname);
      else 
        isf = 0;
    }
  }
  return(isf);

}

/*
   Assign an integer value to the file type.  Use the same code as above.
*/
int file_to_int(xtemp)
  char *xtemp;
{
  int isf;
  isf = 0;
  if (!strcmp("cross.inp",xtemp)     )  isf=31; 
  if (!strcmp("RSplus.inp",xtemp)     )  isf=33; 
  if (!strcmp("qtlcart.sas",xtemp)     )  isf=34; 
  if (!strcmp("plabqtl0.qdt",xtemp)     )  isf=35; 
  if (!strcmp("plabqtl1.qdt",xtemp)     )  isf=36; 
  if (!strcmp("c#t#.?#",xtemp)       )  isf=81;
  if (!strcmp("map.inp",xtemp)       )  isf=11; 
  if (!strcmp("mapmaker.maps",xtemp) )  isf=12; 
  if (!strcmp("mapmaker.raw",xtemp)  )  isf=32;
  if (!strcmp("qtls.inp",xtemp)      )  isf=21; 
  if (!strcmp("qtlcart.hlp",xtemp)   )  isf=3; 
  if (!strcmp("qtlcart.rc",xtemp)    )  isf=2;
  if (!strcmp("Chrom.plt",xtemp)     )  isf=13; 
  if (!strcmp("Eqtl.out",xtemp)      )  isf=70; 
  if (!strcmp("Eqtl.tex",xtemp)      )  isf=71; 
  if (!strcmp("Eqtl.html",xtemp)     )  isf=72;
  if (!strcmp("HKmapqtl.out",xtemp)  )  isf=90; 
  if (!strcmp("HKipermE.out",xtemp)  )  isf=91; 
  if (!strcmp("HKipermC.out",xtemp)  )  isf=92;
  if (!strcmp("HKiboot.out",xtemp)   )  isf=93;
  if (!strcmp("LRmapqtl.out",xtemp)  )  isf=50;  
  if (!strcmp("BTmapqtl.out",xtemp) )  isf=52;
  if (!strcmp("Preplot.plt",xtemp)   )  isf=80;  
  if (!strcmp("QTLcart.log",xtemp)   )  isf=1; 
  if (!strcmp("qtlcart.mcd",xtemp)   )  isf=4; 
  if (!strcmp("#FileID",xtemp)   )  isf=4; /* a trick for mcd files.*/
  if (!strcmp("Qstats.out",xtemp)    )  isf=40; 
  if (!strcmp("Rmap.out",xtemp)      )  isf=10; 
  if (!strcmp("Rqtl.out",xtemp)      )  isf=20; 
  if (!strcmp("Rcross.out",xtemp)    )  isf=30;
  if (!strcmp("SRmapqtl.out",xtemp)  )  isf=51;  
  if (!strcmp("Zmapqtl.out",xtemp)   )  isf=60; 
  if (!strcmp("ZipermE.out",xtemp)   )  isf=61; 
  if (!strcmp("ZipermC.out",xtemp)   )  isf=62;
  if (!strcmp("Ziboot.out",xtemp)    )  isf=63;
  if (!strcmp("Zijack.out",xtemp)    )  isf=64;
  if (!strcmp("JZmapqtl.out",xtemp)    )  isf=65;
  if (!strcmp("JZmapqtl.zr",xtemp)    )  isf=66;
  if (!strcmp("MultiRegress.out",xtemp)    )  isf=67;
  if (!strcmp("JZmapqtl.zp",xtemp)    )  isf=68;
  
  if (!strcmp("Ziboots.out",xtemp)    )  isf=73;
  if (!strcmp("Zijacks.out",xtemp)    )  isf=74;
  if (!strcmp("MImapqtl.out",xtemp)    )  isf=100;
  if (!strcmp("Bmapqtl.out",xtemp)    )  isf=110;
  return(isf);
}

/*
  Write the name of the filetype code.  ft is an integer.
*/
void write_file_type(ft,fptr) 
  int ft;
  FILE *fptr;
{
    fprintf(fptr," -filetype ");
	if ( ft == 0 )   fprintf(fptr,"unknown");  
	if ( ft ==  1)   fprintf(fptr,"QTLcart.log");  
	if ( ft ==  2)   fprintf(fptr,"qtlcart.rc");  
	if ( ft ==  3)   fprintf(fptr,"qtlcart.hlp"); /* Obsolete*/  
	if ( ft ==  4)   fprintf(fptr,"qtlcart.mcd");   
	if ( ft ==  10)  fprintf(fptr,"Rmap.out");  
	if ( ft ==  11)  fprintf(fptr,"map.inp");  
	if ( ft ==  12)  fprintf(fptr,"mapmaker.maps");  
	if ( ft ==  13)  fprintf(fptr,"Chrom.plt");  
	if ( ft ==  20)  fprintf(fptr,"Rqtl.out");  
	if ( ft ==  21)  fprintf(fptr,"qtls.inp");  
	if ( ft ==  30)  fprintf(fptr,"Rcross.out");  
	if ( ft ==  31)  fprintf(fptr,"cross.inp");  
	if ( ft ==  32)  fprintf(fptr,"mapmaker.raw");  
	if ( ft ==  33)  fprintf(fptr,"RSplus.inp");  
	if ( ft ==  34)  fprintf(fptr,"qtlcart.sas");  
	if ( ft ==  35)  fprintf(fptr,"plabqtl0.qdt");  
	if ( ft ==  36)  fprintf(fptr,"plabqtl1.qdt");  
	if ( ft ==  40)  fprintf(fptr,"Qstats.out");  
	if ( ft ==  50)  fprintf(fptr,"LRmapqtl.out");  
	if ( ft ==  51)  fprintf(fptr,"SRmapqtl.out");  
	if ( ft ==  52)  fprintf(fptr,"BTmapqtl.out");
	if ( ft ==  60)  fprintf(fptr,"Zmapqtl.out");  
	if ( ft ==  61)  fprintf(fptr,"ZipermE.out");  
	if ( ft ==  62)  fprintf(fptr,"ZipermC.out");  
	if ( ft ==  63)  fprintf(fptr,"Ziboot.out");  
	if ( ft ==  64)  fprintf(fptr,"Zijack.out");  
	if ( ft ==  65)  fprintf(fptr,"JZmapqtl.out");  
	if ( ft ==  66)  fprintf(fptr,"JZmapqtl.zr");  
	if ( ft ==  67)  fprintf(fptr,"MultiRegress.out");  
	if ( ft ==  68)  fprintf(fptr,"JZmapqtl.zp");  
	if ( ft ==  70)  fprintf(fptr,"Eqtl.out");  
	if ( ft ==  71)  fprintf(fptr,"Eqtl.tex");  
	if ( ft ==  72)  fprintf(fptr,"Eqtl.html");  
	if ( ft ==  73)  fprintf(fptr,"Ziboots.out");  	
	if ( ft ==  74)  fprintf(fptr,"Zijacks.out");  	
	if ( ft ==  80)  fprintf(fptr,"Preplot.plt");  
	if ( ft ==  81)  fprintf(fptr,"c#t#.?#");   
	if ( ft ==  90)  fprintf(fptr,"HKmapqtl.out");  /* The HK stuff is mostly in MultiRegress */
	if ( ft ==  91)  fprintf(fptr,"HKipermE.out");  
	if ( ft ==  92)  fprintf(fptr,"HKipermC.out");  
	if ( ft ==  93)  fprintf(fptr,"HKiboot.out");     
	if ( ft ==  100)  fprintf(fptr,"MImapqtl.out");     
	if ( ft ==  110)  fprintf(fptr,"Bmapqtl.out");     
    fprintf(fptr," ");
}


/*
Print a standard output file header.

outmode = 0 => overwrite
          1 => append         to the output file.
          
chptr is the time


*/
void print_head(char *prog,char *filename,char *chptr,int outmode,int filetype,params *theparams)
{
  FILE *outf;
  outf = NULL;
  if (outmode==1)
    outf = fileopen(filename, "a");
  else
    outf = fileopen(filename, "w");
  if (outf == NULL)
    outf = fileopen(filename, "w");
  if (outf != NULL) {
    fprintf(outf, "#   %12ld  ", theparams->seed );
    write_file_type(filetype,outf);
    fprintf(outf,"\n#\n#\t%s\n#\tThis output file (%s) was created by %s...", VERSION,filename,prog);
    fprintf(outf, "\n#\n#\tIt is %s#\n#", chptr);
    fileclose(filename, outf);
  }

}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file params.c
------------------------------------------------------------------ */


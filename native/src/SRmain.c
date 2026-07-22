/* ------------------------------------------------------ XCutXCodeXSkip
     This file (SRmain.c) is part of QTL Cartographer
         
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



/*  Main function driver and subroutines for SRmapqtl.*/

#include "Main.h"

int main(argc, argv)
int argc;
char *argv[];
{
  FILE *errorf;
  char *chptr,*purpose;
  int nopts,jj,ii,kk,automatic,rows,firsttrait,lasttrait,upper,minsamp,explan;
  params *theparams;
  genome  *agptr   ;
  linpakws *lnpk;

#if defined(MACWARRIOR)  
 /* Simulate the UNIX shell with ccommand */
  argc = ccommand(&argv);
#endif
  whichprogram = 61;
/* just to be on the safe side, initialize pointers to null */

  purpose = cvector(0,MAXNAME);

 /**/
  strcpy(purpose, "Do Stepwise Linear Regression");
  theparams = NULL;
  nopts = 10;
  theparams = create_params(theparams, 1, NULL);
  agptr= theparams->thegenome=NULL;
  theparams->thedata   = NULL;
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);


/* Initailize the data structures...*/
  GetTheMap(theparams, theparams->map);
  GetTheData(theparams, theparams->ifile, 1 );

      
  rows = 1;    /* There will be at least one row for the mean. */
  if ( theparams->themap->otraits > 0 ) {  /*There will also be rows for the Otraits. */
      process_otraits(theparams->thedata,theparams->nn,theparams->themap);
      for ( ii = 1 ; ii <= theparams->themap->otraits ; ii++ )
        if ( theparams->themap->onames[ii][0] == '+' )
          rows = rows +  theparams->themap->otypes[ii] - 1;
  }
  explan = rows;
/* Determine the trait and type of cross, print them to the output file */
  if (theparams->cross==3 || theparams->cross==4)
    rows = rows+2*theparams->themap->ml; /*2 rows per marker if 3 genotypic classes  */
  else
    rows = rows+theparams->themap->ml; /*1 row  per marker if 2 genotypic classes  */
/*
  At this point we should be able to figure out an upper bound to 
  how much workspace we need.  
*/
  lnpk = cd_linpakws(NULL,theparams->nn,rows,1);	  
/*	Create a linked list representing the genome and
	translate the distances in recombination frequencies to those
	of Morgans...  */
  theparams->thegenome = create_genome(theparams->themap);
  agptr = theparams->thegenome;
/*This would be a good time to decide whether we can do forward, backward or forward-backward
   stepwise regression.   
  if ( theparams->cross ==   3 || theparams->cross == 4 )
    j = 2;
  else
    j=1;
  upper = theparams->themap->ml*j + rows;*/  
  if ( theparams->whichtrait > theparams->traits ) {
    firsttrait = 1;
    lasttrait = theparams->traits;
  }
  else if ( theparams->whichtrait > 0 )
    firsttrait = lasttrait = theparams->whichtrait;
  minsamp = 0;
  for ( ii = firsttrait ; ii <= lasttrait ; ii++ ) {
    theparams->whichtrait = ii;
	kk =  calc_samplesize(theparams,theparams->thedata);
	if ( kk > minsamp ) 
	  minsamp = kk;
  }
  upper = minsamp - explan - 2;
  if ( theparams->cross ==   3 || theparams->cross == 4 )
    upper = upper/2  ;  /* maximum number of markers that can be added */
  if ( upper > theparams->themap->ml )
    upper = theparams->themap->ml;  /*Now it's the number of steps.  This means all markers can be added.  Else, upper tells how many. */
  if ( theparams->verbosity == 1 )
    printf("\n\n\tMaximal number of markers that can be added is %d\n",upper);
  if ( upper > theparams->srupper) {
    if ( theparams->verbosity == 1 )
      printf("\tResetting to your specified upper bound of %d\n",theparams->srupper);
      upper = theparams->srupper;
  
  }
  if ( theparams->themap->ml > upper && theparams->srm == 1 ) {
    theparams->srm = 0;
    errorf = fileopen(theparams->error, "a");
    fprintf(errorf, "\n\nToo many variables to do backward stepwise regression...doing forward stepwise regression instead.");
    fprintf(errorf, "\n\tMinimum sample size over selected traits: %10d",minsamp);
    fprintf(errorf, "\n\tMarkers    in first step of regression:   %10d",theparams->themap->ml);
    fprintf(errorf, "\n\tParameters in first step of regression:   %10d",rows);
    fprintf(errorf, "\n\n");
    fileclose(theparams->error, errorf);
    if ( theparams->verbosity == 1 ) {
      printf("\n\tSorry, there are too many parameters to do a backward stepwise regression...");
      printf("\n\tDoing forward stepwise regression instead.\n");
    }
  }
    /*  Write a header for the output file.  */
  print_head(argv[0],theparams->srfile,chptr,1,51,theparams);
  explan = explan-1; /* number of rows due to otraits...the 1 is for the mean.*/
  for (jj = firsttrait; jj <= lasttrait; jj++) {
    theparams->whichtrait = jj;  
    if (theparams->verbosity == 1) {
      if ( theparams->srm == 2 )
	    printf("\n\tFB Stepwise regression analysis for trait %d",jj);
      else if ( theparams->srm == 1 )
        printf("\n\tBackward Stepwise regression analysis for trait %d",jj);
      else  
	    printf("\n\tForward Stepwise regression analysis for trait %d",jj);
    }
    if ( theparams->srm == 2 )
	  agptr = for_back_swr(theparams->thedata,theparams,theparams->themap,theparams->thegenome,lnpk,theparams->themap->ml,explan,upper);
    else if ( theparams->srm == 1 )
      agptr = backward_swr(theparams->thedata,theparams,theparams->themap,theparams->thegenome,lnpk,theparams->themap->ml,explan,upper);
    else  
	  agptr = forward_swr(theparams->thedata,theparams,theparams->themap,theparams->thegenome,lnpk,theparams->themap->ml,explan,upper);
	kk =  calc_samplesize(theparams,theparams->thedata);
    write_sr_results(agptr,theparams,theparams->srfile,kk,explan);
  }
  write_trailer(theparams,chptr,1);
  
  /* Clean up...*/
  if ( purpose != NULL )
	 free_cvector(purpose,0,MAXNAME);
  lnpk = cd_linpakws(lnpk,theparams->nn,rows,0);
  theparams = create_params(theparams, 0, NULL);

  return(0);
}

void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag)
{
  int ii, jj;
  if (flag == 1) {
    strcpy(opt[1], "-i");
    strcpy(opt[2], "-o");
    strcpy(opt[3], "-e");
    strcpy(opt[4], "-m");
    strcpy(opt[5], "-s");
    strcpy(opt[6], "-M");
    strcpy(opt[7], "-t");
    strcpy(opt[8], "-F");
    strcpy(opt[9], "-B");
    strcpy(opt[10], "-u");

    strcpy(opt_e[1], "Input File");
    strcpy(opt_e[2], "Output File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Genetic Linkage Map File");
    strcpy(opt_e[5], "Random Number Seed");
    strcpy(opt_e[6], "FS, BS or FB (0,1,2)?");
    strcpy(opt_e[7], "Trait to analyze");
    strcpy(opt_e[8], "Size: p(Fin) = ");
    strcpy(opt_e[9], "Size: p(Fout) = ");
    strcpy(opt_e[10], "Hard bound on number of forward steps ");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->ifile);
  strcpy(opt_v[2], theparams->srfile);
  strcpy(opt_v[3], theparams->error);
  strcpy(opt_v[4], theparams->map);
  sprintf(opt_v[5], "%ld", theparams->seed);
  sprintf(opt_v[6], "%d", theparams->srm);
  sprintf(opt_v[7], "%d", theparams->whichtrait);
  sprintf(opt_v[8], "%f", theparams->srf1);
  sprintf(opt_v[9], "%f", theparams->srb1);
  sprintf(opt_v[10], "%d", theparams->srupper);
}

void update_params(char **opt_v,  params *theparams)
{
  theparams->srm = atoi(opt_v[6]);

  strcpy(theparams->ifile, opt_v[1]);

  strcpy(theparams->srfile, opt_v[2]);
  strcpy(theparams->error, opt_v[3]);
  strcpy(theparams->map, opt_v[4]);
  theparams->seed = atol(opt_v[5]);
  theparams->whichtrait = atoi(opt_v[7]);
  theparams->srf1 = (FPN)atof(opt_v[8]);
  theparams->srb1 = (FPN)atof(opt_v[9]);
  theparams->srupper = atoi(opt_v[10]);

}



/* ------------------------------------------------------- XCutXCodeXSkip
             End of file SRmain.c
------------------------------------------------------------------ */


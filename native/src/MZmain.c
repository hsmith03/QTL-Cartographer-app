/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MZmain.c) is part of QTL Cartographer
         
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


/*  Main Driver for JZmapqtl*/
#include "Main.h"

int main(argc, argv)
int argc;
char *argv[];
{
  char *chptr,*purpose;
  int nopts,ii,automatic,rows,go_on,do_analysis,t ,marks,ipos,**srranks,srm;
  long fileposition;
  FPN *tmp, dumm;
  params *theparams;
  genome  *agptr,*gptr ,*startptr,*endptr ;
  linpakws *lnpk;
#ifndef WIN32
#if defined(MACWARRIOR)  
 /*  Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
#endif
  whichprogram = 71;
  /* just to be on the safe side, initialize pointers to null */
  tmp =  NULL;
  agptr=gptr=startptr=endptr=NULL;
  theparams = NULL;
  /*   Create space for purpose */
  purpose = cvector(0,MAXNAME);
  strcpy(purpose, "Multitrait CIM");
  /*  Create the parameters structure, get the time and process the command line args. */
  nopts = 14;
  theparams = create_params(theparams,1,NULL);
  chptr = asctime2(); 
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  /*   Walk in units of Morgans... */
  theparams->walk = theparams->walk / (FPN) 100.0;
  srm = theparams->srm;
  if ( theparams->Model == 8 ) { /* 8 is a model to tell JZmapqtl to use srm model 3.  */
    theparams->Model = 6;
    theparams->srm = 3;
  }
  else if ( theparams->Model == 10 ) { /* 10 is a model to tell JZmapqtl to use srm model 4.  */
    theparams->Model = 6;
    theparams->srm = 4;
  }
  /* Initialize the data structures...*/
  GetTheMap(theparams, theparams->map);
  if ( theparams->Model == 7 || theparams->ihypo > 39 )  /*  Model 7 requires the output of Eqtl   */
	GetTheModel(theparams, theparams->eqtl);
  else
    theparams->theqtls = NULL;

  GetTheData(theparams, theparams->ifile, 1);

/*	Create a linked list representing the genome and
	translate the distances in recombination frequencies to those
	of Morgans...  */
  theparams->thegenome = create_genome(theparams->themap);
  agptr = theparams->thegenome;
  if ( theparams->Model == 7 || theparams->ihypo > 39 )
     zplace_qtls(theparams,theparams->theqtls, agptr);
/* Determine the trait and type of cross, print them to the output file */


  if ( theparams->Model == 6 ) {/* Get the results from a previous run of SRmapqtl */
    srranks = jzget_srresults(theparams,theparams->themap);
    marks = set_cofactors(theparams,theparams->themap,theparams->thegenome,srranks);
    if (theparams->verbosity == 1 ) 
      show_cofactors(stdout,theparams,theparams->themap,theparams->thegenome,srranks) ;
    if ( marks < 1 ) {
      theparams->Model = 3;
      free_imatrix(srranks,1,theparams->themap->ml,0,theparams->traits);
    }
  }
  /* Determine the number of rows in the design matrix X */
  rows = how_many_rows(theparams,theparams->themap,theparams->thedata,&ii,marks);
  /* How many traits in a multitrait analysis? */
  t = how_many_traits(theparams,theparams->themap);
  /* At this point we should be able to figure out an upper bound to how much workspace we need. */
  lnpk = cd_jzmapws(NULL,theparams->nn,rows,t,1);	  
  lnpk->pp1 = dvector(1,theparams->nn);
  lnpk->pp2 = dvector(1,theparams->nn);
  lnpk->pv = dvector(1,theparams->nn);
  lnpk->qv = dvector(1,theparams->nn);
  lnpk->bp = imatrix(1,2,1,theparams->themap->ml);
  lnpk->wrsd = dmatrix(0,lnpk->t,1,theparams->nn);
  lnpk->wy = dvector(1,theparams->nn);
  lnpk->k = ii;
  lnpk->samplesize = ivector(0, lnpk->t);
  lnpk->work = dvector( 0, lnpk->t);    
  lnpk->kpvt = ivector(0, lnpk->t);  
  lnpk->ts0 = dvector(0,lnpk->t);  
/*  put phenotypes in lnkp->y.  */
  copy_phenotypes(theparams, theparams->themap, theparams->thedata, lnpk);
  init_xsave(theparams, theparams->themap, theparams->thedata, lnpk);
/*  Write a header for the output file.  */
  theparams->perms = theparams->boots = 0; /*Disable permutations and bootstraps for the time being...*/
  if ( theparams->Model == 9 ) {
    if ( theparams->themap->otraits > 0 )
      process_otraits(theparams->thedata,theparams->nn,theparams->themap);
    fileposition = write_altheader(theparams->mrinput, theparams, argv[0], chptr,theparams->themap,1);
    write_alttraits(theparams->mrinput, theparams, theparams->themap, theparams->thedata, lnpk);  
    go_on = determine_endpoints(theparams,theparams->themap,&startptr,&endptr,theparams->thegenome,&do_analysis);
    ipos = do_jzexpecteds(theparams,startptr,endptr,theparams->themap, theparams->thedata,theparams->mrinput, lnpk);
    insert_positions(theparams->mrinput,ipos,fileposition) ;
  }
  else {
    write_jzheader(theparams->zfile, theparams, argv[0], chptr,theparams->themap,1,theparams->thegenome,srranks);
  /*  Here's where we do the analysis.  */
    go_on = determine_endpoints(theparams,theparams->themap,&startptr,&endptr,theparams->thegenome,&do_analysis);
    if (do_analysis == 1) {
	  dumm = do_jzanalysis(theparams,startptr,endptr,theparams->themap,theparams->theqtls,theparams->thedata,theparams->zfile,agptr,lnpk );
	  write_jzheader(theparams->zfile, theparams, argv[0], chptr,theparams->themap,0,theparams->thegenome,srranks);
    }
  } 
/* Clean up...*/
  theparams->srm = srm;

  theparams->walk = theparams->walk * (FPN) 100.0;
  theparams->Inter = 0;
  theparams->reps = (int) REPS;
  write_trailer(theparams,chptr,1);

  lnpk->k = theparams->themap->ml;  /* Need this for proper deallocation...*/
  if ( theparams->Model == 6 )
      free_imatrix(srranks,0,theparams->themap->ml,0,theparams->traits);

  if ( purpose != NULL )
	 free_cvector(purpose,0,MAXNAME);

  lnpk = cd_jzmapws(lnpk,theparams->nn,rows,lnpk->t,0);
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
    strcpy(opt[5], "-S");
    strcpy(opt[6], "-E");
    strcpy(opt[7], "-s");
    strcpy(opt[8], "-M");
    strcpy(opt[9], "-t");
    strcpy(opt[10], "-c");
    strcpy(opt[11], "-d");
    strcpy(opt[12], "-n");
    strcpy(opt[13], "-w");
    strcpy(opt[14], "-I");


    strcpy(opt_e[1], "Input File");
    strcpy(opt_e[2], "Output File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Genetic Linkage Map File");
    strcpy(opt_e[5], "SRmapqtl results (Model 6)");
    strcpy(opt_e[6], "Eqtl results (Model 7)");
    strcpy(opt_e[7], "Random Number Seed");
    strcpy(opt_e[8], "Model [3,6,7], 3=>IM");
    strcpy(opt_e[9], "Trait to analyze");
    strcpy(opt_e[10], "Chromosome to analyze (0=>all)");
    strcpy(opt_e[11], "Walking speed in cM");
    strcpy(opt_e[12], "Number of Background Parameters (Model 6)");
    strcpy(opt_e[13], "Window Size in cM (Model 6)");
    strcpy(opt_e[14], "Hypothesis test");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->ifile);
  strcpy(opt_v[2], theparams->zfile);
  strcpy(opt_v[3], theparams->error);
  strcpy(opt_v[4], theparams->map);
  strcpy(opt_v[5], theparams->srfile);
  strcpy(opt_v[6], theparams->eqtl);
  sprintf(opt_v[7], "%ld", theparams->seed);
  sprintf(opt_v[8], "%d", theparams->Model);
  sprintf(opt_v[9], "%d", theparams->whichtrait);
  sprintf(opt_v[10], "%d", theparams->wchrom);
  sprintf(opt_v[11], "%f", theparams->walk);
  sprintf(opt_v[12], "%d", theparams->nbp);
  sprintf(opt_v[13], "%f", theparams->window);
  sprintf(opt_v[14], "%d", theparams->ihypo);


}

void update_params(char **opt_v, params *theparams)
{
  strcpy(theparams->ifile, opt_v[1]);
  strcpy(theparams->zfile, opt_v[2]);
  strcpy(theparams->error, opt_v[3]);
  strcpy(theparams->map, opt_v[4]);
  strcpy(theparams->srfile, opt_v[5]);
  strcpy(theparams->eqtl, opt_v[6]);
  theparams->seed = atol(opt_v[7]);
  theparams->Model = atoi(opt_v[8]);
  theparams->whichtrait = atoi(opt_v[9]);
  theparams->wchrom = atoi(opt_v[10]);
  theparams->walk = (FPN) atof(opt_v[11]);
  theparams->nbp = atoi(opt_v[12]);
  theparams->window = (FPN) atof(opt_v[13]);
  if (theparams->nbp < 0 )
    theparams->nbp = (int) NUM_SIG;
  if ( theparams->window < 0.0)
    theparams->window = (FPN) WIN_SIZE;
  theparams->ihypo = atoi(opt_v[14]);

}



/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MZmain.c
------------------------------------------------------------------ */


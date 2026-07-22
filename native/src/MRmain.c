/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MRmain.c) is part of QTL Cartographer
         
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



/* Main function driver and subroutines for MultiRegress. */
#include "Main.h"

int main(argc, argv)
int argc;
char *argv[];
{
  char *chptr,*purpose,**tnames,**otnames;
  int nopts,automatic,error ,i,otraits,doit,offset;
  long *positions,*dpositions;
  FPN walk;
  params *theparams;
  linpakws *lnpk;
#ifndef WIN32
#if defined(MACWARRIOR)  
 /*  Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
#endif
  whichprogram = 73;
  /* just to be on the safe side, initialize pointers to null */
   theparams = NULL;
  /*   Create space for purpose */
  purpose = cvector(0,MAXNAME);
  strcpy(purpose, "Multiple Regression Module");
  /*  Create the parameters structure, get the time and process the command line args. */
  nopts = 10;
  theparams = create_params(theparams,1,NULL);
  theparams->ihypo = 10;
  chptr = asctime2(); 
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  /* Initialize the data structures...*/
  error = check_params(theparams,NULL,73); 
  theparams->cross = get_cross(theparams->mrinput,theparams);
  GetMRParams(theparams->mrinput, &walk,&otraits, &theparams->traits , &theparams->nn);
  positions = GetPositions(theparams->mrinput,1);
  if ( theparams->crosstype == 3 || theparams->crosstype == 4 )
    dpositions = GetPositions(theparams->mrinput,2);
  lnpk = cd_multiregressws(NULL,theparams->nn,theparams->srupper ,theparams->traits,1);	
  lnpk->ipos = (int) positions[0] ;
  lnpk->pcnts = ivector(1,lnpk->ipos);
  lnpk->samplesize =   ivector(1,lnpk->ipos);
  lnpk->lratio = dvector(1,lnpk->ipos);
  GetChromLocales(lnpk->samplesize,lnpk->lratio,positions,theparams->mrinput);
/*  put phenotypes in lnkp->y.  */
  tnames = cmatrix(1,theparams->traits,0,MAXNAME);
  
  lnpk->y = GetTraitData(theparams->mrinput,theparams->traits, theparams->nn,tnames);
  if ( otraits > 0 ) {
    otnames = cmatrix(1,otraits,0,MAXNAME);
    lnpk->k = otraits;
    lnpk->bp = GetOTraitData(theparams->mrinput,otraits, theparams->nn,otnames);
    offset = 0;
    for (i=1; i<=lnpk->k; i++ )
      offset = offset + lnpk->bp[i][0] - 1;
    lnpk->ldx = offset;
  }
  else 
    lnpk->ldx = 0;
/*  Now the analysis section... */
  sprintf(theparams->tfile,"%s.tmp",theparams->stem);
/*  insert_wd(gbuffer,theparams->workdir,theparams->tfile);*/
  print_head(argv[0],theparams->mroutput,chptr,1,21,theparams);

  ShowQTLSinp(NULL, theparams,theparams->mroutput, 0, tnames[1],1);
  print_head(argv[0],theparams->tfile,chptr,1,67,theparams);
  ShowMRParams(theparams,theparams->tfile,lnpk,otnames );
  for ( i=1; i<=theparams->traits; i++ ) {
    if ( theparams->whichtrait < 1 && tnames[i][0] == '+' )
      doit = 1;
    else if ( theparams->whichtrait > theparams->traits && tnames[i][0] != '-' )
      doit = 1;
	else if ( theparams->whichtrait == i )
	  doit = 1;
	else 
	  doit = 0;   
    if ( doit == 1 ) {
      if ( theparams->verbosity == 1 ) 
        printf("\n Trait  %d (%s) and Sample size %d       \n",i, tnames[i], theparams->nn);
      theparams->thegenome = multiregress(theparams->mrinput, theparams, lnpk, positions,dpositions,i);
      EstimateEffects(theparams->thegenome, lnpk,  theparams, i) ;
      if ( theparams->verbosity == 1 )
        ShowMultiRegress(theparams->thegenome,lnpk->pcnts,theparams,"-",i,tnames[i]);
      ShowMultiRegress(theparams->thegenome,lnpk->pcnts,theparams,theparams->tfile,i,tnames[i]);
      ShowQTLSinp(theparams->thegenome, theparams,theparams->mroutput, i, tnames[i],2);
      clear_genome(theparams->thegenome);
      theparams->thegenome = NULL;
    }
    else
      ShowQTLSinp(theparams->thegenome, theparams,theparams->mroutput, 0, tnames[i],2);
  }
  ShowQTLSinp(NULL, theparams,theparams->mroutput, 0, tnames[1],3);
  
  AppendAtoB(theparams->tfile,theparams->mroutput);
  remove(theparams->tfile) ;

/* Clean up...*/
  write_trailer(theparams,chptr,1);
  if ( purpose != NULL )
	 free_cvector(purpose,0,MAXNAME);   
  free_lvector(positions,0,positions[0]);
  if ( theparams->crosstype == 3 || theparams->crosstype == 4 )
    free_lvector(dpositions,0,dpositions[0]);
  free_cmatrix(tnames,1,theparams->traits,0,MAXNAME);
  if ( otraits > 0 ) 
    free_cmatrix(otnames,1,otraits,0,MAXNAME);
  lnpk = cd_multiregressws(lnpk,theparams->nn,theparams->srupper ,theparams->traits,0);
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
    strcpy(opt[4], "-s");
    strcpy(opt[5], "-t");
    strcpy(opt[6], "-c");
    strcpy(opt[7], "-S");
    strcpy(opt[8], "-u");
    strcpy(opt[9], "-w");
    strcpy(opt[10], "-I");


    strcpy(opt_e[1], "Input File");
    strcpy(opt_e[2], "Output File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Random Number Seed");
    strcpy(opt_e[5], "Trait to analyze");
    strcpy(opt_e[6], "Use Categorical Traits [1/0]");
    strcpy(opt_e[7], "Size of the test (alpha)");
    strcpy(opt_e[8], "Hard bound on number of parameters");
    strcpy(opt_e[9], "Window Size in cM ");
    strcpy(opt_e[10], "Hypothesis test [10/30]");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->mrinput);   
  strcpy(opt_v[2], theparams->mroutput);
  strcpy(opt_v[3], theparams->error);
  sprintf(opt_v[4], "%ld", theparams->seed);
  sprintf(opt_v[5], "%d", theparams->whichtrait);
  sprintf(opt_v[6], "%d", theparams->boots);  /* Default is 0, do not use Cat. traits */
  sprintf(opt_v[7], "%f", theparams->size);
  sprintf(opt_v[8], "%d", theparams->srupper);
  sprintf(opt_v[9], "%f", theparams->window);
  sprintf(opt_v[10], "%d", theparams->ihypo);


}

void update_params(char **opt_v,  params *theparams)
{
  strcpy(theparams->mrinput, opt_v[1]);
  strcpy(theparams->mroutput, opt_v[2]);
  strcpy(theparams->error, opt_v[3]);
  theparams->seed = atol(opt_v[4]);
  theparams->whichtrait = atoi(opt_v[5]);
  theparams->boots = atoi(opt_v[6]);
  theparams->size = (FPN) atof(opt_v[7]);
  theparams->srupper = atoi(opt_v[8]);
  theparams->window = (FPN) atof(opt_v[9]);
  theparams->ihypo = atoi(opt_v[10]);

}



/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MRmain.c
------------------------------------------------------------------ */


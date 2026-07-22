/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MImain.c) is part of QTL Cartographer
         
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


/* Main function driver and subroutines for MImapqtl.*/
#include "Main.h"

int main(argc, argv)
int argc;
char *argv[];
{
  FILE *outf;
  char *chptr,*purpose;
  int nopts,automatic,ii,wt,nn,t;
  params *theparams;
  linpakws *lnpk;
#if defined(MACWARRIOR) 
 /*  Simulate the UNIX shell command line with ccommand. */
  argc = ccommand(&argv);
#endif
  whichprogram = 10;
  /* just to be on the safe side, initialize pointers to null */
  theparams = NULL;
  /*   Create space for purpose */
  purpose = cvector(0,MAXNAME);
  strcpy(purpose, "Multiple Interval Mapping");
  /*  Create theparameters, get the time, and process the command line args. */
  nopts = 15;
  theparams = create_params(theparams,1,NULL);
  chptr = asctime2(); 
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
/* Initialize the data structures...
         Start with the map of makers */
  GetTheMap(theparams, theparams->map);

/*       Get the estimates of QTL positions from the Eqtl.out file */
  theparams->themap->traits = theparams->traits;
  if ( theparams->mimwork[1] == 'M' )
    GetTheModel(theparams, theparams->mqtfile);
  else
    theparams->theqtls = NULL;
/*      Now get the data from the Rcross.out file*/
  GetTheData(theparams, theparams->ifile, 1);

  nn = theparams->nn;

/*	Create a linked list representing the genome and
	translate the distances in recombination frequencies to those
	of Morgans...  */
  theparams->thegenome = create_genome(theparams->themap);
/* How many traits in a multitrait analysis? */
  t = how_many_traits(theparams,theparams->themap);
/* At this point we should be able to figure out an upper bound to how much workspace we need. */
  lnpk = cd_mimws(NULL,theparams->nn, (int) MAXGT ,  theparams->traits ,1);	  
/*  Write a header for the output file.  */
  write_mimheader(theparams->mimfile, theparams, argv[0], chptr,theparams->themap,1);
/*  Write a header for the MImapqtl model file that looks like Rqtl.out */
  if ( (ii=(int) strlen(theparams->workdir)) > 0   )
       insert_wd(gname,theparams->workdir, theparams->tfile);
  print_head(argv[0],theparams->tfile,chptr,0,20,theparams);
  outf = fileopen(theparams->tfile, "a");  
  PrintRqtlout(outf,theparams,theparams->themap,NULL,1);
  fileclose(theparams->tfile, outf);

/*  Do the analysis for the requested traits.  */
  wt = theparams->whichtrait;
  for (ii = 1; ii <= theparams->traits; ii++) 
      if ( (wt == 0 && theparams->themap->tnames[ii][0] == '+') || (wt > theparams->traits && theparams->themap->tnames[ii][0] != '-') || (wt == ii) ) {
        theparams->whichtrait = ii;
        theparams->nn = MoveMissPhenotypes(nn,theparams->whichtrait,theparams->thedata);                                      
        copy_traits(theparams, theparams->themap, theparams->thedata, lnpk);  /*  put phenotypes in lnkp->y.  */
        do_mimanalysis(theparams,theparams->themap,theparams->theqtls,theparams->thedata, theparams->thegenome,lnpk);
        theparams->nn = nn;
      }
  theparams->whichtrait = wt;
  outf = fileopen(theparams->tfile, "a");  
  PrintRqtlout(outf,theparams,theparams->themap,NULL,3);
  fileclose(theparams->tfile, outf);
  
  if ( theparams->mimwork[7] == 'R' ) {  /* Write a new data file with residuals replacing traits.  */
/* Translate marker data from -1,0,1  to  0, (1,2), 3 */
    for (ii = 1; ii <= theparams->nn; ii++)
      untrans_data( theparams->thedata[ii].markers, theparams->themap);

    if ( theparams->mimphase > 0 )    
      sprintf(theparams->tfile, "%sPhase%d.res",theparams->stem, theparams->mimphase ); 
    else
      sprintf(theparams->tfile, "%s.res",theparams->stem ); 
    if ( (ii=(int) strlen(theparams->workdir)) > 0   )
       insert_wd(gname,theparams->workdir, theparams->tfile);
    print_head(argv[0],theparams->tfile,chptr,0,30,theparams);
    print_individuals(theparams,theparams->thedata, theparams->nn,  theparams->tfile);
    if ( (ii=(int) strlen(theparams->workdir)) > 0   )
      shift_fn(theparams->tfile);

  }
/* Clean up...*/
  if ( theparams->mimphase > 0 )
    theparams->mimphase +=1;
  write_trailer(theparams,chptr,1);
  if ( purpose != NULL )
	 free_cvector(purpose,0,MAXNAME);
  lnpk = cd_mimws(lnpk,theparams->nn, (int) MAXGT ,t,0);
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
    strcpy(opt[5], "-E");
    strcpy(opt[6], "-O");
    strcpy(opt[7], "-s");
    strcpy(opt[8], "-t");
    strcpy(opt[9], "-q");
    strcpy(opt[10], "-k");
    strcpy(opt[11], "-d");
    strcpy(opt[12], "-S");
    strcpy(opt[13], "-L");
    strcpy(opt[14], "-I");
    strcpy(opt[15], "-p");


    strcpy(opt_e[1], "Data File");
    strcpy(opt_e[2], "Output Analysis File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Genetic Linkage Map File");
    strcpy(opt_e[5], "Input model file ");
    strcpy(opt_e[6], "Output model file ");
    strcpy(opt_e[7], "Random Number Seed");
    strcpy(opt_e[8], "Trait to analyze");
    strcpy(opt_e[9], "Maximum number of QTL to fit");
    strcpy(opt_e[10], "Maximum number of epistatic interactions");
    strcpy(opt_e[11], "Walking speed in cM");
    strcpy(opt_e[12], "Information criterion [1-6]");
    strcpy(opt_e[13], "IC threshold for adding/dropping parameters");
    strcpy(opt_e[14], "Work code ");
    strcpy(opt_e[15], "Phase of analysis ");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->ifile);
  strcpy(opt_v[3], theparams->error);
  strcpy(opt_v[4], theparams->map);
  sprintf(opt_v[7], "%ld", theparams->seed);
  sprintf(opt_v[8], "%d", theparams->whichtrait);
  sprintf(opt_v[9], "%d", theparams->maxqtl);
  sprintf(opt_v[10], "%d", theparams->maxepistatics);
  sprintf(opt_v[11], "%f", theparams->walk);
  sprintf(opt_v[12], "%d", theparams->whoseic);
  sprintf(opt_v[13], "%f", theparams->mimlod);
  sprintf(opt_v[14], "%s", theparams->mimwork);
  sprintf(opt_v[15], "%d", theparams->mimphase);

  if ( theparams->mimphase > 0 ) {
    sprintf(opt_v[2], "%sPhase%d.mim",theparams->stem, theparams->mimphase );
    sprintf(opt_v[5], "%sPhase%d.mqt",theparams->stem, theparams->mimphase-1 );
    sprintf(opt_v[6], "%sPhase%d.mqt",theparams->stem, theparams->mimphase );  
  }
  else {
    strcpy(opt_v[2], theparams->mimfile);  
    strcpy(opt_v[5], theparams->mqtfile);
    strcpy(opt_v[6], theparams->tfile);  
  }

}

void update_params(char **opt_v,  params *theparams)
{
  strcpy(theparams->ifile, opt_v[1]);
  strcpy(theparams->error, opt_v[3]);
  strcpy(theparams->map, opt_v[4]);
  theparams->seed = atol(opt_v[7]);
  theparams->whichtrait = atoi(opt_v[8]);
  theparams->maxqtl = atoi(opt_v[9]);
  theparams->maxepistatics =   atoi(opt_v[10]);
  theparams->walk = (FPN) atof(opt_v[11]);
  theparams->whoseic = atoi(opt_v[12]);
  if ( theparams->whoseic > 6 || theparams->whoseic < 1 ) {
    theparams->whoseic = 1;
    sprintf(opt_v[12], "%d", theparams->whoseic);
  }
  theparams->mimlod = (FPN) atof(opt_v[13]);
  strcpy(theparams->mimwork,opt_v[14]);
  theparams->mimphase = (int) atoi(opt_v[15]);


  if ( theparams->mimphase > 0 ) {
    sprintf(theparams->mimfile, "%sPhase%d.mim",theparams->stem, theparams->mimphase );
    sprintf(theparams->mqtfile, "%sPhase%d.mqt",theparams->stem, theparams->mimphase-1 );
    sprintf(theparams->tfile, "%sPhase%d.mqt",theparams->stem, theparams->mimphase );  
  }
  else {
    strcpy(theparams->mimfile, opt_v[2]);
    strcpy(theparams->mqtfile, opt_v[5]);
    strcpy(theparams->tfile,  opt_v[6]);   /* tfile pointer for filename  */    
  }
  
}



/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MImain.c
------------------------------------------------------------------ */


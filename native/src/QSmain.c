/* ------------------------------------------------------ XCutXCodeXSkip
     This file (QSmain.c) is part of QTL Cartographer
         
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
/* Main function driver for Qstats. */

int main(argc, argv)
int argc;
char *argv[];
{
  FILE  *outf;
  char  *chptr,*purpose;
  int nopts;
  params *theparams;
  int ii,automatic;
  FPN  *trtptr,*thestats;
  int traits,nind;

#if defined(MACWARRIOR)  
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  whichprogram = 5;
  purpose = cvector(0,MAXNAME);
  strcpy(purpose,"Calculate some Basic Statistics for the QTL");
  theparams = NULL;
  nopts = 6;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  print_head(argv[0],theparams->qstat,chptr,1,40,theparams);

/* Initailize the data structures...*/
  theparams->theqtls = NULL;
  GetTheMap(theparams, theparams->map);
  GetTheData(theparams, theparams->ifile, 0 );
  theparams->thegenome = create_genome(theparams->themap);
  trtptr = dvector(1,theparams->nn);
  thestats = dvector(0,15);
  for ( traits = 1;  traits <= theparams->themap->traits; traits++) {
    if ( theparams->verbosity == 1 ) {
      if ( theparams->themap->tnames != NULL )
        printf("\nNow processing trait %s",theparams->themap->tnames[traits]);
      else 
        printf("\nNow processing trait %d",traits);
    }
    nind = move_phenotypes(trtptr,theparams->thedata,traits,theparams->nn);
    calc_qstats(trtptr,nind,thestats,15);
	print_qstats(thestats, theparams,traits,theparams->themap);
	print_histogram(trtptr, nind, 50, theparams->qstat);
	miss_mark_summary(theparams->thedata,theparams->themap,theparams,traits);
/*	otrait_regression(theparams->thedata,theparams->themap,theparams,traits);*/
    if (traits < theparams->themap->traits) {
      outf = fileopen(theparams->qstat, "a");
      fprintf(outf, "\f");	/* Put in a form feed after each trait... */
      fileclose(theparams->qstat, outf);
    }
  }
/* Summarize missing data for each individual... */
  ind_data_summary(theparams->thedata,theparams->themap,theparams);		
/* Translate from 0, (1,2), 3 to -1,0,1 */
  for (ii = 1; ii <= theparams->nn; ii++)
    trans_data(theparams->thedata[ii].markers, theparams->themap);

  marker_segregation(theparams->thedata,theparams->themap,theparams);	 

  if ( theparams->mimwork[0] == 'y' )
    CalcMarkProbs(theparams);
  write_trailer(theparams,chptr,1);
/* Clean up...*/
  free_dvector(trtptr, 1, theparams->nn);
  free_dvector(thestats, 0, 15);
  free_cvector( purpose,0,MAXNAME);
  theparams = create_params(theparams, 0, NULL);

  return(0);
}

void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag)
{
  int ii,jj;
  if ( flag == 1 ) {
    strcpy(opt[1],  "-i"  );
    strcpy(opt[2],  "-o"  );
    strcpy(opt[3],  "-e"  );
    strcpy(opt[4],  "-m"  );
    strcpy(opt[5],  "-s"  );
    strcpy(opt[6],  "-p"  );

    strcpy(opt_e[1],  "Data Input File"  );
    strcpy(opt_e[2],  "Output File"  );
    strcpy(opt_e[3],  "Error File"  );
    strcpy(opt_e[4],  "Genetic Linkage Map File"  );
    strcpy(opt_e[5],  "Random Number Seed"  );
    strcpy(opt_e[6],  "Calc. Marker Probabilities"  );
  }
  for ( ii = 1 ; ii <= nopts ; ii++ ) 
    for ( jj = 0 ; jj <= MAXNAME ; jj++ )
      opt_v[ii][jj] = '\0';

    strcpy(opt_v[1],  theparams->ifile  );
    strcpy(opt_v[2],  theparams->qstat  );
    strcpy(opt_v[3],  theparams->error  );
    strcpy(opt_v[4],  theparams->map  );
    sprintf(opt_v[5],"%ld",theparams->seed );
    sprintf(opt_v[6], "%s", theparams->mimwork);

}  

void update_params(char **opt_v,  params *theparams)
{

    strcpy(theparams->ifile, opt_v[1] );
    strcpy(theparams->qstat, opt_v[2] );
    strcpy(theparams->error, opt_v[3]  );
    strcpy(theparams->map, opt_v[4]  );
    theparams->seed = atol(opt_v[5]);
    strcpy(theparams->mimwork,opt_v[6]);

}  






/* ------------------------------------------------------- XCutXCodeXSkip
             End of file QSmain.c
------------------------------------------------------------------ */


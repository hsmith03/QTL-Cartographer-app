/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Emain.c) is part of QTL Cartographer
         
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
/*  Main function and driver for Emap */

int main(int argc, char *argv[])  {
  char *chptr, *purpose;
  int nopts  ;
  params *theparams;
  int   automatic, skipstage1;
  div_t emethod;


#if defined(MACWARRIOR)  
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  skipstage1 = 0;
  whichprogram = 12;
  purpose = cvector(0,MAXNAME);
  strcpy(purpose, "Estimate a genetic linkage map");
  theparams = NULL;
  automatic = 0;
  nopts = 13;
  theparams = create_params(theparams, 1, NULL);
  theparams->gout = 0;
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);

  emethod = div(theparams->emethod, 10 );

  print_head(argv[0],theparams->map,chptr,0,10,theparams);

/* Initailize the data structures...*/
  if ( isfile(theparams->mapin)  ) 
  	GetTheMap(theparams, theparams->mapin );
  else
    theparams->themap = NULL;
  if ( theparams->themap != NULL ) 
  	skipstage1 = 1;

/* Get the data  */
  if (theparams->iinfile[0] != '\0')  
    GetTheData( theparams, theparams->iinfile,1) ;
  else  /*  If no data, then quit */
    exit(1);
    
  if ( skipstage1 == 0 ) /* Do stage 1.   */
      RCDstage1( theparams, chptr, argv[0]);
    

 if ( skipstage1 == 0 ) {  /* reload the results of stage 1 analysis*/
   strcpy(theparams->tfile,theparams->stem);
   strcat(theparams->tfile,".e1m");
   GetTheMap(theparams, theparams->tfile );

	strcpy(theparams->tfile,theparams->stem);
    strcat(theparams->tfile,".e1c");
    GetTheData( theparams, theparams->tfile,1) ;
  }  
  if ( emethod.quot > 0 && emethod.rem > 0 )
    RCDstage2( theparams,    emethod );
  else
    ReEstimateMap( theparams );

      
/*  print out new map and data file in Rmap.out and Rcross.out formats */
  print_head(argv[0],theparams->map,chptr,0,10,theparams);
  print_map(theparams->themap, theparams->map);
  print_head(argv[0],theparams->ifile,chptr,0,30,theparams);
  print_individuals(theparams,theparams->thedata, theparams->nn, theparams->ifile);
               
/*  print_genome_map( theparams->thegenome,theparams,"-",argv[0], chptr);*/
  write_trailer(theparams,chptr,1);


/* Clean up...*/


  free_cvector(purpose,0,MAXNAME);
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
    strcpy(opt[5], "-l");
    strcpy(opt[6], "-s");
    strcpy(opt[7], "-M");
    strcpy(opt[8], "-S");
    strcpy(opt[9], "-L");
    strcpy(opt[10], "-r");
    strcpy(opt[11], "-f");
    strcpy(opt[12], "-p");
    strcpy(opt[13], "-O");


    strcpy(opt_e[1], "Data Input File");
    strcpy(opt_e[2], "Data Output File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Map Input File");
    strcpy(opt_e[5], "Map Output File");
    strcpy(opt_e[6], "Random Number Seed");
    strcpy(opt_e[7], "Linkage Map method");
    strcpy(opt_e[8], "Segregation test size");
    strcpy(opt_e[9], "Linkage test size");
    strcpy(opt_e[10], "Permutations: not active");
    strcpy(opt_e[11],  "Map Function [1,8], 1 => Haldane"  );
    strcpy(opt_e[12],  "Map Function Parameter"  );
    strcpy(opt_e[13],  "Objective Functions (0=>SAL, 1=>SAR)"  );
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->iinfile);
  strcpy(opt_v[2], theparams->ifile);
  strcpy(opt_v[3], theparams->error);
  strcpy(opt_v[4], theparams->mapin);
  strcpy(opt_v[5], theparams->map);
  sprintf(opt_v[6], "%ld", theparams->seed);
  sprintf(opt_v[7], "%d", theparams->emethod);
  sprintf(opt_v[8], "%f", theparams->segsize);
  sprintf(opt_v[9], "%f", theparams->linksize);
  sprintf(opt_v[10], "%d", theparams->reps);
    sprintf(opt_v[11],"%d",theparams->mapfunc );
    sprintf(opt_v[12],"%f",theparams->mapparam );
    sprintf(opt_v[13],"%d",theparams->emapobj );


}

void update_params(char **opt_v,  params *theparams)
{

  strcpy(theparams->iinfile, opt_v[1]);
  strcpy(theparams->ifile, opt_v[2]);
  strcpy(theparams->error, opt_v[3]);
  strcpy(theparams->mapin, opt_v[4]);
  strcpy(theparams->map, opt_v[5]);
  theparams->seed = atol(opt_v[6]);
  theparams->emethod = atoi(opt_v[7]);
  theparams->segsize = (FPN) atof(opt_v[8]);
  theparams->linksize = (FPN) atof(opt_v[9]);
  theparams->reps = atoi(opt_v[10]);
  whosemf = theparams->mapfunc = atoi( opt_v[11] );
  mapparam = theparams->mapparam = (FPN) atof( opt_v[12] );
  theparams->emapobj = atoi(opt_v[13]);
}




/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Emain.c
------------------------------------------------------------------ */


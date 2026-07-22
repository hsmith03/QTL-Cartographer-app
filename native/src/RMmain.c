/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RMmain.c) is part of QTL Cartographer
         
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
/* Main function driver for Rmap. */

int main(int argc, char *argv[])  {
  
  char  *chptr, *purpose;
  int nopts ;
  params *theparams;
  int automatic ;

#if defined(MACWARRIOR)  
/* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and output. */
  argc = ccommand(&argv);
#endif
  whichprogram = 1;
/**/
  purpose = cvector(0,MAXNAME);
  strcpy(purpose,"Create or translate a genetic linkage map");
  theparams = NULL;
  automatic = 0;
  nopts = 14;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  print_head(argv[0],theparams->map,chptr,0,10,theparams);
  mapparam = theparams->mapparam;
  if (theparams->mapin[0] == '\0') {
    theparams->themap = initailize_markermap(0, theparams);
    create_bogus_markers(theparams, theparams->themap );
  }
  else 
    GetTheMap( theparams, theparams->mapin);
  if (theparams->themap == NULL)
    exit(1);
  if (theparams->gout == 1 || theparams->gout == 3)
    print_map(theparams->themap, theparams->map);
  if (theparams->gout == 2 || theparams->gout == 3)
    gnuplot_map(theparams,theparams->themap);
  if (theparams->gout == 4) 
    print_map_inp(theparams->themap, theparams->map);
  theparams->Rmode = 0;
  write_trailer(theparams,chptr,1);
  free_cvector(purpose,0,MAXNAME);
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
    strcpy(opt[4],  "-s"  );
    strcpy(opt[5],  "-f"  );
    strcpy(opt[6],  "-p"  );
    strcpy(opt[7],  "-g"  );
    strcpy(opt[8],  "-c"  );
    strcpy(opt[9],  "-m"  );
    strcpy(opt[10],  "-vm" );
    strcpy(opt[11], "-d"  );
    strcpy(opt[12], "-vd" );
    strcpy(opt[13], "-t"  );
    strcpy(opt[14], "-M"  );

    strcpy(opt_e[1],  "Input File"  );
    strcpy(opt_e[2],  "Output File"  );
    strcpy(opt_e[3],  "Error File"  );
    strcpy(opt_e[4],  "Random Number Seed"  );
    strcpy(opt_e[5],  "Map Function [1,8], 1 => Haldane"  );
    strcpy(opt_e[6],  "Map Function Parameter"  );
    strcpy(opt_e[7],  "Ouput (1,2,3) => Text, Graphics, Both"  );
    strcpy(opt_e[8],  "Chromosomes"  );
    strcpy(opt_e[9],  "Markers per Chromosome"  );
    strcpy(opt_e[10], "Standard Deviation of Markers per Chromosome" );
    strcpy(opt_e[11], "Intermarker Distance (cM)"  );
    strcpy(opt_e[12], "Standard Deviation of Intermarker Distance (cM)" );
    strcpy(opt_e[13], "Tails (Flanking DNA, in cM)"  );
    strcpy(opt_e[14], "Simulation Mode (0,1)"  );
  }
  for ( ii = 1 ; ii <= nopts ; ii++ ) 
    for ( jj = 0 ; jj <= MAXNAME ; jj++ )
      *(*(opt_v+ii)+jj) = '\0';

    strcpy(opt_v[1],  theparams->mapin  );
    strcpy(opt_v[2],  theparams->map  );
    strcpy(opt_v[3],  theparams->error  );
    sprintf(opt_v[4],"%ld",theparams->seed );
    sprintf(opt_v[5],"%d",theparams->mapfunc );
    sprintf(opt_v[6],"%f",theparams->mapparam );
    sprintf(opt_v[7],"%d",theparams->gout );
    sprintf(opt_v[8],"%d",theparams->chrom );
    sprintf(opt_v[9],"%d",theparams->mark );
    sprintf(opt_v[10],"%f",theparams->vmark );
    sprintf(opt_v[11],"%f",theparams->dist );
    sprintf(opt_v[12],"%f",theparams->vdist );
    sprintf(opt_v[13],"%f",theparams->tail );
    sprintf(opt_v[14],"%d",theparams->Rmode );

}  

void update_params(char **opt_v,  params *theparams)
{
    strcpy(theparams->mapin, opt_v[1] );
    strcpy(theparams->map, opt_v[2] );
    strcpy(theparams->error, opt_v[3]  );
    theparams->seed = atol(opt_v[4]);
    theparams->mapfunc = atoi( opt_v[5] );
    whosemf = theparams->mapfunc;
    theparams->mapparam = (FPN) atof( opt_v[6] );
    mapparam = theparams->mapparam;
    theparams->gout = atoi( opt_v[7] );
    theparams->chrom = atoi( opt_v[8] );
    theparams->mark = atoi( opt_v[9] );
    theparams->vmark = (FPN) atof( opt_v[10] );
    theparams->dist = (FPN) atof( opt_v[11] );
    theparams->vdist = (FPN) atof( opt_v[12] );
    theparams->tail = (FPN) atof( opt_v[13] );
    theparams->Rmode = atoi( opt_v[14] );
}  


 

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RMmain.c
------------------------------------------------------------------ */


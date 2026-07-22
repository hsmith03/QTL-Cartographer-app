/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RQmain.c) is part of QTL Cartographer
         
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
/* Main function driver for Rqtl.*/


int main(argc, argv)
int argc;
char *argv[];
{
  char *chptr,*purpose;
  int nopts,init,r,s;
  params *theparams;
  int automatic;
  int ii, jj,   ttraits,  isinfile, error;
  FPN **episADDA,**episAADD,coinflip;
  genome *first;

#if defined(MACWARRIOR)  
 /*Simulate the UNIX shell with ccommand */
  argc = ccommand(&argv);
#endif
  whichprogram = 2;
  purpose = cvector(0,MAXNAME);
  strcpy(purpose,"Create or translate a genetic model");
  theparams = NULL;
  nopts = 13;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();

  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  if ( theparams->qtlin[0] != '\0' )
  	isinfile = 1;
  else
    isinfile = 0;
   print_head(argv[0],theparams->qtl,chptr,0,20,theparams);

/* Need a map,   if this fails the program will quit. */
  GetTheMap(theparams, theparams->map );

  if (isinfile == 0) { /* Simulate  a model */
    if (theparams->traits > 1) {
      theparams->themap->traits = theparams->traits;
      free_ivector(theparams->themap->knum, 1, theparams->qnum);
      theparams->themap->knum = ivector(1, theparams->traits);
      for (ii = 1; ii <= theparams->traits; ii++) {
	    theparams->themap->knum[ii] = (int) ((FPN) theparams->qnum * ranf(ii)) + 1;
	    if (theparams->themap->knum[ii] < 1)
	      theparams->themap->knum[ii] = 1;
      }
    }
    error = check_params(theparams,theparams->themap,whichprogram);

    theparams->theqtls = qtlvector(theparams->themap);


/* Create a set of ttraits QTLs. */
    ttraits = 0;
    for (ii = 1; ii <= theparams->traits; ii++) {
      first = create_genome(theparams->themap);
      init = ttraits;  
      if ( theparams->themap->knum[ii] > 0 ) {
        episAADD = dmatrix(1,theparams->themap->knum[ii],1,theparams->themap->knum[ii]);
	    episADDA = dmatrix(1,theparams->themap->knum[ii],1,theparams->themap->knum[ii]);
      }
      else {
        episAADD = NULL;
	    episADDA = NULL;
      }
      for (jj = 1; jj <= theparams->themap->knum[ii]; jj++) { /* Additive and dominance effects */
	    ttraits = ttraits + 1;
	    first = create_a_qtl(first, (theparams->theqtls + ttraits), theparams);
	    theparams->theqtls[ttraits].trait = ii;
	    theparams->theqtls[ttraits].episAADD = episAADD;
	    theparams->theqtls[ttraits].episADDA = episADDA;
      }
      for ( r=1; r<theparams->themap->knum[ii] ; r++ ) 
        for ( s=r+1; s<=theparams->themap->knum[ii]; s++ ) {/* Do a random set of the combinations*/
          coinflip = ranf((long) ii);
          if ( coinflip < theparams->null_sse ) 
	          theparams->theqtls[ttraits].episAADD[s][r] = qtl_dom_effects(theparams->beta1, theparams->beta2, 4);
          coinflip = ranf((long) ii);
	      if (  coinflip < theparams->null_sse && theparams->theqtls[init+s].d != 0.0 )
	            theparams->theqtls[ttraits].episADDA[s][r] = qtl_dom_effects(theparams->beta1, theparams->beta2, 4);
          coinflip = ranf((long) ii);
	      if ( coinflip < theparams->null_sse && theparams->theqtls[init+r].d != 0.0 )
	            theparams->theqtls[ttraits].episADDA[r][s] = qtl_dom_effects(theparams->beta1, theparams->beta2, 4);
          coinflip = ranf((long) ii);
	      if ( coinflip < theparams->null_sse &&  theparams->theqtls[init+s].d != 0.0 &&  theparams->theqtls[init+r].d != 0.0 )
	            theparams->theqtls[ttraits].episAADD[r][s] = qtl_dom_effects(theparams->beta1, theparams->beta2, 4);
        }
      
      clear_genome(first);
    }
  }
  else  
    GetTheModel(theparams, theparams->qtlin );

  print_aqtl(theparams, theparams->qtl, 0);
  theparams->Rmode = 0;  /*  convention:  Rmode is always 0 */
  write_trailer(theparams,chptr,1);

/* Clean up...*/
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
    strcpy(opt[4],  "-m"  );
    strcpy(opt[5],  "-s"  );
    strcpy(opt[6],  "-t"  );
    strcpy(opt[7],  "-q"  );
    strcpy(opt[8],  "-b"  );
    strcpy(opt[9],  "-1"  );
    strcpy(opt[10],  "-2"  );
    strcpy(opt[11],  "-d" );
    strcpy(opt[12],  "-E" );
    strcpy(opt[13], "-M"  );

    strcpy(opt_e[1],  "Input File"  );
    strcpy(opt_e[2],  "Output File"  );
    strcpy(opt_e[3],  "Error File"  );
    strcpy(opt_e[4],  "Genetic Linkage Map File"  );
    strcpy(opt_e[5],  "Random Number Seed"  );
    strcpy(opt_e[6],  "Number of Traits"  );
    strcpy(opt_e[7],  "Number of QTL per Trait"  );
    strcpy(opt_e[8],  "Additive effect parameter beta"  );
    strcpy(opt_e[9],  "Dominance effect parameter beta1"  );
    strcpy(opt_e[10],  "Dominance effect parameter beta2"  );
    strcpy(opt_e[11],  "Dominance (1,2,3,4) => None,A,a,Random"  );
    strcpy(opt_e[12],  "Proportion of non-zero epistatic terms"  );
    strcpy(opt_e[13], "Simulation Mode (0,1)"  );
  }
  for ( ii = 1 ; ii <= nopts ; ii++ ) 
    for ( jj = 0 ; jj <= MAXNAME ; jj++ )
      *(*(opt_v+ii)+jj) = '\0';

    strcpy(opt_v[1],  theparams->qtlin  );
    strcpy(opt_v[2],  theparams->qtl  );
    strcpy(opt_v[3],  theparams->error  );
    strcpy(opt_v[4],  theparams->map  );
    sprintf(opt_v[5],"%ld",theparams->seed );
    sprintf(opt_v[6],"%d",theparams->traits );
    sprintf(opt_v[7],"%d",theparams->qnum );
    sprintf(opt_v[8],"%f",theparams->beta );
    sprintf(opt_v[9],"%f",theparams->beta1 );
    sprintf(opt_v[10],"%f",theparams->beta2 );
    sprintf(opt_v[11],"%d",theparams->dom );
    sprintf(opt_v[12],"%f",theparams->null_sse );
    sprintf(opt_v[13],"%d",theparams->Rmode );

}  

void update_params(char **opt_v,  params *theparams)
{

    strcpy(theparams->qtlin, opt_v[1] );
    strcpy(theparams->qtl, opt_v[2] );
    strcpy(theparams->error, opt_v[3]  );
    strcpy(theparams->map, opt_v[4]  );
    theparams->seed = atol(opt_v[5]);
    theparams->traits = atoi( opt_v[6] );
    theparams->qnum = atoi( opt_v[7] );
    theparams->beta = (FPN) atof( opt_v[8] );
    theparams->beta1 = (FPN) atof( opt_v[9] );
    theparams->beta2 = (FPN) atof( opt_v[10] );
    theparams->dom = atoi( opt_v[11] );
    theparams->null_sse = (FPN) atof( opt_v[12] );
    if ( theparams->null_sse < (FPN) 0.0 || theparams->null_sse > (FPN) 1.0 )
      theparams->null_sse = (FPN) 0.0;
    theparams->Rmode = atoi( opt_v[13] );
}  






/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RQmain.c
------------------------------------------------------------------ */


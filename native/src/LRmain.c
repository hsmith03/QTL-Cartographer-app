/* ------------------------------------------------------ XCutXCodeXSkip
     This file (LRmain.c) is part of QTL Cartographer
         
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
/* Main  function driver for LRmapqtl */

int main(argc, argv)
int argc;
char *argv[];
{
  FILE *errorf, *outf;
  char  *chptr,*purpose;
  params *theparams;
  int ii,kk,automatic,jj,nopts,rows,itrait,utrait;
  int error, traits;
  FPN *origtraits,*maxlr,**savefstat,pval;
  int *indx,**pcnts;
  FPN   **beta0, **beta1, **fstat,**lrvalues;
  linpakws *lnpk;

#if defined(MACWARRIOR)
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  whichprogram = 6;
  purpose = cvector(0,MAXNAME);
  strcpy(purpose,"Map QTL using Simple Linear Regression");
  theparams = NULL;
  nopts = 7;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  print_head(argv[0],theparams->lrfile,chptr,1,50,theparams);

/* Initailize the data structures...*/
  GetTheMap(theparams, theparams->map);
  GetTheData(theparams, theparams->ifile, 1 );

  if ( theparams->themap->otraits > 0 )
    process_otraits(theparams->thedata,theparams->nn,theparams->themap);


	if (theparams->perms > 0 && theparams->perms < MAX_REPS) {
	  indx = ivector(1, theparams->nn);
	  savefstat = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
	  origtraits = dvector(1, theparams->nn);
	  pcnts = imatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
	  maxlr = dvector(1, theparams->perms);
	}
	else {
	  maxlr = origtraits = NULL;
	  pcnts = NULL;
	}

    beta0 = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
    beta1 = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
    fstat = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
	if ( theparams->themap->otraits > 0 ) {
	  
      rows = 2;
      for ( ii = 1 ; ii <= theparams->themap->otraits ; ii++ )
        if ( theparams->themap->onames[ii][0] == '+' )
          rows = rows + 2 * theparams->themap->otypes[ii] - 2;
	  if ( rows > 2 ) {
	    lrvalues = dmatrix(1,theparams->themap->m, 1, theparams->themap->maxl);
        lnpk = cd_linpakws(NULL,theparams->nn,rows,1);	
      }
      else
        rows = 1;  
	}
	else
	  rows = 1;  
  if ( theparams->whichtrait <= theparams->themap->traits && theparams->whichtrait > 0 )
    itrait = utrait = theparams->whichtrait;
  else {
    itrait = 1;
    utrait = theparams->themap->traits;
  }
  for (traits = itrait; traits <= utrait; traits++) {
    if ( theparams->verbosity == 1 ) {
      if ( theparams->themap->tnames != NULL )
        printf("\nNow processing trait %s",theparams->themap->tnames[traits]);
      else 
        printf("\nNow processing trait %d",traits);
    }
	outf = fileopen(theparams->lrfile, "a");
    if (traits > 1  && itrait != utrait)
      fprintf(outf, "\n\f");	/* Put in a form feed before each trait... */
    fileclose(theparams->lrfile, outf);
    theparams->thegenome = create_genome( theparams->themap );

/* clean out matrices */
    for (ii = 1; ii <= theparams->themap->m; ii++)
      for (jj = 1; jj <= theparams->themap->maxl; jj++)
	    beta0[ii][jj] = beta1[ii][jj] = fstat[ii][jj] = (FPN) 0.0;
	if ( rows > 2 )
      error = calc_mlrstats(theparams,theparams->themap, theparams->thedata,  traits, beta0, beta1, fstat,lrvalues, theparams->thegenome,lnpk);
    else 
      error = calc_lrstats(theparams,theparams->themap, theparams->thedata,  traits, beta0, beta1, fstat, theparams->thegenome);
    if (error != 0) {/* If there is some errror, */
      errorf = fileopen(theparams->error, "a");
      fprintf(errorf, "\n\nHad trouble doing the linear regression of data from %s...\n", theparams->ifile);
		fprintf(errorf, "errorlevel = %4d\n", error);
      fileclose(theparams->error, errorf);
    }
	if ( rows > 2 )
	 	print_mlrstats(theparams,traits,theparams->themap,beta0,beta1,fstat,lrvalues);
    else
	 	print_lrstats(theparams->nn, traits, theparams->themap, theparams->map, theparams->ifile, theparams->lrfile, beta0, beta1, fstat, theparams->cross);
    if (theparams->perms > 0 && theparams->perms < MAX_REPS) { /* do a shuffle test...*/
      for (ii = 1; ii <= theparams->nn; ii++)  /* indx[ii] = ii */
	    indx[ii] = ii;
      for (ii = 1; ii <= theparams->themap->m; ii++)
        for (jj = 1; jj <= theparams->themap->maxl; jj++) {  /* zero out pcnts, save f stats */
	      pcnts[ii][jj] = 0;
          savefstat[ii][jj] = fstat[ii][jj];
        }
      for ( jj = 1 ; jj <= theparams->perms ; jj++ )
        maxlr[jj] = (FPN)  0.0;
      for ( jj = 1 ; jj<= theparams->nn ; jj++ )  /* save original trait values */
        origtraits[jj] = theparams->thedata[jj].y[traits];
      for (kk = 1; kk <= theparams->perms; kk++) {  /* do permutations */
		if ( theparams->verbosity == 1 )
			 printf("\nNow on Permutation %d",kk);
        for (ii = 1; ii <= theparams->themap->m; ii++) /* clean out matrices */
          for (jj = 1; jj <= theparams->themap->maxl; jj++)
	        beta0[ii][jj] = beta1[ii][jj] = fstat[ii][jj] = (FPN) 0.0;

        shuffle_ivector(indx, 1, theparams->nn);
        for (ii = 1; ii <= theparams->nn; ii++)
          theparams->thedata[ii].y[traits] = origtraits[indx[ii]];
	    if ( rows > 2)
	        error = calc_mlrstats(theparams,theparams->themap, theparams->thedata,traits,beta0,beta1,fstat,lrvalues,theparams->thegenome,lnpk);
        else
	        error = calc_lrstats(theparams,theparams->themap, theparams->thedata,  traits, beta0, beta1, fstat,  theparams->thegenome);
        for (ii = 1; ii <= theparams->themap->m; ii++) {
          for (jj = 1; jj <= theparams->themap->mpc[ii]; jj++) {
            if (  savefstat[ii][jj] < fstat[ii][jj] )
	          pcnts[ii][jj] =  pcnts[ii][jj] + 1;
            if (   fstat[ii][jj] > maxlr[kk] )
	          maxlr[kk] = fstat[ii][jj];
          }
        }
      }
	  if ( theparams->verbosity == 1 )
		printf("\nFinished permutation test for trait %d, now writing results.\n",traits);
      sort(theparams->perms, maxlr);
      outf = fileopen(theparams->lrfile, "a");
      if (outf != NULL) {
        fprintf(outf, "\n#\n# Performed %d permutations of the phenotypes and genotypes", theparams->perms);
        fprintf(outf, "\n# Here are the comparisonwise counts of permuted test statistics\n");
        fprintf(outf, "\n#Chrom Mark MarkerName Cnts Pval");
        for (ii = 1; ii <= theparams->themap->m; ii++) {
          for (jj = 1; jj <= theparams->themap->mpc[ii]; jj++) {
            pval = (FPN) pcnts[ii][jj] / (FPN) theparams->perms ;
            if ( theparams->themap->names != NULL )  
              fprintf(outf, "\n %4d %4d %10s %4d %6.2f",ii,jj,theparams->themap->names[theparams->themap->ttable[ii][jj]],pcnts[ii][jj], pval);
            else
              fprintf(outf, "\n %4d %4d            %4d %6.2f",ii,jj,pcnts[ii][jj], pval);
          }
        }
        fprintf(outf, "\n# Here are the experimentwise significance levels for different sizes");
        fprintf(outf, "\n# Permutation significance level for alpha = 0.1   : %7.4f", maxlr[(int) (0.9 * (FPN) theparams->perms) + 1]);
        fprintf(outf, "\n# Permutation significance level for alpha = 0.05  : %7.4f", maxlr[(int) (0.95 * (FPN) theparams->perms) + 1]);
        fprintf(outf, "\n# Permutation significance level for alpha = 0.025 : %7.4f", maxlr[(int) (0.975 * (FPN) theparams->perms) + 1]);
        fprintf(outf, "\n# Permutation significance level for alpha = 0.01  : %7.4f", maxlr[(int) (0.99 * (FPN) theparams->perms) + 1]);
        fprintf(outf, "\n#end of shuffling results\n#\n");
        fileclose(theparams->lrfile, outf);
      }
      for ( jj = 1 ; jj<= theparams->nn ; jj++ ) /* put traits back */
        theparams->thedata[jj].y[traits] = origtraits[jj];
    }

  }
/* Clean up...*/
  write_trailer(theparams,chptr,1);
  if (theparams->perms > 0 && theparams->perms < MAX_REPS) {
	free_ivector(indx,1, theparams->nn);
    free_dmatrix(savefstat,1, theparams->themap->m, 1, theparams->themap->maxl);
	free_dvector(origtraits,1, theparams->nn);
	free_imatrix(pcnts,1, theparams->themap->m, 1, theparams->themap->maxl);
	free_dvector(maxlr,1, theparams->perms);
  }
  free_dmatrix(beta0, 1, theparams->themap->m, 1, theparams->themap->maxl);
  free_dmatrix(beta1, 1, theparams->themap->m, 1, theparams->themap->maxl);
  free_dmatrix(fstat, 1, theparams->themap->m, 1, theparams->themap->maxl);
  if ( rows > 2 )
    free_dmatrix(lrvalues, 1, theparams->themap->m, 1, theparams->themap->maxl);
  free_cvector( purpose,0,MAXNAME);
  if ( rows > 2 )
    lnpk = cd_linpakws(lnpk,theparams->nn,rows,0);
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
    strcpy(opt[6],  "-r"  );
    strcpy(opt[7],  "-t"  );

    strcpy(opt_e[1],  "Data Input File"  );
    strcpy(opt_e[2],  "Output File"  );
    strcpy(opt_e[3],  "Error File"  );
    strcpy(opt_e[4],  "Genetic Linkage Map File"  );
    strcpy(opt_e[5],  "Random Number Seed"  );
    strcpy(opt_e[6],  "Number of Permutations"  );
    strcpy(opt_e[7],  "Trait to Analyze"  );
  }
  for ( ii = 1 ; ii <= nopts ; ii++ ) 
    for ( jj = 0 ; jj <= MAXNAME ; jj++ )
      opt_v[ii][jj] = '\0';

    strcpy(opt_v[1],  theparams->ifile  );
    strcpy(opt_v[2],  theparams->lrfile  );
    strcpy(opt_v[3],  theparams->error  );
    strcpy(opt_v[4],  theparams->map  );
    sprintf(opt_v[5],"%ld",theparams->seed );
    sprintf(opt_v[6],"%d",theparams->perms );
    sprintf(opt_v[7],"%d",theparams->whichtrait );

}  

void update_params(char **opt_v,  params *theparams)
{
    strcpy(theparams->ifile, opt_v[1] );
    strcpy(theparams->lrfile, opt_v[2] );
    strcpy(theparams->error, opt_v[3]  );
    strcpy(theparams->map, opt_v[4]  );
    theparams->seed = atol(opt_v[5]);
    theparams->perms = atoi(opt_v[6]);
    theparams->whichtrait = atoi(opt_v[7]);
}  

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file LRmain.c
------------------------------------------------------------------ */


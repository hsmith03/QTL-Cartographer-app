/* ------------------------------------------------------ XCutXCodeXSkip
     This file (BTL_LRmain.c) is part of QTL Cartographer
         
    		Copyright (C) 1999-2005
	Lauren McIntyre and Jun Wu

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
    Additional copyright belongs to Lauren McIntyre and Jun Wu 
    of Biological Sciences and Statistics at Purdue University


          Main function for  Binary LRmapqtl.          

       Created by Jun Wu from QTL Cartographer source code.
*/

#include "Main.h"

int main(argc, argv)
int argc;
char *argv[];
{
  FILE *errorf, *outf;
  char  *chptr,*purpose;
  params *theparams;
  int ii,kk,automatic,jj,nopts,rows,itrait,utrait;
  int error, lastmarker, traits;
  FPN *origtraits,*maxlr,**savefstat,pval;
  int *indx,**pcnts;
  
  FPN **beta0, **beta1, **beta2,**fstat,**fstatx1,**fstatx2,**lrvalues, **Rmq, **p3;
  linpakws *lnpk;

#if defined(MACWARRIOR)
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  /*whichprogram 11 is BTL_LRmap function.......................................*/
  whichprogram = 11;
  purpose = cvector(0,MAXNAME);
  strcpy(purpose,"Map Binary QTL using Simple Linear Regression");
  theparams = NULL;
  nopts = 16;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();
  /*initialize the p1,p2 and range of p1,p2*/
  theparams->p1=(FPN)0.5;
  theparams->p2=(FPN)0.5;
  theparams->p1_start=(FPN)0.0;
  theparams->p1_end=(FPN)1.0;
  theparams->p1_step=(FPN)0.1;
  theparams->p2_start=0;
  theparams->p2_end=1;
  theparams->p2_step=(FPN)0.1;
  theparams->Btl_mode=0;

  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  print_head(argv[0],theparams->lrfile,chptr,1,50,theparams);

 
/* Initailize the data structures...*/
  GetTheMap(theparams, theparams->map);
  GetTheData(theparams, theparams->ifile, 0 );

  /* Binary Translate for B1 cross from AA(3,2), Aa(1) to 0,1 */
  /* Binary Translate for B2 cross from aa(0), Aa(1) to 0,1 */
  
  
   if(theparams->cross==1 || theparams->cross==2)/* /means B1,B2  backcross*/
  {
	  for (ii = 1; ii <= theparams->nn; ii++)
          BTL_trans_data(theparams->thedata[ii].markers, theparams->themap);
  }


  if ( theparams->themap->otraits > 0 )
    process_otraits(theparams->thedata,theparams->nn,theparams->themap);


	lastmarker = theparams->themap->maxl;
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
	beta2 = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
    fstat = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
    fstatx1 = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
    fstatx2 = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);

	p3 = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);
	Rmq = dmatrix(1, theparams->themap->m, 1, theparams->themap->maxl);

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
	    beta0[ii][jj] = beta1[ii][jj] = beta2[ii][jj]=p3[ii][jj]=Rmq[ii][jj]=fstat[ii][jj] = (FPN)0.0;
	
	if ( rows > 2 )
      error = calc_mlrstats(theparams,theparams->themap, theparams->thedata,  traits, beta0, beta1, fstat,lrvalues, theparams->thegenome,lnpk);
    else if(theparams->cross ==3)
	    BTL_calc_F2(theparams,   theparams->thedata, traits, beta0,beta1,beta2,fstat,fstatx1,fstatx2, theparams->thegenome);
	else 
      error = BTL_calc_lrstats(theparams,  theparams->thedata,  traits, beta0, beta1, fstat, theparams->thegenome);

    if (error != 0) {/* If there is some errror, */
      errorf = fileopen(theparams->error, "a");
      fprintf(errorf, "\n\nHad trouble doing the linear regression of data from %s...\n", theparams->ifile);
		fprintf(errorf, "errorlevel = %4d\n", error);
      fileclose(theparams->error, errorf);
    }

	if ( rows > 2 )
	 	print_mlrstats(theparams,traits,theparams->themap,beta0,beta1,fstat,lrvalues);
    else
	 	BTL_print_lrstats(theparams->nn, traits, theparams->themap, theparams->map, theparams->ifile, theparams->lrfile);

	
	/*measure recombination rate (Rmq ) and estimated P3.*/
	error=calc_Btl(theparams->nn,theparams->lrfile,theparams->themap,theparams->cross,beta0,beta1,beta2,Rmq,p3,theparams,fstat,fstatx1,fstatx2);
   
   if (theparams->perms > 0 && theparams->perms < MAX_REPS) { /* do a shuffle test...*/
      for (ii = 1; ii <= theparams->nn; ii++)  /* indx[ii] = ii */
	    indx[ii] = ii;
      for (ii = 1; ii <= theparams->themap->m; ii++)
        for (jj = 1; jj <= theparams->themap->maxl; jj++) {  /* zero out pcnts, save f stats */
	      pcnts[ii][jj] = 0;
          savefstat[ii][jj] = fstat[ii][jj];
        }
      for ( jj = 1 ; jj <= theparams->perms ; jj++ )
        maxlr[jj] = (FPN)0.0;
      for ( jj = 1 ; jj<= theparams->nn ; jj++ )  /* save original trait values */
        origtraits[jj] = theparams->thedata[jj].y[traits];
      for (kk = 1; kk <= theparams->perms; kk++) {  /* do permutations */
		if ( theparams->verbosity == 1 )
			 printf("\nNow on Permutation %d",kk);
        for (ii = 1; ii <= theparams->themap->m; ii++) /* clean out matrices */
          for (jj = 1; jj <= theparams->themap->maxl; jj++)
	        beta0[ii][jj] = beta1[ii][jj] = fstat[ii][jj] = (FPN)0.0;

        shuffle_ivector(indx, 1, theparams->nn);
        for (ii = 1; ii <= theparams->nn; ii++)
          theparams->thedata[ii].y[traits] = origtraits[indx[ii]];
	    if ( rows > 2)
	        error = calc_mlrstats(theparams,theparams->themap, theparams->thedata,traits,beta0,beta1,fstat,lrvalues,theparams->thegenome,lnpk);
        else
	        error = calc_lrstats(theparams,theparams->themap, theparams->thedata,  traits, beta0, beta1, fstat,  theparams->thegenome);
        for (ii = 1; ii <= theparams->themap->m; ii++) {
          if ( theparams->themap->sigl > 0 ) 
            lastmarker = theparams->themap->mpc[ii];
          for (jj = 1; jj <= lastmarker; jj++) {
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
          if ( theparams->themap->sigl > 0 ) 
            lastmarker = theparams->themap->mpc[ii];
          for (jj = 1; jj <= lastmarker; jj++) {
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
  free_dmatrix(fstatx1, 1, theparams->themap->m, 1, theparams->themap->maxl);
  free_dmatrix(fstatx2, 1, theparams->themap->m, 1, theparams->themap->maxl);
  free_dmatrix(beta2, 1, theparams->themap->m, 1, theparams->themap->maxl);
  free_dmatrix(p3, 1, theparams->themap->m, 1, theparams->themap->maxl);
  free_dmatrix(Rmq, 1, theparams->themap->m, 1, theparams->themap->maxl);
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
   strcpy(opt[8],  "-p1"  );
   strcpy(opt[9],  "-p2"  );
    strcpy(opt[10],  "-p1s"  );
   strcpy(opt[11],  "-p1e"  );
    strcpy(opt[12],  "-p1t"  );
   strcpy(opt[13],  "-p2s"  );
   strcpy(opt[14],  "-p2e"  );
   strcpy(opt[15],  "-p2t"  );
   strcpy(opt[16],  "-pmode");
   
    strcpy(opt_e[1],  "Data Input File"  );
    strcpy(opt_e[2],  "Output File"  );
    strcpy(opt_e[3],  "Error File"  );
    strcpy(opt_e[4],  "Genetic Linkage Map File"  );
    strcpy(opt_e[5],  "Random Number Seed"  );
    strcpy(opt_e[6],  "Number of Permutations"  );
    strcpy(opt_e[7],  "Trait to Analyze"  );
    strcpy(opt_e[8],  "input probability value of trait for parent 1 " );
    strcpy(opt_e[9],  "input probability value of trait for parent 2 " );
    strcpy(opt_e[10],  "start trait prob value for parent 1 " );
    strcpy(opt_e[11],  "end trait prob value for parent 1 " );
    strcpy(opt_e[12],  "step size for parent 1 trait prob. " );
    strcpy(opt_e[13],  "start trait prob value for parent 2 " );
    strcpy(opt_e[14],  "end trait prob value for parent 2 " );
    strcpy(opt_e[15],  "step size for parent 2 trait prob" );
	strcpy(opt_e[16],  "mode 0: input p1,p2, mode 1: range test");
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
    sprintf(opt_v[8],"%f",theparams->p1  );
    sprintf(opt_v[9],"%f",theparams->p2  );
   sprintf(opt_v[10],"%f",theparams->p1_start  );
   sprintf(opt_v[11],"%f",theparams->p1_end  );
   sprintf(opt_v[12],"%f",theparams->p1_step  );
   sprintf(opt_v[13],"%f",theparams->p2_start  );
   sprintf(opt_v[14],"%f",theparams->p2_end  );
   sprintf(opt_v[15],"%f",theparams->p2_step  );
   sprintf(opt_v[16],"%d",theparams->Btl_mode  );
}  

void update_params(char **opt_v, params *theparams)
{
    strcpy(theparams->ifile, opt_v[1] );
    strcpy(theparams->lrfile, opt_v[2] );
    strcpy(theparams->error, opt_v[3]  );
    strcpy(theparams->map, opt_v[4]  );
    theparams->seed = atol(opt_v[5]);
    theparams->perms = atoi(opt_v[6]);
    theparams->whichtrait = atoi(opt_v[7]);
    theparams->p1=(FPN) atof(opt_v[8]);
    theparams->p2=(FPN) atof(opt_v[9]);
    theparams->p1_start=(FPN) atof(opt_v[10]);
    theparams->p1_end=(FPN) atof(opt_v[11]);
    theparams->p1_step=(FPN) atof(opt_v[12]);
    theparams->p2_start=(FPN) atof(opt_v[13]);
    theparams->p2_end=(FPN) atof(opt_v[14]);
    theparams->p2_step=(FPN) atof(opt_v[15]);
	theparams->Btl_mode= atoi(opt_v[16]);
    
    
}  

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file BTL_LRmain.c
------------------------------------------------------------------ */


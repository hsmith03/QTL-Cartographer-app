/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Zmain.c) is part of QTL Cartographer
         
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


/*  Main function driver and subroutines for Zmapqtl.*/

#include "Main.h"

int main(argc, argv)
int argc;
char *argv[];
{
  char *chptr, *purpose,*tfile, *bootfile,*jackfile;
  FPN  *tmp,**yy, *yyy, **lin_reg,*lratio, maxlr, dumm;
  int nopts,jj, ii, automatic,rows, go_on, do_analysis,noshuffle;
  int *indx,*samplesize,kk,error;
  int *pcnts, ipos,bootnum, startperm,jacknum,wt,nbp,srm,savemodel;
  params *theparams;
  genome  *agptr,*tptr,*gptr,   *startptr, *endptr;
  linpakws *lnpk;

#if defined(MACWARRIOR)  
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif

  whichprogram = 7;
/* just to be on the safe side, initialize pointers to null */
  yyy = NULL;
  yy = lin_reg  = NULL;
  tmp = lratio = NULL;
  indx = samplesize = pcnts= NULL;
  agptr=tptr=gptr= startptr=endptr=NULL;

  purpose = cvector(0,MAXNAME);
  tfile = cvector(0,MAXNAME);
  bootfile = cvector(0,MAXNAME);

 /**/
  strcpy(purpose, "Do (Composite) Interval Mapping");
  theparams = NULL;
  nopts = 15;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);

  free_cvector(theparams->tfile, 0 , MAXNAME);
  theparams->tfile = tfile;
  theparams->walk = theparams->walk / (FPN) 100.0;
  
  srm = theparams->srm;
  if ( theparams->Model == 8 ) { /* 8 is a model to tell Zmapqtl to use srm model 3.  */
    theparams->Model = 6;
    theparams->srm = 3;
  }
/*Are we going to do a basic analysis on all traits?  If so, turn off
  permutations, bootstraps and jackknives.  Also, restrict to models 3 and 6 */
  if (theparams->whichtrait < 1 || theparams->whichtrait > theparams->traits ) {
    theparams->perms = theparams->boots = 0;
    if ( theparams->Model != 3 && theparams->Model != 6 )
      theparams->Model = 3;
  }
/* Initialize the data structures...*/
  GetTheMap(theparams, theparams->map);

 

/*  theparams->themap->knum = NULL;*/
  if ( theparams->Model == 4 || theparams->Model == 5 )
	 lin_reg = get_lin_reg(theparams->lrfile, theparams->themap->m, theparams->themap->maxl, theparams->whichtrait);
  else
	 lin_reg = NULL;
  if ( theparams->Model == 7 )
    GetTheModel(theparams, theparams->eqtl) ; 
  else
    theparams->theqtls = NULL;
  GetTheData(theparams, theparams->ifile, 1 );





  switch (theparams->Model ) {
    default: case 1: case 2: 
      rows = theparams->themap->ml; break;
    case 3:
      rows = 1; break;
    case 4:
      rows = theparams->themap->m; break;
    case 5:
      rows = theparams->themap->maxl + 2*(theparams->themap->m - 1) ; break;
    case 6:
      rows = 2*theparams->nbp + 1; break;
    case 7:
      rows = theparams->themap->knum[theparams->whichtrait] + 1;
    }
  if ( theparams->themap->otraits > 0 ) {
      process_otraits(theparams->thedata,theparams->nn,theparams->themap);
      for ( ii = 1 ; ii <= theparams->themap->otraits ; ii++ )
        if ( theparams->themap->onames[ii][0] == '+' )
          rows = rows +  2*theparams->themap->otypes[ii] - 2;
  }
  nbp = theparams->nbp;
/*
  At this point we should be able to figure out an upper bound to 
  how much workspace we need.  
*/
  lnpk = cd_linpakws(NULL,theparams->nn,rows,1);	
  lnpk->maxlr = NULL;  
  lnpk->pp1 = dvector(1,theparams->nn);
  lnpk->pp2 = dvector(1,theparams->nn);
  lnpk->pv = dvector(1,theparams->nn);
  lnpk->qv = dvector(1,theparams->nn);
  lnpk->bp = imatrix( 1, 2, 1, theparams->themap->ml);
  lnpk->wrsd = dmatrix(1,1,1,theparams->nn);
  lnpk->wy = dvector(1,theparams->nn);
  lnpk->k = theparams->themap->ml;
  yy = dmatrix(1, theparams->traits, 1, theparams->nn);
  samplesize = ivector(1, theparams->traits);
  if (theparams->whichtrait > 0 && theparams->whichtrait <= theparams->traits )
    init_phenotypes(theparams,theparams->thedata,yy,samplesize);
/*	Create a linked list representing the genome and
	translate the distances in recombination frequencies to those
	of Morgans...  */
  theparams->thegenome = create_genome(theparams->themap);
  agptr = theparams->thegenome;
  if ( theparams->Model == 7 )
        zplace_qtls(theparams,theparams->theqtls, theparams->thegenome);

  if ( theparams->Model == 6 ) {/* Get the results from a previous run of SRmapqtl */
    if (theparams->whichtrait > 0 && theparams->whichtrait <= theparams->traits )
      error = get_srresults(theparams,theparams->themap,theparams->thegenome);
  }
/*  Determine the starting and ending points for the analysis. */
  go_on = do_analysis = 1;

  endptr = startptr = theparams->thegenome;
  while (endptr->next != NULL)
    endptr = endptr->next;
  if (theparams->wchrom > 0 && theparams->wchrom <= theparams->themap->m) {
    while (startptr->chrom != theparams->wchrom)
      startptr = startptr->next;
    endptr = startptr;
    while (endptr->next != NULL && endptr->next->chrom == theparams->wchrom)
      endptr = endptr->next;
  }

  if (do_analysis == 1) {
	  if ( theparams->perms > 0 && theparams->perms < MAX_REPS ) { /* Do Permutation test...*/
        lnpk->maxlr = dvector(0,theparams->themap->m);
		indx = ivector(1, theparams->nn);     /*trait index vector, to be permuted*/
		yyy = dvector(1, theparams->nn);      /*permuted trait vector*/
		for (ii = 1; ii <= theparams->nn; ii++)
		  indx[ii] = ii;                 /* start by making indx[i]=i */
		ipos = 0;         /*How many test positions?  This will count them*/
	    for (gptr = startptr; gptr != NULL && gptr != endptr->next; gptr = gptr->next)
	      ipos = ipos + 1 + (int) (gptr->dist / theparams->walk);
	    lratio = dvector(1, ipos);   /* Keep track of the LR of the sample */
	    pcnts = ivector(1, ipos);    /* Count number of shuffled sets with LR > SampleLR */
	    for (ii = 1; ii <= ipos; ii++)
	      pcnts[ii] = 0;    /* Set them equal to zero  */
	
	    sprintf(tfile,"%s.z%dc", theparams->stem,theparams->Model);
	    sprintf(bootfile,"%s.z%de", theparams->stem,theparams->Model);
	    kk = isfile(tfile); /* 1=>file, 0=>none */
	    if ( kk == 1)  {
	      startperm = get_CWTR(theparams,tfile,pcnts,lratio);
	      startperm = startperm+1;
	    }
	    else {
          startperm = 1;
	      get_lratio(theparams,lratio);
	      print_head(argv[0],bootfile,chptr,1,61,theparams);  /* write header for experimentwise results...*/
          write_zheader(bootfile, theparams,   chptr,2,0,theparams->themap);
	    }
	    noshuffle = 0;
	    if ( theparams->perms == 1 ) {
	      noshuffle = 1;
	      theparams->perms = startperm;
	    }
	    
	    
	    
	    for (jj = startperm; jj <= theparams->perms; jj++) {
	      if ( noshuffle != 1 )
	        shuffle_ivector(indx, 1, samplesize[theparams->whichtrait]);
	      for (ii = 1; ii <= samplesize[theparams->whichtrait]; ii++)
	        yyy[ii] = yy[theparams->whichtrait][indx[ii]];
	      
	      maxlr = do_zanalysis(theparams, startptr, endptr, theparams->themap, theparams->theqtls, lin_reg, yyy, theparams->thedata, NULL, lratio, pcnts,agptr,lnpk);
          if ( theparams->rwd == 1 )
            write_maxlr(bootfile,jj,lnpk->maxlr,theparams->themap->m);
          else
            write_maxlr(bootfile,jj,lnpk->maxlr,0);
          print_head(argv[0],tfile,chptr,0,62,theparams); 
          write_zheader(tfile, theparams,   chptr,1,0,theparams->themap);
	      write_CWTR(theparams,startptr,endptr,theparams->themap,theparams->walk,tfile,lratio,pcnts,jj);
	    }
	    free_dvector(yyy, 1, theparams->nn);
	    free_ivector(indx, 1, theparams->nn);    
	    free_dvector(lratio, 1, ipos);
	    free_ivector(pcnts, 1, ipos);
        free_dvector(lnpk->maxlr, 0,theparams->themap->m);
	  }
      else if (theparams->boots == 1 ) { /* Do a bootstrap, assuming that
                                            Prune had been run to create the datafile.  */
	      sprintf(tfile,"%s.z%da", theparams->stem,theparams->Model);/*Old bootstrap file*/
	      sprintf(bootfile,"%s.z%db", theparams->stem,theparams->Model);/*New bootstrap file*/
	      kk = isfile(tfile); /* 1=>file, 0=>none */
	      if ( kk == 0)  {
	        print_head(argv[0],tfile,chptr,0,63,theparams);  /* write header old bootstrap file*/
            write_zheader(tfile, theparams,  chptr,3,0,theparams->themap);
            error = convert_zfile(theparams,tfile);
            bootnum = 0;
	      }
	      else
	        bootnum = get_bootnum(tfile);
	      bootnum = bootnum+1;
          print_head(argv[0],bootfile,chptr,0,63,theparams); 
          write_zheader(bootfile, theparams,  chptr,3,bootnum,theparams->themap);
	      dumm = do_zanalysis(theparams,startptr,endptr,theparams->themap,theparams->theqtls,lin_reg,yy[theparams->whichtrait],theparams->thedata,bootfile,lratio,NULL,agptr,lnpk);
      }
      else if (theparams->boots == 2 ) { /* Do a jackknife.  */
	      sprintf(tfile,"%s.z%di", theparams->stem,theparams->Model);/* Old Jackknife summary statistics */
	      sprintf(bootfile,"%s.z%dj", theparams->stem,theparams->Model);/* Jackknife summary statistics */
	      kk = isfile(tfile); /* 1=>file, 0=>none */
	      if ( kk == 0)  {
	        print_head(argv[0],tfile,chptr,0,63,theparams);  /* write header old bootstrap file*/
            write_zheader(tfile, theparams,  chptr,3,0,theparams->themap);
            error = convert_zfile(theparams,tfile);
            bootnum = 0;
	      }
	      else
	        bootnum = get_bootnum(tfile);
	      jacknum = bootnum = bootnum+1;
	      theparams->tfile = tfile;
	      jackfile = bootfile;
          print_head(argv[0],bootfile,chptr,0,64,theparams); 
          write_zheader(bootfile, theparams,  chptr,3,bootnum,theparams->themap);
          
	      for (kk = bootnum; kk <= theparams->nn; kk++) 
	        if ( theparams->thedata[kk].print_flag == 'y' ) {
	          theparams->thedata[kk].print_flag = 'n';
              print_head(argv[0],jackfile,chptr,0,64,theparams); 
              write_zheader(jackfile,theparams,chptr,3,jacknum,theparams->themap);
	          dumm = do_zanalysis(theparams,startptr,endptr,theparams->themap,theparams->theqtls,lin_reg,yy[theparams->whichtrait],theparams->thedata,jackfile,lratio,NULL,agptr,lnpk);
              theparams->thedata[kk].print_flag = 'y';
	          jacknum = jacknum+1;
	          if ( jackfile == bootfile ) {
	            jackfile = tfile;
	            theparams->tfile = bootfile;
	          }
	          else {
	            jackfile = bootfile;
	            theparams->tfile = tfile;
	          }
	        }
      }
	  else { /* Do basic analysis. */ 
	        /*  Write a header for the output file.  */
	    savemodel = theparams->Model;
		print_head(argv[0],theparams->zfile,chptr,1,60,theparams);
		lratio = NULL;
		pcnts = NULL;
	    wt = theparams->whichtrait;
	    for ( jj = 1; jj <= theparams->traits; jj++ ) 
	      if ( (wt == 0 && theparams->themap->tnames[jj][0] == '+') || (wt > theparams->traits && theparams->themap->tnames[jj][0] != '-') || (wt == jj) ) {
            if ( theparams->verbosity == 1 ) 
              printf("\n\tNow analyzing trait %d named %s\n",jj,theparams->themap->tnames[jj]);
            theparams->nbp = nbp;
            theparams->whichtrait = jj;
            theparams->Model = savemodel;
            init_phenotypes(theparams,theparams->thedata,yy,samplesize);
            if ( theparams->Model == 6 ) /* Get the results from a previous run of SRmapqtl */
              error = get_srresults(theparams,theparams->themap,theparams->thegenome);
		    write_zheader(theparams->zfile,theparams,chptr,0,0,theparams->themap);
	      	dumm = do_zanalysis(theparams,startptr,endptr,theparams->themap,theparams->theqtls,lin_reg,yy[theparams->whichtrait],theparams->thedata,theparams->zfile,lratio,NULL,agptr,lnpk);
	      }
	  }
    }
/*  reset parameters to proper values.*/
  theparams->srm = srm;
  theparams->nbp = nbp;
  theparams->walk = theparams->walk * (FPN) 100.0;
  theparams->Inter = 0;
  write_trailer(theparams,chptr,1);
  
/* Clean up...*/
  if ( purpose != NULL )
	 free_cvector(purpose,0,MAXNAME);
  if ( tfile != NULL )
	 free_cvector(tfile,0,MAXNAME);
  if ( bootfile != NULL )
	 free_cvector(bootfile,0,MAXNAME);
  theparams->tfile = NULL;
  if ( lin_reg != NULL )
	 free_dmatrix(lin_reg, 1, theparams->themap->m, 1, theparams->themap->maxl);
  if ( yy != NULL )
	 free_dmatrix(yy, 1, theparams->traits, 1, theparams->nn);
  if ( samplesize != NULL )
	 free_ivector(samplesize, 1, theparams->traits);
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
    strcpy(opt[5], "-l");
    strcpy(opt[6], "-S");
    strcpy(opt[7], "-s");
    strcpy(opt[8], "-M");
    strcpy(opt[9], "-t");
    strcpy(opt[10], "-c");
    strcpy(opt[11], "-d");
    strcpy(opt[12], "-n");
    strcpy(opt[13], "-w");
    strcpy(opt[14], "-r");
    strcpy(opt[15], "-b");

    strcpy(opt_e[1], "Input File");
    strcpy(opt_e[2], "Output File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Genetic Linkage Map File");
    strcpy(opt_e[5], "LRmapqtl Results file (Models 4&5)");
    strcpy(opt_e[6], "SRmapqtl Results file (Model 6)");
    strcpy(opt_e[7], "Random Number Seed");
    strcpy(opt_e[8], "Model [1-6], 3=>IM");
    strcpy(opt_e[9], "Trait to analyze");
    strcpy(opt_e[10], "Chromosome to analyze (0=>all)");
    strcpy(opt_e[11], "Walking speed in cM");
    strcpy(opt_e[12], "Number of Background Parameters (Model 6)");
    strcpy(opt_e[13], "Window Size in cM (Models 5&6)");
    strcpy(opt_e[14], "Number of Permutations");
    strcpy(opt_e[15], "Number of Bootstraps"); 
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->ifile);
  strcpy(opt_v[2], theparams->zfile);
  strcpy(opt_v[3], theparams->error);
  strcpy(opt_v[4], theparams->map);
  strcpy(opt_v[5], theparams->lrfile);
  strcpy(opt_v[6], theparams->srfile);
  sprintf(opt_v[7], "%ld", theparams->seed);
  sprintf(opt_v[8], "%d", theparams->Model);
  sprintf(opt_v[9], "%d", theparams->whichtrait);
  sprintf(opt_v[10], "%d", theparams->wchrom);
  sprintf(opt_v[11], "%f", theparams->walk);
  sprintf(opt_v[12], "%d", theparams->nbp);
  sprintf(opt_v[13], "%f", theparams->window);
  sprintf(opt_v[14], "%d", theparams->perms);
  sprintf(opt_v[15], "%d", theparams->boots); 


}

void update_params(char **opt_v,  params *theparams)
{
  theparams->nbp = atoi(opt_v[12]);
  theparams->window = (FPN) atof(opt_v[13]);
  if (theparams->nbp < 0 )
    theparams->nbp = (int) NUM_SIG;
  if ( theparams->window < (FPN) 0.0)
    theparams->window = (FPN) WIN_SIZE;
  theparams->Model = atoi(opt_v[8]);

  strcpy(theparams->ifile, opt_v[1]);

  strcpy(theparams->zfile, opt_v[2]);
  strcpy(theparams->error, opt_v[3]);
  strcpy(theparams->map, opt_v[4]);
  strcpy(theparams->lrfile, opt_v[5]);
  strcpy(theparams->srfile, opt_v[6]);
  theparams->seed = atol(opt_v[7]);
  theparams->whichtrait = atoi(opt_v[9]);
  theparams->wchrom = atoi(opt_v[10]);
  theparams->walk = (FPN) atof(opt_v[11]);

  theparams->perms = atoi(opt_v[14]);
  if (theparams->perms >= (int) MAX_REPS)
    theparams->perms = (int) MAX_REPS - 1;
  theparams->boots = atoi(opt_v[15]);
  if (theparams->boots >2)
    theparams->boots = 1;

}



/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Zmain.c
------------------------------------------------------------------ */


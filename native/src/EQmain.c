/* ------------------------------------------------------ XCutXCodeXSkip
     This file (EQmain.c) is part of QTL Cartographer
         
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


/*Main function driver for Eqtl.


Assume that theparams->siglevel is a likelihood ratio threshold.
   if theparams->lodflag is turned on (to 1), then convert LRs to lod scores
   leave theparams->siglevel alone, use a temporary variable for the threshold.

*/


#include "Main.h"
void process_jzoutput( params *theparams,char *bootin,   genome *first);
void process_zoutput(params *theparams,markermap *themap, char *progname, char *chptr);
void  AutoProcessResults(params *theparams);

int main(argc, argv)
int argc;
char *argv[];
{
  char *chptr, *purpose, *bootin,*bootout;
  int nopts,automatic,error,kk,i;
  params *theparams;

#if defined(MACWARRIOR)
 /* 
  Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
 /**/
  whichprogram = 8;
  purpose = cvector(0,MAXNAME);
  bootin = cvector(0,MAXNAME);
  bootout = cvector(0,MAXNAME);
  strcpy(purpose, "Estimate the postions and effects of QTLs");
  theparams = NULL;
  nopts = 11;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();
  automatic =   process_arguments(argc,argv,chptr,purpose,nopts,theparams);
/* 
Initailize the data structures...
*/
  theparams->cross = get_cross(theparams->ifile,theparams);    
  GetTheMap(theparams, theparams->map );
  theparams->thegenome   = create_genome(theparams->themap);

  theparams->chrom = theparams->themap->m;
  error = get_the_nn(&theparams->nn, &theparams->traits, theparams->ifile);
  theparams->themap->traits = theparams->traits;
  get_traitnames(theparams,theparams->themap);
  if (theparams->traits > 1) {
    theparams->themap->traits = theparams->traits;
    free_ivector(theparams->themap->knum, 1, 1);
    theparams->themap->knum = ivector(1, theparams->traits); /*
    for (ii = 1; ii <= theparams->traits; ii++) {
      theparams->themap->knum[ii] = (int) ((FPN) theparams->qnum * ranf(ii)) + 1;
      if (theparams->themap->knum[ii] < 1)
	    theparams->themap->knum[ii] = 1;
    }*/
  }
/* 
Write header for outut file...
*/
  if ( theparams->boot != 0 )
    print_head(argv[0],theparams->eqtl,chptr,1,70,theparams);
  else
    print_head(argv[0],theparams->eqtl,chptr,0,70,theparams);

  for ( i=0; i<=8 && theparams->mimwork[i] != '\0' ; i++ ) {
    if ( theparams->mimwork[i] == 'B' ) {                 
      sprintf(bootin,"%s.z%da",theparams->stem,theparams->Model);
      sprintf(bootout,"%s.z%db",theparams->stem,theparams->Model);
      kk = isfile(bootin);
      if ( kk == 1 ) /*  If there is a bootfile, then process it. */
        kk = process_bootfile(theparams,bootin,bootout,argv[0],chptr);
    }
    if ( theparams->mimwork[i] == 'J' ) { 
      sprintf(bootin,"%s.z%di",theparams->stem,theparams->Model);
      sprintf(bootout,"%s.z%dj",theparams->stem,theparams->Model);
      kk = isfile(bootin);
      if ( kk == 1 ) /*  If there is a jackfile, then process it. */
        kk = process_bootfile(theparams,bootin,bootout,argv[0],chptr);
    }
    if ( theparams->mimwork[i] == 'P' ) { 
      sprintf(bootin,"%s.z%de",theparams->stem,theparams->Model);
      sprintf(bootout,"%s.z%dc",theparams->stem,theparams->Model);
      kk = isfile(bootin);
      if ( kk == 1 ) /*  If there is a permfile, then process it */
        kk = process_permEfile(theparams,bootin,bootout);
    }
    if ( theparams->mimwork[i] == 'Z' )   
      process_zoutput(theparams,theparams->themap,argv[0],chptr);
    if ( theparams->mimwork[i] == 'M' ) 
      process_jzoutput(theparams,bootin, theparams->thegenome); 
    if ( theparams->mimwork[i] == 'A' ) 
      AutoProcessResults(theparams); 
  }
  theparams->lodflag = 0;
  
  write_trailer(theparams,chptr,1);
  free_cvector( purpose,0,MAXNAME);
  free_cvector(bootin,0,MAXNAME);
  free_cvector(bootout,0,MAXNAME);
  theparams = create_params(theparams, 0, NULL);

  return(0);
}


/*
    I'm setting this up as an alternative way to process
    all analytical results.   

*/
void  AutoProcessResults(params *theparams) {
  FILE *out;
  aresults *theZresults, *theJZresults, *thePCLresults, *cnode, *AllResults;
  char infile[MAXNAME];
  int indfile;
/*  open the Zmapqtl output file and process all results */  
  theZresults = theJZresults = thePCLresults = NULL;
  indfile =   isfile(theparams->zfile);
  if ( indfile == 1 ) {
    theZresults = AnalysisFirstPass(theparams->zfile);
    for ( cnode=theZresults; cnode!=NULL; cnode=cnode->next)
      SetAnalysisHeaders(cnode,theparams);
    AnalysisSecondPass(theZresults);  
  }
  AllResults = theZresults;
  
/* open the JZmapqtl file for joint analysis and process the results*/  
  sprintf(infile, "%s0", theparams->zfile);
  indfile =   isfile(infile);
  if ( indfile == 1 ) {
    theJZresults = AnalysisFirstPass(infile);
    for ( cnode=theJZresults; cnode!=NULL; cnode=cnode->next)
      SetAnalysisHeaders(cnode,theparams);
    AnalysisSecondPass(theJZresults);  
  }

/*  append them to the first set if they exists, otherwise make them the first set. */
  if ( AllResults != NULL )
    AppendResults(AllResults, theJZresults);
  else
    AllResults = theJZresults;
/*  check for Pleiotropy v Close linkage results and process them. */
  sprintf(infile, "%sp", theparams->zfile);
  indfile =   isfile(infile);
  if ( indfile == 1 ) {
    thePCLresults = AnalysisFirstPass(infile);
    for ( cnode=thePCLresults; cnode!=NULL; cnode=cnode->next)
      SetAnalysisHeaders(cnode,theparams);
    AnalysisSecondPass(thePCLresults);  
  }

/* Append if they exist.  */
  if ( AllResults != NULL )
    AppendResults(AllResults, thePCLresults);
  else
    AllResults = thePCLresults;

/*  print it all out. */
  if ( AllResults != NULL ) {
    out = fileopen( theparams->eqtl, "a");
    for ( cnode=AllResults; cnode!=NULL; cnode=cnode->next )
      ShowAnalysisResults(out, cnode, theparams);
    fileclose(theparams->eqtl, out);
  }
  
/*  clean up. */  
  if ( AllResults != NULL )
    UnAllocAnalysisResults(AllResults);


}



/*
Process the JZmapqtl output files. 
*/

void process_jzoutput( params *theparams,char *bootin,   genome *first) {
  int kk;
/*  if ( theparams->lodflag == 1 ) 
    theparams->siglevel = lodtolr(theparams->siglevel);*/

  if ( theparams->traits > 1 ) {
/*  If JZmapqtl was run, then process its results. */
    sprintf(bootin,"%s.z0",theparams->stem);
    kk = isfile(bootin);
    if ( kk == 1 ) {
      kk = process_jzfiles(theparams,bootin,first);
      
      if ( theparams->ihypo > 30 && theparams->ihypo < 34 && kk != -1) {
        sprintf(bootin,"%s.z0",theparams->stem);
        kk = process_jzfiles(theparams,bootin,first);
        sprintf(bootin,"%s.z0",theparams->stem);
        kk = process_jzfiles(theparams,bootin,first);      
      }
    }
  }
/*  if ( theparams->lodflag == 1 )
    theparams->siglevel = lrtolod(theparams->siglevel);*/
}

/*
Process the Zmapqtl output file. 

What to display?  For ngt = 3, 
---------------------------------
---------------------------------
Hypo    H     r2/tr2   S   a   d
---------------------------------
10    H1:H0    1:0     1   1   -
20    H2:H0    2:0     2   -   2
30    H3:H0    3:0     3   3   3
31    H3:H1    3:0     3   3   3
32    H3:H2    3:0     3   3   3
---------------------------------

*/
void process_zoutput(params *theparams,markermap *themap, char *progname, char *chptr) {
  int ii,j,k;
  FPN slevel;
  zresult *tzres,*zreslist;
  FILE *outf;
  outf = fileopen(theparams->eqtl, "a");
  fprintf(outf,"\n# The following is for hypothesis test %d and Zmapqtl model %d",theparams->ihypo,theparams->Model);  
  fileclose(theparams->eqtl, outf);
/* 
  if ( theparams->lodflag == 1 ) 
    theparams->siglevel = lodtolr(theparams->siglevel);
*/
  if (theparams->zfile[0] != '\0') {
    zreslist = zmapqtl_list(theparams,231);  
    zreslist = zresult_elim_nonmodel(zreslist,theparams);     
    for (ii = 1; ii <= theparams->traits; ii++) 
      themap->knum[ii] = 0;
    if ( zreslist != NULL ) {
	    for ( tzres = zreslist ; tzres != NULL ; tzres = tzres->next ) 
	      themap->knum[tzres->trait] = themap->knum[tzres->trait] + znode_QTLs(tzres,theparams);
	    if ( theparams->verbosity == 1 ) {
	      if ( theparams->lodflag == 1 ) {
	        slevel = lrtolod(theparams->siglevel);
	        printf("\n\nSignificance level for LOD score is %f",slevel);
	      }
	      else {
	        slevel = theparams->siglevel;
	        printf("\n\nSignificance level for LR statistic is %f",slevel);
	      }
	      for (ii = 1; ii <= theparams->traits; ii++) 
	        printf("\nNumber of QTLs for Trait %d: %d", ii,themap->knum[ii]);
	      printf("\n");
	    }
	    theparams->theqtls  = qtlvector(themap);
	    j=1;	
	    for ( tzres = zreslist ; tzres != NULL ; tzres = tzres->next ) 
	      j = znode_QTLs_estimate(tzres,theparams,theparams->theqtls,themap,j);
	    if ( theparams->lodflag == 1 )
	      print_aqtl(theparams, theparams->eqtl, 2);
	    else
	      print_aqtl(theparams, theparams->eqtl, 1);
	    calc_recomb_dist(theparams->theqtls,themap,j-1);
	    print_aqtl(theparams, theparams->eqtl, 0);
	    
	    if ( themap->ParentalDiff != NULL ) {
          outf = fileopen(theparams->eqtl, "a");
	      strcat(gname, "y\0");
          k = OttoJones(theparams ,outf, themap->ParentalDiff, gname, NULL );
          fileclose(theparams->eqtl, outf);
        }
        print_head(progname,theparams->srfile,chptr,1,51,theparams);
        write_zmapqtl_ranks(theparams,theparams->theqtls,theparams->srfile);

	    zresult_list_abolish(zreslist);
	    free_qtlvector(theparams->theqtls,themap);
	    theparams->theqtls = NULL; 
    }
  }
/*  if ( theparams->lodflag == 1 )
    theparams->siglevel = lrtolod(theparams->siglevel);*/
    
}


void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag)
{
  int ii, jj;
  if (flag == 1) {
    strcpy(opt[1], "-z");
    strcpy(opt[2], "-e");
    strcpy(opt[3], "-o");
    strcpy(opt[4], "-m");
    strcpy(opt[5], "-s");
    strcpy(opt[6], "-M");
    strcpy(opt[7], "-H");
    strcpy(opt[8], "-a");
    strcpy(opt[9], "-S");
    strcpy(opt[10], "-L");
    strcpy(opt[11], "-I");

    strcpy(opt_e[1], "(Composite) Interval Mapping Results");
    strcpy(opt_e[2], "Error File");
    strcpy(opt_e[3], "Output File");
    strcpy(opt_e[4], "Genetic Linkage Map File");
    strcpy(opt_e[5], "Random Number Seed");
    strcpy(opt_e[6], "Model from Zmapqtl/JZmapqtl");
    strcpy(opt_e[7], "Hypothesis Test (10,14,30,31,32,34)");

    strcpy(opt_e[8], "Significance level (alpha)");
    strcpy(opt_e[9], "Significance Threshold");
    strcpy(opt_e[10], "Output LOD scores? (0=no,1=yes)");
    strcpy(opt_e[11], "Work code ");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->zfile);
  strcpy(opt_v[2], theparams->error);
  strcpy(opt_v[3], theparams->eqtl);
  strcpy(opt_v[4], theparams->map);
  sprintf(opt_v[5], "%ld", theparams->seed);
  sprintf(opt_v[6], "%d", theparams->Model);
  sprintf(opt_v[7], "%d", theparams->ihypo);
  sprintf(opt_v[8], "%f", theparams->size);
  sprintf(opt_v[9], "%f", theparams->siglevel);
  sprintf(opt_v[10], "%d", theparams->lodflag);
  sprintf(opt_v[11], "%s", theparams->mimwork);

}

void update_params(char **opt_v,  params *theparams)
{
  strcpy(theparams->zfile, opt_v[1]);
  strcpy(theparams->error, opt_v[2]);
  strcpy(theparams->eqtl, opt_v[3]);
  strcpy(theparams->map, opt_v[4]);
  theparams->seed = atol(opt_v[5]);
  theparams->Model = atoi(opt_v[6]);
  theparams->ihypo = atoi(opt_v[7]);
  theparams->size = (FPN) atof(opt_v[8]);
  theparams->siglevel = (FPN) atof(opt_v[9]);
  if (theparams->siglevel < 0.0)
    theparams->siglevel = (FPN) SIG_LEVEL;
  theparams->lodflag = atoi(opt_v[10]);
  strcpy(theparams->mimwork,opt_v[11]);
}



  
/*  
  if (theparams->zfile[0] != '\0') {
    zreslist = zmapqtl_list(theparams,231);  
    zreslist = zresult_elim_nonmodel(zreslist,theparams);     
    for (ii = 1; ii <= theparams->traits; ii++) 
      themap->knum[ii] = 0;
    if ( zreslist != NULL ) {
	    for ( tzres = zreslist ; tzres != NULL ; tzres = tzres->next ) 
	      themap->knum[tzres->trait] = themap->knum[tzres->trait] + znode_QTLs(tzres,theparams);
	    if ( theparams->verbosity == 1 ) {
	      if ( theparams->lodflag == 1 ) {
	        slevel = lrtolod(theparams->siglevel);
	        printf("\n\nSignificance level for LOD score is %f",slevel);
	      }
	      else {
	        slevel = theparams->siglevel;
	        printf("\n\nSignificance level for LR statistic is %f",slevel);
	      }
	      for (ii = 1; ii <= theparams->traits; ii++) 
	        printf("\nNumber of QTLs for Trait %d: %d", ii,themap->knum[ii]);
	      printf("\n");
	    }
	    theqtls = qtlvector(themap);
	    j=1;	
	    for ( tzres = zreslist ; tzres != NULL ; tzres = tzres->next ) 
	      j = znode_QTLs_estimate(tzres,theparams,theqtls,themap,j);
	    if ( theparams->lodflag == 1 )
	      print_aqtl(theqtls, theparams->eqtl, 2);
	    else
	      print_aqtl(theqtls, theparams->eqtl, 1);
	    calc_recomb_dist(theqtls,themap,j-1);
	    print_aqtl(theqtls, theparams->eqtl, 0);
        
        print_head(argv[0],theparams->srfile,chptr,1,51,theparams);
        write_zmapqtl_ranks(theparams,theqtls,theparams->srfile);

	    zresult_list_abolish(zreslist);
	    free_qtlvector(theqtls,themap);
    }
  }*/
/*  if ( theparams->traits > 1 ) {
    sprintf(bootin,"%s.z0",theparams->stem);
    kk = isfile(bootin);
    if ( kk == 1 ) {
      kk = process_jzfiles(theparams,bootin,argv[0],chptr,1);
      if ( theparams->ihypo > 30 && theparams->ihypo < 34 && kk != -1) {
        sprintf(bootin,"%s.z0",theparams->stem);
        kk = process_jzfiles(theparams,bootin,argv[0],chptr,2);
        sprintf(bootin,"%s.z0",theparams->stem);
        kk = process_jzfiles(theparams,bootin,argv[0],chptr,3);      
      }
    }
  } 
  if ( theparams->lodflag == 1 )
    theparams->siglevel = lrtolod(theparams->siglevel);
*/ 

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file EQmain.c
------------------------------------------------------------------ */


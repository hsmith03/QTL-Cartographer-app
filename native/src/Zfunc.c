/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Zfunc.c) is part of QTL Cartographer
         
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

/*  Subroutines to implement composite interval mapping. */

#include "Main.h"

void showprgenotypes(params *theparams, linpakws *lnpk);

/*
Determine how many bootstraps have already been done.  Read tokens until 
you find  a -boots token.  The token following will be the number of bootstraps already done.
*/
int get_bootnum(char *infile)
{
  int boots;
  FILE *inputf ;
  int ch,start;
  inputf = fileopen(infile,"r");
  ch = 10;
  start = 0;
  do {
    ch = get_next_token(gname, MAXNAME, inputf);
    if ( ch == EOF ) {
      fileclose(infile, inputf);  
      return(-2);
    }
    else if ( !strcmp(gname,"-boots") ) {
        start = 1;
        ch = get_next_token(gname, MAXNAME, inputf);
        boots = atoi(gname);
    }
  } while ( start != 1 && ch != EOF );  /*We should be at the right place... */
  fileclose(infile, inputf);  
  if ( ch == EOF ) 
      return(-3); 
  return(boots);
}

/* 
   This gets the current state of the comparisonwise test results
   for a permutation test.

Fetch the interim comparison wise test results.

Read tokens until you get to a -start token.  
At a -perms token, read the next token and convert to an integer.
It will be the number of permutations already done.

Then read the current state of the comparison wise test results.  

*/
int get_CWTR(params *theparams, char *infile, int *pcnts, FPN *lratio)
{
  FILE *inputf;
  int  ch,  ii, trait,model,start,cntr;
  model = trait = -1;
  ch = 10;
  start = 0;
  inputf = fileopen(infile, "r");
  if ( inputf == NULL )
    return(-1);

  do {
    ch = get_next_token(gname, MAXNAME, inputf);
    if ( ch == EOF ) {
      fileclose(infile, inputf);  
      return(-2);
    }
    else if ( !strcmp(gname,"-rwd") ) {
      theparams->rwd = 1;
    }
    else if ( !strcmp(gname,"-perm") ) {
      ch = get_next_token(gname, MAXNAME, inputf);
      cntr = atoi(gname);
    }
    else if ( !strcmp(gname,"-start") ) 
      start = 1;
  } while ( start != 1 && ch != EOF );

  ch = get_next_token(gname, MAXNAME, inputf); /*row*/
  do {
    ii = atoi(gname);
    ch = get_next_token(gname, MAXNAME, inputf);/*chrom*/
    ch = get_next_token(gname, MAXNAME, inputf);/*mark*/
    ch = get_next_token(gname, MAXNAME, inputf);/*pos*/
    ch = get_next_token(gname, MAXNAME, inputf);/*orig*/
    lratio[ii] = (FPN) atof(gname);
    ch = get_next_token(gname, MAXNAME, inputf);/*p val*/
    ch = get_next_token(gname, MAXNAME, inputf);/*pcnts*/
    pcnts[ii] = atoi(gname);    
    ch = get_next_token(gname, MAXNAME, inputf);/*row or -end */
    if ( !strcmp(gname,"-end") ) 
      start = 0;
  } while ( start == 1 && ch != EOF );
  fileclose(infile, inputf);  
  return(cntr);
}

/*
Fetch the likelihood ratio from the Zmapqtl.out file.  
Since the Zmapqtl.out file may have many sets of results, look for a group
that have the same Model and trait numbers.  Then find the -s token.
The results follow the -s token and stop at the -e token.

Also, we want to read the proper columns.  This will depend on the value of 
theparams->ihypo.  

*/
int get_lratio(params *theparams, FPN *lratio)
{
  FILE *inputf;
  int  thecol,row,col,ch,trait,model,back,start;
  
   thecol = 4;	/* LR H1 v H0 for 2 gt,   H3 v H0 for 1 gt*/
  /* determine the columns to process...*/
  if ( (theparams->cross == 3 || theparams->cross == 4) && (theparams->tcross != 1 || theparams->tcross != 2) ) {
    if (theparams->ihypo == 32  || theparams->ihypo == 3 )  
      thecol = 6;	/* LR H3 v H2 */
    else if (theparams->ihypo == 31 || theparams->ihypo == 2)  
      thecol = 5;	/* LR H3 v H1 */
  }

  back = model = trait = -1;
  ch = 10;
  start = 0;
  inputf = fileopen(theparams->zfile, "r");
  if ( inputf == NULL )
    return(-1);
  do {
    ch = get_next_token(gname, MAXNAME, inputf);
    if ( ch == EOF ) {
      fileclose(theparams->zfile, inputf);  
      return(-2);
    }
    else if ( !strcmp(gname,"-Model") ) {
      ch = get_next_token(gname, MAXNAME, inputf);
      model = atoi(gname);
    }
    else if ( !strcmp(gname,"-background") ) {
      ch = get_next_token(gname, MAXNAME, inputf);
      back = atoi(gname);
    }
    else if ( !strcmp(gname,"-trait") ) {
      ch = get_next_token(gname, MAXNAME, inputf);
      trait = atoi(gname);
    }
    else if ( !strcmp(gname,"-rwd") ) {
      theparams->rwd = 1;
    }
    else if ( !strcmp(gname,"-s") ) {
      if ( model == theparams->Model  && trait == theparams->whichtrait )
        start = 1;
    }
  } while ( start != 1 && ch != EOF );  /*We should be at the right place... */
  if ( ch == EOF ) {
      fileclose(theparams->zfile, inputf);  
      return(-3);
  }
  row = col = 1;
  
  do {
    ch = get_next_token(gname, MAXNAME, inputf);
    if ( col == thecol )
      lratio[row] = (FPN) atof(gname);
    else if ( col == 21) {
      row = row+1;
      col = 0;
    }
    if ( !strcmp(gname,"-e") ) 
      start = 0;
    col = col+1;
  } while ( start == 1 && ch != EOF );
  fileclose(theparams->zfile, inputf);  
  return(row);
}



/*
  Read in the results of the program LRmapqtl to be used
  in picking markers for a background fit.
  
  This routine is fast becoming obsolete.  
*/
FPN **get_lin_reg(char *lrfile, int chroms, int maxl, int t)
{
  int  ch, negctr, row, col, ii, attrait;
  FILE *infile;
  FPN **lin_reg;

  lin_reg = dmatrix(1, chroms, 1, maxl);
  negctr = 0;
  infile = fileopen(lrfile, "r");
  attrait = 0;
  while (attrait == 0) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if (gbuffer[0] == '-' && gbuffer[1] == 't') {
      get_field(2, gname, gbuffer);
      attrait = atoi(gname);
      if (attrait != t)
	attrait = 0;
    }
  }

  while (negctr < 4) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if ((negctr == 3) && (gbuffer[0] != '-')) {
      get_field(1, gname, gbuffer);
      row = atoi(gname);
      get_field(2, gname, gbuffer);
      col = atoi(gname);
      get_field(5, gname, gbuffer);
      lin_reg[row][col] = (FPN) atof(gname);
    }
    if (gbuffer[0] == '-')
      negctr += 1;
    if (ch == EOF)
      negctr = 5;
  }
  fileclose(lrfile, infile);

  return (lin_reg);
}

/*  put phenotypes in yy.  Use only those phenotypes that are > (FPN) MISS_VAL 
    do this for theparams->whichtrait*/
void init_phenotypes(params *theparams,individual *individs,FPN **yy,int *samplesize)
{
  FPN sum,sum2;
  int jj,kk,ii;

  sum=sum2=(FPN) 0.0;
/*  for (jj = 1; jj <= theparams->traits; jj++) {*/
     jj = theparams->whichtrait;
	 kk = 1;
	 for (ii = 1; ii <= theparams->nn; ii++)
		if ( individs[ii].y[jj] > (FPN) MISS_VAL) {
	      yy[jj][kk] = individs[ii].y[jj];
	      kk = kk + 1;
	      if ( jj == theparams->whichtrait ) {
	        sum = sum + individs[ii].y[jj];
	        sum2 = sum2 + individs[ii].y[jj] * individs[ii].y[jj];
	      }
	     individs[ii].print_flag = 'y';
        }
        else if ( jj == theparams->whichtrait) 
          individs[ii].print_flag = 'n';   /*Important:  print_flag is 'y' by default.  
                                                 Setting it to 'n' takes this individual 
                                                 out of the analysis. */
     samplesize[jj] = kk - 1;
/*  }*/
  theparams->total_var =   (sum2 - sum*sum/(FPN) samplesize[theparams->whichtrait])/ (FPN)(samplesize[theparams->whichtrait]);
/* SAMPLINGVARIANCE Should we use the sampling variance instead?  No, for simplicity, use n in the denominator throughout.  For Proper 
estimation of variance components, re-estimate everything in MImapqtl.    */
}


/*
Go into the SRmapqtl results file and get the ranks of the 
markers.

*/
int get_srresults(params *theparams, markermap *themap, genome *first)
{
  FILE *infile;
  int  ch, trait,model,start,mark,chrom,rank,cntr,k;
  genome *tptr;
  k=themap->traits;
  model = trait = -1;
  k=0;
  ch = 10;
  start = 0;
  infile = fileopen(theparams->srfile, "r");
  if ( infile == NULL )
    return(-1);


  for ( tptr = first ; tptr != NULL ; tptr = tptr->next )
    tptr->whichqtl = 0;
  
  do {
    ch = get_next_token(gname, MAXNAME, infile);
    if ( ch == EOF ) {
      fileclose(theparams->srfile, infile);  
      return(-2);
    }
    else if ( gname[0] == '-' ) {
      if ( gname[2] == 'B' ) /*tells model type */
        model =2;
      else if  ( gname[2] == 'o' )
        model = 0;
      else if ( gname[2] == 'a' ) 
        model = 1;
      else if ( gname[2] == 'm' ) 
        model = 3;
      else if ( gname[2] == 'r'  ) { /*tells trait */
        ch = get_next_token(gname, MAXNAME, infile);
        trait = atoi(gname);
      }
      else if ( gname[2] == 't' ) { /*tells start */
        if ( model == theparams->srm && trait == theparams->whichtrait )
          start = 1;
      }    
    }  
  } while ( start != 1 && ch != EOF );
  cntr = 0;
  do {
    ch = get_next_token(gname, MAXNAME, infile);
    if ( ch != EOF && gname[0] != '-' ) {
      chrom = atoi(gname);
      ch = get_next_token(gname, MAXNAME, infile);
      mark = atoi(gname);
      ch = get_next_token(gname, MAXNAME, infile);
      rank = atoi(gname);
      ch = get_next_token(gname, MAXNAME, infile);/*This would be the Fstatistic*/
      ch = get_next_token(gname, MAXNAME, infile);/*This would be the DOF*/
      for ( tptr = first ; tptr != NULL && !( tptr->chrom == chrom && tptr->markr == mark ) ; tptr = tptr->next ) k+=1;
      if ( tptr != NULL ) {
          tptr->whichqtl = rank;
          cntr = cntr+1;
      }
    }
    else
      start = 0;
  } while ( ch != EOF && start == 1 );
  
  fileclose(theparams->srfile, infile);  

    if ( cntr < 1 )
      theparams->Model = 3;
    else if ( cntr < theparams->nbp )
      theparams->nbp = cntr;
     
    if ( theparams->verbosity == 1) {
      if ( cntr == -1 )
        printf("\nNo SRmapqtl output file found...Will do interval mapping\n");
      else if ( cntr == -2 )
        printf("\nSRmapqtl file has no results...Will do interval mapping\n");
      else if ( cntr == 0) 
        printf("\nSRmapqtl file has no ranking markers...Will do interval mapping\n");
      else if ( cntr > 0) 
        printf("\nWe could use up to %d markers, but will use %d.\n",cntr,theparams->nbp);
    }
  return(cntr);
}


/*
  This is used for bootstrapping.  
  Bootstrapping requires two extra files:  The old bootstrap and the
  new one.  The old one holds the current state of the bootstrap experiment.
  The new one will hold the state after the current analysis.  If This is the 
  first repetition of the bootstrap, then there will be no old bootstrap file.
  This routine simply opens the zfile, finds the proper block of results,
  and then reads the chromosome, marker and position colums to be printed into
  a starter old bootstrap file.  The values of the lr, lr2, etc, are all set to 
  0.0.
*/
int convert_zfile(params *theparams, char *outfile)
{
  FILE *inputf,*outputf;
  int col,back,model,trait,ch,start;
  
  
  inputf = fileopen(theparams->zfile,"r");
  outputf = fileopen(outfile,"a");
  
  back = model = trait = -1;
  ch = 10;
  start = 0;
  do {
    ch = get_next_token(gname, MAXNAME, inputf);
    if ( ch == EOF ) {
      fileclose(theparams->zfile, inputf);  
      return(-2);
    }
    else if ( !strcmp(gname,"-Model") ) {
      ch = get_next_token(gname, MAXNAME, inputf);
      model = atoi(gname);
    }
    else if ( !strcmp(gname,"-trait") ) {
      ch = get_next_token(gname, MAXNAME, inputf);
      trait = atoi(gname);
    }
    else if ( !strcmp(gname,"-s") ) {
      if ( model == theparams->Model  && trait == theparams->whichtrait )
        start = 1;
    }
  } while ( start != 1 && ch != EOF );  /*We should be at the right place... */
  if ( ch == EOF ) {
      fileclose(theparams->zfile, inputf);  
      return(-3);
  }
  col = 1;

  do {
    ch = get_next_token(gname, MAXNAME, inputf);
    if ( !strcmp(gname,"-e") ) 
      start = 0;
    else {
	    if ( col == 1 )
	      fprintf(outputf,"\n %5s ",gname);
	    else if ( col == 2 )  
	      fprintf(outputf," %5s ",gname);
	    else if ( col == 3 ) 
	      fprintf(outputf," %10s   0.0   0.0   0.0   0.0   0.0   0.0",gname);
	    else if ( col == 21 )  
	      col = 0;
	}
    col = col+1;
  } while ( start == 1 && ch != EOF );
  fprintf(outputf,"   -end\n " );
  fileclose(theparams->zfile, inputf);
  fileclose(outfile, outputf);
  return(0);
}



/* 
	eliminate node gptr, returning gptr->next if it exists, or gptr->prev.
	return null if it is the only node left;
*/
genome *delete_gnode(genome *gptr)
{
  genome *rptr;
  if ( gptr == NULL )
    return(NULL);

  if ( gptr->prev != NULL && gptr->next != NULL ) {
    rptr = gptr->next;
    rptr->prev = gptr->prev;
    gptr->prev->next = rptr;
  }
  else if ( gptr->prev != NULL && gptr->next == NULL ) {
    rptr = gptr->prev;
    rptr->next = NULL;
  }
  else if ( gptr->prev == NULL && gptr->next != NULL ) {
    rptr = gptr->next;
    rptr->prev = NULL;
  }
  else if ( gptr->prev == NULL && gptr->next == NULL )
    rptr = NULL;
 
  free((char *) gptr);
  return(rptr);
}

/*
  Output the positional results for the current bootstrap rep.
  Need to read the values from the old bootstrap file, and 
  update them with the current results.
*/
long  write_bootline(params *theparams, char *outfile, char *infile, genome *gptr, FPN abspos, FPN *estimates, long int inset)
{
  FILE *outf,*inputf;
  int lri,addi,domi,model,trait,start,ch;
  FPN lr,lr2,add,add2,dom,dom2;
  long offset;
  
  lr=lr2=add=add2=dom=dom2=(FPN) 0.0;
  if ( theparams->cross == 3 ||  theparams->cross == 4 ) {
    lri = 1;
    addi = 5;
    domi = 7;
  }/*  Modify above to take into account ihypo. */
  else {
    lri = 1;
    addi = 4;
    domi = 1;
  }
  
  inputf = fileopen(infile,"r");
  if ( inset == 0L ) {  /* First time to access this file */
	  model = trait = -1;
	  start = 0;
	  do {
	    ch = get_next_token(gname, MAXNAME, inputf);
	    if ( ch == EOF ) {
	      fileclose(infile, inputf);  
	      return(-2L);
	    }
	    else if ( !strcmp(gname,"-Model") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      model = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-trait") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      trait = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-start") ) {
	      if ( model == theparams->Model  && trait == theparams->whichtrait )
	        start = 1;
	    }
	  } while ( start != 1 && ch != EOF );  /*We should be at the right place... */
	  
  
  }
  else
    fseek(inputf, inset, 0L);
    
  ch = get_next_token(gname, MAXNAME, inputf);/*c*/
  ch = get_next_token(gname, MAXNAME, inputf);/*m*/
  ch = get_next_token(gname, MAXNAME, inputf);/*p*/
  ch = get_next_token(gname, MAXNAME, inputf);/*lr*/
  lr = (FPN) atof(gname) + estimates[lri];
  ch = get_next_token(gname, MAXNAME, inputf);/*lr2*/
  lr2 = (FPN) atof(gname) + estimates[lri] * estimates[lri];
  ch = get_next_token(gname, MAXNAME, inputf);/*a*/
  add = (FPN) atof(gname)  + estimates[addi];
  ch = get_next_token(gname, MAXNAME, inputf);/*a2*/
  add2 = (FPN) atof(gname)  + estimates[addi]* estimates[addi];
  ch = get_next_token(gname, MAXNAME, inputf);/*d*/
  dom = (FPN) atof(gname)   + estimates[domi] ;
  ch = get_next_token(gname, MAXNAME, inputf);/*d2*/
  dom2 = (FPN) atof(gname)   + estimates[domi]  * estimates[domi];
  offset = ftell(inputf);
  fileclose(infile, inputf);
  
  
  
  outf = fileopen(outfile, "a");
  if ( theparams->cross == 3 ||  theparams->cross == 4 )
    fprintf(outf, "\n%2d %2d %7.4f %20.5f %20.5f %20.15f %20.15f %20.15f %20.15f", gptr->chrom, gptr->markr, abspos,lr,lr2,add,add2,dom,dom2);
  else
    fprintf(outf, "\n%2d %2d %7.4f %20.5f %20.5f %20.15f %20.15f  0.0 0.0", gptr->chrom, gptr->markr, abspos,lr,lr2,add,add2);
  
  fileclose(outfile, outf);

  return(offset);
}

/*

  Write the current results of the analysis.

Out              Fx        Bi,Ri
Col  estimate  Quant       Quant.
 1              c           c
 2              m           m
 3              position    position
 4   1          H0:H3       H0:H1
 5   2          H1:H3       R2(0:1)
 6   3          H2:H3       TR2(0:1)
 7   4          H1:a        H1:a 
 8   5          H3:a        S
 9   6          H2:d       
10   7          H3:d      
11   8          H0:H1      
12   9          H0:H2
13  10          R2(0:1)      
14  11          R2(0:2)      
15  12          R2(0:3)      
16  13          TR2(0:1)    
17  14          TR2(0:2)     
18  15          TR2(0:3)     
19  16          S1   
20  17          S2     
21  18          S3     


*/
void write_zline(params *theparams, char *outfile, genome *gptr, FPN abspos, FPN *estimates)
{
  FILE *outf;
  int ii;
  ii=theparams->traits;
  outf = fileopen(outfile, "a");
  fprintf(outf, "\n%2d %2d %7.4f ", gptr->chrom, gptr->markr, abspos);
  for (ii = 1; ii <= 18; ii++)
    fprintf(outf, "%14.7f ", estimates[ii]);
  fileclose(outfile, outf);
}
/*
  Write out the current state of the comparisonwise test results.
  Col.    Quant
  1       row
  2       c
  3       m
  4       p
  5       lr(sample)
  6       pvalue
  7       pcnts
*/
void write_CWTR(params *theparams,genome *startptr, genome *endptr, markermap *themap, FPN deltax, char *outfile, FPN *lratio, int *pcnts, int reps)
{
  FPN abspos, relpos, pvalue;
  int  ii;
  genome *gptr;
  FILE *outf;
  abspos = relpos = (FPN) MIN_DIST;
  gptr = startptr;
  ii = 0;
  	outf = fileopen(outfile, "a");
    if ( theparams->rwd == 1) 
    	fprintf(outf, "     -rwd " );

	fprintf(outf, " -perm %d   -start ",reps );
	fileclose(outfile, outf);

  while (gptr <= endptr && gptr != NULL) {
    if ((gptr->markr >= 0 && gptr->markr <= themap->mpc[gptr->chrom]) && gptr->dist > (FPN) 0.0) {
      while (relpos < gptr->dist) {
	ii = ii + 1;
	pvalue = (FPN) pcnts[ii] / (FPN) reps;
	outf = fileopen(outfile, "a");
	fprintf(outf, "\n%5d %3d %3d %8.5f %10.5f %10.6f %10d",ii, gptr->chrom, gptr->markr, abspos, lratio[ii], pvalue, pcnts[ii]);
	fileclose(outfile, outf);
	relpos = relpos + deltax;
	abspos = abspos + deltax;
      }
      abspos = abspos - relpos + gptr->dist + (FPN) MIN_DIST;
      relpos = (FPN) MIN_DIST;
    }
    gptr = gptr->next;
    if (gptr != NULL)
      if (gptr->chrom != gptr->prev->chrom)
	abspos = relpos = (FPN) MIN_DIST;
  }
  
	outf = fileopen(outfile, "a");
	fprintf(outf, "    -end\n" );
	fileclose(outfile, outf);
}

/*
  Write header for the Zmapqtl output files.  
           0    Zmapqtl.out
  eorc =   1    ZpermC.out
           2    ZpermE.out
           3    Zboot.out
           
           
   boots is only used in the eorc==3 case.
*/         
void write_zheader(char *outfile, params *theparams, char *chptr, int eorc, int boots, markermap *themap)
{
  int cross,ii;
  FILE *outf;
  if (themap->traits < 0 )
    printf("%s",chptr);
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  outf = fileopen(outfile, "a");
  if (outf == NULL)
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
    fprintf(outf, "\n#The position is from the left telomere on the chromosome");
    fprintf(outf, "\n-window     %6.2f      Window size for models 5 and 6",theparams->window);
    fprintf(outf, "\n-background %6d      Background parameters in model 6",theparams->nbp);
    fprintf(outf, "\n-Model      %6d      Model number",theparams->Model);
    if ( themap->tnames != NULL )
      fprintf(outf, "\n-trait      %6d      Analyzed trait [%s]", theparams->whichtrait,themap->tnames[theparams->whichtrait]);
    else
      fprintf(outf, "\n-trait      %6d      Trait analyzed", theparams->whichtrait);
    fprintf(outf, "\n-cross      %6s      Cross", theparams->thecross);
    fprintf(outf, "\n#\n#  Note that our Likelihood ratio test statistic compares two nested hypotheses");
    fprintf(outf, "\n#  and is two times the negative natural log of the ratio of the likelihoods.  For example,");
    fprintf(outf, "\n#  assume that  hypothesis H0 is nested within H1 and that they have likelihoods L0 and L1 respectively.");
    fprintf(outf, "\n#  Then, the \"Likelihood Ratio Test Statistic\" is -2ln(L0/L1). \n#");
    if ( eorc == 0 ) {
      fprintf(outf, "\n#  Test Site   *");
      fprintf(outf, " Like. Ratio Test Statistics   *     Additive       *     Dominance        * Misc. HT");
      fprintf(outf, "\n c  m position");
      if ( cross == 1 || cross == 2 || cross == 5 )
        fprintf(outf, "     H0:H1          R2(0:1)        TR2(0:1)        H1:a            S1            ....           ....          ....           ....           ....           ....          ....             .....           ....          .....         .....           ....          .....\n-s");
      else 
        fprintf(outf, "     H0:H3          H1:H3          H2:H3           H1:a           H3:a           H2:d           H3:d          H0:H1          H0:H2         R2(0:1)        R2(0:2)        R2(0:3)         TR2(0:1)       TR2(0:2)       TR2(0:3)         S1             S2             S3 \n-s");
    }   
    else if ( eorc == 1 ) {
      fprintf(outf, "\n#Position is from the left telomere on the chromosome");
      fprintf(outf, "\n#Row Chrom Mark Position  LR(Samp)    P-Val     Count     ");
    }
    else if ( eorc == 2)  {
      if ( theparams->rwd == 1 ) 
        fprintf(outf, "\n#                -rwd          Chromosome");
      fprintf(outf, "\n#  Repetition    GlobalMax ");
      if ( theparams->rwd == 1 ) 
        for ( ii=1; ii<=theparams->chrom; ii++ )
          fprintf(outf, "      %3d      ",ii);
      fprintf(outf, " -start  ");
    }
    else if (eorc == 3) {
      fprintf(outf, "\n#Position is from the left telomere on the chromosome");
      fprintf(outf, "\n#Chrom   Mark   Position    LR     LR2    A    A2     D    D2  -boots %d  -start",boots);
    }
    fileclose(outfile, outf);
  }
}

/*
  This is used for the Permutation test.  It simply adds one line to the
  ZpermE.out file indicating the maximal lr for that permutation.

*/
void write_maxlr(char *outfile, int jj, FPN *maxlr,int chroms)
{
  FILE *outf;
  int ii;
  outf = fileopen(outfile, "a");
  if (outf == NULL)
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
    fprintf(outf, "\n%5d",jj);
    for( ii=0; ii<= chroms; ii++ )
      fprintf(outf, " %15.5f",maxlr[ii]);
    fileclose(outfile, outf);
  }
}

/*
Do the analysis, now that we have the starting and ending points. 
theparams,startptr,endptr,themap,theqtls,lin_reg,yy,individs,outfile,lratio,cnts,agptr,lnpk)
*/
FPN do_zanalysis(params *theparams, genome *startptr, genome *endptr, markermap *themap, aqtl *theqtls, FPN **lin_reg, FPN *yy, individual *individs, char *outfile, FPN *lratio, int *cnts, genome *agptr, linpakws *lnpk)
{
  FPN abspos, relpos,   ts0, thetaL, thetaR, *estimates,   maxe;
  int     ii, jj,last,roffset,savemark;
  long inset;
  genome *gptr,*tgptr;
  FILE *outf;
  inset = 0L;
  roffset = 0;
  for ( ii = 1 ; ii <= themap->otraits ; ii++ )
    if ( themap->onames[ii][0] == '+' )
      roffset = roffset + themap->otypes[ii] - 1;
  estimates = lnpk->estimates[1];
  abspos = relpos = (FPN) MIN_DIST;
  gptr = startptr;
  last = jj = 0;
  maxe = (FPN) 0.0;
  tgptr = create_genome(themap);
  if ( lnpk->maxlr != NULL )
    for ( ii=0; ii<= themap->m; ii++ )
      lnpk->maxlr[ii] = maxe;

  while ( gptr != NULL) {
    savemark = gptr->markr;
  	if ( theparams->verbosity == 1 ) {
  	  if ( gptr->prev == NULL ) 
	    fprintf(stdout,"\nChromosome %d, Marker ",gptr->chrom);      
	  else if ( gptr->prev != NULL && gptr->prev->chrom != gptr->chrom )
	     fprintf(stdout,"R\nChromosome %d, Marker ", gptr->chrom);
	  fflush(stdout);
	}
	if ((gptr->markr >= 0 && gptr->markr <= themap->mpc[gptr->chrom]) && gptr->dist > (FPN) 0.0) {
	   if ( theparams->verbosity == 1 ) {
		 if ( gptr->markr > 0 )
		   fprintf(stdout,"%d.",gptr->markr);
		 else
		   fprintf(stdout, "L.");
		fflush(stdout);
	   }
	   pick_markers(lnpk->bp,themap,theqtls,lin_reg,theparams,gptr,agptr);
       ts0 = fit_back_params(theparams,individs, themap,theqtls, lnpk->xx, yy, lnpk->qraux, lnpk->rsd[1], gptr->chrom, gptr->markr, lnpk->bp, lnpk->bp[1][1],tgptr);
       while (relpos < gptr->dist) {
	     jj = jj + 1;
	     if (gptr->markr == 0)
	       thetaL = -(FPN) 2.0;
	     else
	       thetaL = relpos;
	     if (gptr->markr == themap->mpc[gptr->chrom])
	       thetaR = -(FPN) 2.0;
	     else
	       thetaR = gptr->dist - relpos;
	    zmap(theparams, lnpk,  yy, individs, gptr, lnpk->bp[1][1] , thetaL, thetaR,  ts0,roffset);
	    if ( debugging > 1 && (thetaR < (FPN) 0.0 || thetaL < (FPN) 0.0) ) 
          showprgenotypes( theparams,  lnpk);

	    
	    if (estimates[1] > maxe ) {
	      maxe  =estimates[1];
	      if (lnpk->maxlr != NULL )
	        lnpk->maxlr[0] = maxe;
	    }
	    if ( lnpk->maxlr != NULL &&  estimates[1] > lnpk->maxlr[gptr->chrom])
	      lnpk->maxlr[gptr->chrom] = estimates[1];
	    if (outfile != NULL) {
	      if (lratio != NULL)
	        lratio[jj] = estimates[1];
	    if (theparams->boots == 0 ) 
	      write_zline(theparams,outfile,gptr,abspos,estimates);  
	    else
	      inset = write_bootline(theparams,outfile,theparams->tfile,gptr,abspos,estimates,inset); 
	  }
	  else {
	    if ( estimates[1] >= lratio[jj])
	      cnts[jj] = cnts[jj] + 1;
	  }

	  relpos = relpos + theparams->walk;
	  abspos = abspos + theparams->walk;
    }
    abspos = abspos - relpos + gptr->dist + (FPN) MIN_DIST;
    relpos = (FPN) MIN_DIST;
  }
    gptr = gptr->next;
    if ( last == 1 )
      gptr = NULL;
    else if ( gptr == endptr ) 
      last = 1;
    if (gptr != NULL)
      if (gptr->chrom != gptr->prev->chrom)
	abspos = relpos = (FPN) MIN_DIST;
  }
  if (theparams->verbosity == 1 ) 
    fprintf(stdout,"R\n" );
  if (outfile != NULL) {
    outf = fileopen(outfile, "a");
    if (outf != NULL) {
      fprintf(outf, "\n-e\n");
      if (lratio != NULL)
	fprintf(outf, "\tComparisonwise p values for the shuffles\n c  m   pos      LR       p                n\n-b");
      fileclose(outfile, outf);
    }
  }
  clear_genome(tgptr);

  return (maxe);
}
/*
This drives the em algorithm.  
*/
void     zmap(params *theparams, linpakws *lnpk, FPN *yy, individual *individs, genome *gptr, int nbp, FPN thetaL, FPN thetaR, FPN ts0, int roffset)
{
  FPN ts3, *estimates ;
  int ii, ss,t,nn,cross,bail_status1,bail_status2, bail_status3 ;
  estimates = lnpk->estimates[1];
  t = theparams->whichtrait;
  nn = theparams->nn;
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;

  ss = 0;
  for (ii = 1; ii <= nn; ii++)
    if (individs[ii].print_flag == 'y') 
      ss = ss + 1;
  for (ii = 1; ii <= 18; ii++)
   estimates[ii] = (FPN) 0.0;
  

 
  calc_priors(theparams,lnpk->pp1,lnpk->pp2,nn,individs,gptr,thetaL,thetaR);

  

  if (cross == 1 || cross == 2 || cross == 5 ) {
    bail_status1 = em_solve( lnpk, yy, ss, nbp+roffset, 1,theparams);
    if ( bail_status1 == 0 ) {
      estimates[1] = (FPN) 2.0 * (estimates[1] - ts0);
      estimates[2] = (theparams->null_sse  - estimates[10])/ theparams->total_var ;
      estimates[3] = (theparams->total_var - estimates[10])/ theparams->total_var ;
      estimates[5] = lnpk->thestats[14];

    }
    else {
      estimates[1] = -(FPN) SIGNOFDEVIL;
      estimates[2] = (FPN) 0.0;
      estimates[3] = (FPN) 0.0;
    }
  }

  if (cross == 3 || cross == 4) {
    bail_status1 = em_solve( lnpk, yy, ss, nbp+roffset, 1,theparams);
    estimates[16] = lnpk->thestats[14];
    bail_status2 = em_solve( lnpk, yy, ss, nbp+roffset, 2,theparams);
    estimates[17] = lnpk->thestats[14];

    bail_status3 = em_solve( lnpk, yy, ss, nbp+roffset,  3,theparams);
    estimates[18] = lnpk->thestats[14];


    ts3 = estimates[3];

    if ( bail_status1 == 0 ) 
      estimates[8] = (FPN) 2.0 * (estimates[1] - ts0);
    else 
      estimates[8] = -(FPN) SIGNOFDEVIL;
    if ( bail_status2 == 0 ) 
      estimates[9] = (FPN) 2.0 * (estimates[2] - ts0);
    else 
      estimates[9] = -(FPN) SIGNOFDEVIL;    

    
    if ( bail_status3 == 0 && bail_status2 == 0  ) {
      estimates[3] = (FPN) 2.0 * (ts3 - estimates[2]);
      estimates[11] = (theparams->null_sse - estimates[11])/ theparams->total_var ;
      estimates[14] = (theparams->total_var - estimates[14])/ theparams->total_var ;
    }
    else {
       estimates[3] = -(FPN) SIGNOFDEVIL;
       estimates[11] = (FPN) 0.0;
       estimates[14] = (FPN) 0.0;
    }
    
    
    if ( bail_status3 == 0 && bail_status1 == 0  ) { 
      estimates[2] = (FPN) 2.0 * (ts3 - estimates[1]);
      estimates[10] = (theparams->null_sse - estimates[10])/ theparams->total_var ;
      estimates[13] = (theparams->total_var - estimates[13])/ theparams->total_var ;
    }
    else {
      estimates[2] = -(FPN) SIGNOFDEVIL;
      estimates[10] = (FPN) 0.0;
      estimates[13] = (FPN) 0.0;
    }
    
    
    
   if ( bail_status3 == 0 ) {
      estimates[1] = (FPN) 2.0 * (ts3 - ts0);
      estimates[12] = (theparams->null_sse - estimates[12])/ theparams->total_var ;
      estimates[15] = (theparams->total_var - estimates[15])/ theparams->total_var ;
    }
   else {
      estimates[1] = -(FPN) SIGNOFDEVIL;
      estimates[12] = (FPN) 0.0;
      estimates[15] = (FPN) 0.0;
    }
  }
}

/*
  Pick the background markers.
          1  => use all markers except lfm, rfm
model =   2  => use all unlinked markers (all except for those on chrom)
          3  => use no markers, just the mean (LB)
          4  => use subset of markers:  one marker per unlinked chromosome.
                                        Pick the one with the highest correlation.
          5  => use top two markers from each unlinked chromosome, and all linked markers outside
                the 'window' around the interval in question.
          6  =>   Specify two parameters:  ns and ws.  ns is the number of
                 background markers.  Pick the top ns of them as determined by a stepwise 
                 regression.  ws is a window size that blocks out markers on either side.
                 it is 10 cM by default.
          7  => Use the results of a previous scan with Zmapqtl to pick the background
                markers.  Eqtl should have been run, and the results summarized in the
                ->eqtl file will be used.  They should be in Rqtl.out format.
          8  => Use the results of a previous scan with Zmapqtl to pick the background
                markers.  Eqtl should have been run, and the results summarized in the
                ->srmapqtl file will be used.  They should be in SRmapqtl.out format.
                
*/
void pick_markers(int **bp, markermap *themap, aqtl *theqtls, FPN **lin_reg, params *theparams, genome *tgptr, genome *agptr)
{
  int col,   ii, jj,k, marks, tosub, thischrom;
  int model,chrom,lfm,ns;
  FPN  lpos, rpos,wind,temp;
  genome *gptr, *cgptr;
  k=0;
  model = theparams->Model;
  if ( model == 5 || model == 4 )
    if ( lin_reg == NULL )
      theparams->Model = model = 3;
  chrom = tgptr->chrom;
  lfm = tgptr->markr;
  ns = theparams->nbp;
  for ( ii = 1;ii <= themap->ml;ii++ )
    bp[1][ii]  = bp[2][ii]   = 0;
  
  if ( theparams->window > (FPN) 0.0 )
    wind = theparams->window * (FPN) 0.01;
  else
    wind = (FPN) WIN_SIZE * (FPN) 0.01;

  col = 1;
  if (model == 1) {	/* use all except lfm and rfm. */
    tosub = 1;
    if (lfm == 0 || lfm == themap->mpc[chrom])
      tosub = 0;

    bp[1][1]  = themap->ml - tosub;
    bp[2][1]  = 1;
    for (ii = 1; ii <= themap->m; ii++) {


      for (jj = 1; jj <= themap->mpc[ii]; jj++)
	if (!((ii == chrom && jj == lfm) || (ii == chrom && jj == lfm + 1))) {
	  col += 1;
	  bp[1][col]  = ii;
	  bp[2][col]   = jj;
	}
    }
  }
  else if (model == 2) {	/* use all except those on chrom */

      marks = themap->ml - themap->mpc[chrom] + 1;

    bp[1][1]   = marks;
    bp[2][1]   = 2;
    for (ii = 1; ii <= themap->m; ii++) {

      for (jj = 1; jj <= themap->mpc[ii]; jj++)
	if (ii != chrom) {
	  col += 1;
	  bp[1][col]   = ii;
	  bp[2][col]   = jj;
	}
    }
  }
  else if (model == 3) {	/* use only mean */
    bp[1][1]  = 1;
    bp[2][1]   = 3;
  }
  else if (model == 4) {	/* use the top marker on each other chromosome */
    marks = themap->m;
    bp[1][1]   = marks;
    bp[2][1]   = 4;
    for (ii = 1; ii <= themap->m; ii++) {
      if (ii != chrom) {
	if (ii < chrom)
	  col = ii + 1;
	else
	  col = ii;
	bp[1][col]   = ii;
	bp[2][col]   = isamax(themap->mpc[ii], lin_reg[ii], 1);
      }
    }
  }
  else if (model == 5  ) {
/* Use top two from each unlinked group, plus all outside window of  wind M. */
	 rpos = (FPN) MAXCHROMLEN;
	 lpos = (FPN) 0.0;
	 for ( gptr = agptr ; gptr != NULL ; gptr = gptr->next ) 
		if ( gptr->chrom == chrom && gptr->markr == lfm )
		  lpos = gptr->pos - wind;
		else if ( gptr->chrom == chrom && gptr->markr == lfm+1 )
		  rpos = gptr->pos + wind;
/* Use top two from each unlinked group, plus all outside window of  wind M. */
    for ( cgptr = agptr ; cgptr->chrom != chrom ; cgptr = cgptr->next ) k+=1;
    marks = 1;
    for ( gptr = cgptr ; gptr != NULL && gptr->chrom == chrom ; gptr = gptr->next )
      if ( (gptr->pos < lpos || gptr->pos > rpos) && gptr->markr > 0 ) {
        marks = marks+1;
        bp[1][marks]   = gptr->chrom;
        bp[1][marks]  = gptr->markr;
      }
    bp[1][1]   = marks + (themap->m-1)*2;
    bp[2][1]   = 5;
    thischrom = 0;
    for ( ii = 1 ; ii <= themap->m ; ii++ ) 
      if ( ii != chrom ) {
        bp[1][marks+1+2*thischrom] = bp[1][marks+2+2*thischrom] = ii;
        bp[2][marks+1+2*thischrom] = isamax(themap->mpc[ii], lin_reg[ii], 1);
        temp = lin_reg[ii][ bp[2][marks+1+2*thischrom] ];
        lin_reg[ii][ bp[2][marks+1+2*thischrom] ] = (FPN) 0.0;
        bp[2][marks+2+2*thischrom] = isamax(themap->mpc[ii], lin_reg[ii], 1);
        lin_reg[ii][ bp[2][marks+1+2*thischrom] ] = temp;
        thischrom = thischrom+1;
      }
  }
  else if ( model == 6 ) {
/* Use up to ns background markers   */
	 rpos = (FPN) MAXCHROMLEN;
	 lpos = (FPN) 0.0;
	 for ( gptr = agptr ; gptr != NULL ; gptr = gptr->next ) {
		if ( gptr->chrom == chrom && gptr->markr == lfm )
		  lpos = gptr->pos - wind;
		else if ( gptr->chrom == chrom && gptr->markr == lfm+1 )
		  rpos = gptr->pos + wind;
		gptr->mxo = (FPN) 0.0;
	}
	 ii = 1;
	 gptr = agptr;
	 do {  /* How many will we actually use?  */
	  if ( gptr->markr != 0 && gptr->whichqtl <= ns && gptr->whichqtl > 0 )
		if (gptr->chrom != chrom || (gptr->chrom == chrom && gptr->pos < lpos) || (gptr->chrom == chrom && gptr->pos > rpos) )
          ii = ii+1;

      gptr = gptr->next;
    } while ( gptr != NULL );
	 bp[1][1]   = ii;
	 bp[2][1]   = 6;
    ii = 1;
    gptr = agptr;
    do { /* Now put them in the bp matrix.  */
	  if ( gptr->markr != 0 && gptr->whichqtl <= ns && gptr->whichqtl > 0)
	    if (gptr->chrom != chrom || (gptr->chrom == chrom && gptr->pos < lpos) || (gptr->chrom == chrom && gptr->pos > rpos) ) {
          ii = ii+1;
          bp[1][ii]   = gptr->chrom;
          bp[2][ii]   = gptr->markr;
        }
      gptr = gptr->next;
    } while ( gptr != NULL );
    
  }
  else if ( model == 7 ) {
/* Use the results summarized in the Eqtl.out file.
   Use all that are for this trait. 
   1. get eqtl.out results and put into the standard ... done in Zmain.c
   2. use those results for the background...we assume that the
      data have been run through Zmapqtl for the trait we are now
      analyzing.     
      
      Use virtual markers for the background.
*/
     ns = 0;
     for ( ii = 1 ; ii <= themap->traits ; ii++ )
       ns = ns + themap->knum[ii];
     for ( ii = 1 ; ii <= ns ; ii++ )
       theqtls[ii].d = -(FPN) 1.0;

	 rpos = (FPN) MAXCHROMLEN;
	 lpos = (FPN) 0.0;
	 for ( gptr = agptr ; gptr != NULL ; gptr = gptr->next ) {
		if ( gptr->chrom == chrom && gptr->markr == lfm )
		  lpos = gptr->pos - wind;
		else if ( gptr->chrom == chrom && gptr->markr == lfm+1 )
		  rpos = gptr->pos + wind;
	 }
	 ii = 1;
	 gptr = agptr;
	 do {
	  if ( gptr->whichqtl <= ns && gptr->whichqtl > 0 ) 
		if (gptr->chrom != chrom || ( gptr->pxo < lpos  ||  gptr->pxo > rpos) ) {
          ii = ii+1;
          theqtls[gptr->whichqtl].d = (FPN) 1.0;
        }
      gptr = gptr->next;
    } while ( gptr != NULL );
    
    bp[1][1] = 0;
    for ( ii = 1 ; ii <= ns ; ii++ )
       if ( theqtls[ii].d > (FPN) 0.0) {
         bp[1][1] = bp[1][1] + 1;
         if ( theparams->verbosity == 2 )
           printf("\n Chrom %d, Marker %d, Trait %d",theqtls[ii].chrm,theqtls[ii].mrk,theqtls[ii].trait);
       }
     bp[2][1] = 7;
    
  }
  else {	/* default.  Let's make it Lander Botstein. */
	 bp[1][1] = 1;
	 bp[2][1] = 3;
  }
}

/* Fit the null hypothesis for  the background markers specified in bp. */
FPN   fit_back_params(params *theparams, individual *individs, markermap *themap, aqtl *theqtls, FPN **xw, FPN *yy, FPN *qraux, FPN *rsd, int chrom, int markr, int **bp, int nbp, genome *gptr)
{
  FPN rssp1, s2n, ts0,*stmp;
  int info,  ii, jj, kk,go_on,nn,cross,trait,roffset,row,ind,ns, *tmp;
  genome *tgptr,*lgptr,*rgptr,lg,rg;
  cleanse_genome(gptr);
  ii = chrom;
  jj = markr;
  nn = theparams->nn;
  cross = theparams->cross;
  if ( cross == 4 )
    cross = 3;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  trait = theparams->whichtrait;
  tmp = NULL; stmp=NULL;
  for (ii = 1; ii <= nn; ii++) /* mean in row 1  */
    xw[1][ii] = (FPN) 1.0;
  if (nbp > 1 && theparams->Model != 7 )
    for (ii = 2; ii <= nbp; ii++) {     /* background markers in rows 2 through nbp  */
      go_on = 1;
      tgptr = gptr;
      while ( go_on == 1 ) 
        if ( tgptr==NULL ) 
          go_on = 0;
        else if (tgptr->chrom== bp[1][ii] && tgptr->markr== bp[2][ii] ) 
           go_on = 0;
        else
          tgptr = tgptr->next;   
      jj = 1;
      for (kk = 1; kk <= nn; kk++)
	    if ( individs[kk].print_flag != 'n') {
	      switch ( individs[kk].markers[bp[1][ii]][bp[2][ii]] ) {
	        case 1:
  	        case 0:
	        case -1:
	          xw[ii][jj] = (FPN) individs[kk].markers[bp[1][ii]][bp[2][ii]];
	          break;
	        default: 
              xw[ii][jj] = expect_mark(theparams,individs, kk,tgptr);
	        break;
	      }
	      if (cross == 3)
	        xw[ii][jj] = xw[ii][jj] + (FPN) 1.0;
	      jj = jj + 1;
        }
  }
  if (nbp > 1 && theparams->Model == 7 ) {
    ii = 1;
    ns = 0;
    for ( ii = 1 ; ii <= themap->traits ; ii++ )
      ns = ns + themap->knum[ii];
    for ( tgptr = gptr ; tgptr != NULL ; tgptr = tgptr->next ) {    
     /* background markers in rows 2 through nbp taken from theqtls */
      if ( tgptr->whichqtl > 0 && tgptr->whichqtl <= ns && theqtls[tgptr->whichqtl].d > (FPN) 0.0 ) {
        /* Put a new node in the chain for the virtual marker... */
        lgptr = &lg;
        rgptr = &rg;
        lgptr->prev = tgptr->prev;
	    if ( tgptr->prev != NULL )
	      lgptr->prev->next = lgptr;
   	    rgptr->prev = lgptr;
	    lgptr->next = rgptr;

	    rgptr->next = tgptr->next;
	    if ( tgptr->next != NULL )
	      rgptr->next->prev = rgptr;

    	lgptr->dist = mapfunc((theqtls+tgptr->whichqtl)->c1, 1 );
	    rgptr->dist = mapfunc((theqtls+tgptr->whichqtl)->c2, 1 );
	    lgptr->chrom = tgptr->chrom;
	    lgptr->markr = tgptr->markr;
	    rgptr->chrom = tgptr->chrom;
	    rgptr->markr = -10;      
        ii = ii+1;
        jj = 1;
        for (kk = 1; kk <= nn; kk++)
	      if ( individs[kk].print_flag != 'n') {
            xw[ii][jj] = expect_mark(theparams,individs, kk,rgptr);
 	      
	        if (cross == 3)
	          xw[ii][jj]  = xw[ii][jj] + (FPN) 1.0;
	        jj = jj + 1;
          } /* Remake the chain */
        if ( tgptr->prev != NULL )
  	      tgptr->prev->next = tgptr;
	    if ( tgptr->next != NULL )
	      tgptr->next->prev = tgptr;
      }
    }
  }  
  
/* other traits in the rest of the rows nbp+1 to nbp+roffset  */
  roffset = 0;
  for ( ii = 1 ; ii <= themap->otraits ; ii++ )
    if ( themap->onames[ii][0] == '+' )
      roffset = roffset + themap->otypes[ii] - 1;
  if ( roffset > 0 ) {
    ind = 0;    
    for ( ii = 1 ; ii <= theparams->nn ; ii++ ) {
      row = nbp+1;
      if ( individs[ii].print_flag != 'n') {
        ind = ind+1; /* The rest of the otraits in rows nbp+1p onward Otraits   */
        for ( kk = 1 ; kk <= themap->otraits ; kk++ ) {
          if ( themap->onames != NULL && themap->onames[kk][0] == '+' ) {
	        for ( jj=row; jj < row + themap->otypes[kk]-1 ; jj++ ) 
	          if ( individs[ii].oyt[3] == jj-row+1 ) 
	            xw[jj][ind] = (FPN) 1.0 ;
	          else 
	            xw[jj][ind] =  (FPN) 0.0; 
	          row = jj;
	      }	      
        }
      }
    }
  }  
  
  jj = 0;
  for (kk = 1; kk <= nn; kk++)
    if ( individs[kk].print_flag != 'n') 
      jj = jj + 1;
  for (ii = 1; ii <= nn; ii++)
    rsd[ii] = (FPN) 0.0;
  for (ii = 1; ii <= nbp; ii++)
    qraux[ii] = (FPN) 0.0;

  info = sqrdc(xw, jj, jj, nbp+roffset, qraux, tmp, 0);
  info = sqrsl(xw, jj, jj, nbp+roffset, nbp+roffset, qraux, yy, stmp, rsd, stmp, rsd, stmp, 10L);
  rssp1 = sdot(jj, rsd, 1, rsd, 1);
  s2n = rssp1 / (FPN) (jj);  /*  Should we use the sampling variance here? No...seach for SAMPLINGVARIANCE in this file  */
  theparams->null_sse = s2n;
  if (s2n > (FPN) 0.0) 
    ts0 = -(FPN) jj / (FPN) 2.0 * ((FPN) log(s2n) + (FPN) 1.0);
  else
    ts0 = -(FPN) jj / (FPN) 2.0 * ((FPN) log(-s2n) + (FPN) 1.0);
  return (ts0);
}



/*
  Equation numbers are those in Z.-B. Zeng (1994) Genetics 136:1457-1468.

lnpk->estimates[1]  ests is a pointer that will hold the estimates.
lnpk->xw    xw is the output of sqrdc
            y  is the vector of phenotypes
            ldx = nn is the sample size
            pp = kk is the number of rows in xw
lnpk->rsd   rsd is the residual vector from the null hypothesis, the result of sqrdc and sqrsl
lnpk->qraux qraux is aux. info from sqrdc
            gval is the conditional probabilities of the QTL given flanking markers
            marks is the matrix of flanking markers
            table is a way to determine which value of gval should be used, based on marks
            hypo is the hypothesis to be tested
            cross is the type of cross

lnpk is workspace.  

*/
int     em_solve(linpakws *lnpk, FPN *y, int nn, int pp, int ihypo, params *theparams)
{
   FPN *tmp, cp,cq,ctemp1,ctemp2,ctemp3;
  FPN test_stat,test_stat0,const2, s2, rssp,be,de;
  FPN tc1,tc2,tc3,temp, *estimates;
  int ktime, ii, go_on, info,cross,bailout;
  estimates = lnpk->estimates[1];
  bailout = 0;  /* Just an indicator if one should bail out...0 means no, >0 means yes*/
/*
  With all the new experimental designs, we have to modify this section
  to estimate the effects correctly.  Most notably, the Design III is probably
  not correct.     Otherwise, T(XX)SFx seems to be ok (it's like a backcross).
*/  
  cross = theparams->cross;
  if ( theparams->tcross ==1 || theparams->tcross == 2 ) 
    cross = theparams->tcross;  /* If this is a test cross to P1 or P2, treat it like a backcross.*/
  if ( theparams->tcross == 12 ) /* otherwise, treat it like an SF2 */
    cross = 3;
  if ( cross == 4 ) /* RFx is like SFx as far as the estimates go. */
    cross = 3;
  tmp = NULL;
  de = be = test_stat = (FPN) 0.0;
  scopy(nn, lnpk->pp1, 1, lnpk->pv, 1);  /* copy a priori probabilities into */
  scopy(nn, lnpk->pp2, 1, lnpk->qv, 1);  /* a posteriori probs. initially */
  scopy(nn, lnpk->rsd[1], 1, lnpk->wrsd[1], 1); /* copy residuals into a residual workspace */
  ktime = 0;
  go_on = 1;
  de = be = (FPN) 0.0;
  while (go_on == 1) {
/* This allows you to iconify the program when it runs in MS Windows. 
   I need to have something like this for the Mac versions...*/
#if defined(BORLAND)
    WYield();
#endif
/*  Do the M step first...*/
    cp = ssum(nn, lnpk->pv, 1);	            /* cp is the Sum of the elements of pv */
    ctemp1 = sdot(nn, lnpk->wrsd[1], 1, lnpk->pv, 1);	/* Numerator of Eqn 7 */
    if (cross == 3) {
      cq = ssum(nn, lnpk->qv, 1);          	/* cq is the Sum of the elements of qv */
      ctemp2 = sdot(nn, lnpk->wrsd[1], 1, lnpk->qv, 1);	/* ctemp2 is rsd.qv */
      if (ihypo == 2)
	    be = (FPN) 0.0;
      else if ( cp != (FPN) 0.0 )
	    be = ctemp1 / ((FPN) 2.0 * cp);
	  else /* if the denominator is 0, something's awry. */
	    bailout = bailout+1; 
      if (ihypo == 1)
	    de = (FPN) 0.0;
      else if ( cq != (FPN) 0.0 )
	    de = ctemp2 / cq - be;
	  else /* if the denominator is 0, something's awry. */
	    bailout = bailout+10; 
    }
    else if ( cp != (FPN) 0.0 )
      be = ctemp1 / cp;	/* Eqn 7 */
	else /* if the denominator is 0, something's awry. */
	  bailout = bailout+1;

/* Here is the trick:  Use (Y - XB - Pbe), then add back Pbe. */
    scopy(nn, y, 1, lnpk->wy, 1);
    if ( cross == 1 || cross == 2 || cross == 5)
      saxpy(nn, -be, lnpk->pv, 1, lnpk->wy, 1);
    else if ( cross == 3 ) {
      saxpy(nn, -2*be, lnpk->pv, 1, lnpk->wy, 1);
      saxpy(nn, -(be+de), lnpk->qv, 1, lnpk->wy, 1);
    }


    info = sqrsl(lnpk->xx,nn,nn,pp,pp,lnpk->qraux,lnpk->wy,tmp,lnpk->wrsd[1],tmp,lnpk->wrsd[1],tmp,10L);
    
    scopy(nn, lnpk->wrsd[1], 1, lnpk->wy, 1); /* Save residuals */    
    if ( cross == 1 || cross == 2 || cross == 5)
      saxpy(nn, be, lnpk->pv, 1, lnpk->wrsd[1], 1);
    else if ( cross == 3 ) {
      saxpy(nn, 2*be, lnpk->pv, 1, lnpk->wrsd[1], 1);
      saxpy(nn, (be+de), lnpk->qv, 1, lnpk->wrsd[1], 1);
    }

    rssp = sdot(nn, lnpk->wrsd[1], 1, lnpk->wrsd[1], 1);
    if (cross == 3)
      s2 = (rssp - (FPN)4.0 * cp * be * be - cq * (be + de) * (be + de)) / (FPN) (nn);
    else
      s2 = (rssp - cp * be * be) / (FPN) (nn);	/* This is equation 9 */
    test_stat0 = test_stat;
    if (s2 > (FPN) 0.0)	 
      test_stat = -((FPN) nn * (FPN) log(s2)  + rssp / s2) / (FPN) 2.0;
    else if ( s2 < (FPN) 0.0 ) 
      test_stat =  -((FPN) nn * (FPN) log(-s2)  + rssp / s2) / (FPN) 2.0;
    else /* if s2 is 0, something's awry. */
      bailout = bailout+100;
    if ( bailout == 0 )
      for (ii = 1; ii <= nn; ii++) {
        if (cross == 3) {
	      temp = lnpk->pp1[ii] * (FPN)exp(2.0 * be * (lnpk->wrsd[1][ii] - be) / s2);
	      temp = temp + lnpk->pp2[ii] * (FPN)exp((be + de) * (lnpk->wrsd[1][ii] - (be + de) / 2.0) / s2);
	      temp = temp + (FPN) 1.0 - lnpk->pp1[ii] - lnpk->pp2[ii];
        }
        else if (cross == 5)
	      temp = lnpk->pp1[ii] * (FPN)exp(be * (2.0 * lnpk->wrsd[1][ii] - be) / (2.0 * s2)) + (FPN) 1.0 - lnpk->pp1[ii];
        else
	      temp = lnpk->pp1[ii] * (FPN)exp(be * (2.0 * lnpk->wrsd[1][ii] - be) / (2.0 * s2)) + lnpk->pp2[ii];

        if (temp > (FPN) 0.0)	/*********************Something to check on************************************/
	      test_stat = test_stat + (FPN) log(temp);
      }
    else if ( bailout > 0 ) {
      if ( theparams->verbosity == 1 )
        printf("\nBailing out (%d) of ECM algorithm...",bailout);
      go_on = 0;
    }
    if ( ktime >= (int) M_TIME )
      go_on = 0;
    else if ( ktime >= 2 && bailout == 0 ) {
      if ( test_stat0 > (FPN) 1.0 &&  fabs(test_stat-test_stat0) / test_stat0 < (FPN) STOP_EM  ) 
	    go_on = 0;
      else if ( fabs(test_stat - test_stat0) < (FPN) STOP_EM )  
	    go_on = 0;
    }

    if (go_on == 1) {	/* Do the E step */
      const2 = (FPN) 1.0 / ((FPN) 2.0 * s2);
      ktime = ktime + 1;
      for (ii = 1; ii <= nn; ii++) {
	    tc2 = lnpk->wrsd[1][ii];
	    tc1 = lnpk->wrsd[1][ii] - be;
	    if (cross == 3) {
	      tc3 = lnpk->wrsd[1][ii];
	      tc1 = lnpk->wrsd[1][ii] - (FPN) 2.0 * be;
	      tc2 = lnpk->wrsd[1][ii] - be - de;
	    }
	    ctemp1 = lnpk->pp1[ii] * (FPN)exp(-const2 * tc1 * tc1);
	    if ( cross == 5 )
	      ctemp2 = ((FPN) 1.0- lnpk->pp1[ii]) * (FPN)exp(-const2 * tc2 * tc2);
	    else
	      ctemp2 = lnpk->pp2[ii] * (FPN)exp(-const2 * tc2 * tc2);
	    if (cross == 3)
	      ctemp3 = ((FPN) 1.0 - lnpk->pp1[ii] - lnpk->pp2[ii]) * (FPN)exp(-const2 * tc3 * tc3);
	    else
          ctemp3 = (FPN) 0.0;
        ctemp3 = ctemp1 + ctemp2 + ctemp3;
        if ( ctemp3 == (FPN) 0.0 )
          bailout = bailout+1000;
        else {
	      lnpk->pv[ii] = ctemp1 /  ctemp3;
	      if (cross == 3)
	        lnpk->qv[ii] = ctemp2 /  ctemp3;
	    }
      }
    }
  }
  if ( bailout > 0 ) 
    test_stat = be = de = (FPN) SIGNOFDEVIL;
  calc_qstats(lnpk->wy,nn,lnpk->thestats,15);
    

  if (ihypo == 1) {
    estimates[1] = test_stat;
    estimates[4] = be;
    if ( cross == 5)  /*For Ri lines, we've estimated 2 times the additive effect...*/
      estimates[4] = (FPN)0.5*be;
    else if ( cross == 2 )
      estimates[4] = -be ;   /* For BC2, we have estimated  -1 times the additive effect. */
    estimates[13] = estimates[10] = s2 ;
  }
  else if (ihypo == 2) {
    estimates[2] = test_stat;
    estimates[6] = de;
    estimates[14] = estimates[11] = s2;
  }
  else if (ihypo == 3) {
    estimates[3] = test_stat;
    estimates[5] = be;
    estimates[7] = de;
    estimates[15] = estimates[12] = s2 ;
  }
  return(bailout);
}


/*
Strictly for debugging purposes.   Print the prior and posterior genotypes for
the QQ, Qq and qq genotypes.

How many?  This does them all.  

*/

void showprgenotypes(params *theparams, linpakws *lnpk) {
  FILE *errorf;
  int ii;
  
  	errorf = fileopen(theparams->error, "a" );
  	fprintf(errorf, "\n\n   Prior [p()] and Posterior [P()] QTL genotype probabilities...\n\n");
    putline(errorf,'-',50);
    putline(errorf,'-',50);
	fprintf(errorf, "\n  i  p(QQ)   p(Qq)   p(qq)   P(QQ)   P(Qq)   P(qq)");  	
    putline(errorf,'-',50);
  	for ( ii=1; ii<=lnpk->n; ii++ )
  	  fprintf(errorf, "\n%4d  %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f",ii, lnpk->pp1[ii], lnpk->pp2[ii], (FPN) 1.0-(lnpk->pp1[ii]+lnpk->pp2[ii]), lnpk->pv[ii], lnpk->qv[ii], (FPN) 1.0-(lnpk->pv[ii]+lnpk->qv[ii])); 
    putline(errorf,'-',50);
    putline(errorf,'-',50);
  	fprintf(errorf, "\n\n");
  	fileclose(theparams->error, errorf);
  

}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Zfunc.c
------------------------------------------------------------------ */


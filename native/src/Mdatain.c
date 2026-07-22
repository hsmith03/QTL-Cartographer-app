/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Mdatain.c) is part of QTL Cartographer
         
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



/* Mdatain.c, Subroutines for reading Genetic Linkage Maps*/

#include "Main.h"


/*

Routine to get a map
*/
void GetTheMap(params *theparams, char *infile) {
	int markers,error;
	FILE *outf;
	
    markers = get_file_type(infile);

    if ( markers == 0 )
      markers = det_if_mm(infile);
    
    if ( markers == 10 )
      theparams->themap = get_markermap(theparams->qnum, infile);
    else if ( markers == 11 || markers == 4 )  /*  It is either a map.inp file or a qtlcart.mcd file.   Same routine for the map. */
      theparams->themap = trans_stdm(infile);
    else if ( markers == 12 )
      theparams->themap = trans_maps(infile);
    else {
      LogTheError(theparams->error,"\n  !!! Could not read the map.  Exiting the program.  \n"); 
      exit(1);
      
      
    }  
    error = check_params(theparams,theparams->themap,whichprogram);
    if ( markers == 10 &&  theparams->themap->names == NULL ) { /* give names to standard map if they don't exist  */
      theparams->themap->traits = 0;
      create_bogus_markers(theparams, theparams->themap );
      theparams->themap->traits = 1;
      outf = fileopen(infile, "a" );
      print_marker_names(theparams->themap, outf);
      fileclose(infile,outf);
    }

}

/*
  This takes a genome structure and prints a map in map.inp format.
*/
void print_genome_map(genome *thegenome,params *theparams, char *outfile,char *progname,char *chptr) {
  int chrom, maxmark;
  FPN thedist;
  genome *gptr;
  FILE *fptr;
  
  if ( outfile[0] == '-')
    fptr = stdout;
  else 
    fptr = fileopen(outfile, "w");
  
  chrom = maxmark = 0;
  for ( gptr=thegenome; gptr!=NULL; gptr=gptr->next ) {  
    if ( gptr->chrom > chrom )
      chrom = gptr->chrom;
    if ( gptr->markr > maxmark)
      maxmark = gptr->markr;
   }
  fprintf(fptr,"# %ld    bychromosome  -filetype map.inp",theparams->seed);   
  fprintf(fptr,"\n# This file created by %s ",progname );   
  fprintf(fptr,"\n# on %s#",chptr ); 
  if (theparams->emethod < 20 ) 
     fprintf(fptr,"\n# This map was genereated using the Rapid Chain Delineation Method with " ); 
  else if (theparams->emethod < 30 )  
     fprintf(fptr,"\n# This map was genereated using the Seriation Method with " ); 
  else if (theparams->emethod < 40 )  
     fprintf(fptr,"\n# This map was genereated using the Simulated Annealing Method with " ); 
  else if (theparams->emethod < 50 )  
     fprintf(fptr,"\n# This map was genereated using the Branch and Bound Method with " ); 
  fprintf(fptr,"\n# Segregation test size %f and linkage threshold %f.",theparams->segsize,theparams->linksize); 
  fprintf(fptr,"\n-type intervals ");
  fprintf(fptr,"\n-function     %d  ",theparams->mapfunc); 
  fprintf(fptr,"\n-parameter    %f  ",theparams->mapparam);
  fprintf(fptr,"\n-Units        cM  ");
  fprintf(fptr,"\n-chromosomes  %4d  ",chrom);
  fprintf(fptr,"\n-maximum      %4d ",maxmark);
  fprintf(fptr,"\n-named        yes  ");
  fprintf(fptr,"\n-start ");
  chrom = 0;
  for ( gptr=thegenome; gptr!=NULL; gptr=gptr->next ) 
    if ( gptr->chrom > 0 ) {
      if (  gptr->chrom != chrom ) {
        fprintf(fptr,"\n-Chromosome %d", gptr->chrom);
        chrom = gptr->chrom;
      }
      thedist = mapfunc(gptr->dist,2);
      fprintf(fptr,"\n   %15s  %f",  gptr->markername, thedist);
    }
  fprintf(fptr,"\n-stop\n-end\n");   
  if (outfile[0] != '-')
      fileclose(outfile, fptr);
}


/*
 * This prints the map of markers, either to a file or stdout.  This is the standard
 format for maps in the QTL Cartographer system.
 */
void print_map(markermap *themap,char *outfile)
{
  int ii, jj, maxl, mtot, printspace;
  FPN dist;
  FILE *outf;
  if (*outfile == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  if (outf != NULL) {
      mtot = maxl = 0;
      for (ii = 1; ii <= themap->m; ii++) {
	  if (maxl < themap->mpc[ii])
	    maxl = themap->mpc[ii];
	    mtot = mtot + themap->mpc[ii];
      }

    fprintf(outf, "\n#\n#    Here is the Map of Markers...\n-s");
    fprintf(outf, "\n-f  %2d  Map function [1,8]:  Haldane is 1", whosemf);
    fprintf(outf, "\n-p  %8.4f  Extra parameter for map functions [4-8]", mapparam);
    fprintf(outf, "\n-u  c    The units of measurement is centiMorgans.\n#\n# Markermap parameters\n#");
    fprintf(outf, "\n-c  %8d  Number of chromosomes.", themap->m);
    fprintf(outf, "\n-i  %8d  Total number of markers.", mtot);
    fprintf(outf, "\n-m  %8d  Mean, and...", themap->l);
    fprintf(outf, "\n-vm %8.4f  ...standard deviation for Markers/Chromosome.", themap->sigl);
    fprintf(outf, "\n-d  %8.4f  Mean, and...", themap->s);
    fprintf(outf, "\n-vd %8.4f  ...standard deviation for the intermarker distance.", themap->sigs);
    fprintf(outf, "\n-t  %8.4f  Mean length of material outside the flanking markers.", themap->brdrs);
    fprintf(outf, "\n#\n          |   Chromosome---->\n");
    for (ii = 1; ii <= themap->m * 8 + 12; ii++)
      fprintf(outf, "-");
    fprintf(outf, "\nMarker    |");
    for (ii = 1; ii <= themap->m; ii++)
      fprintf(outf, " %7d", ii);
    fprintf(outf, "\n");
    for (ii = 1; ii <= themap->m * 8 + 12; ii++)
      fprintf(outf, "-");
    for (ii = 0; ii <= maxl; ii++) {
      fprintf(outf, "\n-l %6d | ", ii);
      for (jj = 1; jj <= themap->m; jj++) {
	    printspace = 0;
	  if (ii > themap->mpc[jj])
	    printspace = 1;
	if (printspace == 1)
	  fprintf(outf, "        ");
	else {
	  dist = mapfunc( themap->mrf[jj][ii], 2);
	  fprintf(outf, " %7.4f", dist);
	}
      }
    }
    fprintf(outf, "\n");
    for (ii = 1; ii <= themap->m * 8 + 12; ii++)
      fprintf(outf, "-");
    fprintf(outf, "\n-Number   |");
      for (ii = 1; ii <= themap->m; ii++)
	fprintf(outf, " %7d", themap->mpc[ii]);

    fprintf(outf, "\n\n");

    if (themap->names != NULL)
      print_marker_names(themap, outf);

    if (*outfile != '-')
      fileclose(outfile, outf);
  }
}


/*
 * This prints the map of markers, either to a file or stdout.
      The output is in map.inp format.   It assumes that the first line 
      has already been printed, so that it can be used for map.inp or
      qtlcart.mcd files.   
 */
void print_map_inp(markermap *themap,char *outfile)
{
  int ii, jj, maxl;
  FPN dist;
  FILE *outf;
  if (*outfile == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  if (outf != NULL) {
     maxl = 0;
     for (ii = 1; ii <= themap->m; ii++) 
	   if (maxl < themap->mpc[ii])
	     maxl = themap->mpc[ii];
 
    fprintf(outf, "\n#\n#    Here is the Map of Markers...\n#");
    fprintf(outf, "\n-type       positions");
    fprintf(outf, "\n-function          %2d   Map function [1,8]:  Haldane is 1", whosemf);
    fprintf(outf, "\n-parameter   %8.4f   Extra parameter for map functions [4-8]", mapparam); 
    fprintf(outf, "\n-Units             cM   The units of measurement is centiMorgans. ");
    fprintf(outf, "\n-chromosomes %8d   Number of chromosomes", themap->m);
    fprintf(outf, "\n-maximum     %8d   Max. number of markers per chromosome",  maxl );
    fprintf(outf, "\n-named            yes   Markers have names"   );
    fprintf(outf, "\n-start");
    for (ii = 1; ii <= themap->m  ; ii++) {
      if ( themap->cnames == NULL )
        fprintf(outf, "\n-Chromosome c%d",ii);
      else
        fprintf(outf, "\n-Chromosome %s",themap->cnames[ii]);
      dist = mapfunc( themap->mrf[ii][0], 2);
      for ( jj= 1 ; jj<=themap->mpc[ii]; jj++ ) {
	    fprintf(outf, "\n%s %7.4f", themap->names[ themap->ttable[ii][jj] ], dist);
	    dist = dist + mapfunc( themap->mrf[ii][jj], 2);
      }
      if ( themap->mrf[ii][themap->mpc[ii]] > (FPN) 0.0 )
	    fprintf(outf, "\nTelomere %7.4f", dist);
    }
    fprintf(outf, "\n-stop\n\n---------------------------------------------------");

    if (*outfile != '-')
      fileclose(outfile, outf);
  }
}

/*
 * This prints the map of markers, either to a file or stdout.  This is the standard
 format for maps in the QTL Cartographer system.
 */
void print_plabqtl_map(markermap *themap,FILE *outf)
{
  int ii, jj;
  FPN dist;

    for (ii = 1; ii <= themap->m; ii++) {
      if ( themap->cnames != NULL ) 
        fprintf(outf, "\n*%s",  themap->cnames[ii]);
      else
        fprintf(outf,"\n*chrom%d",ii);


        fprintf(outf," %d",themap->mpc[ii]);
        for (jj = 1; jj < themap->mpc[ii]; jj++) {
          dist = mapfunc( themap->mrf[ii][jj], 2);
          fprintf(outf,"\n%s %f",themap->names[themap->ttable[ii][jj]],dist);
        }
        fprintf(outf,"\n%s ",themap->names[themap->ttable[ii][themap->mpc[ii]]] );



   }
   fprintf(outf,"\n");
}
/*
  These subroutines allow one to read in a map of markers in the format
  that is specified in the output subroutines of Rmap.
*/

/*
  * read_marker_params merely reads in the parameters from the
  * data file. It returns -1 for an error, either because
  * there is no data file or because all the parameters were
  * not found.
  */
int read_marker_params(int *mm, int *ll, FPN *ssigl, FPN *ss, FPN *ssigs, FPN *brdrs,char *inputf)
{
  FILE *infile;
  int  numparams, ch;
  whosemf = 1;
  numparams = 0;
  infile = fileopen(inputf, "r");
  if (infile == NULL) 
    return (-1);
  do {
    ch = get_next_line(gbuffer, MAXLINE, infile);
    if (gbuffer[0] == '-')
      switch (gbuffer[1]) {
       case 'c':
		get_field(2, gname, gbuffer);
		*mm = atoi(gname);
		numparams = numparams + 1;
		break;
       case 'f':
		get_field(2, gname, gbuffer);
		whosemf = atoi(gname);
		break;
       case 'p':
		get_field(2, gname, gbuffer);
		mapparam = (FPN) atof(gname);
		break;
       case 'm':
		get_field(2, gname, gbuffer);
		*ll = atoi(gname);
		numparams = numparams + 1;
		break;
       case 'd':
		get_field(2, gname, gbuffer);
		*ss = (FPN) atof(gname);
		numparams = numparams + 1;
		break;
       case 't':
		get_field(2, gname, gbuffer);
		*brdrs = (FPN) atof(gname);
		numparams = numparams + 1;
		break;
       case 'v':
		get_field(2, gname, gbuffer);
		if (gbuffer[2] == 'm')
		  *ssigl = (FPN) atof(gname);
		else if ( gbuffer[2] == 'd')
		  *ssigs = (FPN) atof(gname);
		numparams = numparams + 1;
		break;
       case 'l':
		fileclose(inputf, infile);
		if (numparams != 6) {
		  printf("\nYou did not specify enough parameters...\n");
		  return (-1);
		}
		else 
		  return (1);
		break;
       default:
		break;
      }
  } while (ch != EOF);
  if (ch == EOF)
    printf("\nTime to close the input file after an EOF...\n");
  fileclose(inputf, infile);
  return (-1);
}

/*
  This is almost a small driver for other routines that read
  the parameters and make the map of markers.
*/
markermap *get_markermap(int knum,char *infile)
{
  markermap *themap;
  int mm, ii, ll, error;
  FPN ssigl, ss, ssigs, bbrdrs;
  error = read_marker_params(&mm, &ll, &ssigl, &ss, &ssigs, &bbrdrs, infile);
  if (error == -1)
    exit(1);
  themap = initialize_markermap(knum, mm, ll, ssigl, ss, ssigs, bbrdrs, infile);
  themap->ml = 0;
    for (ii = 1; ii <= themap->m; ii++)
      themap->ml = themap->ml + themap->mpc[ii];
  do_map_stats(themap);
  get_markernames(themap, infile);
  return (themap);
}

/*
  see if there are any marker names and read them in if so.
*/

void get_markernames(markermap *themap,char *inputf)
{
  FILE *infile;
   int  jj, cc, mm, ll, ch;
  char   *chptr;
  infile = fileopen(inputf, "r");
  ll = 0;
  do {
    ch = get_next_line(gbuffer, MAXLINE, infile);
    if ( gbuffer[0] == '-')
      switch ( gbuffer[1]) {
       case 'b':
        get_field(2, gname, gbuffer);
        chptr = strlwr(gname);
        if ( !strcmp(gname,"markernames" )   ) {
	      themap->names = cmatrix(1, themap->ml, 0, MAXNAME);
	      themap->ttable = imatrix(1, themap->m, 1, themap->maxl);
	      for (jj = 1; jj <= themap->ml; jj++) {
	        ch = get_next_token(gname, MAXNAME, infile);
	        cc = atoi(gname);
	        ch = get_next_token(gname, MAXNAME, infile);
	        mm = atoi(gname);
	        ch = get_next_token(gname, MAXNAME, infile);
	        ll = ll + 1;
	        strcpy(themap->names[ll], gname);
	        themap->ttable[cc][mm] = ll;
          }
        }
        else if ( !strcmp(gname,"chromosomenames" ) ) {
	      themap->cnames = cmatrix(1, themap->m, 0, MAXNAME);
	      for (jj = 1; jj <= themap->m; jj++) {
	        ch = get_next_token(gname, MAXNAME, infile);
	        cc = atoi(gname);
	        ch = get_next_token(gname, MAXNAME, infile);
	        strcpy(themap->cnames[cc], gname);
          }
        }
        else if ( !strcmp(gname,"traitnames" ) ) {
	      themap->tnames = cmatrix(1, themap->traits, 0, MAXNAME);
	      for (jj = 1; jj <= themap->traits; jj++) {
	        ch = get_next_token(gname, MAXNAME, infile);
	        cc = atoi(gname);
	        ch = get_next_token(gname, MAXNAME, infile);
	        strcpy(themap->tnames[cc], gname);
          }
        }
        else if ( !strcmp(gname,"othertraitnames" ) ) {
	      themap->onames = cmatrix(1, themap->otraits, 0, MAXNAME);
	      for (jj = 1; jj <= themap->otraits; jj++) {
	        ch = get_next_token(gname, MAXNAME, infile);
	        cc = atoi(gname);
	        ch = get_next_token(gname, MAXNAME, infile);
	        strcpy(themap->onames[cc], gname);
          }
        }
	    break;
      default:
	    break;
    }
  } while (ch != EOF);  
  fileclose(inputf, infile);
}


/*  Calculate means and variances of intermarker distances, 
    and number of markers per chromosome.  */
void do_map_stats(markermap *themap)
{
  int ii, jj;
  FPN t1, t2, ss, dist;
  t1 = t2 = 0;
  if (  themap->m > 1 ) {
    for (ii = 1; ii <= themap->m; ii++) {
      t1 = t1 + (FPN) themap->mpc[ii];
      t2 = t2 + ((FPN)  themap->mpc[ii]) * ((FPN)  themap->mpc[ii]);
    }
    ss = (FPN) themap->m;
    themap->l = (int) (t1 / ss);
    themap->sigl = t2 / (ss - (FPN) 1.0) - t1 * t1 / (ss * (ss - (FPN) 1.0));
    if ( themap->sigl > (FPN) 0.0 )
      themap->sigl = (FPN) sqrt(themap->sigl);
  }
  else if ( themap->m == 1 ) {    
    themap->l = themap->mpc[1];
    themap->sigl = (FPN) 0.0;
  }  
  t1 = t2 = (FPN) 0.0;
  ss = (FPN) 0.0;
  for (ii = 1; ii <= themap->m; ii++)
    for (jj = 1; jj < themap->mpc[ii]; jj++)
      if ( themap->mrf[ii][jj] > (FPN) 0.0 ) {
        ss +=1;
	    dist = mapfunc( themap->mrf[ii][jj], 2);  /*This is 2.  r -> cM */
	    t1 = t1 + dist;
	    t2 = t2 + dist * dist;
      }
  themap->s = t1 / ss;
  themap->sigs = t2 / (ss - (FPN) 1.0) - t1 * t1 / (ss * (ss - (FPN) 1.0));
  if ( themap->sigs > (FPN) 0.0)
    themap->sigs = (FPN) sqrt(themap->sigs);
    
/*  Do the map stats for the tails (telomeric DNA)*/    
  t1 = t2 = (FPN) 0.0;
  ss = (FPN) 0.0;
  for (ii = 1; ii <= themap->m; ii++) {
      if ( themap->mrf[ii][0] > (FPN) 0.0 ) {
        ss +=1.0;
	    dist = mapfunc( themap->mrf[ii][0], 2);  /*This is 2.  r -> cM */
	    t1 = t1 + dist;
	    t2 = t2 + dist * dist;
      }
      if ( themap->mrf[ii][themap->mpc[ii]] > (FPN) 0.0 ) {
        ss +=1;
	    dist = mapfunc( themap->mrf[ii][themap->mpc[ii]], 2);  /*This is 2.  r -> cM */
	    t1 = t1 + dist;
	    t2 = t2 + dist * dist;
      }
  } 
  if ( ss>0.0 )  
    themap->brdrs = t1 / ss;
  else
    themap->brdrs = 0.0;
    
    
    
    
}



void read_markers(markermap *themap,char *inputf)
{
  FILE *infile;
  int  jj, chromosome, row, number, total, units, ch, linelength;
  char *linebuffer;
  units = 0;	/* 0 => recombination frequencies, 1 => Morgans and 2 => centiMorgans */
  linelength = 10 * (themap->m+3); 
  linebuffer = cvector(0,linelength);
  infile = fileopen(inputf, "r");            /* 0 default is changed in reading file. */
  if (infile == NULL) {
    return;
  }
  ch = 'c';
  do {
    get_next_line(linebuffer, linelength, infile);
    if (linebuffer[0] == '-')
      switch ( linebuffer[1]) {
       case 'u':
	     get_field(2, gname, linebuffer);
	    if ( gname[0] == 'M' || gname[0] == 'm')
	      units = 1;
	    else if (gname[0] == 'c' || gname[0] == 'C')
	      units = 2;
	   break;
       case 'l':
	     get_field(2, gname, linebuffer);
	     row = atoi(gname);
	     total = 0;
	     for (jj = 1; jj <= themap->m; jj++) {
	         number = themap->mpc[jj];
	       if (row <= number)
	         total = total + 1;
	     }
	     if (total > 0) {
	       chromosome = 1;
	       for (jj = 1; jj <= total; jj++) {
	         do {
		         number = themap->mpc[chromosome];
	           if (row > number)
		         chromosome = chromosome + 1;
	         } while (row > number);
	         get_field(jj + 3, gname, linebuffer);
	         themap->mrf[chromosome][row] = (FPN) atof(gname);
	         /* The internal default is recombination frequencies  (r)*/
	         if (units > 1)	/* convert from cM -> M */
	           themap->mrf[chromosome][row] = themap->mrf[chromosome][row] / (FPN) 100.0;
	         if (units > 0)	/* convert from M -> r */
	           themap->mrf[chromosome][row] = mapfunc(themap->mrf[chromosome][row], -1);
	         chromosome = chromosome + 1;
	       }
	    }
	    break;
       case 'N':
	     ch = EOF;
	     break;
       default:
	     break;
      }
  } while (ch != EOF);
  free_cvector(linebuffer,0,linelength);
  fileclose(inputf, infile);
}

/*
  This opens the input file, seeks a line that begins '-N', and then reads the integers on
  that line into the themap->mpc array.


  do {
    get_next_line(gbuffer, MAXLINE, infile);
    if ( gbuffer[0] == '-' )
      switch ( gbuffer[1] ) {
       case 'N':
	     maxl = 0;
	     for (jj = 1; jj <= themap->m; jj++) {
	       get_field(jj + 2, gname, gbuffer);
	       themap->mpc[jj] = atoi(gname);
	       if (themap->mpc[jj] > maxl)
	         maxl = themap->mpc[jj];
	     }
	     ch = EOF;
	     break;
       default:
	     break;
      }
  } while (ch != EOF);


*/
int read_marker_numbers(markermap *themap,char *inputf)
{
  FILE *infile;
  int maxl, ch;
  int  jj;
  infile = fileopen(inputf, "r");
  if (infile == NULL) {
    return (-1);
  }
  maxl = 0;
  do {
    ch = get_next_token(gbuffer, MAXLINE, infile);
    if ( !strcmp( gbuffer, "-Number" ) ) {
      ch = get_next_token(gbuffer, MAXLINE, infile);
      ch = EOF;
    }
  } while (ch != EOF);
  for (jj = 1; jj <= themap->m; jj++) {
    ch = get_next_token(gbuffer, MAXLINE, infile);
	themap->mpc[jj] = atoi(gbuffer);
	if (themap->mpc[jj] > maxl)
	  maxl = themap->mpc[jj];
  }
  
  fileclose(inputf, infile);
  return (maxl);
}


/*
 * Create the space for a map of markers, and place the markers on it.
 */
markermap *initialize_markermap(int kk, int mm,int ll,FPN ssigl,FPN ss,FPN ssigs,FPN bbrdrs,char *infile)
{
  int ii, jj;
  markermap *themap;
  themap = allocate_markermap();
  themap->traits = 1;
  themap->knum = ivector(1, themap->traits);
  themap->knum[1] = kk;
  themap->m = mm;
  themap->l = ll;
  themap->sigl = ssigl;
  themap->s = ss;
  themap->sigs = ssigs;
  themap->brdrs = bbrdrs;
  themap->names = NULL;
  themap->ttable = NULL;
  themap->mpc = ivector(1, mm);

  if (ssigl < (FPN) 0.0) {
    themap->ml = mm * ll;
    themap->maxl = ll;
    for (ii=1; ii<=mm; ii++ )
      themap->mpc[ii] = ll;
  }
  else {
    themap->ml = 0;
    themap->maxl = read_marker_numbers(themap, infile);
  }
  themap->mrf = dmatrix(1, mm, 0, themap->maxl);
  for (ii = 1; ii <= mm; ii++)
	 for (jj = 0; jj <= themap->maxl; jj++)
	   themap->mrf[ii][jj] = (FPN) 0.0;
  themap->types = imatrix(1, mm, 1, themap->maxl);
  for (ii = 1; ii <= mm; ii++)
    for (jj = 1; jj <= themap->maxl; jj++)
      themap->types[ii][jj] = 0;
  read_markers(themap, infile);
  return themap;
}

/* allocate some space for the markermap and initialize the values */
markermap *allocate_markermap(void )
{
   markermap *themap;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  themap = (markermap *) malloc( (size_t) sizeof(markermap));
#else
  themap = (markermap *) malloc( (unsigned) sizeof(markermap));
#endif
 if ( debugging > 2 ) {
        sprintf(gwarn,"In allocate_markermap(), allocated 1 markermap  at %p\n",themap);
        MemoryAccount(gwarn);
 }
/* struct MapofMarkers { */

  themap->m = 0;         /* number of chromosomes */
  themap->l = 0;         /* average number of markers per chromosome */
  themap->maxl = 0;      /* maximum number of markers on any chromosome */
  themap->sigl = (FPN) 0.0;   /* variance in number of markers per chromosome.  <= 0 => fixed */
  themap->s = (FPN) 0.0;      /* average distance between consecutive markers in cM */
  themap->sigs = (FPN) 0.0;   /* variance of s */
  themap->mpc = NULL;   /* = ivector(1,m); number of markers for each chromosome.  = l forall iff sigl <= 0 */
  themap->mrf = NULL;   /* = dvector(1,m,0,max(mpc)+1);  
                           pointer to a matrix of recombination frequencies between markers i and i+1 */ 
  themap->ml = 0;       /* total number of markers */
  themap->brdrs = (FPN) 0.0;  /* Will there be chromosomal material outsite the flanking chromosomal markers? 
                          If this is negative, then there will be no borders.  */
  themap->knum = NULL;     /* = ivector(1,traits) Number of QTLs on this map for each of the traits */
  themap->traits = 0;      /* Number of traits to simulate */
  themap->otraits = 0;     /* Number of other traits */
  themap->otypes = NULL;   /* classes for each other trait*/
  themap->tnames = NULL;   /* Names of the traits  = cmatrix(1,traits,0,MAXNAME) */
  
  themap->ParentalDiff = NULL;  /* differences of Parental means. */
  themap->onames = NULL;   /* Names of other traits  = cmatrix(1,otraits,0,MAXNAME)*/
  themap->names = NULL;    /* Names of the markers  = cmatrix(1,ml,0,MAXNAME)*/
  themap->cnames = NULL;   /* Names of the chromosomes  = cmatrix(1,m,0,MAXNAME)*/
  themap->ttable = NULL;   /* Table to indicate where in names the marker name is  = imatrix(1,m,1,maxl) */
  themap->types = NULL;    /* Indicate the type of marker:  = imatrix(1,m,1,maxl)
                              -1  =>  a-           a dom
                               0  =>  codominant
                               1  =>  A-           A dom           */
  return(themap);
}




/*
 * free up the space used by the map of markers.
 */
void deallocate_markermap(markermap *themap)
{
  int ii, maxl;
  if ( themap->otypes != NULL )
    free_ivector(themap->otypes , 1, themap->otraits);
  if (themap->knum != NULL)
    free_ivector(themap->knum, 1, themap->traits);
  if (themap->names != NULL)
    free_cmatrix(themap->names, 1, themap->ml, 0, MAXNAME);

  if (themap->tnames != NULL)
    free_cmatrix(themap->tnames, 1, themap->traits, 0, MAXNAME);
  if ( themap->ParentalDiff != NULL )
    free_dvector( themap->ParentalDiff,  1, themap->traits);
  if (themap->onames != NULL)
    free_cmatrix(themap->onames, 1, themap->otraits, 0, MAXNAME);

  if (themap->cnames != NULL)
    free_cmatrix(themap->cnames, 1, themap->m, 0, MAXNAME);

    maxl = 0;
    for (ii = 1; ii <= themap->m; ii++)
      if (maxl < themap->mpc[ii])
	    maxl = themap->mpc[ii];
  if (themap->mpc != NULL)
    free_ivector(themap->mpc, 1, themap->m);
  if (themap->types != NULL)
    free_imatrix(themap->types, 1, themap->m, 1, maxl);
  if (themap->ttable != NULL)
    free_imatrix(themap->ttable, 1, themap->m, 1, maxl);
  if (themap->mrf != NULL)
    free_dmatrix(themap->mrf, 1, themap->m, 0, maxl);
 if ( debugging > 2 ) {
        sprintf(gwarn,"In deallocate_markermap(), deallocated 1 markermap  at %p\n",themap);
        MemoryAccount(gwarn);
 }
  free((char *) themap);
}

/*
   determine if the file is mapmaker *.maps file.  
   
   If only two tokens on first line, then it is an Rmap.out file  10
   If second token is 'mapmaker', then it is a mapmaker.maps file 12
   If third token is 'bychromosome', then it is a map.inp file    11
   
*/

int det_if_mm(char *infile)
{
  FILE *fptr;
  int   nummarks, ch;
  char *chptr;
  nummarks = 0;
  fptr = fileopen(infile, "r");
  if (fptr == NULL)  /* No input file...bail out  */
    nummarks = -1;
  if ( nummarks == 0 ) {
    ch = get_next_line(gbuffer, MAXLINE, fptr);
    if ((chptr = strstr(gbuffer, "bychromosome")) != NULL) /* This is a map.inp file*/
      nummarks = 11;
    else if ((chptr = strstr(gbuffer, "mapmaker")) != NULL)  /* This is a mapmaker.maps file*/
      nummarks = 12;
  }
  if ( nummarks == 0 ) {
    ch = get_next_line(gbuffer, MAXLINE, fptr);
    if ((chptr = strstr(gbuffer, "Cartographer")) != NULL)  /* This is an Rmap.out file*/
      nummarks = 10;
  }
  fileclose(infile, fptr);
  return (nummarks);
}

/*
  Translate a mapmaker.maps file.
*/
markermap *trans_maps(char *infile)
{
  FILE *fptr;
  markermap *themap;
  int ch,cc,ii,jj,kk,go_on,markcount;
  char *chptr, **mnames;
  markcount = -1;
  themap = allocate_markermap();
  fptr = fileopen(infile, "r");
  if (fptr == NULL) {
	 deallocate_markermap(themap);
	 return (NULL);
  }
  do { /* get number of chromosomes and total number of markers */
    ch = get_next_token(gname, MAXNAME, fptr);
    if ( ch == EOF )
      themap->m = -1;
    if ( !strncmp(gname , "*Chromosome:",10) ) {
      ch = get_next_token(gname, MAXNAME, fptr);
      themap->m = atoi(gname);
    }
    else if ((chptr = strstr(gname, "*OrderInfo:")) != NULL ) 
      markcount = 0;
    else if ( markcount >= 0 && (chptr = strstr(gname, "*Classes:")) != NULL )
      markcount = -markcount;
    else if ( markcount >= 0 && *(gname+0) == '*' )
      markcount = markcount+1;
    
  } while (themap->m == 0);
  markcount = -markcount; /* This will be the number of markers in the raw file */
  if ( themap->m < 0 ) {
	 fileclose(infile, fptr);
	 deallocate_markermap(themap);
	 return (NULL);
  }
  if ( markcount > 0 )
    mnames =  cmatrix(0,markcount-1,0,MAXNAME);

  themap->mpc = ivector(1,themap->m);
  themap->cnames = cmatrix(1,themap->m,0,MAXNAME);

  cc = 1;
  do {  /* get mpc and chrom names */
    ch = get_next_token(gname, MAXNAME, fptr);
    if ( *(gname+0) == '*' ) {
      strcpy(themap->cnames[cc], gname);
      ch = get_next_token(gname, MAXNAME, fptr);
      themap->mpc[cc] = atoi(gname);
      themap->ml = themap->ml + themap->mpc[cc];
      if ( themap->mpc[cc] > themap->maxl )
        themap->maxl = themap->mpc[cc];
      cc = cc+1;
    }
  } while ( cc <= themap->m );

  themap->names = cmatrix(1,themap->ml,0,MAXNAME);   
  themap->ttable = imatrix(1,themap->m,1,themap->maxl);
  themap->mrf = dmatrix(1,themap->m,0,themap->maxl);

  for ( ii = 1 ; ii <= themap->m ; ii++ )
    for ( jj = 0 ; jj <= themap->maxl ; jj++ )
	  themap->mrf[ii][jj] = (FPN) 0.0;
  fseek(fptr, 0L, SEEK_SET);
  go_on = 1;
  do { 
    ch = get_next_line(gbuffer, MAXLINE, fptr);
    if ((chptr = strstr(gbuffer, "*OrderInfo:")) != NULL)
      go_on = 0;

  } while ( go_on == 1);

  for ( ii = 0 ; ii < markcount ; ii++ ) { /* Put marker names in mnames */
    ch = get_next_line(gbuffer, MAXLINE, fptr);
    get_field(1, gname, gbuffer);
    strcpy( mnames[ii], gname );
  }
  go_on = 1;

  do { 
    ch = get_next_line(gbuffer, MAXLINE, fptr);
    if ((chptr = strstr(gbuffer, "*Chromosomes:")) != NULL)
      go_on = 0;
  } while ( go_on == 1);
  jj = 1;
  for ( cc = 1 ; cc <= themap->m ; cc++ ) {  /* get the recombination frequencies */
    do {  /*get tokens until one starts with an asterisk. */
      ch = get_next_token(gname, MAXNAME, fptr);
    } while ( gname[0] != '*' );
    strcpy(themap->cnames[cc], gname);
    for ( ii = 1 ; ii <= 5 ; ii++ )  /* Ignore 5 numbers */
      ch = get_next_token(gname, MAXNAME, fptr);
    for ( ii = 1 ; ii <= themap->mpc[cc] ; ii++ ) { /* Get integers that tell us the marker names */
      ch = get_next_token(gname, MAXNAME, fptr);
      themap->ttable[cc][ii] = jj;
      jj = jj+1;
      kk = atoi( gname );
      strcpy( themap->names[ themap->ttable[cc][ii] ] , mnames[kk] );
    }
    for ( ii = 1 ; ii< themap->mpc[cc] ; ii++ ) {  /* Get marker recombination frequencies */
      ch = get_next_token(gname, MAXNAME, fptr);
      themap->mrf[cc][ii] = (FPN) atof(gname);
    }
  }

  for ( ii = 1 ; ii <= themap->ml ; ii++ ) 
    for ( jj = 0 ; jj < MAXNAME ; jj++ )  /* Get rid of the * that begins the names. */
      themap->names[ii][jj] = themap->names[ii][jj+1] ;
  for ( ii = 1 ; ii <= themap->m ; ii++ ) 
    for ( jj = 0 ; jj < MAXNAME ; jj++ )
      themap->cnames[ii][jj] = themap->cnames[ii][jj+1];

  fileclose(infile, fptr);
  free_cmatrix(mnames,0,markcount-1,0,MAXNAME);
  themap->sigl = (FPN) 1.0;
  do_map_stats(themap);    
  return (themap);
}

/*
  Translate a map.inp file.
  
Tokens:
  -type     [positions || intervals]        for postions, Telomere can be last marker
                                         for intervals, Telomere can be first marker
  -function  [integer value] which mapping function   values 1-8                                     
  -parameter    [real value]  Extra parameter for map functions [4-8]
  -Units   [cM || M || r]
  -chromosomes    [integer value] # of chromosomes
  -maximum   [integer value] max. markers on any chrom.
  -end  or -quit   stop reading
  -named    [yes || no]       are markers named?
  -start
  -stop    end of map
  -skip   to -unskip    
  -Chromosome   [integer || chromosome name]    
  
  Telomere is a special marker name indicating a non-typed telomere
*/
markermap *trans_stdm(char *infile)
{
  FILE *fptr;

  markermap *themap;
  int  ch;
  int chroms, maxmark,pchrom;
  int whchrom, themark, marker;
  FPN  telomeric;
  int  ii, jj, ist, isf, isu, isc, ism, isn, iss, isC, ise, isskip, go_on;
  whchrom = 0;
  fptr = fileopen(infile, "r");
  if (fptr == NULL) 
    return (NULL);
  go_on = 1;
  whchrom = themark = marker = 0;
  ist = isf = isu = isc = ism = isn = iss = isC = ise = isskip = 0;
  while (go_on == 1) {
    ch = get_next_token(gname, MAXNAME, fptr);

    if ( gname[0] == '-' && isskip == 0)
      switch ( gname[1] ) {
       case 't':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     if (!(strncmp(gname, "position", 8)))
	       ist = -1;
	     else
	       ist = -2;
	   break;
       case 'f':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     whosemf = atoi(gname);
	   break;
       case 'p':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     mapparam = (FPN) atof(gname);
	   break;
       case 'U':
	     ch = get_next_token(gname, MAXNAME, fptr);
		 switch (gname[0]) {
		   case 'M':  case 'm':  isu = -1;  break;
		   case 'c':  case 'C':  isu = -2;  break;
		   case 'r':  case 'R':  isu = 0;   break;
		   default:  isu = -2;  break;
		 }
		 break;
       case 'u': isskip = 0;  break;
       case 'c':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     chroms = atoi(gname);
	   break;
       case 'm':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     maxmark = atoi(gname);
	     break;
       case 'n':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     if (gname[0] == 'n' || gname[0] == 'N')
	       isn = 0;
	     else
	       isn = 1;
	   break;
       case 's':
	     if (gname[3] == 'a')
	       iss = 3;
	     else if ( gname[3] == 'i') {
	       iss = 0;
	       isskip = 1;
	     }
	     else {
	       if (whchrom > 0)
	         themap->mpc[whchrom] = themark;
	       go_on = 0;
	      iss = 0;
	     }
	   break;
       case 'C':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     if (whchrom > 0)
	       themap->mpc[whchrom] = themark;
           pchrom = is_pinteger(gname);
         if ( pchrom > 0 && pchrom <= chroms )
	       whchrom = pchrom;
         else {
           whchrom = whchrom+1;
           if ( themap->cnames == NULL )
             themap->cnames = cmatrix(1,chroms,0,MAXNAME);
           if ( whchrom <= chroms )
             strcpy( themap->cnames[whchrom], gname );
           else {
             go_on = 0;
             deallocate_markermap(themap);
             themap = NULL;
           }
        }
	    themark = 0;
	    ch = get_next_token(gname, MAXNAME, fptr);
	    iss = 1;
	   break;
       case 'e':
	    ise = 1;
	   break;
       default:
	   break;
      }
      else if (!strcmp(gname, "-unskip"))
      iss = isskip = 0;

    if (iss == 1 && isn == 1) {
      if ( !strcmp( "Telomere", gname ) ) {
         ch = get_next_token(gname, MAXNAME, fptr);
         themap->mrf[whchrom][0] = mapfunc( (FPN) atof(gname), isu );
      }
      else {
        themark = themark + 1;
        marker = marker + 1;
        strcpy( themap->names[marker], gname);
        themap->ttable[whchrom][themark] = marker;
        ch = get_next_token(gname, MAXNAME, fptr);
        themap->mrf[whchrom][themark] = mapfunc((FPN) atof(gname), isu );
      }
    }
    else if (iss == 1 && isn == 0) {
      themark = themark + 1;
      marker = marker + 1;
      themap->mrf[whchrom][themark] = mapfunc((FPN) atof(gname), isu );
    }
    else if (iss == 3) {
/* allocate memory, etc */
      themap = allocate_markermap();
      themap->ml = chroms * maxmark;
      themap->maxl = maxmark;
      themap->m = chroms;
		themap->mrf = dmatrix(1, themap->m, 0, themap->maxl);
		for ( ii = 1 ; ii <= themap->m ; ii++ )
		  for ( jj = 0 ; jj <= themap->maxl ; jj++ )
			 themap->mrf[ii][jj] = (FPN) 0.0;
      themap->mpc = ivector(1, chroms);
      themap->brdrs = (FPN) 0.0;
      themap->sigl = (FPN) 1.0;
      themap->s = (FPN) 1.0;
      themap->sigs = (FPN) 1.0;
      themap->knum = ivector(1, 1);
      themap->traits = 1;
      themap->knum[1] = 1;
      if (isn == 1) {
	    themap->names = cmatrix(1, themap->ml, 0, MAXNAME);
	    themap->ttable = imatrix(1, themap->m, 1, themap->maxl);
      }
      iss = 0;
      marker = 0;
    }
  }
  fileclose(infile, fptr);
  if (ist == -1)
    for (ii = 1; ii <= themap->m; ii++) {
      telomeric = (FPN) 0.0;
      if ( themap->mrf[ii][0] > (FPN) 0.0 ) {
        telomeric = themap->mrf[ii][0];
        themap->mrf[ii][0] = (FPN) 0.0;
      }
      for (jj = 0; jj < themap->mpc[ii]; jj++)
	    themap->mrf[ii][jj] = mapfunc(mapfunc(themap->mrf[ii][jj+1], 1) - mapfunc(themap->mrf[ii][jj], 1), -1);
      if ( telomeric > (FPN) 0.0 ) 
        themap->mrf[ii][themap->mpc[ii]] = mapfunc(mapfunc(telomeric, 1) - mapfunc(themap->mrf[ii][themap->mpc[ii]], 1), -1);
      else
        themap->mrf[ii][themap->mpc[ii]] = (FPN) 0.0;
    }

  do_map_stats(themap);
  return (themap);
}



int what_chrom(char *lastline)
{
  int chrom, ii;
  for (ii = 0; lastline[ii] != 'r'; ii++)
    lastline[ii] = ' ';
  lastline[ii] = ' ';
  if ( isdigit(lastline[ii+1]) )
    chrom = atoi(lastline);
  else
    chrom = -1;
  return (chrom);
}


void print_marker_names(markermap *themap,FILE *fptr)
{
  int chrom, mark, marks;
  marks = 0;
  if ( themap->names != NULL ) {
    fprintf(fptr, "\n\nNames and positions of the markers \n\n");
    fprintf(fptr, " Chrom  Mark  Name \n-b MarkerNames    ");
    for (chrom = 1; chrom <= themap->m; chrom++) {
      for (mark = 1; mark <= themap->mpc[chrom]; mark++) {
        if (themap->ttable != NULL)
          marks = themap->ttable[chrom][mark];
        else
	  marks = marks + 1;
        fprintf(fptr, "\n %5d %5d %s", chrom, mark, themap->names[marks]);
      }
    }
    fprintf(fptr, "\n-e MarkerNames\n\n");
  }
  if ( themap->cnames != NULL ) {
     fprintf(fptr,"\n\nNames of the Chromosomes\n\n-b  ChromosomeNames    ");
     for ( chrom = 1 ; chrom <= themap->m ; chrom++ )
      fprintf(fptr, "\n %5d            %s", chrom,  themap->cnames[chrom]);
     fprintf(fptr, "\n-e ChromosomeNames\n");
  }
}


void get_traitnames(params *theparams,markermap *themap)
{
  FILE *infile;
  int ii,ch,go_on;
  go_on=1;
  
  if ( (ch=isfile(theparams->ifile)) == 1 )
    infile = fileopen(theparams->ifile, "r");
  else 
    infile = NULL;
  if ( infile != NULL ) {
    while ( go_on == 1 ) {
      ch=get_next_line(gbuffer,MAXLINE,infile);
      if (gbuffer[1] == 'N') {
          if ( themap->tnames == NULL  )
           themap->tnames = cmatrix(1,  themap->traits, 0, MAXNAME );
	      for ( ii = 1 ; ii <= themap->traits ; ii++ ) {
	        ch = get_next_token(gbuffer, MAXLINE, infile);
	        ch = get_next_token(gbuffer, MAXLINE, infile);
            strcpy( themap->tnames[ii], gbuffer );
          }
          
      }
      else if (gbuffer[1] == 'P') {
          if ( themap->ParentalDiff == NULL  )
           themap->ParentalDiff = dvector(1,  themap->traits  );
	      for ( ii = 1 ; ii <= themap->traits ; ii++ ) {
	        ch = get_next_token(gbuffer, MAXNAME, infile);
            themap->ParentalDiff[ii] = (FPN) atof( gbuffer );
          }
          
      }
      if ( ch == EOF || gbuffer[1] == 's' )
        go_on = 0;
    }
    fileclose(theparams->ifile, infile);
  }
}


/*  The three routines below will create names for markers, traits and categorical traits.
    These names will be required for printouts.   */

/* Create a set of names for the markers*/
void    create_bogus_markers(params *theparams,markermap *themap)
{
  int c,m,i;
  if (themap->names == NULL ) {
    themap->names = cmatrix(1,themap->ml,0,MAXNAME);
    if ( themap->ttable == NULL )
      themap->ttable = imatrix(1,themap->m,1,themap->maxl);
    i = 0;
    for ( c = 1 ; c <= themap->m ; c++ ) { 
	  for (m = 1; m <= themap->mpc[c]; m++) {
	    i = i+1;
        sprintf(themap->names[i], "c%dm%d", c,m);
        themap->ttable[c][m] = i;
      }
     }
   }
   create_bogus_traits(theparams, themap );
   create_bogus_otraits(theparams, themap );
}

/* Create a set of names for the traits.   */
void    create_bogus_traits(params *theparams,markermap *themap)
{
  int  i ;
  if ( theparams->traits != themap->traits )
    themap->traits = theparams->traits;
  if ( themap->tnames == NULL && themap->traits > 0 ) {
     themap->tnames = cmatrix(1,themap->traits,0,MAXNAME);
     for ( i = 1 ; i <= themap->traits ; i++ )
       sprintf(themap->tnames[i], "t%d",i);
  }
}

/*  Create a set of names for the other traits.  */
void    create_bogus_otraits(params *theparams,markermap *themap)
{
  int  i ;
  if ( theparams->traits != themap->traits )
    themap->traits = theparams->traits;
   if ( themap->onames == NULL && themap->otraits > 0 ) {
     themap->onames = cmatrix(1,themap->otraits,0,MAXNAME);
     for ( i = 1 ; i <= themap->otraits ; i++ )
       sprintf(themap->onames[i], "o%d",i);
   }
}


/*
   Create a map of markers from a genome.   This assumes that the marker types
   are encoded in whichqtl.   
*/
markermap *mapfromgenome(genome *thegenome) {
  int  i;
  markermap *themap;
  genome *gptr;
  themap = allocate_markermap();
  
  themap->m = 0;
  for ( gptr=thegenome; gptr != NULL ; gptr=gptr->next ) 
    if (gptr->chrom > themap->m )
      themap->m = gptr->chrom;
  themap->mpc = ivector(1,themap->m); 
  themap->ml = 0;
  for ( gptr=thegenome; gptr != NULL ; gptr=gptr->next ) {
    themap->ml +=1;
    themap->mpc[gptr->chrom] +=1;
  }
  themap->maxl = 0;
  for (i=1; i<=themap->m ; i++ )
    if ( themap->mpc[i] > themap->maxl )
      themap->maxl = themap->mpc[i];  
  themap->names = cmatrix(1, themap->ml, 0, MAXNAME);
  themap->ttable = imatrix(1, themap->m, 1, themap->maxl);

  themap->types = imatrix(1, themap->m, 1, themap->maxl);
  themap->mrf = dmatrix(1, themap->m, 0, themap->maxl);
  i=1;
  for ( gptr=thegenome; gptr != NULL ; gptr=gptr->next ) {
    strcpy(themap->names[i],gptr->markername) ;
    themap->ttable[gptr->chrom][gptr->markr] = i;
    i +=1;
    themap->mrf[gptr->chrom][gptr->markr] = gptr->dist;
    themap->types[gptr->chrom][gptr->markr] = gptr->whichqtl;
  }
  themap->traits = 1;
  themap->knum = NULL;

  do_map_stats(themap);


  return(themap);
}

/*
  This allocates space for a markermap.  
*/
markermap *bogus_markermap(int chrom,int mark,int maxl,int traits,FPN dist)
{
  int ii,jj;
  markermap *themap;

  themap = allocate_markermap();
  themap->traits = traits;
  themap->knum = ivector(1, themap->traits);
  themap->knum[1] = traits;
  themap->m =  chrom;
  themap->l = maxl;
  themap->s = dist;
  themap->sigs =  (FPN) 1.0;
  themap->ml = maxl*chrom;
  themap->maxl = maxl;
  themap->sigl = (FPN) 1.0;
  themap->mpc = ivector(1,chrom);
  for (ii = 1; ii <= chrom; ii++) 
    themap->mpc[ii] = mark;
  themap->mrf = dmatrix(1, chrom, 0, maxl);
  for ( ii = 1 ; ii <= chrom ; ii++ ) {
     themap->mrf[ii][0] =  themap->mrf[ii][maxl] = (FPN) 0.0;
    for ( jj = 1 ; jj < maxl ; jj++ )
      themap->mrf[ii][jj] = dist;
  }
  themap->ttable = imatrix(1,themap->m,1,themap->maxl  );
  themap->names = cmatrix(1, themap->ml, 0, MAXNAME);
  themap->tnames = cmatrix(1, themap->traits, 0, MAXNAME);
  themap->types = imatrix(1, chrom, 1, themap->maxl);
  for (ii = 1; ii <= themap->m; ii++)
    for (jj = 1; jj <= themap->maxl; jj++)
      themap->types[ii][jj] = 0;


  return(themap);

}
/*
   Open a cross.inp file and determine how many markers there are.
   If they have names, then get the names as well.  Return a 
   markermap with the markers on a single chromosome, ordered
   as they appear in the input file and evenly spaced at 10cM.  
*/
markermap *crossbogusmap(char *infile) {
  FILE *fptr;
  markermap *themap; 
  int ml,traits,ch,go_on,nn,named,i,j;
  fptr = fileopen(infile, "r");
  go_on = 1;
  ml = traits = nn = 0;  
  
  
  ch = get_next_token(gname, MAXNAME, fptr);  /*  first token.  */
  if (!strcmp("#FileID",gname)   ) 
    go_on = MoveToToken(fptr, "#bycross", 1 );

  while ( go_on == 1 ) {  /* first pass, get the sample size, number of markers and traits */
     ch = get_next_token(gname, MAXNAME, fptr);
     if (!strcmp(gname, "-traits")) {
            ch = get_next_token(gname, MAXNAME, fptr);
            traits = atoi(gname);
     }
     else if (!strcmp(gname, "-SampleSize")) {
            ch = get_next_token(gname, MAXNAME, fptr);
            nn = atoi(gname);
     }
     else if ( !strcmp(gname, "-start")) {
       named = 1;
       ch = get_next_token(gname, MAXNAME, fptr);
       if  ( !strcmp(gname, "individuals") ) {
         named -=1;
         ch = get_next_token(gname, MAXNAME, fptr);
       }
       if ( !strcmp(gname, "markers") ) {
         do {
           ch = get_next_token(gname, MAXNAME, fptr);
           ml +=1;
           if ( !strcmp(gname, "-stop")) {
             go_on = 0;
             ml = ml / (nn+1);
           }
         } while (go_on == 1 ); 
       }
     }
  }


  fileclose(infile, fptr);
   
  themap = bogus_markermap(1, ml, ml,traits,(FPN) 0.1);
  if ( named == 0 )
    for ( i=1; i<=ml; i++ ) {
      sprintf(themap->names[i],"m%d",i);
      themap->ttable[1][i] = i;
    }
  else {  /*  Go back and get the marker names. */
    fptr = fileopen(infile, "r");
    ch = get_next_token(gname, MAXNAME, fptr);  /*  first token.  */
    if (!strcmp("#FileID",gname)   ) 
      go_on = MoveToToken(fptr, "#bycross", 1 );
    go_on = 1;
    while ( go_on == 1 ) {
      ch = get_next_token(gname, MAXNAME, fptr);
      if ( !strcmp(gname, "-start")) {
       ch = get_next_token(gname, MAXNAME, fptr);
       if ( !strcmp(gname, "markers") ) {
         for ( i=1; i<=ml; i++ ) {
           ch = get_next_token(themap->names[i], MAXNAME, fptr);
           themap->ttable[1][i] = i;

           for ( j=1; j<=nn; j++ )
             ch = get_next_token(gname, MAXNAME, fptr);
          }  
          go_on = 0; 
       }
     }
    }
    fileclose(infile, fptr);
  }
  return(themap);
}




/*
   Open an Rcross.out file, get the number of markers and
   number of traits, and create a bogus map.  Return a 
   markermap with the markers on a single chromosome, ordered
   as they appear in the input file and evenly spaced at 10cM.
*/
markermap *rcrossbogusmap(char *infile) {
  FILE *fptr;
  markermap *themap; 
  int ml,traits,ch,go_on,i;
  fptr = fileopen(infile, "r");
  go_on = 1;
  ml = traits = 0;  
  while ( go_on == 1 ) {  /* first pass, get the sample size, number of markers and traits */
     ch = get_next_token(gname, MAXNAME, fptr);
     if (!strcmp(gname, "-traits")) {
            ch = get_next_token(gname, MAXNAME, fptr);
            traits = atoi(gname);
            go_on = 0;
     }
     else if (!strcmp(gname, "-p")) {
            ch = get_next_token(gname, MAXNAME, fptr);
            ml = atoi(gname);
            ml -=1;
     }
  }


  fileclose(infile, fptr);
   
  themap = bogus_markermap(1, ml, ml,traits,(FPN) 0.1);
  for ( i=1; i<=ml; i++ ) {
    sprintf(themap->names[i],"m%d",i);
    themap->ttable[1][i] = i;
  }
  return(themap);
}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Mdatain.c
------------------------------------------------------------------ */


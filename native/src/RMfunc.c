/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RMfunc.c) is part of QTL Cartographer
         
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
/* Functions to simulate a genetic linkage map.*/

/*
  Create a set of files from the genetic linkage map that would allow
  one to see the map using GNUPLOT.  
  
  As of September, 2002, 
  
  The marker names are placed on the map.
  Chromosomes are labelled 'Chromosome 1', 'Chromosome 2', etc . 
  The plot command and all plot data are in the same file.
  The output file is the filename stem plus "map.plt".  
     For a stem of corn, this means cornmap.plt would have the map.
  There are commented lines in the plot file that allow
     one to specify postscript output to a file.  The default output file
     is stem with "map.ps", i.e.  cornmap.ps for the above example. 
  
  
*/
void gnuplot_map(params *theparams,markermap *themap)
{
  char  *gnufile ;
  FILE *goutput;
  int ii, jj, lastmarker;
  FPN position,maxposition;

  gnufile = cvector(0,MAXNAME);
    sprintf(gnufile, "%s%s", theparams->stem, "map.plt");
  goutput = fileopen(gnufile, "w");
/*  Find the longest chromosome */
  maxposition = (FPN)  0.0;
  for (ii = themap->m; ii > 0; ii--) {
      lastmarker = themap->mpc[ii];
    position = (FPN) 0.0;
    if (themap->brdrs > (FPN) 0.0) 
      position = mapfunc(themap->mrf[ii][0], 1);
    else
      lastmarker = lastmarker - 1;
    for (jj = 1; jj <= lastmarker; jj++)
      position = position + mapfunc(themap->mrf[ii][jj], 1);
    if ( position > maxposition )
      maxposition = position;
  }
  fprintf(goutput, "# For use in gnuplot.");
  fprintf(goutput, "\n#  This will plot the markers at their positions on the chromosomes.");
  fprintf(goutput, "\n#  All the data to be plotted are now in this file.");
  write_file_type(13,goutput);
  fprintf(goutput, "\nset term %s",theparams->term);
  fprintf(goutput, "\n#Delete the above line and uncomment the next two for postscript");
  fprintf(goutput, "\n#set term postscript\n#set output \"qtlcartmap.ps\"" );
  fprintf(goutput, "\nset autoscale");
  fprintf(goutput, "\nset data style linespoints\nset noyzeroaxis");
  fprintf(goutput, "\nset noytics\nset noborder\nset nokey\nset yrange[0:%d]", themap->m+1 );
  fprintf(goutput, "\nset xrange[0:%f]", maxposition*1.15 );
  fprintf(goutput, "\nset title \"Molecular Map\"");
  fprintf(goutput, "\nset xlabel \"Position in Morgans\"");



/*  First, put all the labels where they belong   */
  for (ii = themap->m; ii > 0; ii--) {
      lastmarker = themap->mpc[ii];

    position = (FPN) 0.0;
    if (themap->brdrs > (FPN) 0.0) {
      fprintf(goutput,"\nset label \"Telomere\" at %f,%f rotate", position, 0.1+(FPN) ii);
      position = mapfunc(themap->mrf[ii][0], 1);
    }
    else
      lastmarker = lastmarker - 1;
    for (jj = 1; jj <= lastmarker; jj++) {
      fprintf(goutput,"\nset label \"%s\" at %f,%f rotate",themap->names[themap->ttable[ii][jj]],position, 0.1+(FPN) ii);
      position = position + mapfunc(themap->mrf[ii][jj], 1);
    }
    if (themap->brdrs > (FPN) 0.0)
      fprintf(goutput,"\nset label \"Telomere\" at %f,%f rotate",position, 0.1+(FPN) ii);
    else
      fprintf(goutput,"\nset label \"%s\" at %f,%f rotate",themap->names[themap->ttable[ii][jj]],position, 0.1+(FPN) ii);
    if ( themap->cnames != NULL )
      fprintf(goutput,"\n set label \"%s\" at %f,%f", themap->cnames[ii], position*1.01 , (FPN) ii );    
    else
      fprintf(goutput,"\n set label \"Chromosome %d\" at %f,%f", ii, position*1.01 , (FPN) ii );

  }

/* Now, print the plot line...*/
  fprintf(goutput, "\nplot ");
  for (ii = themap->m; ii > 0; ii--) {
    fprintf(goutput, "\"-\"");  
    if ( ii > 1 )
      fprintf(goutput, ", ");  
  }
/* Now, print data...*/
  for (ii = themap->m; ii > 0; ii--) {
      lastmarker = themap->mpc[ii];
    fprintf(goutput, "#  This is marker position data for chromosome %d\n",  ii);

    position = (FPN) 0.0;
    if (themap->brdrs > (FPN) 0.0) {
      fprintf(goutput, "%9.4f %9.4f\n", position, (FPN) ii);
      position = mapfunc(themap->mrf[ii][0], 1);
    }
    else
      lastmarker = lastmarker - 1;
    for (jj = 1; jj <= lastmarker; jj++) {
      fprintf(goutput, "%9.4f %9.4f\n", position, (FPN) ii);
      position = position + mapfunc(themap->mrf[ii][jj], 1);
    }
    fprintf(goutput, "%9.4f %9.4f\ne\n", position, (FPN) ii);
    
  }

/*  output = fileopen(gnufile, "a");*/
  fprintf(goutput, "\npause -1 \"Hit return to continue\" \n");
  fileclose(gnufile, goutput);
  free_cvector( gnufile,0,MAXNAME);

}






/*
 * Create the space for a map of markers, and place the markers on it.
 *  This version is for simulated maps: The parameters are taken from
 *  theparams. 
 */

markermap *initailize_markermap(int kk,params *theparams)
{
  int ii, maxl;
  markermap *themap;
  themap = allocate_markermap();
  themap->traits = 1;
  themap->knum = ivector(1, themap->traits);
  themap->knum[1] = kk;
  themap->m = theparams->chrom;
  themap->l = theparams->mark;
  themap->sigl = theparams->vmark;
  themap->s = theparams->dist/(FPN) 100.0;
  themap->sigs = theparams->vdist/(FPN) 100.0;
  themap->mpc = ivector(1,theparams->chrom);
  if (theparams->vmark <= (FPN) 0.0) {
    themap->ml = theparams->chrom * theparams->mark;
    maxl = theparams->mark;
    for (ii = 1; ii <= theparams->chrom; ii++)
      themap->mpc[ii] = theparams->mark;
  }
  else {
    themap->ml = 0;
    maxl = 0;
    for (ii = 1; ii <= theparams->chrom; ii++) {
      themap->mpc[ii] = mkr_num(theparams->mark, theparams->vmark);
      if (maxl < themap->mpc[ii])
	    maxl = themap->mpc[ii];
      themap->ml = themap->ml + themap->mpc[ii];
    }
  }
  themap->maxl = maxl;
  themap->mrf = dmatrix(1, theparams->chrom, 0, maxl);
  themap->brdrs = theparams->tail/(FPN) 100.0;
/*  This part needs to be documented.   You can simulate the 
  markermap in two ways:   The first is the old standard, and
  the second is to set the size of the chromosomes first, and 
  put the markers on second.    */
  if ( theparams->Rmode == 0 )
    mrkr_rec(themap);
  else
    mrkr_rec2(themap); /*
  themap->s = theparams->dist*(FPN) 100.0;
  themap->sigs = theparams->vdist*(FPN) 100.0;*/
  themap->brdrs = theparams->tail*(FPN) 100.0;
  do_map_stats(themap);
  return(themap);
}



/*
 * Determine the number of markers for each chromosome. If ssigl <= 0.0 then
 * there will be exactly ll markers on each chromosome.  If not, the number
 * will be N(ll,ssigl).
 */
int mkr_num(int ll,FPN ssigl)
{
  int itemp,k;
  k=0;
  if (ssigl > 0)
    while ((itemp = (int) (ssigl * gasdev(&ll) + ll)) < 1) k+=1;
  else
    itemp = ll;
  return (itemp);
}



/*
 * Determine the recombination frequency between consecutive markers on each
 * chromosome.  If themap->ssigs <= 0.0, then there will be exactly s Morgans
 * between consecutive markers and the recombination frequency will be (1 -
 * exp(-2s))/2.  If themap->sigs > 0.0, then s will have a N(s,sigs)
 * distribution and the recombination frequencies will again be calculated by
 * the Haldane mapping function.
 */

void mrkr_rec(markermap *themap)
{
  int ii, jj;
  FPN ss, scale;
  if (themap->brdrs > (FPN) 0.0)
    scale = (FPN) pow((themap->brdrs / themap->s),  2.0);

    if (themap->sigs > (FPN) 0.0)
      for (ii = 1; ii <= themap->m; ii++) {
	if (themap->brdrs > 0) {
	  ss = scale * themap->sigs * gasdev(&ii) + themap->brdrs;
	  if (ss <= (FPN) 0.0 || ss >= 0.5)
	    ss = themap->brdrs;
	  themap->mrf[ii][0] = mapfunc(ss, -1);
	  ss = scale * themap->sigs * gasdev(&ii) + themap->brdrs;
	  if (ss <= (FPN) 0.0 || ss >= 0.5)
	    ss = themap->brdrs;
	  themap->mrf[ii][themap->mpc[ii]] = mapfunc(ss, -1);
	}
	for (jj = 1; jj < themap->mpc[ii]; jj++) {
	  ss = themap->sigs * gasdev(&ii) + themap->s;
	  if (ss <= (FPN) 0.0 || ss >= 0.5)
	    ss = themap->s;
	  themap->mrf[ii][jj] = mapfunc(ss, -1);
	}
      }
    else
      for (ii = 1; ii <= themap->m; ii++) {
	if (themap->brdrs > 0)
	  themap->mrf[ii][themap->mpc[ii]] = themap->mrf[ii][0] = mapfunc(themap->brdrs, -1);
	for (jj = 1; jj < themap->mpc[ii]; jj++)
	  themap->mrf[ii][jj] = mapfunc(themap->s, -1);
      }
    if (themap->brdrs <= (FPN) 0.0)
      for (ii = 1; ii <= themap->m; ii++)
	    themap->mrf[ii][0] = themap->mrf[ii][themap->mpc[ii]] = (FPN) 0.0;


}


/*
 * Determine the recombination frequency between consecutive markers on each
 * chromosome.  
 *  In this version, the length of the chromosomes will be N(s,ssigs).
 *  The markers will then be placed uniformly on the chromosome.
 */

void mrkr_rec2(markermap *themap)
{
  int ii, jj,nints;
  FPN ss, *intervals;
  
  
  intervals = dvector(0,themap->maxl+1);
  
  
  for ( ii = 1; ii <= themap->m ; ii++ ) {
    intervals[0] = (FPN) 0.0;

      nints = themap->mpc[ii]-1;

    if ( themap->brdrs > (FPN) 0.0 ) 
      nints +=2;
	do {
	  ss = themap->sigs * gasdev(&ii) + themap->s;  /* Length of the chromosome */
	} while ( ss <= (FPN) 0.0 );
	
	for ( jj = 1  ; jj<=nints ; jj++ )
	  intervals[jj] = ss * ranf(jj) ;
	sort(nints,intervals); 
	for ( jj = nints ; jj > 0 ; jj-- )
	  intervals[jj] = intervals[jj] - intervals[jj-1];
	
	
	if ( themap->brdrs > (FPN) 0.0 )
	  for ( jj = 1 ; jj <= nints ; jj++ ) 
	     themap->mrf[ii][jj-1] = mapfunc( intervals[jj], -1);
	else
	  for ( jj = 1 ; jj <= nints ; jj++ ) 
	     themap->mrf[ii][jj] = mapfunc( intervals[jj], -1);
  
  }
  
  free_dvector(intervals,0,themap->maxl+1);
}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RMfunc.c
------------------------------------------------------------------ */


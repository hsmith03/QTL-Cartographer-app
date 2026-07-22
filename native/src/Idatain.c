/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Idatain.c) is part of QTL Cartographer
         
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



/*  
          Idatain.c,  subroutines to read in data files, either from
          MAPMAKER *.raw files, cross.inp files or Rcross.out files.
*/


#include "Main.h"

/*
undo the trans_data function.
   
   GT     Rcross.out Read in as   Analysis format
   -------------------------------------------
   AA        2          3             1
   Aa        1          1             0
   aa        0          0            -1
   A-       12         12             2
   a-       10         10            -2
   --    all else      -3            -3
   -------------------------------------------   
   
*/
void untrans_data(int **dataptr,markermap *themap)
{
  int ii, jj;
  
  for (ii = 1; ii <= themap->m; ii++) {
    for (jj = 1; jj <= themap->mpc[ii]; jj++)
      switch (dataptr[ii][jj]) {
       case -2:   dataptr[ii][jj] = 10;  break;	
       case -1:   dataptr[ii][jj] =  0;  break;
       case  0:   dataptr[ii][jj] =  1;  break;
       case  1:   dataptr[ii][jj] =  3;  break;
       case  2:   dataptr[ii][jj] = 12;  break;
       default:   dataptr[ii][jj] = -3;  break;
    }
  }
}

/*
   Translate the data read in from the Rcross.out formatted data
   file into a format expected by the rest of the programs.  
   
   GT     Rcross.out Read in as  Translated to
   -------------------------------------------
   AA        2          3             1
   Aa        1          1             0
   aa        0          0            -1
   A-       12         12             2
   a-       10         10            -2
   --    all else      -3            -3
   -------------------------------------------   
   
   Note:  the old style had Aa encoded by 1 or 2.  If the data file is in
   the old style, 2's read in by get_the_datapoints are not changed to 3.
   New data files are read in as the above.
   
   Originally, a 1 and a 2 stood for aA and Aa heterozygotes, respectively.
   In the simulation programs, they still do.  Some people thought that
   genotypes would be better encoded via the last column. That's why we have this
   byzantine recoding of data points.
*/
void trans_data(int **dataptr,markermap *themap)
{
  int ii, jj;
  for (ii = 1; ii <= themap->m; ii++) {
    for (jj = 1; jj <= themap->mpc[ii]; jj++)
      switch (dataptr[ii][jj]) {
       case 0:   dataptr[ii][jj] = -1; break;
       case 1: 
       case 2:   dataptr[ii][jj] =  0;  break;
       case 3:   dataptr[ii][jj] =  1;  break;
       case 12:  dataptr[ii][jj] =  2;  break;
       case 10:  dataptr[ii][jj] = -2; break;	
       default:  dataptr[ii][jj] = -3;
    }
  }
}



/*
  Get data from the specified file.
  
  This is general, and allows any program to read any of the appropriate input formats.
*/
void GetTheData(params *theparams, char *infile, int translate) {
  int ii,i,error;
  FILE *errorf;

/* First, determine the type of file... */
    ii = get_file_type(infile);
    if ( ii == 0 || ii == -1) 
      ii = recheck_file_type(infile);
    if ( ii == 30 ) { /*Rcross.out*/
      if ( theparams->themap == NULL ) /* need to create a bogus map */
        theparams->themap = rcrossbogusmap(infile);
      error = get_the_nn(&theparams->nn, &theparams->traits,  infile);
      theparams->cross = get_cross(infile,theparams);
      theparams->themap->traits = theparams->traits;
/*      theparams->themap->knum = NULL;*/
      if ( whichprogram == 10 )
        theparams->thedata = indvector(theparams->nn+1, theparams->themap, theparams->theqtls );
      else
        theparams->thedata = indvector(theparams->nn, theparams->themap, theparams->theqtls );
      error = get_the_datapoints(theparams->thedata, infile, theparams->nn);
    }
    else if ( ii == 32 )   /* mapmaker.raw */
      theparams->thedata = get_mm_data(theparams, infile ); 
    else if ( ii == 31 || ii == 4 ) { /*  This is either a cross.inp file or a qtlcart.mcd file. */ 
      theparams->thedata = get_std_data(theparams,    infile);
    }
    else if ( ii == 35 || ii == 36)   /* PLABQTL file */
      theparams->thedata = convert_plabqtl(theparams);
    else { /*can't figure it out, so bail.*/
      errorf = fileopen(theparams->error, "a");
      if (errorf != NULL) {
        fprintf(errorf, "\n\nBailing out, because filetype of ...%s... can't be determined.\n\n ", infile);      
        fileclose(theparams->error, errorf);
      }
      exit(1);
    }    
    determine_markers(theparams->thedata,theparams->nn);
    for ( i=1; i<=theparams->nn; i++ )/*  This may now be redundant.  */
        theparams->thedata[i].map = theparams->themap; 

	if ( translate == 1 ) {	    /* Translate from 0, (1,2), 3 to -1,0,1 */
      for (ii = 1; ii <= theparams->nn; ii++)
        trans_data(theparams->thedata[ii].markers, theparams->themap);
    }
	theparams->themap->traits = theparams->traits;
}



/* determine_markers    determines whether the markers are dominant, recessive or
   non-existant.  It sets the flags in themap->types = imatrix(1,themap->m , 1, themap->maxl) to

   -99  bad system (not enough data)
    -1  a dominant
     0  no dominance
     1  A dominant                 */

void determine_markers(individual *thedata,int nn)
{
  markermap *themap;
  int ii, jj, kk,   nocntr ;
  themap = thedata[1].map;
  if (themap->types == NULL )
      themap->types = imatrix(1, themap->m, 1, themap->maxl);

  if (themap->types != NULL)
    for (ii = 1; ii <= themap->m; ii++) {
      for (jj = 1; jj <= themap->mpc[ii]; jj++) {
	    themap->types[ii][jj] = 0;
	    nocntr = 0;
	    for (kk = 1; kk <= nn; kk++) {
	    switch ( thedata[kk].markers[ii][jj] ) {
	      case 0: case 1: case 2: case 3:             break;
	      case 10: themap->types[ii][jj] = -1;         break;
	      case 12: case 13: themap->types[ii][jj] = 1; break;
	      default: nocntr = nocntr + 1;                break;
	    }
	  }
	  if (nn - nocntr < MIN_DATA)
	    themap->types[ii][jj] = -99;
      }
    }
}

/*

 Example of a translation table:

-TranslationTable
     AA    3    2
     Aa    1    1
     aa    0    0
     A-   12   12
     a-   10   10
     --   -1   -1

The array will have 5 rows corresponding to 2, 1, 0, 12 and 10, respectively.

*/
void get_TransTable(FILE *fptr,char **TransTable)
{
  int ii, indx, ch;
  for (ii = 1; ii <= 5; ii++) {
    ch = get_next_token(gname, MAXNAME, fptr);
    ch = get_next_token(gname, MAXNAME, fptr);
    if (!strcmp(gname, "2"))
      indx = 1;
    else if (!strcmp(gname, "1"))
      indx = 2;
    else if (!strcmp(gname, "0"))
      indx = 3;
    else if (!strcmp(gname, "12"))
      indx = 4;
    else if (!strcmp(gname, "10"))
      indx = 5;

    ch = get_next_token(gname, MAXNAME, fptr);
    strcpy(TransTable[indx], gname);
  }
}

/*  Get genotypes for a specific marker.  Individuals in a standard order.*/
void get_marki(FILE *fptr, markermap *themap, individual *thedata, int nn, char **TransTable, int iscase)
{
  int ch,ii, go_on, ismark, cc, mm;
  go_on = 1;
  while (go_on == 1) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if (!strcmp(gname, "-stop") || ch == EOF)
      go_on = 0;
    else {
	  ismark = which_mark(gname, themap, &cc, &mm, iscase,0);
      for (ii = 1; ii <= nn; ii++) {
	    ch = get_next_token(gname, MAXNAME, fptr);
	    if (ismark == 1)
	      thedata[ii].markers[cc][mm] = trans_mark(gname, TransTable);
      }
    }
  }
} 

/*Try to match a marker name to one of the markers in the genetic linkage map */
int which_mark(char *xtemp,markermap *themap,int *cc, int *mm, int iscase,int mapmaker)
{
  int ii, jj,tocomp;
  char *ch;
  if ( mapmaker == 1 )
    tocomp = 8;
  else 
    tocomp = (int) MAXNAME;
  if (iscase == 0)
    ch = strlwr(xtemp);
  for (ii = 1; ii <= themap->m; ii++) {
    for (jj = 1; jj <= themap->mpc[ii]; jj++)/* Mapmaker quirk.*/
      if (!strncmp(xtemp, themap->names[ themap->ttable[ii][jj] ], (size_t) tocomp) ) {
	    *cc = ii;
	    *mm = jj;
	    return (1);
      }
  }
  return (0);
}

/*Translate the marker genotype to the standard. */
int trans_mark(char *xtemp, char **TransTable)
{
  int ii, rval, rvals[6] ;
  rvals[0] = -1;
  rvals[1] = 3;
  rvals[2] = 1;
  rvals[3] = 0;
  rvals[4] = 12;
  rvals[5] = 10;
  rval = 0;
  for (ii = 1; ii <= 5; ii++)
    if (!strcmp(xtemp, TransTable[ii]))
      rval = ii;
  return (rvals[rval]);
}

/* get the phenotypes for the specified trait.  Individuals in a standard
order. */
int get_traiti(FILE *fptr,markermap *themap,individual *thedata,int nn, int traitst, int iscase,char *miss)
{
  int ch,go_on, trait, ii;
  ii = iscase;
  go_on = 1;
  trait = traitst - 1;
  while (go_on == 1) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if (!strcmp(gname, "-stop") || ch == EOF)
      go_on = 0;
    else {
      trait = trait + 1;
      strcpy(themap->tnames[trait], gname);
      for (ii = 1; ii <= nn; ii++) {
	    ch = get_next_token(gname, MAXNAME, fptr);
	    if (!strcmp(gname, miss))
	      thedata[ii].y[trait] = (FPN) MISS_VAL - (FPN)1.0 ;
	    else
	      thedata[ii].y[trait] = (FPN) atof(gname);
      }
    }
  }
  return (trait + 1);
}

/*Get the categorical trait values for the specified trait, with the Individuals in
a standard order. */
int get_otraiti(FILE *fptr,markermap *themap,individual *thedata,int nn,int otraitst,int iscase,char *miss)
{
  int ch,go_on,trait,ii;
  ii = iscase;
  go_on = 1;
  trait = otraitst-1;
  while (go_on == 1) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if (!strcmp(gname, "-stop") || ch == EOF)
      go_on = 0;
    else {
      trait = trait + 1;
      strcpy(themap->onames[trait], gname);
      for (ii = 1; ii <= nn; ii++) {
	    ch = get_next_token(gname, MAXNAME, fptr);
	    if (!strcmp(gname, miss))
	      strcpy( thedata[ii].oy[trait], "." );
	    else
	      strcpy( thedata[ii].oy[trait], gname);
      }
    }      
  }
  return (trait + 1);
}

/*
  Data is organized in the nn times repeating sequence,

  Ind.name Mark1.1, Mark1.2, Mark1.3,....MarkC.CM

where Marki.j means marker j on chromosome i.
*/
void get_imark(FILE *fptr,markermap *themap,individual *thedata,int nn, char **TransTable,int iscase)
{
  int cc, mm, ii;
  char  *chptr;
  int ch;
  for (ii = 1; ii <= nn; ii++) {
    ch = get_next_token(gname, MAXNAME, fptr);
    thedata[ii].name = cvector(0, MAXNAME);
    if (iscase == 0)
      chptr = strlwr(gname);
    strcpy( thedata[ii].name, gname);
    for (cc = 1; cc <= themap->m; cc++) {
      for (mm = 1; mm <= themap->mpc[cc]; mm++) {
	    ch = get_next_token(gname, MAXNAME, fptr);
	    thedata[ii].markers[cc][mm] = trans_mark(gname, TransTable);
      }
    }
  }
  ch = get_next_token(gname, MAXNAME, fptr);	/* get the stop token */
}

/*match an individual name from what is known*/
int which_ind(char *xtemp,individual *thedata,int nn)
{
  int ii;
  for (ii = 1; ii <= nn; ii++)
    if (!strcmp(xtemp, thedata[ii].name))
      return (ii);
  return (0);
}


/*
 This subroutine gets the trait values.  If the individuals are
unnamed, then the data should be ordered from individuals 1 to nn.
*/
int get_itrait(FILE *fptr,markermap *themap,individual *thedata,int nn, int traitst,int iscase,char *miss)
{
  int ch,go_on, ntraits, named, whichind, ii;
  char  *chptr;
  whichind = 0;
  ch = get_next_token(gname, MAXNAME, fptr);
  ntraits = atoi(gname);
  for (ii = 1; ii <= ntraits; ii++) {
    ch = get_next_token(gname, MAXNAME, fptr);
    strcpy(themap->tnames[traitst + ii - 1], gname);
  }
  ch = get_next_token(gname, MAXNAME, fptr);
  if (!strcmp(gname, "named"))
    named = 1;
  else
    named = 0;
  go_on = 1;
  while (go_on == 1) {
    if (named == 1) {
      ch = get_next_token(gname, MAXNAME, fptr);
      if (!strcmp(gname, "-stop") || ch == EOF) {
	    return (traitst + ntraits);
      }
      if (iscase == 0)
        chptr = strlwr(gname);
      
      whichind = which_ind(gname, thedata, nn);
    }
    else
      whichind = whichind + 1;
    for (ii = 1; ii <= ntraits; ii++) {
      ch = get_next_token(gname, MAXNAME, fptr);
      if (!strcmp(gname, "-stop") || ch == EOF)
	    return (traitst + ntraits);
      if (!strcmp(gname, miss))
	    thedata[whichind].y[traitst + ii - 1] = (FPN) MISS_VAL - (FPN)1.0;
      else
	    thedata[whichind].y[traitst + ii - 1] = (FPN) atof(gname);
    }
  }
  return (traitst + ntraits);
}

/* get categorical traits for the individual*/
int get_iotrait(FILE *fptr,markermap *themap,individual *thedata,int nn, int otraitst, int iscase,char *miss)
{
  int ch,go_on, ntraits, named, whichind, ii;
  char  *chptr;
  whichind = 0;
  go_on = 1;
  whichind = 0;
  ch = get_next_token(gname, MAXNAME, fptr);
  ntraits = atoi(gname);
  for (ii = 1; ii <= ntraits; ii++) {
    ch = get_next_token(gname, MAXNAME, fptr);
    strcpy(themap->onames[otraitst + ii - 1], gname);
  }
  ch = get_next_token(gname, MAXNAME, fptr);
  if (!strcmp(gname, "named"))
    named = 1;
  else
    named = 0;
  go_on = 1;
  while (go_on == 1) {
    if (named == 1) {
      ch = get_next_token(gname, MAXNAME, fptr);
      if (!strcmp(gname, "-stop") || ch == EOF) {
	return (otraitst + ntraits);
      }
      if (iscase == 0)
        chptr = strlwr(gname);
      whichind = which_ind(gname, thedata, nn);
    }
    else
      whichind = whichind + 1;
    for (ii = 1; ii <= ntraits; ii++) {
      ch = get_next_token(gname, MAXNAME, fptr);
      if (!strcmp(gname, "-stop") || ch == EOF)
	    return (otraitst + ntraits);
      if (!strcmp(gname, miss))
	    strcpy( thedata[whichind].oy[otraitst + ii - 1], ".");
      else
	    strcpy( thedata[whichind].oy[otraitst + ii - 1], gname);
    }
  }
  return (otraitst + ntraits);
}

/*Just in case a -filetype token is not present, give it another check.*/
int recheck_file_type(char *infile)
{
  int ch;
  FILE *fptr;
  fptr = fileopen(infile, "r");
  if (fptr == NULL) 
    return (-1);
  ch = get_next_token(gname, MAXNAME, fptr);
  if ( !strcmp( "data", gname ) ) {
    fileclose(infile, fptr);
    return(32);  
  } 
  else if ( !strcmp("FileID", gname) ) {
    fileclose(infile, fptr);
    return(4);  
  } 
  do {
    if ( gname[0] == '-' && gname[1] == 'n' ) {
      fileclose(infile, fptr);
      return(30);
    }
    else if ( gname[0] == '-' && gname[1] == 'C' ) {
      fileclose(infile, fptr);
      return(31);
    }
    else if ( !strcmp( "mapmaker", gname ) ) {
      fileclose(infile, fptr);
      return(32);  
    }  
  } while (  (ch = get_next_token(gname, MAXNAME, fptr)) != EOF ) ;
  return(-1);
}

/*driver to get cross.inp data format*/
individual *get_std_data(params *theparams, char *infile)
{
  FILE *fptr;
  individual *thedata;
  int ch,go_on, traitst, otraitst, ii, jj, isskip, iscase;
  char   **TransTable, *chptr;
  char  *miss;

    


  if ( theparams->themap == NULL ) {  
    ii = get_file_type(infile);
    if ( ii == 0 || ii == -1) 
      ii = recheck_file_type(infile);
    if ( ii == 4 && whichprogram != 12 )
      theparams->themap = trans_stdm(infile);  
    else
      theparams->themap = crossbogusmap(infile);
  }
  
  miss = cvector(0,MAXNAME);
  thedata = NULL;
  TransTable = cmatrix(0, 6, 0, MAXNAME);
  iscase = traitst = otraitst = go_on = 1;
  isskip = 0;
  fptr = fileopen(infile, "r");
  if (fptr == NULL) {
    free_cvector(miss,0,MAXNAME);
    free_cmatrix( TransTable , 0, 6, 0, MAXNAME);
    return (NULL);
  }
  ch = get_next_token(gname, MAXNAME, fptr);  /*  first token.  */
  if (!strcmp("#FileID",gname)   ) 
    go_on = MoveToToken(fptr, "#bycross", 1 );
    
  if ( go_on == -1 ) {
    fileclose(infile, fptr);
    free_cvector(miss,0,MAXNAME);
    free_cmatrix( TransTable , 0, 6, 0, MAXNAME);
    return (NULL);
  }
  go_on = 1;
  while (go_on == 1 && ch != EOF) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if (gname[0] == '-' && isskip == 0) {
      switch (gname[1]) {
       case 'C':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     theparams->cross = parse_cross(theparams,gname);
	     break;
       case 'S':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     theparams->nn = atoi(gname);
	     break;
       case 't':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     theparams->traits = atoi(gname);
	     if ( theparams->themap->traits != theparams->traits ) {
           theparams->themap->traits = theparams->traits;
           free_ivector(theparams->themap->knum, 1, 1);
           theparams->themap->knum = ivector(1, theparams->themap->traits);
	     } 
	     if ( theparams->themap->traits > 0 )
	       theparams->themap->tnames = cmatrix(1, theparams->themap->traits, 0, MAXNAME);
	     break;
       case 'P':
        theparams->themap->ParentalDiff = dvector(1, theparams->themap->traits);        
	    for ( ii = 1 ; ii <= theparams->themap->traits ; ii++ ) {
	        ch = get_next_token(gname, MAXNAME, fptr);
            theparams->themap->ParentalDiff[ii] = (FPN) atof( gname );
         }
         break;
       case 'o':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     theparams->themap->otraits = atoi(gname);
	     if ( theparams->themap->otraits  > 0 ) 
	       theparams->themap->onames = cmatrix(1, theparams->themap->otraits , 0, MAXNAME);
	     break;
       case 'c':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     if (gname[0] == 'n' || gname[0] == 'N') {
	       iscase = 0;
	       if (theparams->themap->names != NULL)
	         for (ii = 1; ii <= theparams->themap->ml; ii++)
	           chptr = strlwr(theparams->themap->names[ii]);
	     }
	     else if (gname[0] == 'Y' || gname[0] == 'Y')
	     iscase = 1;
	     break;
       case 'm':
	     ch = get_next_token(gname, MAXNAME, fptr);
	     strcpy(miss, gname);
	     break;
       case 'T':
	     get_TransTable(fptr, TransTable);
	     break;
       case 's':
	     switch (gname[3]) {
	       case 'a':
	         if (thedata == NULL) {
	           thedata = indvector(theparams->nn, theparams->themap, theparams->theqtls );
	           if ( theparams->themap->otraits > 0 ) 
	           for ( jj = 1 ; jj <= theparams->nn ; jj++ )
                 if (  thedata[jj].oy == NULL ) 
                   thedata[jj].oy = cmatrix(1,theparams->themap->otraits,0,MAXNAME); 
		     }
	         ch = get_next_token(gname, MAXNAME, fptr);
	         if (gname[0] == 'm')
	           get_marki(fptr, theparams->themap, thedata, theparams->nn, TransTable, iscase);
	         else if (gname[0] == 't')
	           traitst = get_traiti(fptr, theparams->themap, thedata, theparams->nn, traitst, iscase, miss);
	         else if (gname[0] == 'o') {
	           otraitst = get_otraiti(fptr, theparams->themap, thedata, theparams->nn, otraitst, iscase, miss);
             }
	         else if (gname[0] == 'i') {
	           ch = get_next_token(gname, MAXNAME, fptr);
	           if (gname[0] == 'm')
	             get_imark(fptr, theparams->themap, thedata, theparams->nn, TransTable, iscase);
	           else if (gname[0] == 't')
	             traitst = get_itrait(fptr, theparams->themap, thedata, theparams->nn, traitst, iscase, miss);
	           else if (gname[0] == 'o')
	             otraitst = get_iotrait(fptr, theparams->themap, thedata, theparams->nn, otraitst, iscase, miss);
	         }
	         break;
	       case 'i':
	         isskip = 1;
	         break;
	       case 'o':
	       default:
	         break;
	    }
	    break;
        case 'e':
        case 'q':
	      go_on = 0;
	      break;
       default:
	      break;
      }
    }
    else if (!strcmp(gname, "-unskip"))
      isskip = 0;
  }
  for (ii = 1; ii <= theparams->nn; ii++)
    thedata[ii].map = theparams->themap;
  if ( go_on != -98 )
    fileclose(infile, fptr);
  free_cmatrix(TransTable, 0, 6, 0, MAXNAME);
  free_cvector(miss,0,MAXNAME);
  return (thedata);
}

/*  get_the_datapoints gets the values of the markers and quantitative traits.
    There is not yet any
    checking for the correct number of observations.  The observations are numbered,
    and can be in any order.  The -s to start a line must go before the first individual.
    The -e ends the data.  After the -s, there should be an integer identifying which individual,
    a 1, then the markers, then the quantitative trait value.  All of these should be separated by
    whitespace, ie, tabs, newlines or spaces.  There should be a value for all the markers in proper
    order (chrom1 marker1, 2, 3, ..., chrom2 marker1, 2, 3,...).  Put a 0 for no data, but make sure
    that there is a value from 0...9 for each marker defined in the map that was input earlier.
*/
int get_the_datapoints(individual *thedata,char *inputf,int nn)
{
  FILE *infile;
  int ch, ii, ind, chrom, mrkr,  old,otraits,intraits;
  intraits = otraits = 0;
  old = 0;
  infile = fileopen(inputf, "r");
  if (infile == NULL) 
    return (-1);
  do {
    ch = get_next_line(gbuffer, MAXLINE, infile);
    if (gbuffer[0] == '-') {
      if ( gbuffer[1] == 's') {	/* start to read in data, just after a -s */
	    do {
	      ch = get_next_token(gname, MAXNAME, infile);
	      if ( gname[0] == '-')
	        break;
	      else
	        ind = atoi(gname);	/* this is the individual number. */
	      ch = get_next_token(gname, MAXNAME, infile);
	      thedata[ind].bc = atoi(gname);	  
		  for (chrom = 1; chrom <= thedata[ind].map->m; chrom++) {
		    for (mrkr = 1; mrkr <= thedata[ind].map->mpc[chrom]; mrkr++) {
		      ch = get_next_token(gname, MAXNAME, infile);
		      thedata[ind].markers[chrom][mrkr] = atoi(gname);
		      if (old == 0 && thedata[ind].markers[chrom][mrkr] == 2)
			    thedata[ind].markers[chrom][mrkr] = 3;
		    }
		  }
		  for (ii = 1; ii <= thedata[ind].t; ii++) {
		    ch = get_next_token(gname, MAXNAME, infile);
		    if (gname[0] != '.')
		      thedata[ind].y[ii] = (FPN) atof(gname);
		    else
		      thedata[ind].y[ii] = (FPN) MISS_VAL- (FPN)1.0;
	      }
	      if ( (thedata+1)->map->otraits > 0 ) 
		    for (ii = 1; ii <= thedata[1].map->otraits; ii++) {
		      ch = get_next_token(gname, MAXNAME, infile);
		      strcpy( thedata[ind].oy[ii] , gname);
	        }
		} while (ch != EOF);
      }
      else if (gbuffer[1] == 'o' && gbuffer[2] == 't' ) {
	    intraits = 2;
        get_field(2, gname, gbuffer);
		  thedata[1].map->otraits = atoi(gname);
	    if ( thedata[1].map->otraits > 0 )
	      for ( ii = 1 ; ii <= nn ; ii++ )
			 thedata[ii].oy = cmatrix(1,thedata[1].map->otraits,0,MAXNAME);
      }
      else if (gbuffer[1] == 't') 
	   intraits = 1;
      else if (gbuffer[1] == 'P') {
        thedata[1].map->ParentalDiff = dvector(1, thedata[1].map->traits);        
	    for ( ii = 1 ; ii <= thedata[1].map->traits ; ii++ ) {
	        ch = get_next_token(gname, MAXNAME, infile);
            thedata[1].map->ParentalDiff[ii] = (FPN) atof( gname );
         }
      }
      else if (gbuffer[1] == 'N') {
        if ( intraits == 1 ) {
          if ( thedata[1].map->tnames == NULL  )
            thedata[1].map->tnames = cmatrix(1,  thedata[1].map->traits, 0, MAXNAME );
	      for ( ii = 1 ; ii <= thedata[1].map->traits ; ii++ ) {
	        ch = get_next_token(gname, MAXNAME, infile);
	        ch = get_next_token(gname, MAXNAME, infile);
            strcpy( thedata[1].map->tnames[ii], gname );
          }
        }
        else if ( intraits == 2 ) {
          if ( thedata[1].map->onames == NULL  ) {
            thedata[1].map->onames = cmatrix(1,  thedata[1].map->otraits, 0, MAXNAME );
            thedata[1].map->otypes = ivector(1, thedata[1].map->otraits);
          }
	      for ( ii = 1 ; ii <= thedata[1].map->otraits ; ii++ ) {
	        ch = get_next_token(gname, MAXNAME, infile);
	        ch = get_next_token(gname, MAXNAME, infile);
            strcpy( thedata[1].map->onames[ii], gname );
          }
        }
      }
      else if (gbuffer[1] == 'o')
	     old = 1;
/*	  else if ( gbuffer[1] == 'q' )
	     ch = EOF; */
	  else if ( gbuffer[1] == 'I' ) {
	      
	      do {    
	        ch = get_next_token(gname, MAXNAME, infile);
	        if ( gname[0] != '-' ) {
	          ii = atoi(gname);
	          ch = get_next_token(gname, MAXNAME, infile);
	          thedata[ii].name = cvector(0, MAXNAME);
              strcpy( thedata[ii].name, gname );
            }
            else
              ch = EOF;
          } while ( ch != EOF ); 
	  }
    }
  } while (ch != EOF);
  fileclose(inputf, infile);
  if (thedata[1].map->types != NULL)
    determine_markers(thedata, nn);
  return (0);
}



/* 
 get_the_nn merely reads in nn from the data file.
 It returns -1 for an error, either because there is no
 data file or because n and p were not found.  This should
 be used as the first pass through the input file.
*/
int get_the_nn(int *n, int *t, char *inputf)
{
  FILE *infile;
  int ch, ii, jj, nyes;
  nyes = 0;
  infile = fileopen(inputf, "r");
  if (infile == NULL)
    return (-1);
  jj = 0;
  do {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii]  = (char) ch;
    if (gbuffer[0]  == '-')
      switch (gbuffer[1] ) {
       case 'n':	/* sample size */
	     get_field(2, gname, gbuffer);
	     *n = atoi(gname);
	     nyes = nyes + 1;
	     break;
       case 't':	/* number of traits */
	     get_field(2, gname, gbuffer);
	     *t = atoi(gname);
	     nyes = nyes + 1;
	   break;
       case 'q':
	     fileclose(inputf, infile);
	     printf("\nYou did not specify n and t...\n");
	     return (-1);
	     break;
       default:
	     break;
      }
    if (nyes == 2) {
      fileclose(inputf, infile);
      return (0);
    }
  } while (ch != EOF);
  if (ch == EOF)
    printf("\nTime to close the input file after an EOF...\n");
  fileclose(inputf, infile);
  return (-1);
}

/*  
  This reads in MAPMAKER/QTL data files.  

  Two cases are supported.  

  1.  You have already produced a genetic linkage map in MAPMAKER.
      You should have translated that output into the Rmap.out format.
      

  2.  You don't have a genetic linkage map.  In this case, a bogus map
      with one chromosome is created.  The markers are placed on this  
      map in the order read in, in intervals of 10 cM.

  It returns a pointer to the individual data structure.  Note that if
  themap was null upon entry, the only way to access the linkage map will
  be through the indvector->map  pointer.


*/

individual *get_mm_data(params *theparams,char *infile )
{
/* 
Define a pointer to a local map.  If theparams->themap is NULL, then create a bogus map,
put all marker names into it and then assign it to theaparams->themap at the end.  If theparams->themap
exists, then assign lmap to it and use the existing map.
*/
   FILE *fptr;

  char   **TransTable;
  individual *thedata;
  FPN **trait_v;
  int go_on,ii,i,jj,ismark,cc,mm,symbols;
  int ch,nind, isequation,col,iscase,ml,trans[6],bogusmarkermap;
  char    *chptr;
  char   *miss;
  TransTable = cmatrix(0, 6, 0, 4);
  miss = cvector(0,MAXNAME);
  iscase = 0;
  go_on = 1;
  bogusmarkermap = 0;
  fptr = fileopen(infile, "r");
  if (fptr == NULL) {
    free_cvector(miss,0,MAXNAME);
    free_cmatrix( TransTable , 0, 6, 0, MAXNAME);
    return (NULL);  /* Return nothing if there is no file.*/
  }
/* 
  TransTable...
      column      0   1   2   3
-------------------------------
      iscase      0   1   0   1         1 if markers are case sensitive
      symbols     0   0   1   1         1 if different symbols are used

Default Symbols
row Int Gt        MM
             ------------
             BC   F2   RI
-------------------------
0    3  AA   A    A    A
1    1  Aa   H    H
2    0  aa   A    B    B
3   12  A-        D
4   10  a-        C
5   <0  --   -    -    -

*/
  for ( ii = 0 ; ii <= 6 ; ii++ )
    TransTable[ii][3] = TransTable[ii][2] =  TransTable[ii][1] = TransTable[ii][0] = '-';
  TransTable[0][0] = 'A';
  TransTable[1][0] = 'H';
  TransTable[2][0] = 'B';
  TransTable[3][0] = 'D';
  TransTable[4][0] = 'C';
  trans[0] = 3;
  trans[1] = 1;
  trans[2] = 0;
  trans[3] = 12;
  trans[4] = 10;
  trans[5] = -1;

  go_on = 1;
  while ( go_on == 1 ) {
    ch = get_next_token(gname, MAXNAME, fptr);  /* should be `type'  */
    if ( !strcmp(gname, "type") )
      go_on = 0;
    if ( ch == EOF ) {
      fileclose(infile, fptr);
      free_cvector(miss,0,MAXNAME);
      free_cmatrix( TransTable , 0, 6, 0, MAXNAME);
      return (NULL);
    }      
  }
  go_on = 1;
  
  
  ch = get_next_token(gname, MAXNAME, fptr);  /* should be `f2', `f3', or `ri' */
  ch = get_next_token(miss, MAXNAME, fptr);   /* should be `intercross', `backcross', `self' or `sib' */
  chptr = strlwr(gname);
  chptr = strlwr(miss);

  if (  !strcmp(gname, "f2") ) {
    if ( !strcmp(miss, "intercross") ) {
      strcpy(gbuffer, "RF2\0");
    }
    else if ( !strcmp(miss, "backcross") ) {
      TransTable[2][0] = 'A';
      strcpy(gbuffer, "B1\0");
    }
    else if ( !strcmp(miss, "self") ) {
      strcpy(gbuffer, "SF2\0");
    }
  }
  else if ( !strcmp(gname, "f3") ) {
      strcpy(gbuffer, "SF3\0");
  }
  else if ( !strcmp(gname, "ri") ) {
    if ( !strcmp(miss, "self") ) {
      strcpy(gbuffer, "RI1\0");
    }
    else if ( !strcmp(miss, "sib") ) {
      strcpy(gbuffer, "RI2\0");
    }
  }
  else {
    strcpy(gbuffer, "1\0");
    go_on = -1;
  }
  ii = parse_cross(theparams, gbuffer);
  
    
    

  ch = get_next_token(gname, MAXNAME, fptr);  /* should be `nn', sample size */
  theparams->nn = atoi(gname);
  ch = get_next_token(gname, MAXNAME, fptr);  /* should be `ml', number of markers  */
  ml = atoi(gname);
  ch = get_next_token(gname, MAXNAME, fptr);  /* should be `traits', number of traits  */
  theparams->traits  = atoi(gname);  
  for ( ii = 0 ; ii <= 6 ; ii++ ) /* all columns same initially */
    TransTable[ii][3] = TransTable[ii][2] = TransTable[ii][1] = TransTable[ii][0];
  while ( ch != EOF && ch != '\n' ) {  /* check for `case' and `symbols' */
	 ch = get_next_token(gname, MAXNAME, fptr);
    if ( !strcmp(gname, "case") )
      iscase = 1;
    else if ( !strcmp(gname, "symbols") ) {
      while ( ch != EOF && ch != '\n' ) {
        ch = get_next_token(gname, MAXNAME, fptr);
        if ( gname[0] == '*' )
          ch = '\n';
        else
          for ( ii = 0 ; ii <= 6 ; ii++ ) 
            if ( gname[2] == TransTable[ii][0] )
              TransTable[ii][3] = gname[0];
      }
      symbols = 1;
    }
    else if ( gname[0] == '*' )
      ch = '\n';
  }
  for ( ii = 0 ; ii <= 6 ; ii++ )
        TransTable[ii][2] = (char) toupper(TransTable[ii][3]);
  col = 1;
  if ( symbols == 1 ) {
    col = 3;
    if ( iscase == 0 ) 
      col = 2;
  }



  if ( theparams->themap != NULL ) {
    for ( ii = 1 ; ii <= theparams->themap->ml ; ii++ )
      chptr = strlwr(theparams->themap->names[ii]);  /* Change case of markers to lower...*/
	if ( theparams->themap->traits !=  theparams->traits ) {
      theparams->themap->traits =  theparams->traits;
      free_ivector(theparams->themap->knum, 1, 1);
      theparams->themap->knum = ivector(1, theparams->themap->traits);
	} 
    theparams->themap->tnames = cmatrix(1,theparams->traits,0,MAXNAME);
  }
  else  { /* if no map, put all on one chromosome */
    theparams->themap  = bogus_markermap(1, ml, ml,theparams->traits,(FPN) 0.1);
    bogusmarkermap = 1;
  }

  thedata = indvector(theparams->nn, theparams->themap, theparams->theqtls);

  for ( jj = 1 ; jj <= ml ; jj++ ) {
    while (  gname[0] != '*' )  /* go to first locus name */
      ch = get_next_token(gname, MAXNAME, fptr);
    for ( ii = 0 ; ii < MAXNAME ; ii++ )  /* move it over one space */
      gname[ii] = gname[ii+1];
    if ( bogusmarkermap == 0 )
      ismark = which_mark(gname, theparams->themap, &cc, &mm, 0,1);
    else {
	  strcpy( theparams->themap->names[jj] , gname );
	  theparams->themap->ttable[1][jj] = jj;
      ismark = cc = 1;
      mm = jj;
	 }
	 if ( ismark == 1 )
	   nind = get_mm_markers(fptr,thedata, cc,mm,theparams->nn,TransTable,col,trans);
  }

  trait_v = dmatrix(1,theparams->traits,1,theparams->nn);
  ch = 'g';
  for ( jj = 1 ; jj <= theparams->traits ; jj++ ) {
    
    while ( ch != EOF &&  gname[0] != '*' )  /* go to first trait name */
      ch = get_next_token(gname, MAXNAME, fptr);
    for ( ii = 0 ; ii < MAXNAME ; ii++ )  /* move it over one space */
       gname[ii] = gname[ii+1];
	strcpy( theparams->themap->tnames[jj] , gname );
    isequation = 0;
    for ( i = 0 ; *(gname+i) != '\0' ; i++ )
      if ( gname[i] == '=' )
        isequation = 1;

    
    if ( isequation == 0 )
      nind = get_mm_traits(fptr, trait_v, jj, theparams->traits, theparams->nn);
    else
      nind = do_trait_equation(fptr, trait_v, jj, theparams->traits, theparams->nn);

  }

  for ( ii = 1 ; ii <= theparams->nn ; ii++ ) 
    for ( jj = 1 ; jj <= theparams->traits ; jj++ )
      thedata[ii].y[jj] = trait_v[jj][ii];

  fileclose(infile, fptr);


  free_dmatrix(trait_v,1,theparams->traits,1,theparams->nn);
  free_cmatrix(TransTable,0, 6, 0, 4);

  free_cvector(miss,0,MAXNAME);
  return(thedata);
}  

/*
  This will be a driver for the eqn.c and eqn.h files
  from MAPMAKER/QTL.  
  
  This has not been implemented, and probably won't be.  
*/
int  do_trait_equation(FILE *fptr, FPN **trait_v, int wtrait, int traits, int nn)
{
  int ii,  ch;
  ii = traits;
  ch = get_next_line(gbuffer,MAXLINE,fptr);  /* put the RHS of the equation in gbuffer */

  printf("\n%s\n",gbuffer);

  for ( ii = 1; ii <= nn ; ii++ )
     trait_v[wtrait][ii] = (FPN) MISS_VAL - (FPN) 1.0;


  return(ii-1);

}


/*  
  This assumes that you have read in the trait name.  It also
  assumes that this trait is not a function of other traits.
  It reads in trait values one at a time and puts them into
  the trait_v matrix (row is wtrait).  It returns the number
  of trait values successfully read.  If an end of file is
  encountered before all traits are read in, then nn+1 is
  returned. 

  If a # symbol is encountered, then all tokens will be read
  until a '\n' symbol is encountered, then the subroutine will 
  continue to read in traits.
*/
int get_mm_traits(FILE *fptr, FPN **trait_v, int wtrait, int traits, int nn)
{
  int ii,incom,ch;
  ii = traits;
  ii = 1;
  incom = 0;
  while ( ii <= nn ) {
    ch = get_next_token(gname,MAXNAME,fptr);
    if ( ch == EOF )
      ii = nn+2;
    else if ( gname[0] == '#' ) 
      incom = 1;
    else if ( incom == 1 && ch == '\n' )
      incom = 0;
    else if ( gname[0] == '-' && gname[1] == '\0' ) {
      trait_v[wtrait][ii] = (FPN) MISS_VAL- (FPN) 1.0;
      ii = ii + 1;
    }
    else {
      trait_v[wtrait][ii] = (FPN) atof(gname);
      ii = ii + 1;
    }
      
  }

  return(ii-1);
}



/*  
   This assumes that you have just read the name of the locus and are
   ready to read in the nn marker values.  Whenever a # sign is encountered,
   the subroutine reads until a newline character or an end of file.
*/

int get_mm_markers(FILE *fptr,individual *thedata,int cc,int mm,int nn, char **TransTable,int col,int *trans)
{
  int ii,jj,ch;
  for ( ii = 1 ; ii <= nn ; ii++ ) {
    ch = get_next_mark(fptr);
    if ( ch == EOF )
      return(nn+1);
    jj = 0 ;
    while ( jj < 6 ) 
      if ( ch == TransTable[jj][col] ) {
        thedata[ii].markers[cc][mm] = trans[jj];
        jj = 7;
      }
      else
        jj = jj+1;
  }
  return(ii-1);
}

/*  
   This just gets the next character on the file stream.  
   Valid characters are letters, numbers, '+' and '-'.
   If a # symbol is encountered, then characters are read
   until a newline or an end of file.  
*/
int get_next_mark(FILE *fptr)
{
  int ch, go_on,test;
  go_on = 1;
  while ( go_on == 1 ) {
    ch = fgetc(fptr);
    if (  ch == '#' ) /* if a #, then read until newline or end of file */
      while ( ch != EOF || ch != '\n' )
        ch = fgetc(fptr);
    else if ( ch == '"' )  /* if a ", then read until next " or EOF */
      while ( ch != EOF || ch != '"' )
        ch = fgetc(fptr);
    else if ( (test=isalnum(ch) != 0) || ch == '+' || ch == '-' ) /* if alphanumeric, + or - then return it */
      go_on = 0;
    else if (ch == EOF ) /*if end of file, then return end of file */
      go_on = 0;
  }
  return(ch);
}



/*
 
 Read from a PLABQTL file and convert it into themap and individual.
 
*/
individual *convert_plabqtl(params *theparams) {
  int ch,nmark,nmin,nindent,recf,fgen,ralph,i,j,k,l,m;
  FILE *infptr;
  FPN **traits,*recfreq;
  int **markers,*mpc,*transtable,environs,kindent;
  char **traitnames,**markernames,**chromnames,*chptr;
  markermap *themap;
  individual *thedata;
  
  infptr = fileopen(theparams->iinfile, "r");
  ch = get_next_line(gbuffer, (int) MAXLINE, infptr); /* Skip comment line */
  
  chptr = strstr(gbuffer,"-environment");
  if ( chptr != NULL ) {
    get_field(2,gname,chptr);
    environs = atoi(gname);
  }
  else
    environs = 1;

  ch = get_next_token(gname, MAXNAME, infptr);
  theparams->nn = atoi(gname);
  
  ch = get_next_token(gname, MAXNAME, infptr);
  nmark  = atoi(gname);
  ch = get_next_token(gname, MAXNAME, infptr);
  theparams->traits = atoi(gname);

  ch = get_next_token(gname, MAXNAME, infptr);
  theparams->chrom = atoi(gname);
  ch = get_next_token(gname, MAXNAME, infptr);
  nmin  = atoi(gname);
  nmin = nmin + 35;
  ch = get_next_token(gname, MAXNAME, infptr);
  nindent  = atoi(gname);
  ch = get_next_token(gname, MAXNAME, infptr);
  recf  = atoi(gname);
  ch = get_next_token(gname, MAXNAME, infptr);
  fgen  = atoi(gname);
  ch = get_next_token(gname, MAXNAME, infptr);
  ralph  = atoi(gname);
  if ( fgen == 1 ) 
    strcat(gname,"Ri0");
  else if ( fgen < 34 )  
    sprintf(gname,"SF%d",fgen);
  else 
    strcat(gname, "Ri1");
  theparams->cross = parse_cross(theparams,gname);
  
  
  traits = dmatrix(1,theparams->nn,1,theparams->traits);
  traitnames = cmatrix(1,(environs+1)*theparams->traits,0,MAXNAME);
  markers = imatrix(1,theparams->nn,1,nmark);
  markernames = cmatrix(1,nmark,0,MAXNAME);
  mpc = ivector(1,theparams->chrom);
  recfreq = dvector(1,nmark);
  transtable = ivector(1,nmark);
  chromnames = cmatrix(1,theparams->chrom,0,MAXNAME);
  
  for ( i=1; i<=theparams->traits; i++ )
    ch = get_next_token(gname, MAXNAME, infptr);/*   Skip heritabilities   */

  
  if ( nmin == 35 ) {  /*  Get the trait and marker names from matrix formatted file. */  
    for ( i=1; i<=theparams->traits; i++ )
      ch = get_next_token( traitnames[i], MAXNAME, infptr);
    for ( i=1; i<=nmark; i++ )
      ch = get_next_token( markernames[i], MAXNAME, infptr);

    for ( i=1; i<=theparams->nn; i++ ) { 
      for ( j=1; j<=nindent; j++ )
        ch = get_next_token(gname, MAXNAME, infptr);
      for ( j=1; j<=theparams->traits; j++ ) {
        ch = get_next_token(gname, MAXNAME, infptr);
        traits[i][j] = (FPN) atof(gname);
      }
      for ( j=1; j<=nmark; j++ ) {
        ch =  get_next_mark(infptr);
        markers[i][j] = convert_mark(ch);
      }        
    }
  }
  else {  /* Get vector form data*/
    for ( j=1; j<=nmark;j++ ) {
      ch = get_next_token(markernames[j], MAXNAME, infptr);
      for ( i=1; i<=theparams->nn; i++ ) {
        ch =  get_next_mark(infptr);
        markers[i][j] = convert_mark(ch);
      }
    }
    for ( j=1 ; j<=theparams->traits; j++ ) {
      ch = get_next_token(traitnames[j], MAXNAME, infptr);
      for ( i=1; i<=theparams->nn; i++ ) {  
        ch = get_next_token(gname, MAXNAME, infptr);
        traits[i][j] = (FPN) atof(gname);
      }    
    }  
  }
  k=1;
  for ( i=1; i<=theparams->chrom; i++ ) {
    ch = get_next_token(chromnames[i], MAXNAME, infptr);
    ch = get_next_token(gname, MAXNAME, infptr);
    mpc[i] = atoi(gname);
    for ( j=1 ; j< mpc[i] ; j++ ) {
      ch = get_next_token(gname, MAXNAME, infptr);
      m = is_pinteger(gname);
      if ( m <0 )  
        m = find_token(nmark, 7, gname,markernames);
      transtable[k] = m;
      ch = get_next_token(gname, MAXNAME, infptr);
      recfreq[k] = (FPN) atof(gname);
      if ( recf == 1 && recfreq[k] ) 
        recfreq[k] = mapfunc(recfreq[k], 2);
      k = k+1;
    }
    ch = get_next_token(gname, MAXNAME, infptr);
      m = is_pinteger(gname);
      if ( m <0 )  
        m = find_token(nmark, 7, gname,markernames);

      transtable[k] = m;
    recfreq[k] = (FPN) 0.0;
    k=k+1;
  
  }
  if ( recf == 2 )
    for ( i=nmark; i>1; i-- )
      if ( recfreq[i] > 0.0 && recfreq[i-1] > 0.0 )
        recfreq[i] = recfreq[i] - recfreq[i-1] ;
        
  for ( i=nmark; i>0; i-- )      
    if ( recfreq[i] > 0.0 )
      recfreq[i] = mapfunc(recfreq[i], -2);
      
      
  if ( environs > 1 ) 
    for ( i=1; i<=environs; i++ )
      for ( j=1; j<=theparams->traits; j++ ) 
        sprintf(traitnames[i*theparams->traits+j], "Tr%dEnv%d", j ,i);

  themap = allocate_markermap( );  
  themap->traits = theparams->traits;
  themap->knum = ivector(1,  themap->traits);
 
  themap->m = theparams->chrom;         /* number of chromosomes */
  for ( i=1; i<=themap->m; i++ )
    if ( themap->maxl < mpc[i] )
      themap->maxl = mpc[i];
  themap->mpc = mpc;   /* = ivector(1,m); number of markers for each chromosome.  = l forall iff sigl <= 0 */
  themap->mrf = dmatrix(1,themap->m,0,themap->maxl+1);   /* = dvector(1,m,0,max(mpc)+1); If sigl <= 0, then l is max(mpc), and mpc isn't used.
                           pointer to a matrix of recombination frequencies between markers i and i+1 */ 
  k = 0;
  themap->ttable = imatrix(1,themap->m,1,themap->maxl);   /* Table to indicate where in names the marker name is  = imatrix(1,m,1,maxl) */
  for ( i=1; i<=themap->m; i++ )
    for ( j=1; j<=mpc[i] ; j++ ) {
      k = k+1;
      themap->mrf[i][j] = recfreq[ k ];
      themap->ttable[i][j] = transtable[k];
    }
  themap->ml = nmark ;       /* total number of markers, = m*l iff sigl <= 0 */
  themap->traits = theparams->traits*(environs+1);      /* Number of traits to simulate */
  themap->tnames = traitnames;   /* Names of the traits  = cmatrix(1,traits,0,MAXNAME) */

  themap->names = markernames;    /* Names of the markers  = cmatrix(1,ml,0,MAXNAME)*/
  themap->cnames = chromnames;   /* Names of the chromosomes  = cmatrix(1,m,0,MAXNAME)*/
   /* themap->types = NULL;    Indicate the type of marker:  = imatrix(1,m,1,maxl)*/
  themap->sigl = (FPN) 0.1;
  do_map_stats(themap);
  free_dvector(recfreq,1,nmark);
  
  thedata = indvector(theparams->nn,themap,NULL);
  for ( i=1; i<=theparams->nn; i++ )  {
    for ( j=1; j<=theparams->traits; j++ )
      thedata[i].y[j] = traits[i][j];
    l = 1;
    for ( j=1; j<=themap->m; j++ )
      for (k=1; k<=themap->mpc[j]; k++ ) {
        thedata[i].markers[j][k] = markers[i][transtable[l] ];
        l = l+1;
      }
  }
  if ( environs > 1 ) {
    ch = get_next_line(gbuffer, (int) MAXLINE, infptr); /* Skip comment line */
    get_field(2,gname,gbuffer);
    kindent = atoi(gname);
    for ( i=1; i<=environs; i++ ) 
      for ( j=1; j<=theparams->nn; j++ ) {
        ch = get_next_line(gbuffer, (int) MAXLINE, infptr);
        for ( k=1; k<=theparams->traits; k++ ) {
          get_field(kindent+k,gname,gbuffer);
          thedata[j].y[i*theparams->traits+k] = (FPN) atof(gname);
        }
      }
  }  
  free_ivector(transtable,1,nmark);
  free_dmatrix(traits,1,theparams->nn,1,theparams->traits);
  free_imatrix(markers,1,theparams->nn,1,nmark);

      
  fileclose(theparams->iinfile, infptr);
  theparams->themap = themap;
  
  return(thedata);
}

/*

  Translate from the PLABQTL standard to ours.
*/
int convert_mark(int ch) {
  int mark;
	    switch (ch) {
	      case '0':  case 'A': case 'a':  mark  = 3 ;     break; 
	      case '1':  case 'B': case 'b':  mark  = 0 ;     break; 
	      case '2':  case 'H': case 'h':  mark  = 1 ;     break; 
	      case '3':  case 'C': case 'c':  mark  = 10 ;     break; 
	      case '4':  case 'D': case 'd':  mark  = 12  ;    break; 
	      default:  mark  = -4  ;              break;
	    }
  return(mark);
}

/*
 * Create a vector of individuals, allocating space, etc.
 */
individual *indvector(int n,markermap *themap,aqtl *qtlptr)
{
  individual *thedata;
  int k, kk, maxl, ii, jj, jjj;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  thedata = (individual *) malloc((size_t) n * sizeof(individual));
#else
  thedata = (individual *) malloc((unsigned) n * sizeof(individual));
#endif
    if ( debugging > 2 ) {
/*    char gname[MAXNAME];*/
        sprintf(gwarn,"In indvector(), allocated %d individuals  at %p\n",n,thedata);
        MemoryAccount(gwarn);
    }
  for (ii = 0; ii < n; ii++) {
    thedata[ii].print_flag = 'y';
    thedata[ii].name = NULL;
    thedata[ii].map = themap;
    thedata[ii].qtls = qtlptr;
    thedata[ii].t = themap->traits;
    if ( themap->traits > 0 )
      thedata[ii].y = dvector(1, themap->traits);
    else
      thedata[ii].y = NULL;
	thedata[ii].g = NULL;
	thedata[ii].oyt = NULL;
	thedata[ii].oy = NULL;
    thedata[ii].vqtls = NULL;
	thedata[ii].bc = 1;
  }
  maxl = 0;
  for (ii = 1; ii <= themap->m; ii++)
    if (maxl < themap->mpc[ii])
	  maxl =  themap->mpc[ii];

 /* allocate space for vqtls and markers, and set all
    elements = -1 */
  k = 0;
  if (themap->knum != NULL && qtlptr != NULL) {
    for (ii = 1; ii <= themap->traits; ii++)
      if (k < themap->knum[ii])
	    k = themap->knum[ii];
  }
  for (ii = 0; ii < n; ii++) {
    if (k > 0) {
      thedata[ii].vqtls = imatrix(1, themap->traits, 1, k);
      for (kk = 1; kk <= themap->traits; kk++)
	    for (jj = 1; jj <= k; jj++)
	      thedata[ii].vqtls[kk][jj] = -1;
    }

    thedata[ii].markers = imatrix(1, thedata[ii].map->m, 0, maxl);
    for (jj = 1; jj <= themap->m ; jj++)
      for (jjj = 0; jjj <= maxl; jjj++)
	    thedata[ii].markers[jj][jjj] = -2;
  }

  return(thedata - 1);
}

/*
Copy the individual from one structure to another.

*/
void cp_individ(individual *ind1,individual  *ind2)
{
  int k, kk, jj, jjj, ii;
  ind2->map = ind1->map;
  ind2->qtls = ind1->qtls;
  ind2->t = ind1->t;
  ind2->bc = ind1->bc;
  ind2->print_flag = ind1->print_flag;
  if ( ind1->name != NULL && ind2->name != NULL )
    strcpy(ind2->name, ind1->name);
  
  if ( ind1->oy != NULL && ind2->oy != NULL) 
    for ( ii=1; ii<= ind1->map->otraits ; ii++ ) 
      strcpy( ind2->oy[ii], ind1->oy[ii] );

  if ( ind1->oyt != NULL && ind2->oyt != NULL) 
    for ( ii=1; ii<= ind1->map->otraits ; ii++ ) 
      ind2->oyt[ii] = ind1->oyt[ii];    
      
  k = 0;
  if (ind1->map->knum != NULL) {
    for (ii = 1; ii <= ind1->map->traits; ii++)
      if (k < ind1->map->knum[ii])
            k = ind1->map->knum[ii];
  }

  for (kk = 1; kk <= ind1->map->traits; kk++) {
    ind2->y[kk] = ind1->y[kk];
    if ( ind1->g != NULL  &&  ind2->g != NULL )
      ind2->g[kk] = ind1->g[kk];
    if (k > 0 && ind1->vqtls != NULL)
      for (jj = 1; jj <= k; jj++)
            ind2->vqtls[kk][jj] = ind1->vqtls[kk][jj];
  }

  for (jj = 1; jj <= ind1->map->m; jj++)
    for (jjj = 1; jjj <= ind1->map->maxl; jjj++)
      ind2->markers[jj][jjj] = ind1->markers[jj][jjj];

}

/*  This assumes that the vector ind1 goes from 1 to nn+1 and that the data are in 1, nn.  
    It moves individuals with missing phenotypes to the end of the
    vector ind1.    It uses position nn+1 as a temporary swap space.  */
int MoveMissPhenotypes(int nn,int wt, individual *thedata) {
  int last,i,k;
  k=0;
  /*make sure that the last node is not missing. */
  for ( last = nn; last > 0 && thedata[last].y[wt]  <= (FPN) MISS_VAL ; last-- ) k+=1; 
  for ( i=1; i<=last; i++ ) 
    if ( thedata[i].y[wt]  <= (FPN) MISS_VAL  )  {
      cp_individ(&thedata[i], &thedata[nn+1]);
      cp_individ(&thedata[last], &thedata[i]);
      cp_individ(&thedata[nn+1], &thedata[last]);
      i = i-1;  /* want to recheck the new individual */
      last = last-1;  
    }
  return(last);

}

/*
 * free up the space allocated to the individuals in the sample
 */
void free_indvector(individual *thedata,int n)
{
  int k, ii, maxl = 0;
  int otraits, traits,chroms;
  otraits = thedata[1].map->otraits;
  traits =  thedata[1].map->traits;
  chroms =  thedata[1].map->m;
  k = 0;
  if ( thedata[1].map->knum != NULL )
    for (ii = 1; ii <= thedata[1].map->traits; ii++)
      if (k < thedata[1].map->knum[ii])
	    k = thedata[1].map->knum[ii];
  
    for (ii = 1; ii <= thedata[1].map->m; ii++)
      if (maxl < thedata[1].map->mpc[ii])
	    maxl = thedata[1].map->mpc[ii];


  for (ii = 1; ii <= n; ii++) {
    if ( thedata[ii].name != NULL)
      free_cvector( thedata[ii].name, 0, MAXNAME);
	 if (thedata[ii].y != NULL )
		free_dvector(thedata[ii].y, 1, traits);
    if (thedata[ii].g != NULL)
		free_dvector(thedata[ii].g, 1, traits);
    if (k > 0 && thedata[ii].vqtls != NULL)
		free_imatrix(thedata[ii].vqtls, 1, traits, 1, k);
	 if (thedata[ii].markers != NULL )
		free_imatrix(thedata[ii].markers, 1, chroms, 0, maxl);
	 if (thedata[ii].oy != NULL )
		free_cmatrix(thedata[ii].oy, 1, otraits, 0, MAXNAME);
	 if (thedata[ii].oyt != NULL )
		free_ivector(thedata[ii].oyt, 1, otraits);
  }
    if ( debugging > 2 ) {
/*    char gname[MAXNAME];*/
        sprintf(gwarn,"In free_indvector(), deallocated %d individuals  at %p\n",n,thedata+1);
        MemoryAccount(gwarn);
    }
  free((char *) (thedata + 1));
}

/*
 * print out the information on all the individuals in the sample. This will
 * include the values of the markers, the QTLs and the phenotypes.  Output to
 * a file or stdout.
 *   Markers should be untranslated.  See trans_data and untrans_data
 */
void print_individuals(params *theparams,individual *indptr,int n, char *outfile)
{
  int ii, jj, kk,  indcntr,gt;
  FILE *outf;
  indcntr = 0;
  for (ii = 1; ii <= n; ii++)
    if ( indptr[ii].print_flag == 'y')
      indcntr = indcntr + 1;
  if (*outfile == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  if (outf != NULL) {
    fprintf(outf, "\n-n       %4d   is the sample size", indcntr);
    fprintf(outf, "\n-p       %4d   is one more than the number of markers", theparams->themap->ml + 1);
    fprintf(outf, "\n-cross   %s   is the type of cross", theparams->thecross);
    fprintf(outf, "\n-traits  %4d   is the number of traits", theparams->themap->traits);
    if ( theparams->themap->tnames != NULL ) {
      fprintf(outf,"\n-Names of the traits...");
      for ( ii = 1 ; ii <= theparams->themap->traits ; ii++ )
        fprintf(outf, "\n %2d %s",ii,theparams->themap->tnames[ii] );
    }
    if ( theparams->themap->ParentalDiff != NULL ) {
      fprintf(outf,"\n-ParentalDifferences \n");
      for ( ii = 1 ; ii <= theparams->themap->traits ; ii++ )
        fprintf(outf, " %f ", theparams->themap->ParentalDiff[ii]);
    
    
    }
    fprintf(outf, "\n-otraits %4d   is the number of other traits", theparams->themap->otraits);
    if ( theparams->themap->onames != NULL ) {
      fprintf(outf,"\n-Names of the other traits...");
      for ( ii = 1 ; ii <= theparams->themap->otraits ; ii++ )
        fprintf(outf, "\n %2d %s",ii,theparams->themap->onames[ii]);
    }
    fprintf(outf, "\n#From here, the first number is the individual\n");
    fprintf(outf, "# the second is a 1 or 2 (for BC1, BC2 in Design III),\n");
    fprintf(outf, "# and then come the %3d marker values and finally the trait values.\n-s", theparams->themap->ml);
    indcntr = 0;
    for (ii = 1; ii <= n; ii++)
      if ( indptr[ii].print_flag == 'y') {
	indcntr = indcntr + 1;
	fprintf(outf, "\n%3d %2d ", indcntr,indptr[ii].bc );
	for (jj = 1; jj <= theparams->themap->m; jj++) {
	  fprintf(outf, "\n    ");

	  for (kk = 1; kk <= theparams->themap->mpc[jj]; kk++) 
	    if ( theparams->cross == 6 ) {
	      if ( theparams->tcross == 1 ) 
	        gt = first_gam( indptr[ii].markers[jj][kk] );
	      else  if ( theparams->tcross == 2 ) 
	        gt = second_gam( indptr[ii].markers[jj][kk] ) - 3;
	      fprintf(outf, "%3d ", gt);	    
	    }
	    else {
	      if (indptr[ii].markers[jj][kk] == 3 )
	        fprintf(outf, "%3d ", indptr[ii].markers[jj][kk] - 1);
	      else if (indptr[ii].markers[jj][kk] == 2  )
	        fprintf(outf, "%3d ", indptr[ii].markers[jj][kk] - 1);
	      else
	        fprintf(outf, "%3d ", indptr[ii].markers[jj][kk] );
	    }
	}
	for (kk = 1; kk <= indptr[ii].t; kk++)
	  if ( indptr[ii].y[kk] > (FPN) MISS_VAL)
	    fprintf(outf, "\n %24.12f", indptr[ii].y[kk]);
	  else
	    fprintf(outf, "\n  .  ");
        if ( theparams->themap->otraits > 0 ) 
	for (kk = 1; kk <= theparams->themap->otraits; kk++)
	    fprintf(outf, "\n           %s", indptr[ii].oy[kk]);
      }
    fprintf(outf, "\n-e\n-q\n\n");
    if ( indptr[1].name != NULL && theparams->boot == 0 ) {
      fprintf(outf, "\n-IndividualNames ");
      indcntr = 0;
      for (ii = 1; ii <= n; ii++) /*  Give a table of names at the end, but only if original data (no bootstrapped, permuted, etc, data)*/
        if ( indptr[ii].print_flag == 'y' ) {
	      indcntr = indcntr + 1;
	      if ( indptr[ii].name != NULL)
	        fprintf(outf, "\n%3d   %s ", indcntr,indptr[ii].name );
	    }
      fprintf(outf, "\n-unIndividualNames\n-q\n");
	}
    if (outfile[0] != '-')
      fileclose(outfile, outf);
  }
}


/*   Print the data in the qtlcart.mcd format for use with the windows version.*/
void print_individuals_mcd(params *theparams,individual *indptr,int n, char *outfile) {

  FILE *outf; 
  
  if (*outfile == '-')
    outf = stdout;
  else {  /*  destroy the old output file  */
    outf = fileopen(outfile, "w");
  }
  if (outf != NULL)  
    fprintf(outf, "#FileID  %ld\n#bychromosome", theparams->seed);
  if (*outfile != '-')
    fileclose(outfile, outf);

  print_map_inp(theparams->themap, outfile);
  print_individuals_std(theparams,indptr,n,outfile);


}

/*
 * print out the information on all the individuals in the sample. This will
 * include the values of the markers, the QTLs and the phenotypes. 
 * The format will be that of cross.inp.   Output to
 * a file or stdout.
 */
void print_individuals_std(params *theparams,individual *indptr,int n, char *outfile)
{
  int ii, jj, kk,  indcntr;
  FILE *outf;
  div_t xx;
  
  
    /*if there are no marker names, create bogus names and print them
    to the map file. */
  if ( theparams->themap->names == NULL ) {
    create_bogus_markers(theparams, theparams->themap );
    outf = fileopen(theparams->map, "a");    
    print_marker_names(theparams->themap, outf);
    fileclose(theparams->map, outf);    
  }
  else if ( theparams->themap->tnames == NULL ||  theparams->themap->onames == NULL ) 
    create_bogus_markers(theparams, theparams->themap );


  
  indcntr = 0;
  for (ii = 1; ii <= n; ii++)
    if ( indptr[ii].print_flag == 'y')
      indcntr = indcntr + 1;
  if (*outfile == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  if (outf != NULL) {
    fprintf(outf, "\n#bycross\n#\tThis is a data set in the standard QTL Cartographer input format.");
    fprintf(outf, "\n#\tIt is defined in the manual and in the file cross.inp.");
    fprintf(outf, "\n#\tThis particular file was translated by Rcross.\n\n");
    fprintf(outf, "\n-SampleSize       %4d   is the sample size", indcntr);
    fprintf(outf, "\n-Cross   %s   is the type of cross", theparams->thecross);
    fprintf(outf, "\n-traits  %4d   is the number of traits",theparams->themap->traits);
    fprintf(outf, "\n-otraits %4d   is the number of other traits", theparams->themap->otraits);
    fprintf(outf, "\n-case  yes ");
    fprintf(outf, "\n-TranslationTable    ");
    fprintf(outf, "\n     AA    2    AA  ");
    fprintf(outf, "\n     Aa    1    Aa  ");
    fprintf(outf, "\n     aa    0    aa  ");
    fprintf(outf, "\n     A-   12    A-  ");
    fprintf(outf, "\n     a-   10    a-  ");
    fprintf(outf, "\n     --   -1    --  ");
    fprintf(outf, "\n-missingtrait .   A missing trait is encoded by a solitary period.");

/*  Write out the marker values...  */
    fprintf(outf, "\n-start markers ");
	for (jj = 1; jj <= theparams->themap->m; jj++) {
	  for (kk = 1; kk <= theparams->themap->mpc[jj]; kk++) {
	    if ( theparams->themap->names != NULL ) 
	      fprintf(outf, "\n%s  ", theparams->themap->names[ theparams->themap->ttable[jj][kk] ]);
	    else
	      fprintf(outf, "\n");
	    for ( ii = 1 ; ii <= n; ii++ ) {
	      xx = div(ii,20);
	      if ( xx.quot > 0 && xx.rem == 0 )
	        fprintf(outf,"\n        ");

	      if ( indptr[ii].print_flag == 'y') 
		    switch(indptr[ii].markers[jj][kk]) {
		      case 10 :           fprintf(outf,"a- ");    break;
		      case 0 :            fprintf(outf,"aa ");	  break;
		      case 2 :  case 1 :  fprintf(outf,"Aa ");    break;
		      case 3 :            fprintf(outf,"AA ");    break;
		      case 12 :           fprintf(outf,"A- ");    break;
		      default :           fprintf(outf,"-- ");    break;
		    }
		  }
	    }	      
	 }
     fprintf(outf, "\n-stop markers ");

/*  Write out the trait values...  */
    fprintf(outf, "\n-start traits ");
	for (jj = 1; jj <= theparams->themap->traits; jj++) {
	    if ( theparams->themap->tnames != NULL ) 
	      fprintf(outf, "\n%s  ", theparams->themap->tnames[jj]);
	    else
	      fprintf(outf, "\n");
	    for ( ii = 1 ; ii <= n; ii++ ){
	      xx = div(ii,10);
	      if ( xx.quot > 0 && xx.rem == 0 )
	        fprintf(outf,"\n        ");
	      if ( indptr[ii].print_flag == 'y') {
	        if (indptr[ii].y[jj]  > (FPN) MISS_VAL)
	          fprintf(outf, "%f ", indptr[ii].y[jj] );
	        else
	          fprintf(outf, ". ");
	      }
	    }
	 }
     fprintf(outf, "\n-stop traits ");
/*  Write out the differences in the parental line means for each trait.  
    They should be in the same order as the traits above.   */

    if ( theparams->themap->ParentalDiff != NULL ) {
      fprintf(outf,"\n-ParentalDifferences ");
      for ( ii = 1 ; ii <= theparams->themap->traits ; ii++ )
        fprintf(outf, " %f ", theparams->themap->ParentalDiff[ii]);
    
    
    }



/*  Write out the categorical trait values...  */
    if (  theparams->themap->otraits  > 0 ) {
	    fprintf(outf, "\n-start otraits ");
		for (jj = 1; jj <= theparams->themap->otraits ; jj++) {
		    if ( theparams->themap->onames != NULL ) 
		      fprintf(outf, "\n%s  ", theparams->themap->onames[jj]);
		    else
		      fprintf(outf, "\n");
		    for ( ii = 1 ; ii <= n; ii++ )	{
	          xx = div(ii,20);
	          if ( xx.quot > 0 && xx.rem == 0 )
	            fprintf(outf,"\n        ");	      
		      if ( indptr[ii].print_flag == 'y')
		          fprintf(outf, "%s ", indptr[ii].oy[jj] );	
		    }	     
		 }
		 if ( theparams->tcross == 12 ) {
		   fprintf(outf,"\nBackcross ");
		   for ( ii = 1 ; ii <= n; ii++ )	{
	         xx = div(ii,20);
	         if ( xx.quot > 0 && xx.rem == 0 )
	          fprintf(outf,"\n        ");	      
		       if ( indptr[ii].print_flag == 'y')
		          fprintf(outf, "%d ", indptr[ii].bc);
		    }
		 }
	     fprintf(outf, "\n-stop otraits ");
     }

     if (*outfile != '-')
       fileclose(outfile, outf);
  }
}



/*
 * print out the information on all the individuals in the sample. This will
 * include the values of the markers, the QTLs and the phenotypes. 
 * The format will be that of a MAPMAKER raw file.   Output to
 * a file or stdout.   This will also print a file suitable for Plabqtl.   
 */
void print_individuals_mm(params *theparams,individual *indptr,int n, char *outfile)
{
  int ii, jj, kk,  indcntr,fgen,ralph;
  FILE *outf;
  div_t xx;
  ralph = 1;
  fgen = 2;
    
  
  
  /*if there are no marker names, create bogus names and print them
    to the map file. */
  if ( theparams->themap->names == NULL ) {
    create_bogus_markers(theparams, theparams->themap );
    outf = fileopen(theparams->map, "a");    
    print_marker_names(theparams->themap, outf);
    fileclose(theparams->map, outf);    
  }
  else if ( theparams->themap->tnames == NULL ||  theparams->themap->onames == NULL ) 
    create_bogus_markers(theparams, theparams->themap );

  indcntr = 0;
  for (ii = 1; ii <= n; ii++)
    if (indptr[ii].print_flag == 'y')
      indcntr = indcntr + 1;
  if (*outfile == '-')
    outf = stdout;
  else 
    outf = fileopen(outfile, "w");
  if ( theparams->gout == 2 )
    fprintf(outf,"# Rcross output file  -filetype mapmaker.raw\n");
  if (outf != NULL) {
    fprintf(outf, "data type ");
    switch (theparams->cross) {
      case 1:  case 2:  default : fprintf(outf, "f2 backcross");  break;
      case 3:  case 4:
        fgen = theparams->crosst;
        if ( theparams->crosst == 3 )
          fprintf(outf, "f3 self");
        else
          fprintf(outf, "f2 intercross");
        break;
      case 5:
        if ( theparams->crosst == 1 )
          fprintf(outf, "ri self");
        else
          fprintf(outf, "ri sib");
        if ( theparams->crosst == 0 )
          fgen = 1;
        else 
          ralph = 2;
      break;
    }  
    if ( theparams->gout == 6 )
      fprintf(outf, "  Rcross output file -filetype plabqtl1.qdt  ");  
    fprintf(outf, "\n%d %d %d ", indcntr, theparams->themap->ml, theparams->themap->traits);
    if ( theparams->gout == 6 ) {
      fprintf(outf, "  %d 1 0 0 %d %d \n",theparams->themap->m, fgen, ralph );  
      for ( jj=1; jj<= theparams->themap->traits; jj++ )
        fprintf(outf,"0.0 ");
    }

/*  Write out the marker values...  */

	for (jj = 1; jj <= theparams->themap->m; jj++) {

	  for (kk = 1; kk <= theparams->themap->mpc[jj]; kk++) {
	    if ( theparams->themap->names != NULL ) 
	      fprintf(outf, "\n*%s  ", theparams->themap->names[ theparams->themap->ttable[jj][kk] ]);
	    else
	      fprintf(outf, "\n");
	    for ( ii = 1 ; ii <= n; ii++ ) {
	      xx = div(ii,60);
	      if ( xx.quot > 0 && xx.rem == 0 )
	        fprintf(outf,"\n        ");
	      if ( indptr[ii].print_flag == 'y') 
		    switch(indptr[ii].markers[jj][kk]) {
		      case 10 :          fprintf(outf,"C");	 break;
		      case 0 :           fprintf(outf,"B");	 break;
		      case 1 : case 2 :  fprintf(outf,"H");	 break;
		      case 3 :           fprintf(outf,"A");  break;
		      case 12 :          fprintf(outf,"D");  break;
		      default :          fprintf(outf,"-");  break;
		    }
	    }
	  }	      
	 }


/*  Write out the trait values...  */

	for (jj = 1; jj <= indptr[1].t; jj++) {
	  if ( theparams->themap->tnames != NULL ) 
	    fprintf(outf, "\n*%s  ", theparams->themap->tnames[jj ]);
	  else
	    fprintf(outf, "\n        ");
	  for ( ii = 1 ; ii <= n; ii++ ) {
	    xx = div(ii,10);
	    if ( xx.quot > 0 && xx.rem == 0 )
	      fprintf(outf,"\n");
	    if ( (indptr + ii)->print_flag == 'y') {
	      if (indptr[ii].y[jj] > (FPN) MISS_VAL)
	        fprintf(outf, "%f ", indptr[ii].y[jj]);
	      else
	        fprintf(outf, " - ");
	    }
	  }
	}
	
	if ( theparams->gout == 6 )
	  print_plabqtl_map(theparams->themap,outf);
	
	
    if (*outfile != '-')
      fileclose(outfile, outf);
  }
}

/*
 *   This will create a  Plabqtl input file in matrix format.   
 */
void print_plabqtl(params *theparams,individual *indptr,int n, char *outfile)
{
  int ii, jj, kk,  indcntr,fgen,ralph;
  FILE *outf;
 /* div_t xx;*/
  ralph = 1;
  fgen = 2;
     
  /*if there are no marker names, create bogus names and print them
    to the map file. */
  if ( indptr[1].map->names == NULL ) {
    create_bogus_markers(theparams, indptr[1].map );
    outf = fileopen(theparams->map, "a");    
    print_marker_names(indptr[1].map, outf);
    fileclose(theparams->map, outf);    
  }
  else if ( indptr[1].map->tnames == NULL ||  indptr[1].map->onames == NULL ) 
    create_bogus_markers(theparams, indptr[1].map );

  indcntr = 0;
  for (ii = 1; ii <= n; ii++)
    if (indptr[ii].print_flag == 'y')
      indcntr = indcntr + 1;
  if (*outfile == '-')
    outf = stdout;
  else 
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
    fprintf(outf, "data type ");
    switch (theparams->cross) {
      case 1:  case 2:  default : fprintf(outf, "f2 backcross");  break;
      case 3:  case 4:
        fgen = theparams->crosst;
        if ( theparams->crosst == 3 )
          fprintf(outf, "f3 self");
        else
          fprintf(outf, "f2 intercross");
        break;
      case 5:
        if ( theparams->crosst == 1 )
          fprintf(outf, "ri self");
        else
          fprintf(outf, "ri sib");
        if ( theparams->crosst == 0 )
          fgen = 1;
        else 
          ralph = 2;
      break;
    }  
    fprintf(outf, "  Rcross output file -filetype plabqtl0.qdt  ");  
    fprintf(outf, "\n%d %d %d ", indcntr, indptr[1].map->ml, indptr[1].t);
    fprintf(outf, "  %d 0 0 0 %d %d \n",indptr[1].map->m, fgen, ralph );  
    for ( jj=1; jj<= indptr[1].t; jj++ )
        fprintf(outf,"0.0 ");
    fprintf(outf,"\n");
    for ( jj=1; jj<= indptr[1].t; jj++ )
        fprintf(outf,"%s ",indptr[1].map->tnames[jj] );
	for (jj = 1; jj <= indptr[1].map->m; jj++) {
      fprintf(outf,"\n");
	  for (kk = 1; kk <=  indptr[1].map->mpc[jj]; kk++) 
	      fprintf(outf, "%s ", indptr[1].map->names[ indptr[1].map->ttable[jj][kk] ]);
    }
	
	for ( ii = 1 ; ii <= n; ii++ ) 
	  if ( indptr[ii].print_flag == 'y') {
        fprintf(outf,"\n");
	    for (jj = 1; jj <= indptr[1].t; jj++) 
 	      if (indptr[ii].y[jj] > (FPN) MISS_VAL)
	        fprintf(outf, "%f ", indptr[ii].y[jj]);
	      else
	        fprintf(outf, " - ");
	    for (jj = 1; jj <= indptr[1].map->m; jj++) {
	      fprintf(outf,"\n");
	      for (kk = 1; kk <=  indptr[1].map->mpc[jj]; kk++) 
		    switch(indptr[ii].markers[jj][kk]) {
		      case 10 :          fprintf(outf,"C ");	 break;
		      case 0 :           fprintf(outf,"B ");	 break;
		      case 1 : case 2 :  fprintf(outf,"H ");	 break;
		      case 3 :           fprintf(outf,"A ");  break;
		      case 12 :          fprintf(outf,"D ");  break;
		      default :          fprintf(outf,"- ");  break;
		    }
	    }
	  }	      

      
      

	
	  print_plabqtl_map(indptr[1].map,outf);
	
	
    if (*outfile != '-')
      fileclose(outfile, outf);
  }
}
/*
 *       This is an way to move your data into R or Splus for any analyses you 
 *    might want to do there.   
 *
 * Create a file with all the trait, marker and categorical data that is suitable
 * for loading into Splus or R.   In addition, create data frames for the different
 * types of data and set up a set of commands to do One Factor ANOVA of each trait
 * on each marker, and each trait on each categorical trait.   From there, the user
 * can do higher level ANOVA analyses in R or Splus.   The output data frames will
 * be Markers$....,  Traits$...., and Otraits$....   If each marker has a name, then
 * data for the marker 'markerXY' is in Markers$markerXY .   If the marker has no name,
 * then the data for marker Y on chromosome X is in Markers$cXmY (note that X and Y will 
 * are integer variables).   A similar convention holds for traits and categorical traits:
 * A trait named 'traitZ' will be in Traits$traitZ, while a categorical trait called 'otraitW'
 * will be in OTraits$otraitW .   If unamed, then the trait numbered Z will go into Traits$tZ and
 * categorical trait W will go into OTraits$oW .   
 */
void print_individuals_R(params *theparams,individual *indptr,int n ,char *outfile)
{
  int ii, jj, kk,  indcntr;
  FILE *outf;
  div_t xx;
  ii = theparams->traits;
  indcntr = 0;
  for (ii = 1; ii <= n; ii++)
    if (indptr[ii].print_flag == 'y')
      indcntr = indcntr + 1;
  if (*outfile == '-')
    outf = stdout;
  else 
    outf = fileopen(outfile, "a");
  if (outf != NULL) {


/*  Write out the marker values in list form...  */

	for (jj = 1; jj <= indptr[1].map->m; jj++) {
	  for (kk = 1; kk <= indptr[1].map->mpc[jj]; kk++) {
	    if ( indptr[1].map->names != NULL ) 
	      fprintf(outf, "\n%s <- c( ", indptr[1].map->names[ indptr[1].map->ttable[jj][kk] ]);
	    else
	      fprintf(outf, "\nc%dm%d <- c( ", jj,kk);
	    
	    for ( ii = 1 ; ii <= n; ii++ ) {
	      xx = div(ii,10);
	      if ( xx.quot > 0 && xx.rem == 0 )
	        fprintf(outf,"\n               ");
	      if ( indptr[ii].print_flag == 'y') 
		    switch(indptr[ii].markers[jj][kk]) {
		      case 10 :          fprintf(outf,"\"a-\"");	 break;
		      case 0 :           fprintf(outf,"\"aa\"");	 break;
		      case 1 : case 2 :  fprintf(outf,"\"Aa\"");	 break;
		      case 3 :           fprintf(outf,"\"AA\"");  break;
		      case 12 :          fprintf(outf,"\"A-\"");  break;
		      default :          fprintf(outf,"NA");  break;
		    }
		  if ( ii < n )
		    fprintf(outf,", ");
	    }
	    fprintf(outf, ") \n" );

	  }	      
	 }

/*   Create data frames for the marker values.  The data frames allow the markers to be used in ANOVA  */
	fprintf(outf, "\nMarkers <- data.frame(" );
	for (jj = 1; jj <= indptr[1].map->m; jj++) {
	  for (kk = 1; kk <= indptr[1].map->mpc[jj]; kk++) {
	      xx = div(kk,8);
	      if ( xx.quot > 0 && xx.rem == 0 )
	        fprintf(outf,"\n     ");
	    if ( indptr[1].map->names != NULL ) 
	      fprintf(outf, "%s", indptr[1].map->names[ indptr[1].map->ttable[jj][kk] ]);
	    else
	      fprintf(outf, "c%dm%d", jj,kk);
	    if ( kk < indptr[1].map->mpc[jj]) 
	      fprintf(outf, "," );
	      

	  }	
	  if ( jj < indptr[1].map->m )
	      fprintf(outf, ",\n     " );
	  else
	      fprintf(outf, ")\n");
	        
	 }



/*  Write out the categorical trait values.  */
    if (  indptr[1].map->otraits  > 0 ) {
		for (jj = 1; jj <= indptr[1].map->otraits ; jj++) {
		    if ( indptr[1].map->onames != NULL ) 
		      fprintf(outf, "\n%s <- c(  ", indptr[1].map->onames[jj]);
		    else
		      fprintf(outf, "\no%d <- c(  ", jj);
		    for ( ii = 1 ; ii <= n; ii++ )	{
	          xx = div(ii,20);
	          if ( xx.quot > 0 && xx.rem == 0 )
	            fprintf(outf,"\n        ");	      
		      if ( indptr[ii].print_flag == 'y')
		          fprintf(outf, "\"%s\"", indptr[ii].oy[jj] );	
		      if ( ii < n )
		          fprintf(outf, ", ");
		    }	     
	     fprintf(outf, ")\n");
		 }
     }

/*  and create data frames for the categorical traits.  Again, these can be use in ANOVA now that they are in frames. */
    if (  indptr[1].map->otraits  > 0 ) {
	  fprintf(outf, "\nOTraits <- data.frame(" );
	  for (jj = 1; jj <= indptr[1].map->otraits ; jj++) {
	    xx = div(jj,8);
	    if ( xx.quot > 0 && xx.rem == 0 )
	      fprintf(outf,"\n     ");	      
		if ( indptr[1].map->onames != NULL ) 
		  fprintf(outf, "%s", indptr[1].map->onames[jj]);
		else
		  fprintf(outf, "o%d", jj);
		if ( jj < indptr[1].map->otraits )
		  fprintf(outf, ", ");
	    else
	      fprintf(outf, ")\n");
	  }
     }
/*  Write out the trait values.    */
	for (jj = 1; jj <= indptr[1].t; jj++) {
	  if ( indptr[1].map->tnames != NULL ) 
	    fprintf(outf, "\n%s <- c( ",indptr[1].map->tnames[jj]);
	  else
	    fprintf(outf, "\nt%d <- c( ", jj);
	  for ( ii = 1 ; ii <= n; ii++ ) {
	    xx = div(ii,6);
	    if ( xx.quot > 0 && xx.rem == 0 )
	      fprintf(outf,"\n     ");
	    if ( (indptr + ii)->print_flag == 'y') {
	      if (indptr[ii].y[jj] > (FPN) MISS_VAL)
	        fprintf(outf, "%f", indptr[ii].y[jj]);
	      else
	        fprintf(outf, "NA");
	    }
		  if ( ii < n )
		    fprintf(outf,", ");
	  }
	  fprintf(outf, ")\n");
	}


/*  Create data frames for the Traits.   We do this to be consistent.     */
	fprintf(outf, "\nTraits <- data.frame(" );
	for (jj = 1; jj <= indptr[1].t; jj++) {
	  xx = div(jj,6);
	  if ( xx.quot > 0 && xx.rem == 0 )
	    fprintf(outf,"\n     ");
	  if ( indptr[1].map->tnames != NULL ) 
	    fprintf(outf, "%s",indptr[1].map->tnames[jj]);
	  else
	    fprintf(outf, "t%d", jj);
	  if ( jj < indptr[1].t )
		fprintf(outf,", ");
	  else
	    fprintf(outf, ")\n");
	}

/*  Do commands for the analysis .   These commands will do 1 way ANOVA for each trait on 
    each categorical trait, and each trait on each marker.    The data frames are called
    Traits, OTraits and Markers .  */

	for (jj = 1; jj <= indptr[1].t; jj++) {

    if (  indptr[1].map->otraits  > 0 ) {
		for (ii = 1; ii <= indptr[1].map->otraits ; ii++) {
	      if ( indptr[1].map->tnames != NULL )
	        fprintf(outf,"\nout <- lm( Traits$%s",indptr[1].map->tnames[jj]);
	      else
	        fprintf(outf,"\nout <- lm( Traits$t%d",jj);

		    if ( indptr[1].map->onames != NULL ) 
		      fprintf(outf, " ~ OTraits$%s", indptr[1].map->onames[ii]);
		    else
		      fprintf(outf, " ~ OTraits$o%d", ii);

	      fprintf(outf," )\nprint(result.t%do%d   <- anova(out))\n",jj,ii);


		 }
     }

	  for (ii = 1; ii <= indptr[1].map->m; ii++) {
	    for (kk = 1; kk <= indptr[1].map->mpc[ii]; kk++) {
	      
	      if ( indptr[1].map->tnames != NULL )
	        fprintf(outf,"\nout <- lm( Traits$%s",indptr[1].map->tnames[jj]);
	      else
	        fprintf(outf,"\nout <- lm( Trait$t%d",jj);
	      if ( indptr[ii].map->names != NULL )
	        fprintf(outf," ~ Markers$%s",indptr[1].map->names[ indptr[1].map->ttable[ii][kk] ]);
	      else
	        fprintf(outf," ~ Markers$c%dm%d",ii,kk);
	       
	      fprintf(outf," )\nprint(result.t%dc%dm%d   <- anova(out))\n",jj,ii,kk);
        }

      }
	}
    if (*outfile != '-')
      fileclose(outfile, outf);
  }
}
/*
 *       This is a way to move your data into SAS.
 */
void print_individuals_SAS(params *theparams,individual *indptr,int n ,char *outfile)
{
  int ii, jj, kk,  indcntr;
  FILE *outf;

  indcntr = 0;
  for (ii = 1; ii <= n; ii++)
    if (indptr[ii].print_flag == 'y')
      indcntr = indcntr + 1;
  if (*outfile == '-')
    outf = stdout;
  else 
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
    shift_fn(theparams->stem);
    fprintf(outf, "DATA %s;\nINPUT ",theparams->stem);
    insert_wd(gname,theparams->workdir, theparams->stem);
    if (  indptr[1].map->otraits  > 0 ) /* Categorical traits first */
		for (jj = 1; jj <= indptr[1].map->otraits ; jj++) 
		    if ( indptr[1].map->onames != NULL ) 
		      fprintf(outf, "%s $ ", indptr[1].map->onames[jj]);
		    else
		      fprintf(outf, "o%d $ ", jj);
	for (jj = 1; jj <= (indptr + 1)->map->m; jj++) { /*  markers next */
	  for (kk = 1; kk <= indptr[1].map->mpc[jj]; kk++) 
	    if ( indptr[1].map->names != NULL ) 
	      fprintf(outf, "%s ", indptr[1].map->names[ indptr[1].map->ttable[jj][kk] ]);
	    else
	      fprintf(outf, "c%dm%d ", jj,kk);
	  fprintf(outf,"\n");
    }
	for (jj = 1; jj <= indptr[1].t; jj++)   /* Now the traits */
	  if ( indptr[1].map->tnames != NULL ) 
	    fprintf(outf, "%s ",indptr[1].map->tnames[jj]);
	  else
	    fprintf(outf, "t%d ", jj);
	fprintf(outf,";\nCARDS;");
/*  Write out the data first. */
    for ( ii = 1 ; ii <= n; ii++ ) 
      if (indptr[ii].print_flag == 'y' ) {
        fprintf(outf,"\n");
        if (  indptr[1].map->otraits  > 0 ) /* Categorical traits first */
		  for (jj = 1; jj <= indptr[1].map->otraits ; jj++) 
		      fprintf(outf, "%2s ",  indptr[ii].oy[jj]  );
	    for (jj = 1; jj <= indptr[1].map->m; jj++) {
	      for (kk = 1; kk <= indptr[1].map->mpc[jj]; kk++) 
		    switch(indptr[ii].markers[jj][kk]) {
		      case 10 :          fprintf(outf,"1 ");	 break;
		      case 0 :           fprintf(outf,"0 ");	 break;
		      case 1 : case 2 :  fprintf(outf,"2 ");	 break;
		      case 3 :           fprintf(outf,"4 ");  break;
		      case 12 :          fprintf(outf,"3 ");  break;
		      default :          fprintf(outf,". ");  break;
		    }
        } 

/*  Write out the trait values.    */
	   for (jj = 1; jj <= indptr[1].t; jj++) {
 	      if (indptr[ii].y[jj] > (FPN) MISS_VAL)
	        fprintf(outf, "%f ", indptr[ii].y[jj]);
	      else
	        fprintf(outf, ". ");
	    }
 	}

	fprintf(outf,"\n;\n");

	for (jj = 1; jj <= indptr[1].t; jj++)  { /* Now write out the proc anova statements */
	  if ( indptr[1].map->tnames != NULL ) 
	    sprintf(gname, "%s\0",indptr[1].map->tnames[jj]);
	  else
	    sprintf(gname, "t%d\0", jj);
      if (  indptr[1].map->otraits  > 0 ) /* Categorical traits first */
		for (ii = 1;ii <= indptr[1].map->otraits ; ii++) {
		    if ( indptr[1].map->onames != NULL ) 
		      sprintf(gbuffer, "%s\0", indptr[1].map->onames[ii]);
		    else
		      sprintf(gbuffer, "o%d\0", ii);
            if ( theparams->gout == 4 )
              write_proc(gbuffer,gname,outf);
            else 
              write_procANOVA(gbuffer,gname,outf);
         }
         
	  for (ii = 1; ii <= (indptr + 1)->map->m; ii++) { /*  markers next */
	    for (kk = 1; kk <=  indptr[1].map->mpc[ii]; kk++) {
	      if ( indptr[1].map->names != NULL ) 
	        sprintf(gbuffer, "%s\0", indptr[1].map->names[ indptr[1].map->ttable[ii][kk] ]);
	      else
	        sprintf(gbuffer, "c%dm%d\0", ii,kk);
            if ( theparams->gout == 4 )
              write_proc(gbuffer,gname,outf);
            else 
              write_procANOVA(gbuffer,gname,outf);
	    }
      }
         
         
    }
    fprintf(outf,"\n\nRUN;\n\n");
    if (*outfile != '-')
      fileclose(outfile, outf);
  }
}

/*
   Write the proc commands for SAS analysis
*/
void write_procANOVA(char *x, char *y,FILE *outf) {
  fprintf(outf,"\n\nPROC ANOVA;");
  fprintf(outf,"\n\tCLASS %s;",x);
  fprintf(outf,"\n\tMODEL %s=%s;",y,x);
  fprintf(outf,"\n\tMEANS %s;",x);
}

/*
   Write the proc commands for SAS analysis
*/
void write_proc(char *x, char *y,FILE *outf) {
  fprintf(outf,"\n\nPROC GLM;");
  fprintf(outf,"\n\tCLASS %s;",x);
  fprintf(outf,"\n\tMODEL %s=%s;",y,x);
  fprintf(outf,"\n\tMEANS %s;",x);
}


/*
 *
 * The following are assumed: Value  Genotype  Explanation
 *
 *
gt          1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
alleles    11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44

11 = A1 A1  (ordered)


 *
 *
 * mgt is the maternal genotype, mgamete tells which gamete is used pgt is the
 * paternal genotype, pgamete tells which gamete is used
 */
int first_gam(int gt)
{
int gam;
  switch (gt) {
    case 1: case 2: case 3: case 4:  gam = 1; break;
    case 5: case 6: case 7: case 8: gam = 2; break;
    case 9: case 10: case 11: case 12: gam = 3; break;
    default:  gam = 4; break;
  } 
  return(gam);
}

int second_gam(int gt)
{
int gam;
  switch (gt) {
    case 1: case 5: case 9: case 13: gam = 1; break;
    case 2: case 6: case 10: case 14: gam = 2; break;
    case 3: case 7: case 11: case 15:  gam = 3; break;
    default:  gam = 4; break;
  } 
  return(gam);
}


int sc_genotype(int mgt,int pgt,FPN mgamete,FPN pgamete)
{
  int ogt, mgam, pgam;
  
  if ( mgamete <=0.5 )
    mgam = first_gam(mgt);
  else 
    mgam = second_gam(mgt);
  if ( pgamete <=0.5 )
    pgam = first_gam(pgt);
  else 
    pgam = second_gam(pgt);
  ogt = 4*(pgam-1) + mgam;
  return (ogt);
}



/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Idatain.c
------------------------------------------------------------------ */


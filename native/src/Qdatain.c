/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Qdatain.c) is part of QTL Cartographer
         
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
          Qdatain.c,  subroutines to read in the genetic model from
          an Rqtl.out file.
*/

#include "Main.h"


/*
    Open the file 'infile' and read in a Model:  a set of 
    QTLs with positions and effects.   theparams->theqtls will 
    point to it.  

*/
void GetTheModel(params *theparams, char *infile) {
	int ii, knum, error;


    ii = get_file_type(infile);
    if ( ii == 21 )   
      theparams->theqtls   = get_std_qtls(theparams->themap, infile);
    else if (  ii == 20 || ii == 70 || ii == 100 ) {
      knum = get_knum(infile);
      if (knum > 1) {
        theparams->themap->traits = knum;
        free_ivector(theparams->themap->knum, 1, 1);
        theparams->themap->knum = ivector(1, theparams->themap->traits);
      }
      get_knums(infile, theparams->themap);
      theparams->theqtls = qtlvector(theparams->themap);
      error = get_qtls(theparams->theqtls, infile);
      if (error != 0) {
        sprintf(gwarn, "Had trouble reading the QTL data from %s\n", infile);
        LogTheError(theparams->error, gwarn);
        exit(1);
      }
   }
   else
     theparams->theqtls = NULL ;
}


/*
*     This reads in the number of loci for each trait.
*     It moves to the first model available and reads that
*     information.  If Rqtl or Eqtl output, then it will get
*     the number of loci for each trait.  If MImapqtl output, then
*     it will only get the first available model.
*/

void get_knums(char *inputf,markermap *themap)
{
  FILE *infile;
  int  knum,  trait, ch,go_on,i;
/*  filetype = get_file_type(inputf);*/
  for ( i=1; i<= themap->traits; i++ )
    themap->knum[i] = 0;
  infile = fileopen(inputf, "r");
  
  go_on = movetoqtls(infile);

  do {
    ch =  get_next_token(gbuffer, MAXLINE, infile);
    if (gbuffer[0] == '-')
      switch (gbuffer[1]) {
       case 'k':
         ch =  get_next_token(gbuffer, MAXLINE, infile);
	     knum = atoi(gbuffer);
	     break;
       case 'n':
         if ( gbuffer[2] == 'u' ) {
           ch =  get_next_token(gbuffer, MAXLINE, infile);
	       trait = atoi(gbuffer);
	       themap->knum[trait] = knum;
	     }
	     break;
       case 'e':
         ch =  EOF;
	     break;
       default:
	     break;
      }
  } while (ch != EOF);  
  fileclose(inputf, infile);
   
}

/*
*  Reads the number of traits in the first encountered model.
*/

int get_knum(char *inputf)
{
  FILE *infile;
  int  knum, ch,go_on;
/*  filetype = get_file_type(inputf);*/
  
  infile = fileopen(inputf, "r");
  if (infile == NULL) {
    return (-1);
  }
  
  go_on = movetoqtls(infile);
  
  do {
    ch =  get_next_line(gbuffer, MAXLINE, infile);
    if (gbuffer[0] == '-')
      switch ( gbuffer[1]) {
        case 't':
		  get_field(2, gname, gbuffer);
		  knum = atoi(gname);
		  fileclose(inputf, infile);
		  return (knum);
		  break;
        default:
	      break;
      }
  } while (ch != EOF);
  if (ch == EOF)
    printf("\nTime to close the input file after an EOF...No QTLS? \n");
  fileclose(inputf, infile);
  return (-1);
}

/*
Open the Rqtl.out formatted file and get the genetic model.  This only reads
the first model encountered.


episAADD = dmatrix(1,themap->knum[trait],1,themap->knum[trait])
episADDA = dmatrix(1,themap->knum[trait],1,themap->knum[trait])

episAADD has AA interactions in lower half of matrix ( r > s and episAADD[r][s] )
episAADD has DD interactions in upper half of matrix ( r < s and episAADD[r][s] )
 
episADDA has AD interactions in lower half of matrix ( r > s and episADDA[r][s] )
episADDA has DA interactions in upper half of matrix ( r < s and episADDA[r][s] )
 
r must always be less than  s in the input 
-i  r  s  abd   value
 
 r = 1,knum[trait]-1
 s = r+1,knum[trait]
 abd = AA, AD, DA, DD
 value is a FPN
*/
int get_qtls(qtlptr, inputf)
aqtl *qtlptr;
char *inputf;
{
  FILE *infile;
  int  row, trait, ch,go_on,r,s,abd, nqtl;
  FPN **episAADD,**episADDA,value;
  
/*  filetype = get_file_type(inputf);*/
  
  infile = fileopen(inputf, "r");
  if (infile == NULL) {
    return (-1);
  }
  go_on = movetoqtls(infile);
  if (go_on == -1) {
    return (-1);
  }

  trait = row = 0;
      
  do {
    ch =  get_next_token(gbuffer, MAXLINE, infile);
    if (gbuffer[0] == '-')
      switch (gbuffer[1]) {
       case 'k':
         ch =  get_next_token(gbuffer, MAXLINE, infile);
		 nqtl = atoi(gbuffer);
		 if ( nqtl > 0 ) {
	       episAADD = dmatrix(1,nqtl,1,nqtl);
	       episADDA = dmatrix(1,nqtl,1,nqtl);
	     }
	     else {
	       episAADD = NULL;
	       episADDA = NULL;
	     }
	     break;
       case 'l':
		row = row + 1;
        ch =  get_next_token(gbuffer, MAXLINE, infile);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
		qtlptr[row].chrm = atoi(gbuffer);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
		qtlptr[row].mrk = atoi(gbuffer);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
		qtlptr[row].c1 = (FPN) atof(gbuffer);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
		qtlptr[row].c2 = (FPN) atof(gbuffer);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
		qtlptr[row].a = (FPN) atof(gbuffer);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
		qtlptr[row].d = (FPN) atof(gbuffer);
		qtlptr[row].trait = trait;
		qtlptr[row].episAADD = episAADD;
		qtlptr[row].episADDA = episADDA;
		break;
       case 'n':
         if ( gbuffer[2] == 'u' ) {
           ch =  get_next_token(gbuffer, MAXLINE, infile);
		   trait = atoi(gbuffer);          
         }
       
         break;
       case 'i':
        ch =  get_next_token(gbuffer, MAXLINE, infile);
        r = atoi(gbuffer);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
        s = atoi(gbuffer);
        ch =  get_next_token(gbuffer, MAXLINE, infile);
        if ( gbuffer[0] == 'A' && gbuffer[1] == 'A' )
          abd = 3;
        else if ( gbuffer[0] == 'A' && gbuffer[1] == 'D' )
          abd = 4;
        else if ( gbuffer[0] == 'D' && gbuffer[1] == 'A' )
          abd = 5;
        else if ( gbuffer[0] == 'D' && gbuffer[1] == 'D' )
          abd = 6;
        ch =  get_next_token(gbuffer, MAXLINE, infile);
        value = (FPN) atof(gbuffer);
        if ( abd == 3)
          qtlptr[row].episAADD[s][r] = value;
        else if ( abd == 4 )
          qtlptr[row].episADDA[s][r] = value;
        else if ( abd == 5 )
          qtlptr[row].episADDA[r][s] = value;
        else if ( abd == 6 )
          qtlptr[row].episAADD[r][s] = value;        
		break;
       default:
	    break;
      }
  } while (ch != EOF);
  fileclose(inputf, infile);
  return (0);
}





/*Get the qtls in the qtl.inp formatted file. */
aqtl *get_std_qtls(markermap *themap,char *infile)
{
  FILE *fptr;
  aqtl *theqtls;

  int ch;
  int   isu, isn, ise, isskip, go_on;

  fptr = fileopen(infile, "r");
  if (fptr == NULL) {
    return (NULL);
  }
  go_on = 1;
  isu = isn = ise = isskip = 0;
  while (go_on == 1) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if ( gname[0] == '-' && isskip == 0)
      switch ( gname[1]) {
        case 'U':
	      switch (gname[0]) {
	        case 'M':  case 'm': isu = -1;  break;
	        case 'c':  case 'C':  isu = -2;  break;
	        default:  isu = -2;  break;
	      }
	      break;
        case 'u': isskip = 0; break;
        case 'i': read_interactions(fptr,themap,theqtls); break;
        case 'n':
	      ch = get_next_token(gname, MAXNAME, fptr);
	      if (gname[0] == 'n' || gname[0] == 'N')
	        isn = 0;
	      else
	        isn = 1;
	      break;
        case 's':
	      if (gname[3] == 'a')
	        theqtls = read_qtls(fptr, isu, themap, isn);
	      else if (gname[3] == 'i')
	        isskip = 1;
	      break;
        case 'e':  case 'q': go_on = 0; break;
        default: break;
      }
    else if (!strcmp(gname, "-unskip"))
      isskip = 0;
  }
  fileclose(infile, fptr);
  return (theqtls);
}

aqtl *read_qtls(FILE *fptr,int isu,markermap *themap, int isn)
{
  int ch;
  genome *first, *gptr;
  FPN position, tpos,**episAADD,**episADDA;
  long filepos;
  int error, total, loci, jj, ii, traits, atqtl, chrom;
  aqtl *theqtls;
  first = create_genome(themap);
  atqtl = 0;
/* Determine how many total loci */
  ch = get_next_token(gname, MAXNAME, fptr);
  ch = get_next_token(gname, MAXNAME, fptr);
  filepos = ftell(fptr);
  traits = atoi(gname);
  themap->traits = traits;
  free_ivector(themap->knum, 1, 1);
  themap->knum = ivector(1, themap->traits);
  if (isn == 1)
    themap->tnames = cmatrix(1, themap->traits, 0, MAXNAME);
  total = 0;
  for (ii = 1; ii <= traits; ii++) {
    if (isn == 1) {
      ch = get_next_token(gname, MAXNAME, fptr);
      strcpy(themap->tnames[ii], gname);
    }
    ch = get_next_token(gname, MAXNAME, fptr);
    themap->knum[ii] = atoi(gname);
    total = total + themap->knum[ii];
    for (jj = 1; jj <= 4 *  themap->knum[ii]; jj++)
      ch = get_next_token(gname, MAXNAME, fptr);
  }
  theqtls = qtlvector(themap);
  error = fseek(fptr, filepos, 0);
  for (ii = 1; ii <= traits; ii++) {
    if (isn == 1)
      ch = get_next_token(gname, MAXNAME, fptr);

    ch = get_next_token(gname, MAXNAME, fptr);
    loci = atoi(gname);
    if ( loci > 0 ) {
	  episAADD = dmatrix(1,loci,1,loci);
	  episADDA = dmatrix(1,loci,1,loci);
	}
	else {
	  episAADD = NULL;
	  episADDA = NULL;
    }
    for (jj = 1; jj <= loci; jj++) {
      atqtl = atqtl + 1;
      theqtls[atqtl].episAADD = episAADD;
      theqtls[atqtl].episADDA = episADDA;
      theqtls[atqtl].map = themap;
      theqtls[atqtl].trait = ii;
      ch = get_next_token(gname, MAXNAME, fptr);
      chrom = atoi(gname);
      theqtls[atqtl].chrm = chrom;
      ch = get_next_token(gname, MAXNAME, fptr);
      position = (FPN) atof(gname);
      if (isu == -2)
	    position = position / (FPN) 100.0;

      gptr = first;

      while (gptr != NULL && gptr->chrom != chrom)
	    gptr = gptr->next;
      tpos = 0;
      while (gptr != NULL && tpos < position) {
	    tpos = tpos + gptr->dist;
	    gptr = gptr->next;
      }
      if (gptr != NULL) {
	    if (gptr->prev != NULL)
	      gptr = gptr->prev;
	    theqtls[atqtl].mrk = gptr->markr;
	    theqtls[atqtl].c2 = mapfunc(tpos - position, -1);
	    theqtls[atqtl].c1 = mapfunc(position - tpos + gptr->dist, -1);
      }
      else {
	    theqtls[atqtl].mrk =  themap->mpc[chrom];
	    theqtls[atqtl].c2 = (FPN) 0.00001;
	    theqtls[atqtl].c1 = (FPN) 0.00001;
      }

      ch = get_next_token(gname, MAXNAME, fptr);
      theqtls[atqtl].a = (FPN) atof(gname);
      ch = get_next_token(gname, MAXNAME, fptr);
      theqtls[atqtl].d = (FPN) atof(gname);

    }
  }
  clear_genome(first);
  return (theqtls);
}

/*
Read the interaction effects from a QTLs.inp file
*/
void read_interactions(FILE *fptr, markermap *themap, aqtl *qtlptr)
{
  int row,ch,go_on,trait,abd,k,kk,r,s;
  FPN value;
  trait = kk = 0;
  for ( k=1; k<=themap->traits; k++ )
    kk = kk + themap->knum[k];  /* Size of qtlptr*/
  go_on = 1;
  do {
    ch =  get_next_line(gbuffer, MAXLINE, fptr);
    get_field(1, gname, gbuffer);
    if ( !strncmp(gname,"-trait",6) ) {
      trait = trait+1;
      for ( k=1; k<=kk; k++ )
        if ( qtlptr[k].trait == trait )
          row = k;  /* Last qtlptr that is of this trait.  Since all for the same trait point to the same matrices, this will work */
    }
    else if ( !strncmp(gname,"-stop",5) ) 
      go_on = 0;
    else {
        get_field(1, gname, gbuffer);
        r = atoi(gname);
        get_field(2, gname, gbuffer);
        s = atoi(gname);
        get_field(3, gname, gbuffer);
        if ( gname[0] == 'A' && gname[1] == 'A' )
          abd = 3;
        else if ( gname[0] == 'A' && gname[1] == 'D' )
          abd = 4;
        else if ( gname[0] == 'D' && gname[1] == 'A' )
          abd = 5;
        else if ( gname[0] == 'D' && gname[1] == 'D' )
          abd = 6;
        get_field(4, gname, gbuffer);
        value = (FPN) atof(gname);
        if ( abd == 3 )
          qtlptr[row].episAADD[s][r] = value;
        else if ( abd == 4 )
          qtlptr[row].episADDA[s][r] = value;
        else if ( abd == 5 )
          qtlptr[row].episADDA[r][s] = value;
        else if ( abd == 6 )
          qtlptr[row].episAADD[r][s] = value;  
     }
   } while (go_on == 1 );      

}

/*

Start at current place in infile stream.
Go to the next QTL model.
If EOF, return -1, else return 0 for a new model.

*/
int movetoqtls(FILE *infile)
{
  int go_on,ch;
/*If this is an Eqtl.out or MImapqtl.out file, go to the first -start token */
    go_on = 1;
    do {
      ch = get_next_token(gname,MAXNAME,infile);
      if ( !strcmp(gname,"-start") ) {
        ch = get_next_token(gname,MAXNAME,infile);
        if ( !strcmp(gname,"model") )
          go_on = 0;
      }
      if (ch == EOF )
        return( -1);
    } while (go_on == 1);
    return(go_on);
}

/*
 * print out the information on all the QTLs, either to a file or stdout.
 */
void print_aqtl(params *theparams, char *outfile,int boot)
{
  int ii, k, kk, wo,r,s;
  aqtl *tptr;
  FILE *outf;
  if (*outfile == '-')
    outf = stdout;
  else {
    outf = fileopen(outfile, "a");
    if (outf == NULL)
      outf = fileopen(outfile, "w");
  }
  if (outf != NULL) {
    if ( boot == 0 ) {
      fprintf(outf,"\n#  The following section is formatted just like an Rqtl.out file.");
      fprintf(outf,"\n#  If this is a part of an Eqtl.out file, it is here just in case");
      fprintf(outf,"\n#  you want to use it as input for Rcross or Zmapqtl.   Otherwise,");
      fprintf(outf,"\n#  the information between the start and end tokens is identical to");
      fprintf(outf,"\n#  what has been written above.  If this is an Rqtl.out file,");
      fprintf(outf,"\n#  disregard this note.");
      fprintf(outf,"\n# ");
      fprintf(outf,"\n#  Below, those lines that start with a minus sign and an l (-l)");
      fprintf(outf,"\n#  define loci, while those that start with a minus sign and an i (-i)");
      fprintf(outf,"\n#  define epistatic interactions between the loci previously defined.");
      fprintf(outf,"\n#  The loci must be defined before the interactions can be defined.");
      fprintf(outf,"\n#  For some experimental designs, some interactions will be ignored.");
      fprintf(outf,"\n#   ");
      fprintf(outf,"\n-start model");
    }
    fprintf(outf, "\n#\n# Here is the data on the QTLS...\n#\n-t %5d   is the number of traits\n#\n#", theparams->themap->traits);
    wo = 0;
    for (kk = 1; kk <= theparams->themap->traits; kk++) {
      k =  theparams->themap->knum[kk];
      if ( theparams->themap->tnames != NULL )
        fprintf(outf, "\n-k %5d\n# for trait -named %s which is -number %5d\n#",k,theparams->themap->tnames[kk],  kk);
      else 
        fprintf(outf, "\n-k %5d\n# for trait -number %5d\n#", k, kk);
      if (boot == 0)
	    fprintf(outf, "\n#      #  ..Chrom..Markr.    .RecombiL.    .RecombiR.    .Additive.    .Dominance");
      else if ( boot == 1 )
	    fprintf(outf, "\n#      #  ..Chrom..Markr.    .Position.   .Test Stat.    .Additive.    .Dominance.      .R2.          .TR2.        .S.");
      else
	    fprintf(outf, "\n#      #  ..Chrom..Markr.    .Position.   .LOD Score.    .Additive.    .Dominance.      .R2.          .TR2.        .S.");
      for (ii = 1; ii <= k; ii++) {
	    wo = wo + 1;
/*	    if ( boot == 2 ) 
	      (qtlptr + wo)->c2 = lrtolod((qtlptr + wo)->c2);*/
	    if ( boot == 0 ) 
	      fprintf(outf, "\n-l %5d  %5d  %5d   %12.4f  %12.4f  %12.4f  %12.4f", ii, theparams->theqtls[wo].chrm, theparams->theqtls[wo].mrk, theparams->theqtls[wo].c1, theparams->theqtls[wo].c2, theparams->theqtls[wo].a, theparams->theqtls[wo].d);
        if ( boot > 0 ) {
	      fprintf(outf, "\n   %5d  %5d  %5d   %12.4f  %12.4f  %12.4f  %12.4f", ii, theparams->theqtls[wo].chrm, theparams->theqtls[wo].mrk, theparams->theqtls[wo].c1, theparams->theqtls[wo].c2, theparams->theqtls[wo].a, theparams->theqtls[wo].d);
	      fprintf(outf, "  %12.4f  %12.4f  %12.4f",  theparams->theqtls[wo].r2, theparams->theqtls[wo].tr2, theparams->theqtls[wo].s );
        }
        tptr = &theparams->theqtls[wo] ;
      }
      fprintf(outf, "\n#");
      fprintf(outf, "\n#   QTL1   QTL2   Type       Value  \n#");
      for ( r = 1; r < k ; r++ )
        for ( s = r+1; s <= k ; s++ ) {
          if ( tptr->episAADD != NULL && tptr->episAADD[s][r] != 0.0 )
	        fprintf(outf, "\n-i %5d  %5d    AA   %12.4f ",  r, s, tptr->episAADD[s][r] );    
          if ( tptr->episADDA != NULL && tptr->episADDA[s][r] != 0.0 )
	        fprintf(outf, "\n-i %5d  %5d    AD   %12.4f ",  r, s, tptr->episADDA[s][r] );    
          if ( tptr->episADDA != NULL && tptr->episADDA[r][s] != 0.0 )
	        fprintf(outf, "\n-i %5d  %5d    DA   %12.4f ",  r, s, tptr->episADDA[r][s] );    
          if ( tptr->episAADD != NULL && tptr->episAADD[r][s] != 0.0 )
	        fprintf(outf, "\n-i %5d  %5d    DD   %12.4f ",  r, s, tptr->episAADD[r][s] );         
        }      
      fprintf(outf, "\n#\n#");
    }
    fprintf(outf, "\n#End of this block of information\n#");
    if ( boot == 0 ) 
      fprintf(outf,"\n-end model\n#");
    if (*outfile != '-')
      fileclose(outfile, outf);
  }
}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Qdatain.c
------------------------------------------------------------------ */


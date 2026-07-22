/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Preplot.c) is part of QTL Cartographer
         
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


/*  Main driver for Preplot, along with subroutines to reformat results
to be plotted with GNUPLOT.  */

#include "Main.h"

void write_header_plt(params *theparams, char *filename, char *term);
void add_pause(params *theparams);
void add_comma(params *theparams);
int do_qfile(params *theparams, char *inputfile, int chrom, genome *gptr,int whichone);
void do_lrfile(params *theparams, int chrom, genome *gptr);
void do_zfile(params *theparams, int chrom, zresult *zresnode);
void do_siglevel(params *theparams, int chrom, genome *gptr);
FPN get_pos(genome *gptr, int row, int col);
void write_preamble(params *theparams, char *filename );
void do_markernames(FILE *outf, params *theparams );


int main(argc, argv)
int argc;
char *argv[];
{
  char  *chptr,*purpose ;
  int nopts,chrom,trait,thestep;
  params *theparams;
  int ii,automatic,error,datacntr;
  char *datafile;
  genome *gptr ;
  zresult *tzres,*zreslist;

#if defined(MACWARRIOR) 
 /*   Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  whichprogram = 9;
 /**/
  purpose = cvector(0,MAXNAME);
  datafile = cvector(0,MAXNAME);
  strcpy(purpose,"Reformat Files for Plotting by Gnuplot");
  theparams = NULL;
  nopts = 12;
  theparams = create_params(theparams, 1, NULL);
  chptr = asctime2();  
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
/*  if ( theparams->lodflag == 1 )
    theparams->siglevel = lodtolr(theparams->siglevel);*/
  for (ii = 0; ii < MAXNAME; ii++) 
    datafile[ii]  =   '\0';
  strcat(datafile, theparams->stem);
  strcat(theparams->stem, ".plt");

  write_header_plt(theparams,datafile,theparams->term); /* put in a header */
  
/*
  Now we have the basic information...Time to use it.  
  First, create the genome... 
*/
    GetTheMap(theparams, theparams->map);

    error = get_the_nn(&theparams->nn, &theparams->traits, theparams->ifile);
    theparams->themap->traits = theparams->traits;
	theparams->thegenome = create_genome(theparams->themap);
	gptr = theparams->thegenome->next;
	while (gptr != NULL) { /* Move distances up one node.  */
	  if (gptr->chrom == gptr->prev->chrom)
	    gptr->dist = gptr->dist + gptr->prev->dist;
	  gptr = gptr->next;
	}    
/* This code automagically processes all results... */
  if ( (ii=isfile(theparams->zfile))==1 ) {
    zreslist = zmapqtl_list(theparams,312);  
    trait = chrom = 0;
    if (theparams->verbosity == 1 )
      thestep = Rotator(0);
    for ( tzres = zreslist ; tzres != NULL ; tzres = tzres->next ) {
      if (theparams->verbosity == 1 )
        thestep = Rotator(thestep);
      theparams->whichtrait = tzres->trait;
      theparams->Model = tzres->model;
      if ( tzres->chrom != chrom || tzres->trait != trait) {
	    if ( tzres != zreslist ) 
          add_pause(theparams);
        chrom = theparams->wchrom = tzres->chrom;
        trait = theparams->whichtrait = tzres->trait;
        write_preamble(theparams,theparams->stem );
        do_siglevel(theparams,chrom,theparams->thegenome);
        add_comma(theparams); 
        if ( (ii=isfile(theparams->qtl))==1 ) {
          datacntr = do_qfile(theparams,theparams->qtl, chrom, theparams->thegenome,'q');	          
          if ( datacntr > 0 )
            add_comma(theparams);
        }
        if ( (ii=isfile(theparams->eqtl))==1 ) {
          datacntr = do_qfile(theparams, theparams->eqtl, chrom, theparams->thegenome, 'e');	          
          if ( datacntr > 0 )
            add_comma(theparams);
        }
        if ( (ii=isfile(theparams->mqtfile))==1 ) {
          datacntr = do_qfile(theparams,theparams->mqtfile, chrom, theparams->thegenome,'m');	          
          if ( datacntr > 0 )
            add_comma(theparams);
        }
        if ( (ii=isfile(theparams->lrfile))==1 ) {
          do_lrfile(theparams,  chrom,   theparams->thegenome);	          
          add_comma(theparams);
        }
      }
      else  
        add_comma(theparams);
      do_zfile(theparams,  chrom, tzres);	    
    }    
    add_pause(theparams);
    zresult_list_abolish(zreslist);
  }  
  if (theparams->verbosity == 1 )
      thestep = Rotator(1);
  write_trailer(theparams,chptr,0);

/* Time to clean up. */
  free_cvector( purpose,0,MAXNAME);
  free_cvector( datafile,0,MAXNAME);
  theparams = create_params(theparams, 0, NULL);

  return(0);
}

void write_header_plt(params *theparams, char *filename, char *term)
{
  FILE *fptr;
   if ( theparams->Inter == 1 )
    fptr = fileopen(theparams->stem, "a");
  else
    fptr = NULL;
  if (fptr == NULL) 
    fptr = fileopen(theparams->stem, "w");
  fprintf(fptr, "#\n# Gnuplot control file created by Preplot\n#\n");
  fprintf(fptr, "#The x axis is automatically scaled to the chromosome size.\n#\n");
  fprintf(fptr, "set autoscale x\n");
  fprintf(fptr, "#If you ran Eqtl, then the y axis is between 0.0 and the maximum LR.\n#\n");
  if ( theparams->maxlr > (FPN) 0.0 ) 
    fprintf(fptr,"set yrange [%f:%f]\n",-0.2*theparams->maxlr,theparams->maxlr);
  else
  fprintf(fptr, "set autoscale y\n");
/*  
  if ( theparams->maxlr > 500.0 )
    yy = 1000.0;
  else if ( theparams->maxlr > 400.0 )
    yy = 500.0;
  else if ( theparams->maxlr > 300.0 )
    yy = 400.0;
  else if ( theparams->maxlr > 200.0 )
    yy = 300.0;
  else if ( theparams->maxlr > 100.0 )
    yy = 200.0;
  else if ( theparams->maxlr > 50.0 )
    yy = 100.0;
  else  
    yy = 50.0;

  fprintf(fptr, "set yrange [%f:%f]\n",-yy,yy);
*/
  if (term[0] != '\0') 
    fprintf(fptr, "set term %s\n", term);
  if ( !strcmp(term,"postscript") ) {
    fprintf(fptr, "set output \"%s.ps\"\n",filename);
    theparams->perms = -1;
  }
  else if (  !strcmp(term,"hpljii") ) {
    fprintf(fptr, "set output \"%s.hp\"\n",filename);
    theparams->perms = -1;
  }
  fprintf(fptr, "set lmargin 8\nset ytics nomirror\nset border 2\nset noxtics\nset xzeroaxis");
  fileclose(theparams->stem, fptr);
}

void write_preamble(params *theparams, char *filename )
{
  FILE *fptr;
  int ii;
  /*
  int last,i,k ;
  k=0;
  last = (int) strlen(filename);
  for ( i = last ; i >= 0 && filename[i] != (char) FILESEP ; i-- ) k+=1;  
  i = i+1;*/
  fptr = fileopen(filename, "a");
  if (fptr == NULL) 
    fptr = fileopen(filename, "w");

  fprintf(fptr, "\nset title \"Results for Trait %d, Chromosome %d\" ",theparams->whichtrait,theparams->wchrom);
  fprintf(fptr, "\nset xlabel \"Position in Morgans\" ");
  if ( theparams->lodflag == 1 )
    fprintf(fptr, "\nset ylabel \"LOD Score\" ");
  else
    fprintf(fptr, "\nset ylabel \"LR Statistic\"");
  fprintf(fptr, "\n");
  if ( (ii=isfile(theparams->lrfile))==1 )  
          do_markernames(fptr,theparams);	          

  
  fprintf(fptr, "\nplot ");
  fileclose(filename, fptr);
}

void do_siglevel(params *theparams, int chrom, genome *gptr)
{
  char *datafile;
  FILE  *outf, *dataf;
  FPN xx,yy;
  genome *tgptr;
  datafile = cvector(0,MAXNAME);
    if ( theparams->workdir != NULL )
      sprintf(datafile, "%sc%dt%d.s", theparams->workdir, chrom,theparams->whichtrait);
    else
      sprintf(datafile, "c%dt%d.s", chrom,theparams->whichtrait);
  dataf = fileopen(datafile, "w");
  if ( theparams->lodflag == 1 )
    yy = lrtolod(theparams->siglevel);
  else
    yy = theparams->siglevel;
  xx = (FPN) 0.0;
  fprintf(dataf, "%f %f\n", xx, yy);
  for ( tgptr = gptr ; tgptr != NULL && tgptr->chrom <= chrom ; tgptr = tgptr->next )
    xx = tgptr->dist;
  fprintf(dataf, "%f %f\n", xx, yy);

  fileclose(datafile, dataf);
  outf = fileopen(theparams->stem, "a");
  if ( theparams->workdir != NULL ) 
    shift_fn(datafile);
  fprintf(outf, "\'%s\' with lines", datafile);
  fileclose(theparams->stem, outf);
  free_cvector(datafile,0,MAXNAME);
}


void do_zfile(params *theparams, int chrom, zresult *zresnode)
{
  int lrcol,acol,dcol,r2col,tr2col,scol;
  int  xcol, ycol, ii, start, chr, attrait, atmodel, themodel, thetrait, ch;
  char   *datafile, *input;
  FILE *infile, *outf, *dataf;
  FPN xx, yy;
 
  datafile = cvector(0,MAXNAME);
  assign_zfilecols(theparams,&lrcol,&acol,&dcol,&r2col,&tr2col,&scol) ;

  attrait = atmodel = 0;
  xcol = 3;
  ycol = lrcol;

  input = theparams->zfile;
  if ( theparams->workdir != NULL )
    sprintf(datafile, "%sc%dt%d.z%d", theparams->workdir, chrom,theparams->whichtrait,theparams->Model);
  else
    sprintf(datafile, "c%dt%d.z%d", chrom,theparams->whichtrait,theparams->Model);
  
  infile = fileopen(input, "r");
  if ( zresnode == NULL ) {
    start = 0;
    while (start == 0) {
      for (ii = 0; ii < MAXLINE; ii++)
        gbuffer[ii] = '\0';
      for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
        gbuffer[ii] = (char) ch;
      if (gbuffer[0] == '-')
        switch (gbuffer[1]) {
          case 't':
	        if (theparams->whichtrait != 0) {
	          get_field(2, gname, gbuffer);
	          thetrait = atoi(gname);
	          if (thetrait == theparams->whichtrait)
	            attrait = 1;
	          else {
	            attrait = 0;
	            if (theparams->Model != 0)
	              atmodel = 0;
	          }
	        }
	      break;
          case 'M':
	        if (theparams->Model != 0) {
	          get_field(2, gname, gbuffer);
	          themodel = atoi(gname);
	          if (themodel == theparams->Model)
	            atmodel = 1;
	          else {
	            atmodel = 0;
	            if (theparams->whichtrait != 0)
	              attrait = 0;
	          }
	        }
	      break;
          case 's':
	        if (attrait == 1 && atmodel == 1)
	          start = 1;
	        break;
          default:
	      break;
        }
    }
  }
  else {
    fseek(infile,zresnode->start_offset,SEEK_SET);
    start = 1;
  }
  dataf = fileopen(datafile, "w");

  while (start == 1) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if (gbuffer[0] == '-' && gbuffer[1] == 'e')
      start = 0;
    else {
      get_field(1, gname, gbuffer);
      chr = atoi(gname);
      if (chr == chrom) {
	    get_field(xcol, gname, gbuffer);
	    xx = (FPN) atof(gname);
	    get_field(ycol, gname, gbuffer);
        if ( !strncmp(gname,"NaN",3) || !strncmp(gname,"-NaN",4) ||  !strncmp(gname,"Inf",3) )
          yy = (FPN) 0.0;
        else
	      yy = (FPN) atof(gname);
	    if ( yy < (FPN) 0.0 )
	      yy = (FPN) 0.0;
	    if ( theparams->lodflag == 1 )
	      yy = lrtolod(yy);
	    fprintf(dataf, "%f %f\n", xx, yy);
      }
    }
    if (ch == EOF)
      start = 0;
  }

  fileclose(datafile, dataf);
  if ( theparams->workdir != NULL ) 
    shift_fn(datafile);
  outf = fileopen(theparams->stem, "a");
  fprintf(outf, "\'%s\' with lines", datafile);
  fileclose(theparams->stem, outf);
  if ( theparams->workdir != NULL ) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    insert_wd(gbuffer,theparams->workdir, datafile);
  }
  fileclose(input, infile);
  free_cvector(datafile,0,MAXNAME);
}

/*
  Put markername info into the qtlcart.plt file
  
  This will allow up to 20 marker names to be on the graph.  
  I put in the restriction because too many marker names clutter up the
  plots.   
  
  
*/
void do_markernames(FILE *outf, params *theparams )
{
  int negctr, row, col, ii, attrait, ch,chrom,nmarkers,factor,cnt;
  char  *datafile;
  FILE *infile ;
  genome *gptr;
  markermap *themap;
  FPN xx, yy;
  gptr = theparams->thegenome;
  themap = theparams->themap;
  chrom = theparams->wchrom;
  nmarkers = themap->mpc[chrom];

  factor = (int) ceil(  (FPN) nmarkers / 20.0 );
  cnt = 0;
  
  datafile = cvector(0,MAXNAME);

  negctr = 0;
  infile = fileopen(theparams->lrfile, "r");
  attrait = 0;
  while (attrait == 0) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if (gbuffer[0] == '-' && gbuffer[1] == 't') {
      get_field(2, gname, gbuffer);
      attrait = atoi(gname);
      if (attrait != theparams->whichtrait)
	    attrait = 0;
    }
  }

  

  while (negctr < 4) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if ((negctr == 3) && (*(gbuffer + 0) != '-')) {
      get_field(1, gname, gbuffer);
      row = atoi(gname);
      if (row == chrom) {
	    get_field(2, gname, gbuffer);
	    col = atoi(gname);
	    get_field(5, gname, gbuffer);
	    yy = (FPN) atof(gname);
	    if ( theparams->lodflag == 1 )
	      yy = lrtolod(yy);
	    if ( theparams->maxlr > (FPN) 0.0 )
	      yy =  (FPN) 0.0; /* - 0.2 * theparams->maxlr ; */
	    else
	      yy = (FPN) 1.1 * yy; 
	    xx = get_pos(gptr, row, col);
	    cnt +=1;
	    if ( cnt == factor ) {
	      if ( theparams->maxlr > (FPN) 0.0 )
	        fprintf(outf, "set label \"%s<%4.2f>\" at %f,%f right rotate\n",themap->names[themap->ttable[chrom][col]],xx, xx, yy);	    
	      else
	        fprintf(outf, "set label \"<%4.2f>%s\" at %f,%f left rotate\n",xx,themap->names[themap->ttable[chrom][col]], xx, yy);
          cnt = 0;
	    }
      }
    }
    if (gbuffer[0] == '-')
      negctr += 1;
    if (ch == EOF)
      negctr = 5;
  }
  fileclose(theparams->lrfile, infile);
}


void do_lrfile(params *theparams, int chrom, genome *gptr)
{
  int negctr, row, col, ii, attrait, ch;
  char  *datafile;
  FILE *infile, *outf, *dataf;
  FPN xx, yy;

  datafile = cvector(0,MAXNAME);

  if ( theparams->workdir != NULL )
    sprintf(datafile, "%sc%dt%d.lr", theparams->workdir, chrom,theparams->whichtrait);
  else
    sprintf(datafile, "c%dt%d.lr", chrom,theparams->whichtrait);
  negctr = 0;
  infile = fileopen(theparams->lrfile, "r");
  attrait = 0;
  while (attrait == 0) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if (gbuffer[0] == '-' && gbuffer[1] == 't') {
      get_field(2, gname, gbuffer);
      attrait = atoi(gname);
      if (attrait != theparams->whichtrait)
	    attrait = 0;
    }
  }
  dataf = fileopen(datafile, "w");


  while (negctr < 4) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if ((negctr == 3) && (*(gbuffer + 0) != '-')) {
      get_field(1, gname, gbuffer);
      row = atoi(gname);
      if (row == chrom) {
	    get_field(2, gname, gbuffer);
	    col = atoi(gname);
	    get_field(5, gname, gbuffer);
	    yy = (FPN) atof(gname);
	    if ( theparams->lodflag == 1 )
	      yy = lrtolod(yy);
	    xx = get_pos(gptr, row, col);
	    fprintf(dataf, "%f %f\n", xx, yy);
      }
    }
    if (gbuffer[0] == '-')
      negctr += 1;
    if (ch == EOF)
      negctr = 5;
  }
  fileclose(theparams->lrfile, infile);
  fileclose(datafile, dataf);
  if ( theparams->workdir != NULL ) 
      shift_fn(datafile);
  outf = fileopen(theparams->stem, "a");
  fprintf(outf, "\'%s\' with points pt 1", datafile);
  fileclose(theparams->stem, outf);
  free_cvector(datafile,0,MAXNAME);
}

/*
  This processes files in the format of Rqtl.out.   These files are produced by
  Rqtl, Eqtl and MImapqtl.   
*/
int do_qfile(params *theparams, char *inputfile, int chrom, genome *gptr,int whichone)
{
  int negctr, marker, chromosome, ii, attrait, qnum, ch,go_on,datacntr;
  char  *datafile;
  FILE *infile, *outf, *dataf;
  FPN xx, yy;

  datafile = cvector(0,MAXNAME);
  if ( theparams->workdir != NULL )
    sprintf(datafile, "%sc%dt%d.%c", theparams->workdir, chrom,theparams->whichtrait,whichone);
  else
    sprintf(datafile, "c%dt%d.%c", chrom,theparams->whichtrait,whichone);
  negctr = 0;
  infile = fileopen(inputfile, "r");
  ii = movetoqtls(infile);
  
    go_on = 1; 
    attrait = 0;
    do {
      ch = get_next_token(gname,MAXNAME,infile);
      if ( !strcmp(gname,"-k") ) {
        ch = get_next_token(gname,MAXNAME,infile);
        qnum = atoi(gname);
      }
      else if  ( !strcmp(gname,"-number") ){
        ch = get_next_token(gname,MAXNAME,infile);
        attrait = atoi(gname);
      }
      else if (ch == EOF ) {
        fileclose(inputfile, infile);
        free_cvector(datafile,0,MAXNAME);
        return(-1) ;
      }
      if (attrait == theparams->whichtrait)  
        go_on = 0;     
    } while (go_on == 1);
  dataf = fileopen(datafile, "w");

  datacntr = 0;
  while (negctr <= qnum) {
    for (ii = 0; ii < MAXLINE; ii++)
      gbuffer[ii] = '\0';
    for (ii = 0; ((ch = fgetc(infile)) != EOF) && (ch != '\n'); ii++)
      gbuffer[ii] = (char) ch;
    if (gbuffer[0] == '-' && gbuffer[1] == 'l') {
      negctr = negctr + 1;
      get_field(3, gname, gbuffer);
      chromosome = atoi(gname);
      if (chromosome == chrom) {
	    get_field(4, gname, gbuffer);
	    marker = atoi(gname);

	    get_field(5, gname, gbuffer);	/* get left recombination */
	    yy = (FPN) atof(gname);

	    if (yy != 0.5)
	      xx = mapfunc(yy, -1) + get_pos(gptr, chromosome, marker);
	    else {
	      get_field(6, gname, gbuffer);	/* get right recombination */
	      yy = (FPN) atof(gname);
	      xx = get_pos(gptr, chromosome, marker + 1) - mapfunc(yy, -1);
	    }
	    get_field(7, gname, gbuffer);	/* additive effect */
	    if ( theparams->maxlr > (FPN) 0.0 && theparams->maxeffect > (FPN) 0.0 )
	      yy = theparams->maxlr * (FPN) atof(gname) / theparams->maxeffect ;	    
	    else
	      yy = (FPN) 100.0 * (FPN) atof(gname);
	    fprintf(dataf, "%f %f\n", xx, yy);
	    datacntr +=1;
      }
    }
    if (ch == EOF)
      negctr = qnum + 1;
  }
  fileclose(inputfile, infile);
  fileclose(datafile, dataf);
    if ( theparams->workdir != NULL ) 
      shift_fn(datafile);
  if ( datacntr > 0 ) {
	  outf = fileopen(theparams->stem, "a");
	  if ( whichone == 'q' ) 
	    fprintf(outf, "\'%s\' with impulses", datafile);
	  else if (whichone == 'm' ) 
	    fprintf(outf, "\'%s\' with points pt 3", datafile);
	  else if (whichone == 'e' ) 
	    fprintf(outf, "\'%s\' with points pt 2", datafile);
	  fileclose(theparams->stem, outf);
  }
  free_cvector(datafile,0,MAXNAME);
  return(datacntr);
}

/* this gives the position of marker col on chromosome row*/
FPN get_pos(genome *gptr, int row, int col)
{
  genome *lgptr;
  FPN pos;
  pos = (FPN) 0.0;
  lgptr = gptr;
  while (lgptr != NULL)
    if (lgptr->chrom == row && lgptr->markr == col) {
      if (lgptr->prev != NULL && lgptr->prev->chrom == lgptr->chrom)
	pos = lgptr->prev->dist;
      else
	pos = (FPN) 0.0;
      return (pos);
    }
    else
      lgptr = lgptr->next;
  return (pos);
}


void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag)
{
  int ii,jj;
  if ( flag == 1 ) {
    strcpy(opt[1],  "-o"  );
    strcpy(opt[2],  "-e"  );
    strcpy(opt[3],  "-m"  );
    strcpy(opt[4],  "-q"  );
    strcpy(opt[5],  "-E"  );
    strcpy(opt[6],  "-M"  );
    strcpy(opt[7],  "-l"  );
    strcpy(opt[8],  "-z"  );
    strcpy(opt[9],  "-S"  );
    strcpy(opt[10],  "-T"  );
    strcpy(opt[11], "-H"  );
    strcpy(opt[12], "-L"  );
 

    strcpy(opt_e[1],  "Gnuplot Control File Name"  );
    strcpy(opt_e[2],  "Error File"  );
    strcpy(opt_e[3],  "Genetic Linkage Map File"  );
    strcpy(opt_e[4],  "Rqtl.out  file"  );
    strcpy(opt_e[5],  "Eqtl.out  file"  );
    strcpy(opt_e[6],  "MImapqtl.mqt file"  );
    strcpy(opt_e[7],  "LRmapqtl Output File"  );
    strcpy(opt_e[8],  "Zmapqtl Output File"  );
    strcpy(opt_e[9],  "Significance Threshold"  );
    strcpy(opt_e[10],  "Terminal"  );
    strcpy(opt_e[11],  "Hypothesis (30,31,32) for Fx design"  );
    strcpy(opt_e[12],  "Output LOD scores? (0=no,1=yes)"  );
  }
  for ( ii = 1 ; ii <= nopts ; ii++ ) 
    for ( jj = 0 ; jj <= MAXNAME ; jj++ )
      opt_v[ii][jj] = '\0';

    strcpy(opt_v[1],  theparams->stem  );
    strcpy(opt_v[2],  theparams->error  );
    strcpy(opt_v[3],  theparams->map  );
    strcpy(opt_v[4],  theparams->qtl  );
    strcpy(opt_v[5],  theparams->eqtl  );
    strcpy(opt_v[6],  theparams->mimfile  );
    strcpy(opt_v[7],  theparams->lrfile );
    strcpy(opt_v[8],  theparams->zfile  );

    sprintf(opt_v[9],"%f",theparams->siglevel );
    strcpy(opt_v[10], theparams->term );
    sprintf(opt_v[11],"%d",theparams->ihypo );
    sprintf(opt_v[12],"%d",theparams->lodflag );



}  

void update_params(char **opt_v,  params *theparams)
{

    strcpy(theparams->stem, opt_v[1] );
    strcpy(theparams->error, opt_v[2] );
    strcpy(theparams->map, opt_v[3]  );
    strcpy(theparams->qtl, opt_v[4]  );
    strcpy(theparams->eqtl, opt_v[5]  );
    strcpy(theparams->mimfile, opt_v[6]  );
    strcpy(theparams->lrfile, opt_v[7]  );
    strcpy(theparams->zfile, opt_v[8] );
    theparams->siglevel = (FPN) atof(opt_v[9]);
    strcpy(theparams->term,opt_v[10]);
    theparams->ihypo = atoi( opt_v[11] );
    theparams->lodflag = atoi( opt_v[12] );
}
  
void add_comma(params *theparams)
{	  
    FILE *outf;     
    outf = fileopen(theparams->stem, "a");
    fprintf(outf, ", ");
    fileclose(theparams->stem, outf);
}

void add_pause(params *theparams)
{
  FILE *outf;
  outf = fileopen(theparams->stem, "a");
  if (theparams->perms != -1 )
    fprintf(outf, "\npause -1 \"Hit return to continue\"\nset nolabel\n");
  else
    fprintf(outf,"\n#\n");
  fileclose(theparams->stem, outf);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Preplot.c
------------------------------------------------------------------ */


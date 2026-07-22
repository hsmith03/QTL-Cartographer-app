/* ------------------------------------------------------ XCutXCodeXSkip
     This file (RCmain.c) is part of QTL Cartographer
         
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
/*    Main function  for Rcross.   */

void print_format_data(params *theparams, individual *indptr, int n,  char *outfile);


int main(int argc, char *argv[])  {
  FILE  *outf;
  char *chptr, *purpose;
  int nopts, j;
  params *theparams;
  int jj, ii, automatic, newmap;
  int error,   iline;
  individual *p1p2f1;
  thelines *lines, *lptr, *p1ptr, *p2ptr, *f1ptr;

#if defined(MACWARRIOR)  
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and  output. */
  argc = ccommand(&argv);
#endif
  whichprogram = 3;
  lines = NULL;
  purpose = cvector(0,MAXNAME);
  strcpy(purpose, "Create or translate a data set");
  theparams = NULL;
  automatic = 0;
  nopts = 13;
  theparams = create_params(theparams, 1, NULL);
  theparams->theqtls = NULL;
  theparams->gout = 0;
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);

  print_head(argv[0],theparams->ifile,chptr,0,(30+theparams->gout),theparams);

/* Initailize the data structures...*/
  if ( isfile(theparams->map)  )
    GetTheMap(theparams, theparams->map); 
  else
    theparams->themap = NULL;

  if (theparams->iinfile[0] == '\0') {
	GetTheModel(theparams, theparams->qtl );

    calc_parental_diffs(theparams->themap,theparams->theqtls);

    error = check_params(theparams,theparams->themap,whichprogram);
    create_bogus_markers(theparams,theparams->themap);
    p1p2f1 = indvector(4, theparams->themap, theparams->theqtls);



/* Create initial populations P1, P2 and the F1 */
    init_pop(p1p2f1,theparams);
    calc_phenotypes(p1p2f1, 4, theparams);

    theparams->thegenome = create_genome(theparams->themap);
    place_qtls(theparams, theparams->thegenome);

    lines = alloc_aline("Parental Group 1 (P1)", &p1p2f1[1], 1, NULL, NULL, 1);
    p1ptr = lptr = lines;
    lptr->next = alloc_aline("Parental Group 2 (P2)", &p1p2f1[2], 1, lptr, NULL, 2);
    p2ptr = lptr = lptr->next;
    lptr->next = alloc_aline("F1 population", &p1p2f1[3], 2, lptr, NULL, 3);
    f1ptr = lptr = lptr->next;
    lptr->next = alloc_aline("", NULL, theparams->nn, lptr, NULL, 4);
    lptr = lptr->next;
    lptr->iptr = indvector(lptr->nn, theparams->themap, theparams->theqtls);
    lptr->iptr = lptr->iptr + 1;

/* Create a new generation of nn individuals...*/
    switch (theparams->cross) {
     default:
     case 1:	/* Backcross to the high line... */
      strcat(lptr->name, "F1 x P1 backcross");
      do_a_cross(theparams,p1ptr, f1ptr, lptr, theparams->thegenome ,0);
      break;
     case 2:	/* Backcross to the low line... */
      strcat(lptr->name, "F1 x P2 backcross");
      do_a_cross(theparams,p2ptr, f1ptr, lptr, theparams->thegenome ,0);
      break;
     case 3:	/* self F1 x F1 to produce an SF2... */
      strcat(lptr->name, "F2 (F1 x F1)  selfed");
      do_a_cross(theparams,f1ptr, f1ptr, lptr, theparams->thegenome ,1);
      break;
     case 4:	/* F1 x F1 to produce an RF2... */
      strcat(lptr->name, "F2 (F1 x F1) intercross");
      do_a_cross(theparams,f1ptr, f1ptr, lptr, theparams->thegenome ,0);
      break;
     case 5:	/* Recombinant inbred lines */
      sprintf(lptr->name, "Ri%d lines",theparams->crosst);
      if ( theparams->crosst > 0 ) {
        expand_map(theparams,theparams->thegenome,theparams->theqtls,theparams->themap);
        theparams->tcross = 0;
        theparams->tcrosst = 0;        
      }
      if ( theparams->crosst == 2 )
        do_a_cross(theparams,f1ptr, f1ptr, lptr, theparams->thegenome ,0);
      else
        do_a_cross(theparams,f1ptr, f1ptr, lptr, theparams->thegenome ,1);
      break;
     case 6:	/* Sue Carson Cross */
      strcat(lptr->name, "Sue Carson Cross");
      do_a_cross(theparams,p1ptr, p2ptr, lptr, theparams->thegenome ,0);
      break;
    }
    calc_phenotypes(lptr->iptr - 1, lptr->nn, theparams);
    
    /* lptr now points to the base cross:  B1, B2, SF2, RF2 or Ri.   Now we complete the cross. */
    if ( theparams->crosst > 2   )  /*  This means that the base cross is SFx or RFx, where x > 2.   Cross to the x'th generation. */
      for ( ii = 3 ; ii<= theparams->crosst ; ii++ ) {
        p1ptr = lptr;
      	lptr->next = alloc_aline("", NULL, theparams->nn, lptr, NULL, lptr->which + 1);
	    lptr = lptr->next;
	    lptr->iptr = indvector(lptr->nn, theparams->themap, theparams->theqtls);
	    lptr->iptr = lptr->iptr + 1;
	    if ( theparams->cross == 3 ) {
	      sprintf(lptr->name,"SF%d line",ii);
	      do_a_cross(theparams,p1ptr, p1ptr, lptr, theparams->thegenome ,1);
	    }
	    else if ( theparams->cross == 4 ) {
	      sprintf(lptr->name,"RF%d line",ii);
	      do_a_cross(theparams,p1ptr, p1ptr, lptr, theparams->thegenome ,0);
	    }
	    calc_phenotypes( lptr->iptr - 1, lptr->nn, theparams);
      }  
/*  Now, we ask whether this is a test cross.   0 means no, greater than 0 means yes.   

    We'll have 
    
      1.  T(B1)SFx,   T(B1)RFx
      2.  T(B2)SFx,   T(B2)RFx
      3.  T(SFy)SFx   y must be x+1
      4.  T(D3)SFx    x > 1


*/      
    if ( theparams->tcross != 0 ) {
      switch (theparams->tcross) {
        case 1:
          p1ptr = which_line(1, lines); break;
        case 2:
          p1ptr = which_line(2, lines); break;
        case 3:      
        case 0:
        default: 
          p1ptr = lptr;   break;
      }  
      p2ptr = lptr;
      lptr->next = alloc_aline("", NULL, theparams->nn, lptr, NULL, lptr->which + 1);
	  lptr = lptr->next;
	  lptr->iptr = indvector(lptr->nn, theparams->themap, theparams->theqtls);
	  lptr->iptr = lptr->iptr + 1;
/*      if ( theparams->tcross == 3 )*/
        sprintf(lptr->name,"%s line",theparams->thecross);

	  if ( theparams->tcross == 12 ) {
	    p1ptr = which_line(1,lines);
	    do_a_cross(theparams,p1ptr, p2ptr, lptr, theparams->thegenome ,-1);
	    p1ptr = which_line(2,lines);
	    do_a_cross(theparams,p1ptr, p2ptr, lptr, theparams->thegenome ,-2);	  
	  }
	  else 
	    do_a_cross(theparams,p1ptr, p2ptr, lptr, theparams->thegenome ,0);
	  calc_phenotypes(lptr->iptr - 1, lptr->nn, theparams);
	    
    }
    if ( theparams->cross == 6 ) { 
      if ( theparams->tcross == 1 )
        f1ptr = lines;
      else
        f1ptr = lines->next;
      outf = fileopen(theparams->ifile, "w");
      if (outf != NULL) {
        if (theparams->iinfile[0] == '\0')
          fprintf(outf, "# %12ld \n\t%s\n\tThis output file (%s) is for %s...\n", theparams->seed, VERSION,theparams->ifile, argv[0]);
        else 
          fprintf(outf, "# %12ld \n\t%s\n\tThis output file (%s) is for %s...\n", theparams->seed, VERSION,theparams->ifile, argv[0]);
        fprintf(outf, "\tIt is %s\n", chptr);
        jj=0;
        if ( theparams->tcross == 1 )
          fprintf(outf,"#\n#\n#Paternally derived markers\n#");
        else
          fprintf(outf,"#\n#\n#Maternally derived markers\n#");
        for (ii = 1; ii <= theparams->themap->traits; ii++) {
          fprintf(outf,"\n\n# Trait %2d\n-Values Chrom Mark     C1       C2      Add.      Dom.  ",ii);
          for (j = 1; j <= theparams->themap->knum[ii]; j++) {
            jj = jj+1;
            fprintf(outf,"\n  %2d  ",f1ptr->iptr->vqtls[ii][j]);
            fprintf(outf,"  %3d  %3d  %8.4f %8.4f %8.4f %8.4f", f1ptr->iptr->qtls[jj].chrm ,f1ptr->iptr->qtls[jj].mrk,  f1ptr->iptr->qtls[jj].c1,  f1ptr->iptr->qtls[jj].c2, f1ptr->iptr->qtls[jj].a,  f1ptr->iptr->qtls[jj].d )	;
          }
        }
        fileclose(theparams->ifile, outf);
      }
      print_format_data(theparams,lptr->iptr - 1, lptr->nn, theparams->ifile);
    }
    else {
      if (theparams->Inter == 0) 
        print_format_data(theparams,lptr->iptr - 1, lptr->nn,  theparams->ifile);
   }
    while (theparams->Inter != 0) {
      printf("\n\n\tWhich of these should be the first line? ");
      show_lines(lines);
      iline = get_int();
      p1ptr = which_line(iline, lines);
      if (p1ptr == NULL)
	    theparams->Inter = 0;
      else {
	    printf("\n\n\tWhich of these should be the second line? ");
	    show_lines(lines);
	    iline = get_int();
	    if ( iline < 0 )
	      p2ptr = p1ptr;
	    else
	      p2ptr = which_line(iline, lines);
	    if (p2ptr == NULL)
	      theparams->Inter = 0;
	    else {
	      printf("\n\n\tHow many offspring from this cross?\n  ");
	      theparams->nn = get_int();
	      if (theparams->nn > 0) {
	        lptr->next = alloc_aline("", NULL, theparams->nn, lptr, NULL, lptr->which + 1);
	        lptr = lptr->next;
	        lptr->iptr = indvector(lptr->nn, theparams->themap, theparams->theqtls);
	        lptr->iptr = lptr->iptr + 1;
	        printf("\n\n\tWhat would you like to call the offspring line?\n ");
	        jj = myfgets(lptr->name, MAXNAME, stdin);
	        if ( iline < 0 ) 
	          do_a_cross(theparams,p1ptr, p2ptr, lptr, theparams->thegenome ,1);
	        else
	          do_a_cross(theparams,p1ptr, p2ptr, lptr, theparams->thegenome ,0);
	       calc_phenotypes(lptr->iptr - 1, lptr->nn, theparams);
	     }
	     else
	       theparams->Inter = 0;
	  }
    }
    printf("\n\n\tSelect a line to print out (or 0 to continue doing crosses)... ");
    show_lines(lines);
    iline = get_int();
    if (iline > 0) {
      lptr = which_line(iline, lines);
      print_format_data(theparams,lptr->iptr - 1, lptr->nn,  theparams->ifile);
    }
    else if (iline < 0)
	    theparams->Inter = 0;
    }

  }

  if ( theparams->iinfile[0] != '\0') {  /* If an input file is specified, translate it... */
    if ( theparams->themap == NULL )
      newmap = 1;
    else
      newmap = 0;
	GetTheData(theparams, theparams->iinfile, 0 );
	if ( newmap == 1 ) {
      print_head(argv[0],theparams->map,chptr,0,10,theparams);
      print_map(theparams->themap, theparams->map);
	}
	print_format_data(theparams,theparams->thedata, theparams->nn,   theparams->ifile);
  }
  
  
  
  
  
  theparams->Inter = 0;
  theparams->reps = (int) REPS;
  theparams->gout = 0; 
  write_trailer(theparams,chptr,1);


/* Clean up...*/
  if (lines != NULL) {
    lptr = lines;	/* P1 */
    lptr->nn = 4;
    lines = lines->next;	/* P2 */
    lines = lines->next;	/* F1 */
    lines = lines->next;	/* subsequent lines */
    free_aline(lptr);
    lptr = lines;
    while (lptr->next != NULL) {
      lptr = lptr->next;
      free_aline(lptr->prev);
    }
	 if (lptr != NULL ) {
	   free_aline(lptr);
	 }
  }

  free_cvector(purpose,0,MAXNAME);
  theparams = create_params(theparams, 0, NULL);

  return(0);
}

void print_format_data(params *theparams, individual *indptr, int n,  char *outfile) {
  if ( theparams->gout == 0 )
      print_individuals(theparams,indptr,  n,  outfile);
  else if ( theparams->gout == 1 )
      print_individuals_std(theparams,indptr,  n,  outfile);
  else if ( theparams->gout == 2 || theparams->gout == 6)
      print_individuals_mm(theparams,indptr,  n, outfile);
  else if ( theparams->gout == 3 )
	  print_individuals_R(theparams,indptr,  n,  outfile);
  else if ( theparams->gout == 4 )
	  print_individuals_SAS(theparams,indptr,  n,  outfile);
  else if ( theparams->gout == 5 )
	  print_plabqtl(theparams,indptr,  n,   outfile);
  else if ( theparams->gout == 7 )
      print_individuals_mcd( theparams,indptr,  n,   outfile);
  else if ( theparams->gout == 8 )
	  print_individuals_SAS(theparams,indptr,  n,  outfile);
}

void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag)
{
  int ii, jj;
  if (flag == 1) {
    strcpy(opt[1], "-i");
    strcpy(opt[2], "-o");
    strcpy(opt[3], "-e");
    strcpy(opt[4], "-m");
    strcpy(opt[5], "-q");
    strcpy(opt[6], "-s");
    strcpy(opt[7], "-g");
    strcpy(opt[8], "-c");
    strcpy(opt[9], "-H");
    strcpy(opt[10], "-r");
    strcpy(opt[11], "-I");
    strcpy(opt[12], "-E");
    strcpy(opt[13], "-n");


    strcpy(opt_e[1], "Input File");
    strcpy(opt_e[2], "Output File");
    strcpy(opt_e[3], "Error File");
    strcpy(opt_e[4], "Genetic Linkage Map File");
    strcpy(opt_e[5], "QTL Data File");
    strcpy(opt_e[6], "Random Number Seed");
    strcpy(opt_e[7], "Output format [0,1,2,...,7]");
    strcpy(opt_e[8], "Cross (1,2,3) => (B1,B2,F2)");
    strcpy(opt_e[9], "Heritability");
    strcpy(opt_e[10], "Replications (Not yet active)");
    strcpy(opt_e[11], "Interactive Crosses? (0,1) => (no,yes)");
    strcpy(opt_e[12], "Environmental Variance (used if > 0)");
    strcpy(opt_e[13], "Sample Size");
  }
  for (ii = 1; ii <= nopts; ii++)
    for (jj = 0; jj <= MAXNAME; jj++)
      opt_v[ii][jj] = '\0';

  strcpy(opt_v[1], theparams->iinfile);
  strcpy(opt_v[2], theparams->ifile);
  strcpy(opt_v[3], theparams->error);
  strcpy(opt_v[4], theparams->map);
  strcpy(opt_v[5], theparams->qtl);
  sprintf(opt_v[6], "%ld", theparams->seed);
  sprintf(opt_v[7], "%d", theparams->gout);
  strcpy(opt_v[8],   theparams->thecross);
  sprintf(opt_v[9], "%f", theparams->Herit);
  sprintf(opt_v[10], "%d", theparams->reps);
  sprintf(opt_v[11], "%d", theparams->Inter);
  sprintf(opt_v[12], "%f", theparams->Environ);
  sprintf(opt_v[13], "%d", theparams->nn);


}

void update_params(char **opt_v,  params *theparams)
{

  strcpy(theparams->iinfile, opt_v[1]);
  strcpy(theparams->ifile, opt_v[2]);
  strcpy(theparams->error, opt_v[3]);
  strcpy(theparams->map, opt_v[4]);
  strcpy(theparams->qtl, opt_v[5]);
  theparams->seed = atol(opt_v[6]);
  theparams->gout = atoi(opt_v[7]);
  theparams->cross = parse_cross(theparams,opt_v[8]);
  theparams->Herit = (FPN) atof(opt_v[9]);
  theparams->reps = atoi(opt_v[10]);
  theparams->Inter = atoi(opt_v[11]);
  theparams->Environ = (FPN) atof(opt_v[12]);
  theparams->nn = atoi(opt_v[13]);

}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file RCmain.c
------------------------------------------------------------------ */


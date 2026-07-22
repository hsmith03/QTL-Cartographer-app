/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Prune.c) is part of QTL Cartographer
         
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


/* Main driver and subroutines for Prune,  which modifies data sets.  */

#include "Main.h"


void pruneotraits(params *theparams, individual *individs, markermap *themap);
void prunetraits(params *theparams, individual *individs, markermap *themap);
void prune2onetrait(params *theparams, individual *individs, markermap *themap);
void check_missmark(markermap *themap, individual *individs, int nn, FPN missmark);
void check_misstrait(markermap *themap, individual *individs, params *theparams);
int do_int_options(  params *theparams, individual *individs, markermap *themap);
void elim_marker(markermap *themap, individual *individs, int nn, int cc, int mm);
void elim_trait(markermap *themap, individual *individs, int nn, int trait);
void sim_missing_markers(params *theparams, individual *individs, markermap *themap);
void sim_dominant_markers(params *theparams, individual *individs, markermap *themap);
void do_a_bootstrap(params *theparams, individual *individs, markermap *themap);
void do_a_permutation(params *theparams, individual *individs, markermap *themap);
void sim_selective_genotyping(params *theparams, individual *individs, markermap *themap);
void indxsort(int n, FPN *ra, int *rb);
void collapsemap(params *theparams, individual *individs, markermap *themap) ;
void ImportMarkerInfo(int nn, int chrom, individual  *newdata,int nmark, individual *olddata, int omark);

int main(argc, argv)
int argc;
char *argv[];
{
  char  *chptr,*purpose;
  int nopts,printmap;
  params *theparams;
  int ii,automatic ;
  int  cc, mm;
  int   savetraits;


#if defined(MACWARRIOR)  
 /* Simulate the UNIX shell with ccommand, this should also
    be the place to specify redirection of input and
    output. */
  argc = ccommand(&argv);
#endif
  whichprogram = 4;
  purpose = cvector(0,MAXNAME);
  strcpy(purpose,"Prune, resample or shuffle data");
  theparams = NULL;
  nopts = 9;
  theparams = create_params(theparams, 1, NULL);
  theparams->p1 = 0.0;
  theparams->boot = 0;
  theparams->Inter = 1;
  theparams->whichtrait = 0;
  chptr = asctime2();
  automatic = process_arguments(argc,argv,chptr,purpose,nopts,theparams);
  
  if ( theparams->boot == 6 && ( theparams->whichtrait < 1 || theparams->whichtrait > theparams->traits ) )
    theparams->whichtrait = 1;
  if ( theparams->boot > 0 )
    theparams->Inter = 0;
/* Initailize the data structures...
  knum = 0;
  theparams->themap = get_markermap(knum, theparams->map);*/
  GetTheMap(theparams, theparams->map );
  GetTheData(theparams, theparams->ifile, 0 );

  savetraits = theparams->themap->traits;

  printmap = 0;
  if ( theparams->Inter == 1 ) 
    printmap = do_int_options( theparams,theparams->thedata,theparams->themap);
  else {
    if ( theparams->boot == 1 )  /*  Bootstrap the data */
      do_a_bootstrap(theparams,theparams->thedata,theparams->themap);  
    else if ( theparams->boot == 2 || theparams->boot == 12 || theparams->boot == 6)  /* Permute the data */ 
      do_a_permutation(theparams,theparams->thedata,theparams->themap);  
    else if ( theparams->boot == 3 )  /*  Simulate missing markers */
      sim_missing_markers(theparams,theparams->thedata,theparams->themap);    
    else if ( theparams->boot == 4 ) /*  Simulate dominant markers */
      sim_dominant_markers(theparams,theparams->thedata,theparams->themap);    
    else if ( theparams->boot == 5 )/*  Simulate selective genotyping */
      sim_selective_genotyping(theparams,theparams->thedata,theparams->themap);    
    else if ( theparams->boot == 7 )/*  Use only one trait */
      prune2onetrait(theparams,theparams->thedata,theparams->themap);    
    else if ( theparams->boot == 8 )/*  use only selected traits */
      prunetraits(theparams,theparams->thedata,theparams->themap);    
    else if ( theparams->boot == 9 )/*  use only selected traits */
      pruneotraits(theparams,theparams->thedata,theparams->themap);    
    else if ( theparams->boot == 10 ) {/* collapse map and data */
      printmap = 1;
      collapsemap(theparams,theparams->thedata,theparams->themap);    
    }
  }


  for (ii = 0; ii < MAXNAME; ii++)
    theparams->map[ii] = theparams->ifile[ii] = '\0';
  strcpy(theparams->map, theparams->stem);
  strcpy(theparams->ifile, theparams->stem);

  strcat(theparams->map, ".mpb");
  if ( theparams->boot == 7 ) 
    sprintf(theparams->ifile, "%s%d.crb",theparams->stem,theparams->whichtrait);    
  else
    strcat(theparams->ifile, ".crb");
  if ( printmap == 1 ) {
      print_head(argv[0],theparams->map,chptr,0,10,theparams);
	  print_map(theparams->themap, theparams->map);
  }
  if ( theparams->boot == 7 )
    theparams->themap->traits = 1;
  for (ii = 1; ii <= theparams->nn; ii++)
    theparams->thedata[ii].t = theparams->themap->traits;
  cc = 0;
  for (ii = 1; ii <= theparams->themap->m; ii++)
    cc = cc + theparams->themap->mpc[ii];
  mm = theparams->themap->ml;
  theparams->themap->ml = cc;
  if ( printmap != -1 ) {
      print_head(argv[0],theparams->ifile,chptr,0,30,theparams);
      print_individuals(theparams,theparams->thedata, theparams->nn,   theparams->ifile);
  }
  for (ii = 1; ii <= theparams->nn; ii++)
    theparams->thedata[ii].t = savetraits;    
  theparams->themap->traits = savetraits;
  theparams->themap->ml = mm;
  write_trailer(theparams,chptr,0);
/* Clean up...*/
  free_cvector( purpose,0,MAXNAME);
  theparams = create_params(theparams, 0, NULL);
  return(0);
}

/*
  For each individual, if the proportion of missing marker data
  is greater than some level (missmark, in %), then eliminate that
  individual from the data set.

*/
void check_missmark(markermap *themap, individual *individs, int nn, FPN missmark)
{
  int ii, jj, kk, ismiss;
  FPN misspercent;
  misspercent = missmark / (FPN)100.0;
  for (ii = 1; ii <= nn; ii++) {
    ismiss = 0;
    for (jj = 1; jj <= themap->m; jj++) {
      for (kk = 1; kk <= themap->mpc[jj] ; kk++)
	    if (individs[ii].markers[jj][kk] < 0)
	      ismiss = ismiss + 1;
    }
    if ((FPN) ismiss / (FPN) themap->ml >= misspercent)
	  individs[ii].print_flag = 'n';
  }
}

/*
   If an individual has no trait values, get rid of it.
*/
void check_misstrait(markermap *themap, individual *individs, params *theparams)
{
  int ii,nn;
  ii = themap->traits;
  nn = theparams->nn;
  for (ii = 1; ii <= nn; ii++) 
    if ( individs[ii].y[theparams->whichtrait] < (FPN) MISS_VAL)
      individs[ii].print_flag = 'n';
}


/*
  Remove marker mm on chromosome cc, and remove the data from
  the inivids array.
*/
void elim_marker(markermap *themap, individual *individs, int nn, int cc, int mm)
{
  int ii, jj;

  themap->mrf[cc][mm-1] = mapfunc(themap->mrf[cc][mm-1], 2) + mapfunc(themap->mrf[cc][mm], 2);
  themap->mrf[cc][mm-1] = mapfunc(themap->mrf[cc][mm-1], -2);

  for (ii = mm; ii <= themap->mpc[cc] - 1 ; ii++) {
    themap->mrf[cc][ii] = themap->mrf[cc][ii+1];
    themap->types[cc][ii] = themap->types[cc][ii+1];

    if (themap->ttable != NULL)
      themap->ttable[cc][ii] = themap->ttable[cc][ii+1];

    for (jj = 1; jj <= nn; jj++)
      individs[jj].markers[cc][ii] = individs[jj].markers[cc][ii+1];
  }
  themap->mpc[cc] = themap->mpc[cc] - 1;
}

/*
  Eliminate the specified trait.
*/
void elim_trait(markermap *themap, individual *individs, int nn, int trait)
{
  int ii, tt;
  if (  themap->tnames != NULL )
    printf("\nEliminating trait %d  named |%s| ",trait,themap->tnames[trait]);
  for (ii = 1; ii <= nn; ii++)
    for (tt = trait; tt <= themap->traits - 1; tt++)  
      individs[ii].y[tt] = individs[ii].y[tt+1];

  if ( themap->tnames != NULL )
    for (tt = trait; tt <= themap->traits - 1; tt++) 
        strcpy( themap->tnames[tt], themap->tnames[tt+1] );

  themap->traits = themap->traits - 1;
}

void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag)
{
  int ii,jj;
  if ( flag == 1 ) {
    strcpy(opt[1],  "-o"  );
    strcpy(opt[2],  "-e"  );
    strcpy(opt[3],  "-m"  );
    strcpy(opt[4],  "-i"  );
    strcpy(opt[5],  "-s"  );
    strcpy(opt[6],  "-I"  );
    strcpy(opt[7],  "-b"  );
    strcpy(opt[8],  "-M"  );
    strcpy(opt[9],  "-t"  );
 
    strcpy(opt_e[1],  "Output Filename Stem"  );
    strcpy(opt_e[2],  "Error File"  );
    strcpy(opt_e[3],  "Genetic Linkage Map File"  );
    strcpy(opt_e[4],  "Data File"  );
    strcpy(opt_e[5],  "Random Number Seed"  );
    strcpy(opt_e[6],  "Interactive mode (0,1)=>(no,yes)"  );
    strcpy(opt_e[7],  "B (1), P (2), M (3) or D (4)"  );
    strcpy(opt_e[8],  "Percent missing data to simulate"  );
    strcpy(opt_e[9],  "Trait to process"  );
  }
  for ( ii = 1 ; ii <= nopts ; ii++ ) 
    for ( jj = 0 ; jj <= MAXNAME ; jj++ )
      opt_v[ii][jj] = '\0';

    strcpy(opt_v[1],  theparams->stem  );
    strcpy(opt_v[2],  theparams->error  );
    strcpy(opt_v[3],  theparams->map  );
    strcpy(opt_v[4],  theparams->ifile  );
    sprintf(opt_v[5],"%ld",theparams->seed );
    sprintf(opt_v[6],"%d",theparams->Inter );
    sprintf(opt_v[7],"%d",theparams->boot );
    sprintf(opt_v[8],"%f",theparams->p1 );
    sprintf(opt_v[9],"%d",theparams->whichtrait );
}  

void update_params(char **opt_v,  params *theparams)
{
  strcpy(theparams->stem, opt_v[1] );
  strcpy(theparams->error, opt_v[2] );
  strcpy(theparams->map, opt_v[3]  );
  strcpy(theparams->ifile, opt_v[4]  );
  theparams->seed = atol(opt_v[5]);
  theparams->Inter = atoi(opt_v[6]);
  theparams->boot = atoi( opt_v[7] );
  theparams->p1 = (FPN) atof( opt_v[8] );
  theparams->whichtrait = atoi( opt_v[9] );
}
    
/*  
  This is an interactive way of eliminating traits, markers,individuals, etc. 
*/
int do_int_options( params *theparams, individual *individs, markermap *themap)
{
  int cc,mm,printmap;
  int go_on,ans,simless;
  printmap = 0;
  simless = 0;
  go_on = 1;
  fprintf(stdout,"\n\nYou can loop through items 1-6, but 7-11 terminate.\n\n");
  while ( go_on == 1 ) {
	  fprintf(stdout,"\n No.                Action  ");
      fprintf(stdout,"\n  1. Eliminate A- marker systems (P1 dominant)");
      fprintf(stdout,"\n  2. Eliminate a- marker systems (P2 dominant)");
      fprintf(stdout,"\n  3. Eliminate marker m on chromosome c");
      fprintf(stdout,"\n  4. Eliminate trait t");
      fprintf(stdout,"\n  5. Eliminate individuals with missing phenotypes for trait t");
      fprintf(stdout,"\n  6. Eliminate individuals with more than m%% missing markers");
      fprintf(stdout,"\n  7. Bootstrap the data");
      fprintf(stdout,"\n  8. Permute the traits in the data");
      fprintf(stdout,"\n  9. Simulate m%% missing markers");
      fprintf(stdout,"\n 10. Simulate m%% dominant markers");
      fprintf(stdout,"\n 11. Simulate m%% selective genotyping");
      fprintf(stdout,"\n 12. Write modified dataset and exit");
      fprintf(stdout,"\n 13. Exit without writing anything");
      fprintf(stdout,"\n\n\tPick a number to do an action... ");
      ans = get_int();
      if ( ans < 7 )
        simless = 1;
      switch (ans ) {
        default :  break;
        case 1: /* Get rid of markers that are dominant A- */
          printmap = 1;
	      for (cc = 1; cc <= themap->m; cc++) {
	        for (mm = themap->mpc[cc]; mm > 0; mm--)
	          if (themap->types[cc][mm] == 1 || themap->types[cc][mm] == -99)
	            elim_marker(themap, individs, theparams->nn, cc, mm);
	      }
	      break;
        case 2: /* Get rid of markers that are dominant a- */
          printmap = 1;
	      for (cc = 1; cc <= themap->m; cc++) {
	        for (mm = themap->mpc[cc]; mm > 0; mm--)
	          if ( themap->types[cc][mm] == -1 || themap->types[cc][mm] == -99)
	            elim_marker(themap, individs, theparams->nn, cc, mm);
	      }
	      break;
        case 3:  /* Select chromosome and marker to eliminate */
          printmap = 1;
          fprintf(stdout,"\nPick a chromosome from 1 to %d: ",themap->m);
          theparams->wchrom = get_int();
          if ( theparams->wchrom > 0 && theparams->wchrom <= themap->m ) {
            fprintf(stdout,"\nPick a marker from 1 to %d: ",themap->mpc[theparams->wchrom]);
            theparams->mark = get_int();
            if ( theparams->mark > 0 && theparams->mark <= themap->mpc[theparams->wchrom]) 
               elim_marker(themap, individs, theparams->nn, theparams->wchrom, theparams->mark);
          }
          break;
        case 4: /*  Select a trait to eliminate. */
          fprintf(stdout,"\nPick a trait from 1 to %d: ",theparams->traits);
   	      theparams->whichtrait = get_int();
	      if ( theparams->whichtrait > 0 && theparams->whichtrait <= theparams->traits ) 
	        elim_trait(themap, individs, theparams->nn, theparams->whichtrait);
	      break;
        case 5: /* Select a trait, and eliminate individuals that don't have data for that trait. */
          fprintf(stdout,"\nPick a trait from 1 to %d: ",theparams->traits);
   	      theparams->whichtrait = get_int();
          check_misstrait(themap, individs, theparams);
          break;
        case 6: /* Specify the level of missing data to be tolerated. */  
          fprintf(stdout,"\n Specify the percent of missing marker data you will tolerate [0-100]: ");
          get_next_token(gname,MAXNAME,stdin);
          theparams->p1 = (FPN) atof(gname);
	      check_missmark(themap, individs, theparams->nn, theparams->p1);
          break;
        case 7: 
          if ( simless == 1  ) 
            break;
          do_a_bootstrap(theparams,individs,themap);  
          go_on = 0;  
          break;
        case 8: 
          if ( simless == 1  ) 
            break;
          do_a_permutation(theparams,individs,themap);  
          go_on = 0;  
          break;
        case 9: 
          if ( simless == 1  ) 
            break;
          fprintf(stdout,"\n Specify the percent of missing marker data you  want [0-100]: ");
          get_next_token(gname,MAXNAME,stdin);
          theparams->p1 = (FPN) atof(gname);
          sim_missing_markers(theparams,individs,themap);  go_on = 0;  break;
        case 10: 
          if ( simless == 1  ) 
            break;
          fprintf(stdout,"\n Specify the percent of dominant marker data you  want [0-100]: ");
          get_next_token(gname,MAXNAME,stdin);
          theparams->p1 = (FPN) atof(gname);
          sim_dominant_markers(theparams,individs,themap);  go_on = 0;  break;
        case 11: 
          if ( simless == 1  ) 
            break;
          fprintf(stdout,"\n Specify the percent proportion selectively genotyped [0-100]: ");
          get_next_token(gname,MAXNAME,stdin);
          theparams->p1 = (FPN) atof(gname);
          sim_selective_genotyping(theparams,individs,themap);  go_on = 0;  break;
        case 12: /* Print out the data in the main section. Get out of here */   
          go_on = 0; 
          break;
        case 13:
          printmap = -1;
          go_on = 0;  
          break;
      }
  }
  return(printmap);
}

/*
  Simulate missing marker data.  
*/
void sim_missing_markers(params *theparams, individual *individs, markermap *themap)
{
  int ii,cc,mm;
  writeseed=1;
  for ( ii = 1 ; ii <= theparams->nn ; ii++ )
    for (cc = 1; cc <= themap->m; cc++) {
      for (mm = themap->mpc[cc]; mm > 0; mm--)
        if ( 100.0*ranf(cc) < theparams->p1  )
          individs[ii].markers[cc][mm] = -10;
    }
}

/*
  Simulate dominant marker data.  
*/
void sim_dominant_markers(params *theparams, individual *individs, markermap *themap)
{
  int ii,cc,mm;
   FPN rval;
  individs[1].map = themap;
  determine_markers(individs,theparams->nn);
  for (cc = 1; cc <= themap->m; cc++) {
    for (mm = themap->mpc[cc]; mm > 0; mm--) 
      if ( themap->types[cc][mm] == 0 ) {
	      rval = (FPN)100.0*ranf(cc);
	      if ( rval < theparams->p1/2.0  )  {
	        for ( ii = 1 ; ii <= theparams->nn ; ii++ )
	          if ( individs[ii].markers[cc][mm] == 1 || individs[ii].markers[cc][mm] == 0)
	             individs[ii].markers[cc][mm] = 10;
	      }
	      else if ( rval < theparams->p1   )  {
	        for ( ii = 1 ; ii <= theparams->nn ; ii++ )
	          if ( individs[ii].markers[cc][mm] == 1 || individs[ii].markers[cc][mm] == 3)
	             individs[ii].markers[cc][mm] = 12;
          }
      }
  }
  append_seed(theparams, theparams->resource);
}
/*
  Simulate selectively genotyped data.    Rank each individual by trait and use
  only the tails accounting for theparams->p1 proportion of the total.
  
  do it on theparams->whichtrait
*/
void sim_selective_genotyping(params *theparams, individual *individs, markermap *themap)
{
  int ii,*indx,lb,ub,ss;
  FPN *yy;
  ii=themap->traits;
  if ( theparams->whichtrait < 1 || theparams->whichtrait > theparams->traits ) 
    theparams->whichtrait = 1;
  indx = ivector(1,theparams->nn);
  yy = dvector(1,theparams->nn);
  ss = 0;
  for (ii = 1; ii <= theparams->nn; ii++) {
    individs[ii].print_flag = 'n';
    if (individs[ii].y[theparams->whichtrait] > (FPN) MISS_VAL ) {
      ss = ss+1;
	  indx[ii] = ii;                 /* start by making indx[i]=i */
	  yy[ii]= individs[ii].y[theparams->whichtrait] ;
    }
  }
/*  Now, sort the yy vector and bring along the indx vector for the ride. */
  indxsort(ss,yy,indx);
  lb = (int) ( (FPN) ss  * theparams->p1 * 0.005 ) ;
  ub = ss - lb ;

  for ( ii = 1 ; ii <= lb ; ii++ ) 
    individs[indx[ii]].print_flag = 'y';
  for ( ii = ub+1 ; ii <= ss ; ii++ ) 
    individs[indx[ii]].print_flag = 'y';

  free_ivector(indx, 1, theparams->nn);
  free_dvector(yy, 1, theparams->nn);
  append_seed(theparams, theparams->resource);
  
}
 


/*   indxsort    Sort a vector of FPNs and move their indices at the same time.  
      modification of sort2 from Num. Rec. C,  using Heapsort  */
      
void indxsort(int n, FPN *ra, int *rb) {
  int l,j,ir,i,rrb;
  FPN rra;
        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1) {
                        rra=ra[--l];
                        rrb=rb[l];
                } else {
                        rra=ra[ir];
                        rrb=rb[ir];
                        ra[ir]=ra[1];
                        rb[ir]=rb[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                rb[1]=rrb;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                rb[i]=rb[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
                rb[i]=rrb;
        }}


      
      
/*
Do a single bootstrap resampling of individuals in the data set.
*/
void do_a_bootstrap(params *theparams, individual *individs, markermap *themap)
{
  int kk;
  individual *individs2;
   individs2 = indvector(theparams->nn, themap, NULL);
  for (kk = 1; kk <= theparams->nn; kk++)
    cp_individ( &individs[(int) ((FPN) theparams->nn * ranf(kk)) + 1], &individs2[kk]);
    /*
  for (kk = 1; kk <= theparams->nn; kk++)
    cp_individ( &individs2[kk], &individs[kk]);*/
  theparams->thedata = individs2;
  free_indvector(individs, theparams->nn);
  append_seed(theparams, theparams->resource);
}

/*
Permute the trait values.
*/
void do_a_permutation(params *theparams, individual *individs, markermap *themap)
{
  int kk,ii,*indx;
  FPN *yy;
  ii=themap->traits;

  indx = ivector(1,theparams->nn);
  yy = dvector(1,theparams->nn);
  for (ii = 1; ii <= theparams->nn; ii++)
	indx[ii] = ii;                 /* start by making indx[i]=i */
  if ( theparams->boot == 12 || theparams->boot == 6 ) /* Shuffle the index vector once...preserve trait arrays*/
    shuffle_ivector(indx, 1,theparams->nn);
  for ( ii = 1 ; ii <= theparams->traits ; ii++ ) {
    if ( theparams->boot == 2 ) /* Shuffle the index vector for each trait */
      shuffle_ivector(indx, 1,theparams->nn);
    for (kk = 1; kk <= theparams->nn; kk++)
      yy[kk] = individs[kk].y[ii];
    for (kk = 1; kk <= theparams->nn; kk++)
      individs[kk].y[ii] =  yy[ indx[kk] ]  ;
  }
  free_ivector(indx, 1, theparams->nn);
  free_dvector(yy, 1, theparams->nn);
  append_seed(theparams, theparams->resource);
  
}

/*
Create a single trait data set with no missing phenotypes
*/
void      prune2onetrait(params *theparams, individual *individs, markermap *themap)
{
  int ii,trait;
  if ( theparams->whichtrait == 0 )
    trait = 1;
  else if ( theparams->whichtrait <= themap->traits )
    trait = theparams->whichtrait;
  else
    trait = 1;
    
  if ( themap->tnames != NULL )
    strcpy( themap->tnames[1], themap->tnames[trait] );
  for ( ii=1; ii<=theparams->nn; ii++ ) {
    individs[ii].y[1]  = individs[ii].y[trait] ;
    if ( individs[ii].y[1] <= (FPN) MISS_VAL )
      individs[ii].print_flag = 'n';
    else
      individs[ii].print_flag = 'y';
   
  }
}

/*
  Remove unwanted traits from the data set.  
  
  If there are t traits, then the value of theparams->whichtrait will determine which traits to keep.
  
  
  if theparams->whichtrait <= 0, then only traits beginning with a '+' sign will be kept
  if theparams->whichtrait > t,  then all traits beginning with a '-' sign will be deleted
  if 0 < theparams->whichtrait <= t, then only theparams->whichtrait will be kept.
  
  
*/
void      prunetraits(params *theparams, individual *individs, markermap *themap)
{
  int ii,jj;
  
  for ( ii = theparams->traits; ii >= 1 ; ii-- ) {/* Get rid of unwanted traits */
    if ( theparams->whichtrait < 1 && themap->tnames[ii][0] != '+' )
      elim_trait(themap, individs, theparams->nn, ii);
    if ( theparams->whichtrait > theparams->traits &&  themap->tnames[ii][0] == '-' ) 
      elim_trait(themap, individs, theparams->nn, ii);
      
    if (  theparams->whichtrait>0 &&  theparams->whichtrait<= theparams->traits &&  theparams->whichtrait != ii)  
       elim_trait(themap, individs, theparams->nn, ii);
  }

  for ( ii=1; ii<=theparams->nn; ii++ ) {
    individs[ii].print_flag = 'y';
    for ( jj = themap->traits; jj >= 1 ; jj-- )
      if ( individs[ii].y[jj] <= (FPN) MISS_VAL )
        individs[ii].print_flag = 'n';
   
  }
}

/*
Added on 8 June 2004 for Brian Yandell.  R/QTL could  not handle otraits, so here is an option
to get rid of them.   May be unecessary now. 


  Delete otrait information from the data structure.   


from themap:
  int otraits;    Number of other traits 
  int *otypes;    Number of classes for each of the other traits 
  char **onames;   Names of other traits 


from an individ:
  char **oy;     other traits = cmatrix(1,ot,0,MAXNAME) 
  int *oyt;       oyt = ivector(1,otraits)
*/
void pruneotraits(params *theparams, individual *individs, markermap *themap) 
{

  int ii;
  for ( ii=1; ii<=theparams->nn; ii++ ) {
    if ( individs[ii].oy != NULL )
      free_cmatrix( individs[ii].oy,1,themap->otraits,0,MAXNAME);
    individs[ii].oy = NULL;
    if ( individs[ii].oyt != NULL )
      free_ivector( individs[ii].oyt,1,themap->otraits);
    individs[ii].oyt = NULL;
    
   
  }

  if ( themap->otypes != NULL )
    free_ivector(themap->otypes , 1, themap->otraits);
  themap->otypes = NULL;
  if (themap->onames != NULL)
    free_cmatrix(themap->onames, 1, themap->otraits, 0, MAXNAME);
  themap->onames = NULL;
  themap->otraits = 0; 


}


/*
Added in December 2004 for Jessica Maia.  If you have a lot of markers and a small sample size,
then many intervals will be 0.0.   This collapses those intervals, and combines marker data 
to minimize missing markers.   


A marker interval of 0.0 means that the array of marker data is the same for two markers, and will differ
only for missing markers.   This means that we loose no info by combining them.   

*/
void collapsemap(params *theparams, individual *individs, markermap *themap) 
{
  int i,j,wmark,ttable;
  individual *newdata;
  markermap *newmap;   
  FILE *LOGFILE; 
  
  
  newmap = allocate_markermap();
  newmap->m = themap->m;
  newmap->mpc = ivector(1,newmap->m);
  newmap->traits = themap->traits;
  newmap->otraits = themap->otraits;
  newmap->ml = 0;
  
  newmap->maxl = 0;
  for ( i=1; i<=newmap->m; i++ ) {
    newmap->mpc[i] = 1;
    for (j=1; j<=themap->mpc[i]; j++ )
      if ( themap->mrf[i][j] > 0.0 )  /*  decide how many nonzero intervals exist on each chromosome. */
        newmap->mpc[i] +=1;           
    if ( newmap->maxl < newmap->mpc[i] )
       newmap->maxl = newmap->mpc[i] ;
    newmap->ml = newmap->ml + newmap->mpc[i] ;
  }
  
  newmap->mrf = dmatrix(1,newmap->m,0,newmap->maxl+1);
  if ( newmap->traits > 0 ) {  /*  transfer trait names */
    newmap->tnames = cmatrix(1,newmap->traits,0,(int) MAXNAME);
    for ( i=1; i<=newmap->traits; i++ )
      strcpy( newmap->tnames[i], themap->tnames[i] );
  }
  if ( newmap->otraits > 0 ) {  /* transfer otrait names, if they exist */
    newmap->onames = cmatrix(1,newmap->otraits,0,(int) MAXNAME);
    for ( i=1; i<=newmap->otraits; i++ )
      strcpy( newmap->onames[i], themap->onames[i] );
  }
  if ( themap->cnames != NULL ) {  /*  transfer chromosome names, if they exist */
    newmap->cnames = cmatrix(1,newmap->m,0,(int) MAXNAME);
    for ( i=1; i<=newmap->m; i++ ) 
      strcpy( newmap->cnames[i], themap->cnames[i]);  
  }
  newmap->names = cmatrix(1,newmap->ml,0,(int) MAXNAME);
  newmap->ttable = imatrix(1,newmap->m,1,newmap->maxl);
  newmap->types = imatrix(1,newmap->m,1,newmap->maxl);
  newdata = indvector( theparams->nn, newmap, NULL );
  if ( newmap->traits > 0 ) {  /*  transfer trait values */
    for ( i=1; i<=theparams->nn; i++ )
      for ( j=1; j<=newmap->traits ; j++ )
        newdata[i].y[j] =  individs[i].y[j]; 
  }
  
  if ( newmap->otraits > 0 ) {  /*  transfer otrait values */
    for ( i=1; i<=theparams->nn; i++ ) {
      newdata[i].oy = cmatrix(1, newmap->otraits, 0, MAXNAME);
      for ( j=1; j<=newmap->otraits ; j++ )
        strcpy( newdata[i].oy[j] ,  individs[i].oy[j] );
        
    } 
  }
  LOGFILE = fileopen(theparams->error, "a");
  if (LOGFILE == NULL)
    LOGFILE = fileopen(theparams->error, "w");
  
  fprintf(LOGFILE, "\n#\n#     Here are the collapsed marker systems\n#\nChrom  Mark   Marker Names");
  ttable=0;
  for ( i=1; i<=newmap->m; i++ ) {  /* go through map and transfer markers with non-zero intervals.  */
    
    wmark=1;
    fprintf(LOGFILE,"\n %3d %3d   ", i,wmark);
    ttable +=1;
    for (j=1; j<=themap->mpc[i]; j++ ) {
      fprintf(LOGFILE," %s ",themap->names[ themap->ttable[i][j]] );
      strcpy(newmap->names[ttable], themap->names[ themap->ttable[i][j] ] );
      newmap->ttable[i][wmark] = ttable; 
      ImportMarkerInfo(theparams->nn, i, newdata, wmark, individs, j); /* add in data. */
      if ( themap->mrf[i][j] > 0.0 ) { 
        newmap->mrf[i][wmark] = themap->mrf[i][j];
        wmark+=1; ttable+=1;
        fprintf(LOGFILE, "\n     %3d   ",wmark);
      }
    }
  }
  fprintf(LOGFILE, "\n#\n#");
  fileclose(theparams->error, LOGFILE);
  
  if ( theparams->thedata != NULL ) 
      free_indvector(theparams->thedata, theparams->nn);
  if (theparams->theqtls != NULL)
      free_qtlvector(theparams->theqtls,theparams->themap);
  if (theparams->thegenome != NULL)
      clear_genome(theparams->thegenome);
  if ( theparams->themap != NULL )
        deallocate_markermap(theparams->themap);
  theparams->theqtls = NULL;
  theparams->thegenome = NULL;
  theparams->themap = newmap;
  theparams->thedata = newdata; 
  do_map_stats( newmap );
}

/*
   We overlay marker information.  Here, we assume that the previous marker and the current marker 
   are the same locus, or contain exactly the same information.   Thus, we add any positive information
   to the marker.   
*/
void ImportMarkerInfo(int nn, int chrom, individual  *newdata,int nmark, individual *olddata, int omark) {
  int j;
  
  for ( j=1; j<=nn; j++ ) 
    if ( newdata[j].markers[chrom][nmark] > 3 || newdata[j].markers[chrom][nmark] < 0 ) /* we have missing or dom data */
      if ( olddata[j].markers[chrom][omark] >= 0 )   /*  old data adds something.   */
        newdata[j].markers[chrom][nmark] = olddata[j].markers[chrom][omark];
  

}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Prune.c
------------------------------------------------------------------ */


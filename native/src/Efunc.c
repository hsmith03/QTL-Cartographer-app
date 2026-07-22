/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Efunc.c) is part of QTL Cartographer
         
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
/* Subroutines estimate genetic linkage maps */
void ReEstimateMap(params *theparams ) {
  int i,glen,chrom;
  genome  **gvect;
/* 

    Stage 2 :  Refine marker order and output Rmap.out and Rcross.out files

*/
  theparams->thegenome = create_genome(theparams->themap);
  addmarkernames(theparams->themap,theparams->thegenome);
  gvect = vectorizegenome(theparams->thegenome);
  glen = genomelength(theparams->thegenome);
  for ( i=1; i<=glen; i++ )  
    gvect[i]->whichqtl = i;  /* make sure that nodes know where they are */
  
/*  Refine the order.   First, determine the estimated recombination fractions and likelihoods */ 
  for ( chrom=1; chrom<=theparams->themap->m; chrom++ )
        update_chromosome_rec(theparams, gvect,glen,theparams->thedata,chrom);
 for ( i=1; i<=glen; i++ )  /*  Reset the recomb. values in theparams->themap */
   if ( gvect[i]->pxo == 0.5 )
     theparams->themap->mrf[ gvect[i]->chrom ][ gvect[i]->markr ] =(FPN) 0.0;
   else if (gvect[i]->pxo > 0.5 )
     theparams->themap->mrf[ gvect[i]->chrom ][ gvect[i]->markr ] =(FPN) 0.499;
   else
     theparams->themap->mrf[ gvect[i]->chrom ][ gvect[i]->markr ] = gvect[i]->pxo ; 
 for ( i=1; i<=theparams->themap->m; i++ )  /*  zero out marker zero */
   theparams->themap->mrf[ i ][ 0 ] = 0.0;

    glen = genomelength(theparams->thegenome);
    if (gvect != NULL )
      free_genomevector(gvect,1,glen);
    gvect = NULL;

  for (i=1; i<=theparams->nn; i++ )  /*  get data back into what the print function wants. */
    untrans_data(theparams->thedata[i].markers,theparams->themap);



}

void RCDstage2(params *theparams,  div_t emethod ) {
  int i,glen,chrom,changes;
  genome  **gvect;
  FPN thisobj;
/* 

    Stage 2 :  Refine marker order and output Rmap.out and Rcross.out files

*/
  theparams->thegenome = create_genome(theparams->themap);
  addmarkernames(theparams->themap,theparams->thegenome);
  gvect = vectorizegenome(theparams->thegenome);
  glen = genomelength(theparams->thegenome);
  for ( i=1; i<=glen; i++ )  
    gvect[i]->whichqtl = i;  /* make sure that nodes know where they are */
  
/*  Refine the order.   First, determine the estimated recombination fractions and likelihoods */ 
  for ( chrom=1; chrom<=theparams->themap->m; chrom++ )
        update_chromosome_rec(theparams, gvect,glen,theparams->thedata,chrom);



 for ( chrom=1; chrom<=theparams->themap->m; chrom++ ) {
  thisobj = Mapobj(theparams,glen,gvect);
  if ( debugging > 0 ) {
    sprintf(gbuffer,"\nChromosome %d, prior to modifications  objective function value %f",chrom,thisobj);
    LogTheError(theparams->error,gbuffer);
    show_genome(theparams->thegenome, theparams->error);
  }
 
   if (theparams->verbosity == 1)
     printf("\n\n    Modifying chromosome %d",chrom);
   if ( emethod.rem > 0 ) {

     do {
       changes = FirstLevelMod(theparams,theparams->themap,gvect,glen,theparams->thedata,chrom);   
       if (changes > 0 && theparams->thegenome->prev !=NULL)
         theparams->thegenome = firstone( gvect[1] );
     }  while (changes > 0 ); 
     thisobj = Mapobj(theparams,glen,gvect);
     if ( debugging > 0 ) {
       sprintf(gbuffer,"\nChromosome %d, step 1:  objective function value %f",chrom,thisobj);
       LogTheError(theparams->error,gbuffer);
       show_genome(theparams->thegenome, theparams->error);
     }
   }

   if ( emethod.rem > 1 ) {
     do {
       changes = SecondLevelMod(theparams,theparams->themap,gvect,glen,theparams->thedata,chrom);   
       if (changes > 0  && theparams->thegenome->prev !=NULL)
         theparams->thegenome = firstone( gvect[1] );
      }  while (changes > 0 ); 
     thisobj = Mapobj(theparams,glen,gvect);
     if ( debugging > 0 ) {
       sprintf(gbuffer,"\nChromosome %d, step 2:  objective function value %f",chrom,thisobj);
       LogTheError(theparams->error,gbuffer);
       show_genome(theparams->thegenome, theparams->error);
     }
  }
   if ( emethod.rem > 2 ) {

     do {
       changes = ThirdLevelMod(theparams,theparams->themap,gvect,glen,theparams->thedata,chrom);   
       if (changes > 0  && theparams->thegenome->prev !=NULL )
         theparams->thegenome = firstone( gvect[1] );
     }  while (changes > 0 ); 
     thisobj = Mapobj(theparams,glen,gvect);
     if ( debugging > 0 ) {
       sprintf(gbuffer,"\nChromosome %d, step 3:  objective function value %f",chrom,thisobj);
       LogTheError(theparams->error,gbuffer);
       show_genome(theparams->thegenome, theparams->error);
     }
  }
 }   

 for ( i=1; i<=glen; i++ )  /*  Reset the recomb. values in theparams->themap */
   if ( gvect[i]->pxo == 0.5 )
     theparams->themap->mrf[ gvect[i]->chrom ][ gvect[i]->markr ] =(FPN) 0.0;
   else
     theparams->themap->mrf[ gvect[i]->chrom ][ gvect[i]->markr ] = gvect[i]->pxo ; 

  thisobj = Mapobj(theparams,glen,gvect);
  if ( debugging > 0 ) {
    sprintf(gbuffer,"\nAfter modifications,  objective function value %f",thisobj);
    LogTheError(theparams->error,gbuffer);
    show_genome(theparams->thegenome, theparams->error);

  }  
    glen = genomelength(theparams->thegenome);
    if (gvect != NULL )
      free_genomevector(gvect,1,glen);
    gvect = NULL;

  for (i=1; i<=theparams->nn; i++ )  /*  get data back into what the print function wants. */
    untrans_data(theparams->thedata[i].markers,theparams->themap);

}



void RCDstage1(params *theparams,  char *chptr, char *progname) {
  FILE *fptr;
  int i,j,k,glen,chrom;
  genome  **gvect;
  FPN lr;
/*   

      Stage I :   

   1.  Read in data and any old map
   2.  Check for segregation, and remove markers that fail
   3.  Create a set of linkage groups with an initial order
   4.  Remove markers that aren't linked to anything
   5.  Reorder the markers based on anchoring
   6.  print out the map and data set in map.inp, cross.inp formats



*/    
  theparams->thegenome = create_genome(theparams->themap);
  addmarkernames(theparams->themap,theparams->thegenome);
  gvect = vectorizegenome(theparams->thegenome);
  glen = genomelength(theparams->thegenome);
  for ( i=1; i<=glen; i++ )  
    gvect[i]->whichqtl = theparams->themap->types[gvect[i]->chrom][gvect[i]->markr] ;
  
  check_segregation(theparams->thedata,theparams->themap,theparams,theparams->thegenome);	 
  fptr = fileopen(theparams->error, "a");
  fprintf(fptr,"\n   Emap: Marker segregation tests");
  putline(fptr,'-',50);
  putline(fptr,'-',50);
  fprintf(fptr, "\nNumber MarkerName       X^2      Pr(X^2)   Ind.");
  putline(fptr,'-',50);

  for ( i=1; i<=glen; i++ )
    fprintf(fptr,"\n   %3i %-15s %8.4f  %8.4f %4d",i, theparams->themap->names[ theparams->themap->ttable[ gvect[i]->chrom][ gvect[i]->markr]], gvect[i]->pxo, gvect[i]->mxo, gvect[i]->whichqtl);
  putline(fptr,'-',50);
  putline(fptr,'-',50);
  fileclose(theparams->error, fptr);
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n  Finished segregation analysis"); 
  
  for (i=1; i<=glen; i++ )  /*  move all uninformative, unclassified or poorly segregating markers   to the front */
    if ( gvect[i]->whichqtl == -99 )
      theparams->thegenome = mv_genome_front(gvect[i]);
  while ( theparams->thegenome != NULL && theparams->thegenome->whichqtl == -99 )  /*  delete them */
    theparams->thegenome = del_genome_node(theparams->thegenome);
  if ( theparams->thegenome == NULL )
    exit(1);
  free_genomevector(gvect,1,glen);  /*  get rid of the old vector*/
  gvect =  vectorizegenome(theparams->thegenome); /* create a new vector */
  glen = genomelength(theparams->thegenome);
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n  Got rid of old gvect, created new one. "); 
  for ( i=1; i<=glen; i++ ) {
    gvect[i]->values = dvector(1,glen);
    gvect[i]->n = glen;
    gvect[i]->whichqtl = 0;  /* 0 now means that it is not on a chromosome*/
  }
  for (i=1; i<glen; i++)
    for (j=i+1; j<=glen; j++ ) 
        gvect[i]->values[j] = gvect[j]->values[i] = EstimateRecomb(i,j,theparams->thedata,theparams,gvect,&lr);
/*  for (i=1; i<glen; i++)
    for (j=i+1; j<=glen; j++ ) 
      printf("\n  %d  %d  %f",i,j,gvect[i]->values[j]);*/
  chrom = 1;
  for ( i=1; i<=glen; i++ ) {
    gvect[i]->chrom = gvect[i]->markr = 0;  /* 0 now means that it is not on a chromosome*/
    gvect[i]->whichqtl = i;  /* make sure that nodes know where they are */
  }
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n  Ready to do linkage maps"); 
  do {
    chrom = dolinkgroup(chrom,glen,gvect,theparams);
  } while ( chrom != 0 );
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n  Did initial linkage groups...now get rid of singlets"); 
/*
    Need to add some code that gets rid of singlets...they are of no use.  
    Or should we just put them on their own chromosome?
*/
  theparams->thegenome  = firstone(gvect[1]);
/*  for ( gptr=theparams->thegenome; gptr!=NULL; gptr=gptr->next )  
    printf("\n%d %d %d",gptr->chrom, gptr->markr, gptr->whichqtl);*/

  for (i=1; i<=glen; i++ )  /*  move all uninformative, unclassified or poorly segregating markers   to the front */
    if ( gvect[i]->chrom == 0 )
      theparams->thegenome = mv_genome_front(gvect[i]);
  while ( theparams->thegenome != NULL && theparams->thegenome->chrom == 0 )  /*  delete them */
    theparams->thegenome = del_genome_node(theparams->thegenome);
  if ( theparams->thegenome == NULL )
    exit(1);
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n  Got rid of singlets...now reorder genome"); 
/*  Rearrange the map taking into account the 'Anchors' */  

  theparams->thegenome = reorder_genome(theparams->thegenome,theparams,glen,gvect);

  if ( debugging > 0 )
  	LogTheError(theparams->error,"\n  Reordered genome...now print map"); 
/*  
    Print the map to stem.e1m and the data to stem.e1c
    in  the map.inp and cross.inp formats, respectively.   Once done, 
    read in a new map from stem.e1m and the data from stem.e1c and assume 
    that the number markers will not change.  
*/  
  strcpy(theparams->tfile,theparams->stem);
  strcat(theparams->tfile,".e1m");
  print_genome_map(theparams->thegenome,theparams,theparams->tfile,progname, chptr);
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n printed map...now get rid of anchor"); 
  strcpy(theparams->tfile,theparams->stem);
  strcat(theparams->tfile,".e1c");
  for (i=1; i<=theparams->themap->ml; i++ ) {/*  Get rid of Anchor|p|c| at beginning */
    strcpy(gname,theparams->themap->names[i]);
    strlwr(gname);
    if ( !strncmp(gname,"anchor",6) ) {
      for ( j=9; j<MAXNAME && gname[j] != '|'; j++ ) k+=1;
      for ( k=j+1; k< (MAXNAME - j - 1) ; k++ )
        gname[k-j-1] = gname[k];
      strcpy(theparams->themap->names[i],gname);
    }  
  }
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n Gid rid of Anchor...now translate data and print it..."); 
  for (i=1; i<=theparams->nn; i++ )  /*  get data back into what the print function wants. */
    untrans_data(theparams->thedata[i].markers,theparams->themap);
  print_head(progname,theparams->tfile,chptr,0,31,theparams);
  print_individuals_std(theparams,theparams->thedata,theparams->nn,theparams->tfile);
  if ( debugging > 0 )
    LogTheError(theparams->error,"\n Printed data...now get rid of structures");
/*  get rid of the old structures*/
  if (theparams->thegenome != NULL) {
    glen = genomelength(theparams->thegenome);
    if (gvect != NULL )
      free_genomevector(gvect,1,glen);
    gvect = NULL;
    clear_genome(theparams->thegenome);  theparams->thegenome = NULL;
  }
  if ( theparams->thedata != NULL )
    free_indvector(theparams->thedata, theparams->nn);
  theparams->thedata = NULL;
  deallocate_markermap(theparams->themap);  theparams->themap = NULL;
}



/*
  For each marker in turn, test whether it conforms to Hardy-Weinberg
  proportions.
*/
void  check_segregation(individual *individs,markermap *themap,params *theparams,genome *gptr) {
  FPN freqs[3],test_stat,opAA,opAa,opaa,epAA,epAa,epaa ,dof;
  int ii, chrom, mark, nm, nmt,counts[3];
  genome *tgptr;
  assign_q(freqs,theparams);


  
  for (tgptr=gptr; tgptr!= NULL ; tgptr = tgptr->next ) { 
      chrom = tgptr->chrom;
      mark = tgptr->markr;
      counts[0]=counts[1]=counts[2]=0;
      tgptr->whichqtl = themap->types[chrom][mark] ;

      nm = nmt = 0;
      for ( ii = 1 ; ii <= theparams->nn ; ii++ )       
        switch ( individs[ii].markers[chrom][mark] ) {
          case -2 :
          case -1 :
            counts[0] +=1;   
            break;  
          case  0 :
            counts[1] +=1;   
            break;  
          case  1 :
          case  2 :
            counts[2] +=1;   
            break;  
        }
      nm = counts[0] + counts[1] + counts[2];
      opAA = (FPN) counts[2]  ;
      opAa = (FPN) counts[1]  ;
      opaa = (FPN) counts[0]  ;
      epAA = freqs[0] * (FPN) nm;
      epAa = freqs[1] * (FPN) nm;
      epaa = freqs[2] * (FPN) nm;
      test_stat = (FPN)0.0;
      if ( themap->types != NULL )
        switch ( themap->types[chrom][mark] ) {
          case -1 :  /* a dominant system*/
            if ( epAa+epaa > 0.0 )
              test_stat = (opAa+opaa-epAa-epaa)*(opAa+opaa-epAa-epaa)/(epAa+epaa);
            if ( epAA > 0.0 )
              test_stat = test_stat + (opAA-epAA)*(opAA-epAA)/epAA ;   
            break;  
          case  0 :   /* codominant system*/
            if ( epAA > 0.0 )
              test_stat =  (opAA-epAA)*(opAA-epAA)/epAA ;   
            if ( epAa > 0.0 )
              test_stat = test_stat + (opAa-epAa)*(opAa-epAa)/epAa ;   
            if ( epaa > 0.0 )
              test_stat = test_stat  + (opaa-epaa)*(opaa-epaa)/epaa ;   
            break;  
          case  1 : /* A dominant  system*/
            if ( epAa+epAA > 0.0 )
              test_stat = (opAa+opAA-epAa-epAA)*(opAa+opAA-epAa-epAA)/(epAa+epAA);
            if ( epaa > 0.0 )
              test_stat = test_stat  + (opaa-epaa)*(opaa-epaa)/epaa ;   
            break;  
          default :
            test_stat = -(FPN)1.0;   
            break;  
        }
      else
        test_stat = (FPN)0.0;
     if ( theparams->ngt == 3 && themap->types[chrom][mark] == 0) 
       dof = (FPN)2.0;
     else
       dof = (FPN)1.0;
     tgptr->pxo = test_stat; 
     tgptr->mxo =   chiprob( (int) dof, test_stat);
     if ( tgptr->mxo < theparams->segsize )
       tgptr->whichqtl = -99;
     if ( tgptr->whichqtl == 1 ) {
       if ( theparams->cross == 1 && theparams->tcross == 0 )
         tgptr->whichqtl = -99;
       if ( theparams->cross == 4 && theparams->tcross == 1 )
         tgptr->whichqtl = -99;
       if ( theparams->cross == 5 && theparams->tcross == 1 )
         tgptr->whichqtl = -99;
     }
     if ( tgptr->whichqtl == -1 ) {
       if ( theparams->cross == 2 && theparams->tcross == 0 )
         tgptr->whichqtl = -99;
       if ( theparams->cross == 4 && theparams->tcross == 2 )
         tgptr->whichqtl = -99;
       if ( theparams->cross == 5 && theparams->tcross == 2 )
         tgptr->whichqtl = -99;
     }
       
  }


}





/*
    Assume that chrom does not yet exist.  Look at pool of markers and picked two
    most closely linked.   If their rec. fraction beats the criteria, then
    put them on the chromosome at the end of the chain.  
    
    Else, return  0 to tell main to stop mapping. 
*/
int dolinkgroup(int chrom,int glen,genome **gvect,params *theparams) {

  int  l1,l2,i,j,left,right,mmin;
  FPN rmin,lmin,min;
/* First, initialize the linkage group */
  min = (FPN)0.5;
  for (i=1; i<glen; i++ )
    for ( j=i+1; j<=glen ; j++ ) 
      if ( gvect[i]->chrom == 0 && gvect[j]->chrom == 0 && gvect[i]->values[j] < min) {
          min = gvect[i]->values[j] ;
          l1 = i;
          l2 = j;           
      }
   if ( min < theparams->linksize ) {
     gvect[l1]->chrom = gvect[l2]->chrom = chrom;
     gvect[l1]->markr = 1;
     gvect[l2]->markr = 2;
     mv_genome_node(gvect[l1], NULL);
     gvect[l1]->dist = gvect[l1]->values[l2];
     mv_genome_node(gvect[l2], NULL);
     gvect[l2]->dist = (FPN)0.0;
   }
   else 
     return(0);

/* now extend the chain in either direction from l1 and l2 */
  do {
	  lmin = rmin = (FPN)0.5;
	  for ( i=1; i<=glen; i++ ) {
	    if ( gvect[i]->chrom == 0 && gvect[i]->values[l1] < lmin ) {
	          lmin = gvect[i]->values[l1] ;
	          left = i;
	        }
	    if ( gvect[i]->chrom == 0 && gvect[i]->values[l2] < rmin ) {
	          rmin = gvect[i]->values[l2] ;
	          right = i;
	        }
	   }

	   if ( lmin < rmin  && lmin < theparams->linksize ) {
	     gvect[left]->chrom = chrom;
	     gvect[left]->markr = gvect[l1]->markr - 1;
	     mv_genome_node(gvect[left],gvect[l1]);
	     gvect[left]->dist = gvect[left]->values[l1] ;
	     l1 = left;
	   }
	   else if ( rmin <= lmin && rmin < theparams->linksize ) {
	     gvect[right]->chrom = chrom;
	     gvect[right]->markr = gvect[l2]->markr + 1;
	     mv_genome_node(gvect[right],NULL);
	     gvect[l2]->dist = gvect[l2]->values[right];
	     gvect[right]->dist = (FPN)0.0 ;
	     l2 = right;
	   }
	   else
	     l1 = 0;
	   
   } while (l1 > 0 );
   
   mmin = 2;
   for ( i=1; i<=glen; i++ )
     if ( gvect[i]->chrom == chrom && gvect[i]->markr < mmin )
       mmin = gvect[i]->markr;
   if ( mmin <= 0 )
     mmin = -1 * mmin  + 1;
   else
     mmin = mmin - 1;
   for ( i=1; i<=glen; i++ )
     if ( gvect[i]->chrom == chrom  )
        gvect[i]->markr = gvect[i]->markr + mmin;
   
   return(chrom+1);
}

/*

Now, reorder the whole map based on whether there are Anchor nodes.   

the marker name 

    Anchor|p|c|name
    
where    p indicates a postion and takes on values l (left telomere), r (right telomere) or c (centromere)
                p is not yet active
         c indicates the chromosome number
         name indicates the marker name
         
         
The name Anchor|p|c|name will be replaced by 'n' in the output

*/
genome *reorder_genome(genome *thegenome,params *theparams,int glen,genome **gvecto) {
  int oldchrom, newchrom,i,j,k,l,nchrom,mmax,mmin;
  char *chptr,lcr;
  genome *gptr,*tgptr,*first,*last,*right,*left,**gvect;
  i = theparams->traits;
  nchrom = 0;  /*  Determine the number of chromosomes */
  for ( gptr = thegenome; gptr != NULL ; gptr = gptr->next )
    if (gptr->chrom > nchrom)
      nchrom = gptr->chrom;
/*  Check for Anchor Nodes and reassign chromosome numbers and indicate if a chromosome should be reversed */  
  for ( gptr = thegenome; gptr != NULL ; gptr = gptr->next ) {    
    strcpy(gbuffer, gptr->markername); 
    chptr = strlwr(gbuffer);
    if ( !strncmp(gbuffer,"anchor",6 )   ) {  /* This is an anchor node. */
      oldchrom = gptr->chrom;
      lcr = gbuffer[7];
      for ( i=9; gbuffer[i] != '|' ; i++ )
        gname[i-9] = gbuffer[i];
      gname[i-9] = '\0';
      j=i+1;
      newchrom = atoi(gname);
      for ( i=j; gbuffer[i] != '\0' ; i++ )
        gname[i-j] = gbuffer[i];
      gname[i-j] = '\0';
      for ( tgptr=thegenome; tgptr!=NULL; tgptr=tgptr->next )
        if ( tgptr->chrom == oldchrom )  {  /*  swap chromosome assignments */
          tgptr->chrom = newchrom;
          
        }
        else if ( tgptr->chrom == newchrom )
          tgptr->chrom = oldchrom;
      strcpy(gptr->markername, gname);
      if ( lcr == 'l' || lcr == 'r' ) {  /* if the anchor is a left or right telomere, then reorder the chromosome */

        k=0;
        for ( l=1; l<=glen;l++)
          if (gvecto[l]->chrom == gptr->chrom )  
            k +=1 ; /* how many markers are on this chromosome?*/

        if ( (lcr=='r' && gptr->markr<k/2) || (lcr=='l' && gptr->markr>k/2) )   /* switch only if marker is on wrong side*/
          gptr->markr = - gptr->markr;
      }
    }
  }
  
  for ( i=1; i<=glen; i++ ) {   
    if ( gvecto[i]->markr < 0 )  {  /* reverse the order of this chromosome if markr is negative*/
      gvecto[i]->markr = -gvecto[i]->markr;
      mmax = glen;
      mmin = 0;
      k=0;
      for ( j=1; j<=glen; j++ ) /*  find the first and last nodes of the chromosome */
        if ( gvecto[j]->chrom == gvecto[i]->chrom ) {
          k+=1; /* how many markers are on this chromosome?*/
          if ( gvecto[j]->markr < mmax ) {
            mmax = gvecto[j]->markr;
            first = gvecto[j];
          }
          if ( gvecto[j]->markr > mmin ) {
            mmin = gvecto[j]->markr;
            last = gvecto[j];
          }
        
        
        }
      right = last->next;
      left = first->prev;
      gvect = genomevector(1,k); 
      l=0;
      for ( tgptr=first; tgptr!=last->next; tgptr=tgptr->next)  {
        l +=1 ;
        gvect[l] = tgptr;
      }
      for (l=k;l>0;l--) {
        gvect[l]->markr = k-l+1;
        if ( l>1 ) {
          gvect[l]->dist = gvect[l-1]->dist;
          gvect[l]->next = gvect[l-1];
          if ( l<k )
            gvect[l]->prev = gvect[l+1];
          else
            gvect[l]->prev = left;
        }
        else {
          gvect[l]->dist = (FPN)0.0;
          gvect[l]->next = right;
          gvect[l]->prev = gvect[l+1];
        }
        if (right!=NULL)
          right->prev = gvect[1];
        if (left!=NULL )
          left->next = gvect[k];
            
      }
      free_genomevector(gvect,1,k);
      gvect = NULL;
    }
  }
  left =  firstone( gvecto[1] );/* Find the beginning of the genome. */


/*Reorder the chromosomes in the genome to reflect their ordinal */
  gvect = genomevector(1,2*nchrom);    
  for ( tgptr=left; tgptr!=NULL; tgptr=tgptr->next ) {
    if ( tgptr->prev == NULL || (tgptr->prev != NULL && tgptr->prev->chrom != tgptr->chrom) )
      gvect[tgptr->chrom] = tgptr;
    if (  tgptr->next == NULL || (tgptr->next != NULL && tgptr->next->chrom != tgptr->chrom) )
      gvect[nchrom+tgptr->chrom] = tgptr;
  }  
  
  gvect[1]->prev = NULL;
  for ( i=2; i<=nchrom; i++ ) {
    gvect[i]->prev = gvect[nchrom+i-1];
    gvect[i]->prev->next = gvect[i];    
  }  
  gvect[nchrom+nchrom]->next = NULL; 
  left = gvect[1];
  free_genomevector(gvect,1,2*nchrom); gvect = NULL;
/*  Return a pointer to the first node in the genome*/
  return(left);

}

/*
    The first level modification for a chromosome chrom will perform a switch on each of
    the m(m-1)/2 pairs of loci.   The objective function will be calculated for each switch,
    and the best switch will be compared to the original objective function.  If better, 
    a permanant switch is made and a 1 is returned.  If it is not a better order, a zero is returned.
*/
int FirstLevelMod(params *theparams,markermap *themap,genome **gvect,int glen,individual *thedata,int chrom) {
  FPN iobj, minobj, thisobj; 
  int i,j,l1,l2,error;
  genome *g1ptr,*g2ptr,*bg1,*bg2;
  update_chromosome_rec(theparams,gvect,glen,theparams->thedata, chrom);
  iobj = minobj = obj_func(theparams,chrom,glen,gvect);

  bg1=bg2=NULL;

  
  for ( l1=1; l1<themap->mpc[chrom]; l1++ )
    for ( l2=l1+1; l2<=themap->mpc[chrom]; l2++ ) {
      i = GvectNode(gvect,glen,chrom,l1);
      j = GvectNode(gvect,glen,chrom,l2);
      g1ptr = gvect[i];
      g2ptr = gvect[j];
      error = swap_gnodes(g1ptr,g2ptr); 
      if ( error == 1 ) {
        update_recombs( g1ptr, g2ptr, thedata, theparams, gvect,  chrom ) ;
        thisobj  = obj_func(theparams,chrom,glen,gvect);

        if ( thisobj < minobj ) {
          minobj = thisobj;
          bg1 = g1ptr;
          bg2 = g2ptr;
        }
        error = swap_gnodes(g1ptr,g2ptr);  /* swap back*/
        update_recombs( g1ptr, g2ptr, thedata, theparams, gvect,  chrom ) ;


      
      }
    
    }
  if ( bg1 == NULL ) {
    /*update_chromosome_rec(theparams,gvect,glen,theparams->thedata, chrom);*/
    return(0);
  }
  else {
    error = swap_gnodes(bg1,bg2);
    /*update_chromosome_rec(theparams,gvect,glen,theparams->thedata, chrom);*/
    update_recombs( bg1, bg2, thedata, theparams, gvect,  chrom ) ;
    SwapNodesData(thedata,theparams,bg1,bg2,themap) ;   
    if (theparams->verbosity == 1 )
          printf("\n 1st Mod: Swap markers %d <-> %d. Chrom %d, minobj = %f, init_obj = %f",bg1->markr,bg2->markr,chrom,minobj,iobj);
    return(1);
  }
}

/*
  go through a vector of nodes, and find the one the corresponds to 'marker' on 'chrom'
*/
int GvectNode(genome **gvect,int glen,int chrom, int marker) {

  int i;
      for ( i=1; i<=glen; i++ )  
        if ( gvect[i]->chrom == chrom && gvect[i]->markr == marker )
          return(i);
  return(0);     
}  

/*  for a pair of swapped genome nodes, update the recombination fractions */
void update_recombs(genome *g1ptr, genome *g2ptr,individual *thedata, params *theparams, genome **gvect, int chrom ) {


  update_interval(g1ptr,thedata,theparams,gvect,chrom);
  if ( g1ptr !=NULL )
    update_interval(g1ptr->prev,thedata,theparams,gvect,chrom);
  update_interval(g2ptr,thedata,theparams,gvect,chrom);
  if (g2ptr != NULL )
    update_interval(g2ptr->prev,thedata,theparams,gvect,chrom);

}

/*
  update the recombination fraction for an interval
*/
void update_interval(genome *lptr, individual *thedata, params *theparams, genome **gvect, int chrom) {
  
    if ( lptr != NULL && lptr->chrom == chrom ) { 
      if ( lptr->next != NULL && lptr->next->chrom == chrom ) {
        lptr->pxo = EstimateRecomb(lptr->whichqtl,lptr->next->whichqtl,thedata,theparams,gvect,&lptr->mxo);  
        lptr->dist = mapfunc(lptr->pxo, 1);
      }  
      else   {
        lptr->pxo=(FPN)0.5; lptr->mxo=(FPN)0.0;
      }   
    }         
}


/*
    The second level modification for a chromosome chrom will perform a permutation on each of
    the (m-2) consecutive triplets of loci.   The objective function will be calculated for each permutation,
    and the best permuted order will be compared to the original objective function.  If better, 
    a permanant switch is made and a 1 is returned.  If it is not a better order, a zero is returned.
    
    Note that there are 6 ordered permuations.   One of these is the original order, and three have already been
    tried by the level 1 modification.  Thus, there are two new permutations to try against the original. 
*/
int SecondLevelMod(params *theparams,markermap *themap,genome **gvect,int glen,individual *thedata,int chrom) {
  FPN iobj, minobj, thisobj; 
  int i,l1, right;
  genome *aptr,*bptr,*cptr,*dptr,*bg1,*bg2;
  
  iobj = minobj = obj_func(theparams,chrom,glen,gvect);

  bg1=bg2=NULL;

  right = -1;
  for ( l1=1; l1 < themap->mpc[chrom]-1; l1++ )  {
      i = GvectNode( gvect, glen, chrom, l1);
      aptr = gvect[i];

      bptr = aptr->next;
      cptr = bptr->next;
      dptr = cptr->next;
      SwapNodeDataUp(thedata,theparams,aptr,cptr,themap);
      update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
      thisobj  = obj_func(theparams,chrom,glen,gvect);
      if ( thisobj < minobj ) {
          minobj = thisobj;
          bg1 = aptr;
          bg2 = cptr;
          right = 1;  
      }
      SwapNodeDataUp(thedata,theparams,bptr,aptr,themap);
      update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
      thisobj  = obj_func(theparams,chrom,glen,gvect);
      if ( thisobj < minobj ) {
          minobj = thisobj;
          bg1 = cptr;
          bg2 = aptr;
          right = 0; 
      }
      SwapNodeDataUp(thedata,theparams,cptr,bptr,themap);

      update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
  }
  if ( right == 1 ) {
        if (theparams->verbosity == 1 )
          printf("\n 2nd Mod: Marker %d -> after %d. Chrom %d, minobj = %f, init_obj = %f",bg1->markr,bg2->markr,chrom,minobj,iobj);
        SwapNodeDataUp(thedata,theparams,bg1,bg2,themap);
        update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
        return(1);
   }
   else if ( right == 0 ) {
        if (theparams->verbosity == 1 )
          printf("\n 2nd Mod: Marker %d <- before %d. Chrom %d, minobj = %f, init_obj = %f",bg1->markr,bg2->markr,chrom,minobj,iobj);
        SwapNodeDataDown(thedata,theparams,bg1,bg2,themap);
        update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
        return(1);

    }
    else
      return(0);
  
}

/*  
   Estimate the recombination fractions between consecutive markers on a chromosome.   Set
   the last marker to r=0.5 and lr=0.0.  
*/
void update_chromosome_rec(params *theparams,genome **gvect,int glen,individual *thedata,int chrom) {
  int i;
  for ( i=1; i<=glen; i++ ) 
    update_interval(gvect[i],thedata,theparams,gvect,chrom);
    
}


/*  We assume that bg1 and bg2 are on the same chromosome, and that bg1 comes prior to
    bg2 in the chain.  We are moving  bg1  to a point after bg2.   As a result,  all data
    are moved down.  

*/
void SwapNodeDataUp(individual *thedata,params *theparams, genome *bg1, genome *bg2, markermap *themap) {
  genome *gptr;
  int error;
  gptr = bg1;
  while ( gptr->prev != bg2 ) {  /*  Move it to just after bg2 */
    error = swap_gnodes(gptr,gptr->next);
    SwapNodesData(thedata,theparams,gptr,gptr->prev,themap); 

  }
}

/*  We assume that bg1 and bg2 are on the same chromosome, and that bg1 comes after
    bg2 in the chain.  We are moving  bg1  to a point before bg2.   As a result,  all data
    are moved up.  */
void SwapNodeDataDown(individual *thedata,params *theparams, genome *bg1, genome *bg2, markermap *themap) {
  genome *gptr;
  int error;
  gptr = bg1;
  while ( gptr->next != bg2 ) {  /*  Move it to just before bg2 */
    error = swap_gnodes(gptr,gptr->prev);
    SwapNodesData(thedata,theparams,gptr,gptr->next,themap);  

  }
}


/*  We assume that bg1 and bg2 are on the same chromosome  We are switching the data
    between their two positions.  This doesn't rely on whether they have been swapped in the chain. */
void SwapNodesData(individual *thedata,params *theparams, genome *bg1, genome *bg2, markermap *themap) {
  int i,savemark,savename;
    /*swap data  This is a big swap...need to move everything down...*/
	for (i=1; i<= theparams->nn; i++) {
	  savemark = thedata[i].markers[bg1->chrom][bg1->markr] ;
 
	  thedata[i].markers[bg1->chrom][bg1->markr] = thedata[i].markers[bg2->chrom][bg2->markr];
	  thedata[i].markers[bg2->chrom][bg2->markr] = savemark;
	}
	/*swap pointer to names and move down marker numbers*/
	savemark = bg1->markr;
	savename = themap->ttable[bg1->chrom][bg1->markr]; 
	themap->ttable[bg1->chrom][bg1->markr] = themap->ttable[bg2->chrom][bg2->markr];
	bg1->markr = bg2->markr;
	themap->ttable[bg2->chrom][bg2->markr] = savename;
	bg2->markr = savemark;
}

/*
    The third level modification will try each marker in turn in each interval to test if the order is better.
    If a better order is found, a 1 is returned, otherwise a zero is returned.  
*/
int ThirdLevelMod(params *theparams,markermap *themap,genome **gvect,int glen,individual *thedata,int chrom) {
  FPN iobj, minobj, thisobj; 
  int l1,il1,il2;
  genome *aptr,*bptr,*cptr,*dptr;
  iobj = minobj = obj_func(theparams,chrom,glen,gvect);
  bptr=cptr=NULL;


  for ( l1=1; l1<=themap->mpc[chrom]; l1++ )  { 
    il1 = GvectNode( gvect, glen, chrom, l1);
    dptr = gvect[il1]->next;
    il2 =  GvectNode( gvect, glen, chrom, 1);
    if ( il1 >  il2 )
      SwapNodeDataDown(thedata,theparams,gvect[il1],gvect[il2],themap);
    aptr = gvect[il1];
    while ( aptr->next != NULL && aptr->next->chrom == chrom ) {
      update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
          thisobj  = obj_func(theparams,chrom,glen,gvect);
      if ( thisobj < minobj ) {
          minobj = thisobj;
        bptr=aptr;
        if ( aptr->next == NULL || aptr->next->chrom != chrom)
          cptr=aptr->prev;
        else
          cptr=aptr->next;
	  }
	  if (aptr->next != NULL )
	    SwapNodeDataUp(thedata,theparams,aptr,aptr->next,themap);

     
    }   
      SwapNodeDataDown(thedata,theparams,gvect[il1],dptr,themap);
    
  }    

  if ( bptr == NULL ) {
    update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
    return(0);
  }
  else {
    if (   cptr->markr > bptr->markr ) {
         if (theparams->verbosity == 1 )
          printf("\n 3rd Mod: Marker %d -> after %d. Chrom %d, minobj = %f, init_obj = %f",bptr->markr,cptr->markr,chrom,minobj,iobj);
        SwapNodeDataUp(thedata,theparams,bptr,cptr,themap);
        update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
        return(1);    
    }
    else {
        if (theparams->verbosity == 1 )
          printf("\n 3rd Mod: Marker %d <- before %d. Chrom %d, minobj = %f, init_obj = %f",bptr->markr,cptr->markr,chrom,minobj,iobj);
        SwapNodeDataDown(thedata,theparams,bptr,cptr,themap);
        update_chromosome_rec(theparams, gvect,glen,thedata,chrom);
        return(1);
    }
  }
}


/*
   Objective function SAR:  sum of recombination frequencies.   These will be stored in gptr->pxo.
   We want to minimize SAR.  Note that if mxo == 0.0, that the node is on the end of the chromosome
   and shouldn't contribute to the SAR. 
*/
FPN obfuncsar(int chrom,int glen,genome **gvect) {
  FPN total;
  int i;
  total=(FPN)0.0;
  for (i=1; i<=glen; i++)
    if (gvect[i]->chrom == chrom  && gvect[i]->mxo < 0.0 )
      total = total+gvect[i]->pxo;

  return(total);
}

/*
   Objective function SAL:  sum of likelihoods.   These will be stored in gptr->mxo.
   We want to maximize SAL, but to make things easier in the calling routine, we will return
   -SAL and allow it to minimize the quantity.
*/
FPN obfuncsal(int chrom,int glen,genome **gvect) {
  FPN total;
  int i;
  total=(FPN)0.0;
  for (i=1; i<=glen; i++)
    if (gvect[i]->chrom == chrom)
      total = total+gvect[i]->mxo;

  return( -total);
}

/*
  Calculate the objective function base on the value of theparams->emapobj
  obj_func(theparams,chrom,glen,gvect);
*/
FPN obj_func(params *theparams, int chrom,int glen,genome **gvect) {
  FPN iobj;
  
  if ( theparams->emapobj == 0 )
    iobj =  obfuncsal( chrom, glen, gvect);
  else
    iobj =   obfuncsar( chrom, glen, gvect);
  return(iobj);
}


/*

                                      Counts
---------------------------------------------------------------------------------  
---------------------------------------------------------------------------------  
n_i           co/co   A-/co  a-/co    co/B-   co/b-   A-/B-  A-/b-  a-/B-  a-/b-
---------------------------------------------------------------------------------  
  1          AA/BB   A-/BB   AA/BB    AA/B-   AA/BB   A-/B-  A-/BB   AA/B-  AA/BB
  2          AA/Bb   A-/Bb   AA/Bb
  3          AA/bb   A-/bb   AA/bb    AA/bb   AA/b-   A-/bb  A-/b-   AA/bb  AA/b-
  4          Aa/BB                    Aa/B-   Aa/BB   
  5          Aa/Bb
  6          Aa/bb                    Aa/bb   Aa/b-
  7          aa/BB   aa/BB   a-/BB    aa/B-   aa/BB   aa/B-  aa/BB   a-/B-  a-/BB
  8          aa/Bb   aa/Bb   a-/Bb
  9          aa/bb   aa/bb   a-/bb    aa/bb   aa/b-   aa/bb  aa/b-   a-/bb  a-/b-
---------------------------------------------------------------------------------  
---------------------------------------------------------------------------------  
Classes that are collapsed:
                     4->1    4->7     2->1    2->3  2,4,5->1  8->9    6->9  5,6,8->9
                     5->2    5->8     5->4    5->6    6->3  2,5,6->3  2->1    2->3
                     6->3    6->9     8->7    7->8    8->7    4->1  4,5,8->7  4->7

n = sum(ni)   cnts[0] = sum(cnts[i]), i=1,9

Estimates of recombination

B1:  neither locus A-        r = (n2+n4)/n
B2:  neither locus a-        r = (n6+n8)/n

SFx, RFx:  done iteratively in estFXrec using fxlike

-----------------------------------------------------------------------------------------------------------------------------------------    
-----------------------------------------------------------------------------------------------------------------------------------------    
Line          p(AABB),p(aabb)   p(AAbb),p(aaBB)         Est(r)                          l(r) = ln(L(r))
-----------------------------------------------------------------------------------------------------------------------------------------    
RI0 (DH):      (1-r)/2                r/2             (n3+n7)/n                 (n1+n9)ln((1-r)/2) + (n3+n7)ln(r/2)
RI1:          1/2* 1/(1+2r)           r/(1+2r)      0.5 * (n3+n7)/(n1+n9)     (n1+n9)ln(1/2* 1/(1+2r) ) + (n3+n7)ln(r/(1+2r))
RI2:          1/2* (1+2r)/(1+6r)      2r/(1+6r)  (n3+n7)/(4(n1+n9) - 2(n3+n7))    (n1+n9)ln(1/2* (1+2r)/(1+6r)) + (n3+n7)ln(2r/(1+6r))
-----------------------------------------------------------------------------------------------------------------------------------------    
-----------------------------------------------------------------------------------------------------------------------------------------    

Test crosses?  They should just be the first part of the cross.


*/
FPN EstimateRecomb(int l1, int l2,individual *individs,params *theparams,genome **gvect,FPN *lr) {
  FPN rec;
  int cnts[10],    i,c1,m1,c2,m2;
/*  First, count the number of different observed genotypes */
  if ( gvect[l1]->whichqtl != -99 && gvect[l2]->whichqtl != -99 ) {
      for ( i=0; i<10; i++ )
        cnts[i] = 0;
	  c1 = gvect[l1]->chrom;
	  c2 = gvect[l2]->chrom;
	  m1 = gvect[l1]->markr;
	  m2 = gvect[l2]->markr;
	  for ( i=1; i<=theparams->nn; i++ ) {
	     if ( individs[i].markers[c1][m1] > -3   && individs[i].markers[c2][m2] > -3 ) {  /* if there is data */
	       cnts[0] +=1;
	       if ( individs[i].markers[c1][m1] > 0  ) {
		       if ( individs[i].markers[c2][m2] >0  )          cnts[1] +=1;		       
		       else if ( individs[i].markers[c2][m2] == 0 )    cnts[2] +=1;		       
		       else if ( individs[i].markers[c2][m2] <0 )      cnts[3] +=1;		       
	       }
	       else if ( individs[i].markers[c1][m1] == 0 ) {
		       if ( individs[i].markers[c2][m2] >0  )          cnts[4] +=1;		       
		       else if ( individs[i].markers[c2][m2] == 0 )    cnts[5] +=1;		       
		       else if ( individs[i].markers[c2][m2] <0 )      cnts[6] +=1;		       
	       }
	       else if ( individs[i].markers[c1][m1] <0 ) {
		       if ( individs[i].markers[c2][m2] >0  )          cnts[7] +=1;		       
		       else if ( individs[i].markers[c2][m2] == 0 )    cnts[8] +=1;		       
		       else if ( individs[i].markers[c2][m2] <0 )      cnts[9] +=1;		       
	       }	       
	     } 	  
	  }
  }
  if ( cnts[0] == 0 ) {
    rec = (FPN)0.5;
    *lr = ((FPN)theparams->nn) * (FPN)log(rec);
    return(rec);  
  }
/*  Next, calculate the recombination frequency */
  if ( theparams->cross == 1    ) {
    rec = (FPN) (cnts[2]+cnts[4]) / (FPN) cnts[0] ; 
    *lr = ((FPN) (cnts[1]+cnts[5]))*(FPN)log(1.0-rec) + ((FPN)(cnts[2]+cnts[4])*(FPN)log(rec));  
  }
  else if ( theparams->cross == 2  ) {
    rec = (FPN) (cnts[6]+cnts[8]) / (FPN) cnts[0] ;   
    *lr = ((FPN) (cnts[9]+cnts[5]))*(FPN)log(1.0-rec) + ((FPN)(cnts[6]+cnts[8])*(FPN)log(rec));  
  }
  else if ( theparams->cross == 5  ) {
    if ( theparams->crosst == 0 ) {
      rec = (FPN) (cnts[3]+cnts[7]) / (FPN) cnts[0] ;   
      *lr = ((FPN) (cnts[9]+cnts[1]))*(FPN)log(1.0-rec) + ((FPN)(cnts[3]+cnts[7]))*(FPN)log(rec);  
    }
    else if ( theparams->crosst == 1 ) {
      rec = (FPN)0.5 * (FPN) (cnts[3]+cnts[7]) / (FPN) (cnts[1]+cnts[9]) ;   
      *lr = ((FPN) (cnts[9]+cnts[1]))*(FPN)log(0.5*(1.0/(1.0+2.0*rec))) + ((FPN) (cnts[3]+cnts[7]))*(FPN)log(rec/(1.0+2.0*rec));  
    }
    else if ( theparams->crosst == 2 ) {  
      rec = (FPN) (cnts[3]+cnts[7]) / ((FPN) (4 * (cnts[1]+cnts[9]) - 2 * (cnts[3]+cnts[7]))) ;   
      *lr = ((FPN) (cnts[9]+cnts[1]))*(FPN)log(0.5*((1.0+2.0*rec)/(1.0+6.0*rec))) + ((FPN) (cnts[3]+cnts[7]))*(FPN)log(2.0*rec/(1.0+6.0*rec));  
    }
  }
  else if ( theparams->cross == 3 || theparams->cross == 4 ) {
    *lr = estFXrec(theparams->cross,theparams->crosst,cnts,&rec, theparams->themap->types[c1][m1],  theparams->themap->types[c2][m2]);
    if ( debugging > 1 && rec < 0.00006 ) {
      sprintf(gbuffer,"\nNew: %s %s %3d %3d |%3d %3d %3d %3d %3d %3d %3d %3d %3d| %3d %3d %8.6f %8.2f",individs[1].map->names[individs[1].map->ttable[c1][m1]],individs[1].map->names[individs[1].map->ttable[c2][m2]],l1,l2,cnts[1],cnts[2],cnts[3],cnts[4],cnts[5],cnts[6],cnts[7],cnts[8],cnts[9],individs[1].map->types[c1][m1],individs[1].map->types[c2][m2],rec,*lr);
      LogTheError(theparams->error,gbuffer);
    }

  }
  else {
    rec = (FPN)0.5;
    *lr = ((FPN) cnts[0]) * (FPN)log(rec);
  }

  return(rec);
/* Return the recombination fraction  */
}


/*
   Find the mle of r for an F2 via iteration.   The likelihood function is different
   depending upon whether there are dominant markers.  
   
   This will return a value in (0.0,0.5).   It will never return 0.0 or 1/2.  
*/
FPN estFXrec(int cross,int generation,int *cnts ,FPN *r,int l1, int l2) {
  int i,go_on;
  FPN dcnts[10],  delta,low,high,rold,rnew,lnlnew,lnl;
  for ( i=1; i<10; i++ )
    dcnts[i] = (FPN) cnts[i]; 

  delta =(FPN) 0.25;
  rold = (FPN)0.25;   
  for ( i=1; i<=5; i++ ) {    
	  low = rold - delta;
	  if ( low <(FPN) 0.0 )
	    low = (FPN)0.0;
	  high = rold + delta;
	  if ( high > (FPN)0.5 )
	    high =(FPN) 0.5;
	  delta = delta *(FPN) 0.1;
	  rnew = rold = low + delta;
	  lnl = fxlike(cross,generation,l1,l2,dcnts,rold);
	  go_on = 1; 
	  while ( go_on == 1 ) {
	    rnew = rnew + delta;
	    if ( rnew < high) { 
	      lnlnew =  fxlike(cross,generation,l1,l2,dcnts,rnew);
	      if ( lnlnew > lnl ) {
	        lnl = lnlnew;
	        rold = rnew;
	      }
	      else
	        go_on = 0;
	    }
	    else 
	      go_on = 0;
	  } 
  }
  
  *r = rold; 
  

  return(lnl);
}

/*
  This calculates the log of the likelihood for the recombination fraction 

let   p=(1-r)/2
      q=r/2
----------------------------
----------------------------      
          n(i)     P(i)
GT      count    P(GT|r)
----------------------------
AABB    in1       p^2 
AABb    in2      2pq
AAbb    in3       q^2
AaBB    in4      2pq
AaBb    in5      (1-2r+2r^2)
Aabb    in6      2pq
aaBB    in7       q^2
aaBb    in8      2pq
aabb    in9       p^2
----------------------------
----------------------------

L(r) ~ Prod( P(i) ^ n(i) )

l(r) = ln(L(r)) = Sum n(i) * ln(P(i))

   
*/
FPN fxlike(int cross,int generation, int l1, int l2, FPN *dcnts , FPN inr) {
  FPN lnl;
  FPN   *m1[11],*m2[11],*m3[11]; 
  
  FPN thematrix[40][11];
  
  FPN r;
  int i,j;
  
  for ( i=0; i<=10; i++ ) {
    m1[i] = thematrix[i];
    m2[i] = thematrix[i+11];
    m3[i] = thematrix[i+22];
 } 
  
  if ( cross == 3 || generation == 2 )
    r = inr;
  else
    r = (FPN) 0.5 - (FPN) 0.5 * (FPN) pow( (1.0-inr), (FPN) (generation-2)) *((FPN)1.0-(FPN)2.0*inr) ;
/*  m1 = dmatrix(0,10,0,10);
  m2 = dmatrix(0,10,0,10);
  m3 = dmatrix(0,10,0,10);*/
  for ( i=0; i<=10; i++ ) 
    for ( j=0; j<=10; j++ )
      m1[i][j] = m2[i][j] =m3[i][j] =(FPN)0.0;
  m1[0][0] = m1[9][9] = m1[2][2] = m1[7][7] =(FPN)1.0;
  m1[1][0] = m1[1][2] = m1[3][0] = m1[3][7] = m1[6][2] = m1[6][9] = m1[8][7] = m1[8][9] = (FPN)0.25;
  m1[1][1] = m1[3][3] = m1[6][6] = m1[8][8] =(FPN)0.5;
  m1[4][1] = m1[4][3] = m1[5][1] = m1[5][3] = m1[4][6] = m1[4][8] = m1[5][6] = m1[5][8] = (FPN)0.5*((FPN)1.0-r)*r;
  m1[4][0] = m1[4][9] = m1[5][2] = m1[5][7] =(FPN) 0.25*((FPN)1.0-r)*((FPN)1.0-r);
  m1[4][2] = m1[4][7] = m1[5][0] = m1[5][9] = (FPN)0.25*r*r;
  m1[4][5] = m1[5][4] = (FPN)0.5*r*r;
  m1[4][4] = m1[5][5] = (FPN)0.5*((FPN)1.0-r)*((FPN)1.0-r);

  for ( i=10; i>0; i-- ) 
    for ( j=10; j>0; j-- )
      m1[i][j] = m1[i-1][j-1] ;

  
    
  i=dtranspose(  m1 , m2  , 1,1,10,10);
  for ( i=2; i< generation; i++ ) {
    m1dotm2(m1, m2, m3, 10);
    j=dtranspose(  m3, m2, 1,1,10,10);  
  }
   for ( i=0; i<=10; i++ ) 
    m3[0][i] = m3[1][i] = (FPN)0.0;
  m3[0][5] =(FPN) 1.0;
  
  for ( i=1; i<=10; i++ )
    m3[1][i] = sdot(10, m3[0],1,m2[i],1 );
  
  m3[1][5] = m3[1][5]+m3[1][6];
  for (i=6; i<10; i++ )
    m3[1][i] = m3[1][i+1];    
  if ( l1==0 && l2==0 )
    lnl = dcnts[1]*(FPN)log(m3[1][1]) + dcnts[2]*(FPN)log(m3[1][2]) + dcnts[3]*(FPN)log(m3[1][3]) + dcnts[4]*(FPN)log(m3[1][4]) + dcnts[5]*(FPN)log(m3[1][5]) + dcnts[6]*(FPN)log(m3[1][6]) + dcnts[7]*(FPN)log(m3[1][7]) +dcnts[8]*(FPN)log(m3[1][8]) + dcnts[9]*(FPN)log(m3[1][9]);
  else if ( l1==1 && l2==0   ) 
    lnl = (dcnts[1]+dcnts[4])*(FPN) log(m3[1][1]+m3[1][4]) + (dcnts[2]+dcnts[5])*(FPN) log(m3[1][2]+m3[1][5]) + (dcnts[3]+dcnts[6])*(FPN) log(m3[1][3]+m3[1][6]) + dcnts[7]*(FPN) log(m3[1][7]) +dcnts[8]*(FPN) log(m3[1][8]) + dcnts[9]*(FPN) log(m3[1][9]);
  else if ( l1==0 && l2==1   ) 
    lnl = (dcnts[1]+dcnts[2])*(FPN) log(m3[1][1]+m3[1][2]) + dcnts[3]*(FPN) log(m3[1][3]) + (dcnts[4]+dcnts[5])*(FPN) log(m3[1][4]+m3[1][5]) + dcnts[6]*(FPN) log(m3[1][6]) + (dcnts[7]+dcnts[8])*(FPN) log(m3[1][7]+m3[1][8]) + dcnts[9]*(FPN) log(m3[1][9]);
  else if ( l1==-1 && l2==0  ) 
    lnl = dcnts[1]*(FPN) log(m3[1][1]) + dcnts[2]*(FPN) log(m3[1][2]) + dcnts[3]*(FPN) log(m3[1][3]) + (dcnts[4]+dcnts[7])*(FPN) log(m3[1][4]+m3[1][7]) + (dcnts[5]+dcnts[8])*(FPN) log(m3[1][5]+m3[1][8]) + (dcnts[6]+dcnts[9])*(FPN) log(m3[1][6]+m3[1][9]) ;
  else if ( l1==0 && l2==-1  ) 
    lnl = dcnts[1]*(FPN) log(m3[1][1]) + (dcnts[2]+dcnts[3])*(FPN) log(m3[1][2]+m3[1][3]) + dcnts[4]*(FPN) log(m3[1][4]) + (dcnts[5]+dcnts[6])*(FPN) log(m3[1][5]+m3[1][6]) +(dcnts[7]+dcnts[8])*(FPN) log(m3[1][7]+m3[1][8])+ dcnts[9]*(FPN) log(m3[1][9]) ;
  else if ( l1==1 && l2==1   ) 
    lnl = (dcnts[1]+dcnts[2]+dcnts[4]+dcnts[5])*(FPN) log(m3[1][1]+m3[1][2]+m3[1][4]+m3[1][5]) + (dcnts[3]+dcnts[6])*(FPN) log(m3[1][3]+m3[1][6]) + (dcnts[7]+dcnts[8])*(FPN) log(m3[1][7]+m3[1][8]) + dcnts[9]*(FPN) log(m3[1][9]);
  else if ( l1==-1 && l2==-1 ) 
    lnl = dcnts[1]*(FPN) log(m3[1][1]) + (dcnts[2]+dcnts[3])*(FPN) log(m3[1][2]+m3[1][3]) + (dcnts[4]+dcnts[7])*(FPN) log(m3[1][4]+m3[1][7]) + (dcnts[5]+dcnts[6]+dcnts[8]+dcnts[9])*(FPN) log(m3[1][5]+m3[1][6]+m3[1][8]+m3[1][9])  ;
  else if ( l1==1 && l2==-1  ) 
    lnl = (dcnts[1]+dcnts[4])*(FPN) log(m3[1][1]+m3[1][4]) + (dcnts[2]+dcnts[3]+dcnts[5]+dcnts[6])*(FPN) log(m3[1][2]+m3[1][3]+m3[1][5]+m3[1][6])+ dcnts[7]*(FPN) log(m3[1][7]) +(dcnts[8]+dcnts[9])*(FPN) log(m3[1][8]+m3[1][9]);
  else if ( l1==-1 && l2==1  ) 
    lnl = (dcnts[1]+dcnts[2])*(FPN) log(m3[1][1]+m3[1][2]) + dcnts[3]*(FPN) log(m3[1][3]) + (dcnts[4]+dcnts[5]+dcnts[7]+dcnts[8])*(FPN) log(m3[1][4]+m3[1][5]+m3[1][7]+m3[1][8]) + (dcnts[6]+dcnts[9])*(FPN) log(m3[1][6]+m3[1][9]);

  return(lnl);
}

/*
  This calculates the (FPN) log of the likelihood for the recombination fraction for advanced backcrosses
  
  
  DOES THIS WORK?????   It is not yet being used.   

            
   generation is the number of repeated backcrosses, and should be greater than 1.  
   
   If BC1x, then   input (n1, n2, n3, n4) as ( nAABB, nAABb, nAaBB, nAaBb )
      BC2x, then   input (n1, n2, n3, n4) as ( naabb, nAabb, naaBb, nAaBb ) 
      
   
      
*/
   
FPN bxlike(int generation, FPN *dcnts,  FPN  r) {
  FPN lnl ;
  int i,j;


  FPN   *m1[11],*m2[11],*m3[11]; 
  
  FPN thematrix[40][11];
  
   
  for ( i=0; i<=10; i++ ) {
    m1[i] = thematrix[i];
    m2[i] = thematrix[i+11];
    m3[i] = thematrix[i+22];
 } 

 
   for ( i=0; i<=10; i++ ) 
    for ( j=0; j<=10; j++ )
      m1[i][j] = m2[i][j] =m3[i][j] =(FPN)0.0;
    m1[1][1] = (FPN)1.0;
    m1[1][2] = m1[1][3] = (FPN)0.5*((FPN)1.0+r);
    m1[1][4] = m1[4][4] = (FPN)0.5*((FPN)1.0-r);
    m1[2][2] = m1[3][3] = (FPN)0.5;
    m1[2][4] = m1[3][4] = (FPN)0.5 * r;
 
  
    
  i=dtranspose( m1, m2, 1,1,4,4);
  for ( i=2; i< generation; i++ ) {
    m1dotm2(m1, m2, m3, 4);
    j=dtranspose(  m3, m2, 1,1,4,4);  
  }
  for ( i=0; i<=4; i++ ) 
    m3[0][i] = m3[1][i] = (FPN)0.0;
  m3[0][4] = (FPN)1.0;
  
  for ( i=1; i<=4; i++ )
    m3[1][i] = sdot(4, m3[0],1,m2[i],1 );
  
  lnl = dcnts[1]*(FPN) log(m3[1][1]) + dcnts[2]*(FPN) log(m3[1][2]) + dcnts[3]*(FPN) log(m3[1][3]) + dcnts[4]*(FPN) log(m3[1][4]);

  return(lnl);
}

/*
    For matrices m1 and m2, this does m3 = m1.m2'  
    If you want m3 = m1.m2, transpose m2 first.   
*/
void m1dotm2(FPN **m1,FPN **m2,FPN **m3,  int ub) {
  int i,j;
  for ( i=1; i<=ub; i++ )
    for ( j=1; j<=ub; j++ )
      m3[i][j] = sdot(ub , m1[i], 1, m2[j], 1 );

}


/*
Calculate the objective function for the entire map.
*/

FPN Mapobj(params *theparams, int glen,genome **gvect) {
  int i;
  FPN thisobj;


 thisobj =(FPN) 0.0;
 for ( i=1; i<=theparams->themap->m; i++ )
   thisobj = thisobj + obj_func(theparams,i,glen,gvect);

 return(thisobj);
}


/*
Let's see the recombination estimates.

*/


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Efunc.c
------------------------------------------------------------------ */


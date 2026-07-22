/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Genome.c) is part of QTL Cartographer
         
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


/* Genome.c,  subroutines to manage a linked list that is an
          abstraction of the genome.*/

#include "Main.h"

/*
  return a pointer to the first node in a genome chain, starting from
  any point in the chain.
*/
genome *firstone(genome *gptr) {
  genome *thegenome;
  int k;
  k=0;
  for ( thegenome=gptr; thegenome->prev != NULL ; thegenome=thegenome->prev) 
    k+=1;     
  return(thegenome);
}



/*
  This takes a genome structure shows what is in it.  Used mainly for debugging.
*/
void show_genome(genome *thegenome,  char *outfile ) {
  genome *gptr;
  FILE *fptr;
  
  if ( outfile[0] == '-')
    fptr = stdout;
  else 
    fptr = fileopen(outfile, "a");
  fprintf(fptr,"\n\n\tHere is the genome: ");   
  putline(fptr,'-',50);
  putline(fptr,'-',50);
  fprintf(fptr,"\nMarkername Chr Mrk Distan Positi QTL mxo    pxo ");   
  putline(fptr,'-',50);

  for (gptr=thegenome; gptr!=NULL; gptr=gptr->next )
    fprintf(fptr,"\n%10s %3d %3d %6.3f %6.3f %3d %6.3f %6.3f",gptr->markername,gptr->chrom, gptr->markr,gptr->dist, gptr->pos,gptr->whichqtl, gptr->pxo, gptr->mxo);
  putline(fptr,'-',50);
  putline(fptr,'-',50);
  
  fprintf(fptr,"\n\n");   
  if (outfile[0] != '-')
      fileclose(outfile, fptr);
}



/* allocate space for a genome node and set elements to zero or null */
genome *alloc_genome(void) 
{
  genome *gptr;
#if defined(MACWARRIOR) || defined(WINWARRIOR)  
  gptr = (genome *) malloc( (size_t) sizeof(genome));
#else
  gptr = (genome *) malloc( (unsigned) sizeof(genome));
#endif
 if ( debugging > 2 ) {
        sprintf(gwarn,"In alloc_genome(), allocated 1 genome node at %p\n",gptr);
        MemoryAccount(gwarn);
 }

  gptr->chrom=0;     /* Chromosome of the marker interval */
  gptr->markr=0;     /* The marker that precedes the interval, can be 0 if borders are allowed */
  gptr->dist=(FPN)0.0;    /* This distance, from marker markr to markr+1 on chromosome chrm, is in M */
  gptr->pos=(FPN)0.0;     /* Postion from left telomere in Morgans */
  gptr->whichqtl=0;  /* Indicator of whether there is a qtl on the interval 1 if yes, 0 if no */
  gptr->mxo=(FPN)0.0;     /* xo position of maternal    xo = 0.0 if no crossover,   */
  gptr->pxo=(FPN)0.0;     /*             of paternal    xo > 0.0 the distance (M) from the marker to the xo */
  gptr->markername = NULL; /*Only used in MultiRegress*/
  gptr->n = 0;            /*sample size for use in MultiRegress, number of markers in Emap*/
  gptr->values = NULL; /*  Only used in MultiRegress and Emap*/
  gptr->prev=NULL;   /* Pointer to previous marker interval */
  gptr->next=NULL;   /* Pointer to next marker interval */
  return(gptr);
}

void cleanse_genome(genome *gptr)
{
  genome *tgptr;
  
  for (tgptr = gptr ; tgptr != NULL ; tgptr = tgptr->next ) {
    tgptr->whichqtl = 0;
  
  }

}

/*
 * This 'creates' a genome from the linkage map.  Each interval is a node in
 * a chain.  It has a field for its chromosome and the marker that it
 * follows.  A marker of 0 indicates that it is the genetic material before
 * the first marker on a chromosome, ie, tail material.  The data structure
 * also has the length of the interval in Morgans (M).
 */

genome *create_genome( markermap *themap)
{
  genome *first, *last, *gptr;
  int tail, done = 0, chr, mrk, endofchrom;
  first = alloc_genome();
  if (themap->brdrs > 0.0)
    tail = 0;
  else
    tail = 1;
  chr = 1;
  mrk = tail;
  first->chrom = chr;
  first->markr = mrk;
  first->pos = first->dist = mapfunc(themap->mrf[chr][ mrk], 1);
  gptr = first;
    endofchrom = themap->mpc[chr];
  while (done == 0) {
    mrk = mrk + 1;
    if (mrk > endofchrom) {
      chr = chr + 1;
      if (  chr <= themap->m)
	    endofchrom = themap->mpc[chr];
      mrk = tail;
    }
    if (chr > themap->m) {
      done = 1;
      gptr->next = NULL;
    }
    else {
      gptr->next = alloc_genome();
      gptr->next->prev = gptr;
      gptr = gptr->next;
      gptr->chrom = chr;
      gptr->markr = mrk;
	  gptr->pos = gptr->dist = mapfunc(themap->mrf[chr][mrk], 1);
    }
  }

  gptr = first;
    do {    /* Put nextmarker's pos.in gptr->pos */
      if (gptr->prev != NULL && gptr->prev->chrom == gptr->chrom)
	    gptr->pos = gptr->pos + gptr->prev->pos;
      if (gptr->next == NULL)
	    last = gptr;
      gptr = gptr->next;
    } while (gptr != NULL);
    gptr = last;
    do {	/* put marker's position in gptr->pos */
      if (gptr->prev == NULL)
	    gptr->pos = (FPN)0.0;
      else {
	    gptr->pos = gptr->prev->pos;
	    if (gptr->chrom < last->chrom )
	      if ( gptr->chrom != gptr->next->chrom )
		     gptr->next->pos = (FPN)0.0;
      }
      gptr = gptr->prev;
	} while (gptr != NULL);
  return (first);	/* return a pointer to the first node in the chain. */
}


/*Delete a genome node...move it to the front, delete it, return the next node.*/
genome *del_genome_node(genome *gptr) {
  genome *lgptr, *rgptr ;
  if ( gptr == NULL )
    return(NULL);
  lgptr = mv_genome_front(gptr);
  rgptr = lgptr->next;
  if (rgptr != NULL )
    rgptr->prev = NULL;
  
     if ( debugging > 2 ) {
        sprintf(gwarn,"In del_genome_node(), deallocated 1 genome node at %p\n",lgptr);
        MemoryAccount(gwarn);
     }
     if ( lgptr->markername != NULL )
       free_cvector(lgptr->markername,0,MAXNAME);
     if ( lgptr->values != NULL )
       free_dvector(lgptr->values,1,lgptr->n);
	 free((char *) lgptr);
  return(rgptr);

}

/*Put marker names onto a genoome*/
void  addmarkernames(markermap *themap,genome *thegenome) {
  genome *gptr;
  for (gptr=thegenome; gptr!=NULL; gptr=gptr->next) {
    gptr->markername = cvector(0,MAXNAME);
    if ( gptr->markr > 0 ) 
      strcpy(gptr->markername, themap->names[ themap->ttable[ gptr->chrom][gptr->markr] ] );
    else
      sprintf(gptr->markername, "LeftTelomere%d",gptr->chrom);
  }
      

}


/* move gptr to the beginning */
genome *mv_genome_front(genome *gptr)
{
  genome *lgptr;
  if ( gptr->prev == NULL )/* If at the front, leave it. */
    return(gptr);
  lgptr=gptr;
  while ( lgptr->prev != NULL )  /* find first node */
    lgptr = lgptr->prev;
  gptr->prev->next = gptr->next;
  if ( gptr->next != NULL )
    gptr->next->prev = gptr->prev;
  lgptr->prev = gptr;
  gptr->next = lgptr;
  gptr->prev = NULL;
  return(gptr);
}

/* move gptr to the left of rptr.
   if rptr == NULL, then put gptr at the end
 */
void mv_genome_node(genome *gptr, genome *rptr)
{
  genome *first, *last;
  int k;
  k=0;
  for ( first=gptr; first->prev != NULL ; first = first->prev ) k +=1;
  for ( last=gptr; last->next != NULL ; last = last->next ) k -=1;
  
  
  
  if ( gptr->prev == NULL ) {  /* If at the front, make next pointer anchor */
    gptr->next->prev = NULL;
    first = gptr->next;
  }
  else if ( gptr->next == NULL ) { /*  if last, make prev. pointer terminus */  
     gptr->prev->next = NULL;
     last = gptr->prev;  /*  This had been pointing at NULL.*/
  }
  else {  /*  pull gptr and close chain */
    gptr->next->prev = gptr->prev;
    gptr->prev->next = gptr->next;
  }
  if ( rptr == NULL ) {/*  rptr == NULL means put it on the end */
    gptr->prev = last;
    gptr->prev->next = gptr;
    gptr->next = NULL;
  }    
  else {
    gptr->next = rptr;
    gptr->prev = rptr->prev;
    rptr->prev = gptr;
    if ( gptr->prev != NULL )
      gptr->prev->next = gptr;  
  }

}


/*
 * this goes down the chain of nodes that makes up the genome, and frees up
 * the memory.
 */

void clear_genome(genome *gptr)
{
  genome *gtptr;
  int k;
  k=0;
  gtptr = gptr;
  if ( gtptr != NULL ) 
    while ( (gtptr = del_genome_node(gtptr)) != NULL ) k+=1;
  
  gptr = NULL;
}


/*  eliminate one node in the genome.  */
genome *elim_gnode(genome *gptr)
{
  genome *rptr;

  if (gptr->prev == NULL && gptr->next == NULL)
    rptr = gptr;	/* gptr is only node so return it... */
  else if (gptr->prev == NULL) {
    rptr = gptr->next;	/* gptr is the left terminator, so return next node */
    rptr->prev = NULL;
    free((char *) gptr);
  }
  else if (gptr->next == NULL) {
    rptr = gptr->prev;	/* gptr is the right terminator, so return previous node */
    rptr->next = NULL;
    free((char *) gptr);
  }
  else {
    rptr = gptr->next;	/* gptr is in the middle, so return the next node */
    gptr->next->prev = gptr->prev;
    gptr->prev->next = gptr->next;
    free((char *) gptr);
  }
  return (rptr);
}





/*
 * Create a vector of aqtls that will hold the necessary information on the
 * QTLs
 */

aqtl *qtlvector(markermap *themap)
{
  int ii, kk;
  aqtl *qtlptr;
  kk = 0;
  for (ii = 1; ii <= themap->traits; ii++)
    kk = kk + themap->knum[ii];
#if defined(MACWARRIOR) || defined(WINWARRIOR)  
  qtlptr = (aqtl *) malloc((size_t) kk * sizeof(aqtl));
#else
  qtlptr = (aqtl *) malloc((unsigned) kk * sizeof(aqtl));
#endif
  if ( debugging > 2 ) {
    sprintf(gwarn,"In qtlvector(), allocated %d qlts at %p\n",kk,qtlptr);
    MemoryAccount(gwarn);
  }
  for (ii = 0; ii < kk; ii++) {
    qtlptr[ii].episAADD = NULL;
    qtlptr[ii].episADDA = NULL;
    CleanQTL( &qtlptr[ii], kk, themap);
    qtlptr[ii].nptrs = 1;
  }
  
  return qtlptr - 1;
}

/*
  How many nodes in chain?
*/
int genomelength(genome *gptr) {
  genome *tgptr ;
  int counter;
  
  counter = 0;
  for ( tgptr=gptr; tgptr !=NULL; tgptr=tgptr->next)
    counter +=1;
   return(counter);

}
/*
  Given a genome, send back a vector pointing to the nodes.
*/
genome **vectorizegenome(genome *gptr) {
  genome *tgptr,**retgptr;
  int counter;
  
  counter = genomelength(gptr);

  retgptr = genomevector(1,counter);
  counter = 0;
  for ( tgptr=gptr; tgptr !=NULL; tgptr=tgptr->next) {
    counter +=1;
    retgptr[counter] = tgptr;
  }
  return(retgptr);

}

/*   Create a vector of pointers to genome nodes*/
genome **genomevector(int nrl, int nrh )
{
  int i;
  genome **m;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  m = (genome **) malloc((size_t) (nrh - nrl + 1) * sizeof(genome *));
#else
  m = (genome **) malloc((unsigned) (nrh - nrl + 1) * sizeof(genome *));
#endif
  
  if (!m)
    nrerror("allocation failure 1 in genomevector()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)  
    m[i] = NULL;
  if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In genomevector(), allocated %d genome vector of pointers  at %p\n",nrh-nrl+1, m+nrl);
    MemoryAccount(gwarn);
  } 
  return m;
}
void free_genomevector(genome **m, int nrl, int nrh )
{
  
  if ( debugging > 2 ) {
/*  */
    sprintf(gwarn,"In free_lsvector(), deallocated %d long vector of pointers at %p\n",nrh-nrl+1,m+nrl);
    MemoryAccount(gwarn);
  } 
  free((char *) (m + nrl));
}



void CleanQTL(aqtl *qtlptr,int kk, markermap *themap) {

    qtlptr->chrm = 0;
    qtlptr->mrk = 0;
    qtlptr->c1 = (FPN)0.0;
    qtlptr->c2 = (FPN)0.0;
    qtlptr->a =(FPN) 1.0 / (FPN) kk;
    qtlptr->d = (FPN)0.0;
    qtlptr->map = themap;
    qtlptr->r2 = (FPN)0.0;
    qtlptr->tr2= (FPN)0.0;
    qtlptr->s = (FPN)0.0;

}

/*
 * Free up the space used by the qtl vector
 */

void free_qtlvector(aqtl *qtlptr, markermap *themap)

{
  int ii,total,prev;
  total = 0;
  for (ii = 1; ii <= themap->traits; ii++) 
    total = total + themap->knum[ii];
  prev = 0; 
  for (ii = 1; ii <= total; ii++) 
    if ( qtlptr[ii].trait != prev ) {
      if ( qtlptr[ii].episAADD != NULL ) {
        free_dmatrix(qtlptr[ii].episAADD,1,themap->knum[ qtlptr[ii].trait],1,themap->knum[ qtlptr[ii].trait]);
		qtlptr[ii].episAADD = NULL;
	  }
      if ( qtlptr[ii].episADDA != NULL ) {
         free_dmatrix(qtlptr[ii].episADDA,1,themap->knum[ qtlptr[ii].trait],1,themap->knum[ qtlptr[ii].trait]);
		 qtlptr[ii].episADDA = NULL;
	  }
      prev =  qtlptr[ii].trait;   
    }
  if ( debugging > 2 ) {
     sprintf(gwarn,"In free_qtlvector(), deallocated  qlts at %p\n", qtlptr+1);
     MemoryAccount(gwarn);
  }
  free((char *) &qtlptr[1]);
  qtlptr = NULL;
}



/*
 * Indicate whether 1 or more QTLs reside in each interval.  Use this only after
 * QTLs have been created, and a new genome linked list has been created.
 */

void place_qtls(params *theparams,genome *first)
{
  int kk, knum;
  genome *gptr;
  knum = 0;
  for (kk = 1; kk <= theparams->themap->traits; kk++)
    knum = knum + theparams->themap->knum[kk];
  for (kk = 1; kk <= knum; kk++) {
    gptr = first;
    while (gptr != NULL)
      if (gptr->chrom == theparams->theqtls[kk].chrm && gptr->markr == theparams->theqtls[kk].mrk) {
	    gptr->whichqtl = 1;
	    gptr = NULL;
      }
      else
	gptr = gptr->next;
  }
}

/* For a QTL in qtlptr, go to the place in the genome where it resides.  */
genome *aqtl_genomenode( genome *gptr, aqtl *qtlptr) {
  genome *tgptr; 
  int go_on;
  
  go_on = 1;
  tgptr=gptr;
  while ( go_on == 1 ) {  /*go to the correct place in the genome list */
    if ( tgptr->chrom == qtlptr->chrm && tgptr->markr == qtlptr->mrk ) 
      go_on = 0;
    else
       tgptr = tgptr->next;
    if (tgptr == NULL)
        go_on = 0;
  }
  return(tgptr); 
}

/*
   Swap g1 and g2
*/
int swap_gnodes( genome *g1ptr, genome *g2ptr) {
  genome *g1r,*g1l,*g2r,*g2l;
  int left,right;
  
  if ( g1ptr == NULL || g2ptr == NULL ) /*  don't want to swap one or more NULLs */
    return(0);
  if ( g1ptr == g2ptr )  /*  nothing to do if they are the same*/
    return(0);
  
  right = left = 0;
  
    
  if ( g2ptr ==  g1ptr->next )
     left = 1;
  if ( g1ptr == g2ptr->next )
     right = 1; 
     
  g1l = g1ptr->prev;
  g1r = g1ptr->next;
  g2l = g2ptr->prev;
  g2r = g2ptr->next;

  if ( left == 1 ) {  /* g1ptr is to the immediate left of g1ptr */
    g1ptr->prev = g2ptr;
    g2ptr->next = g1ptr;
    g2ptr->prev = g1l;
    g1ptr->next = g2r;
    if ( g1l != NULL )
      g1l->next = g2ptr;
    if ( g2r != NULL )
      g2r->prev = g1ptr;
  }
  else if ( right == 1 ) { /* g1ptr is to the immediate right of g1ptr */
    g1ptr->next = g2ptr;
    g2ptr->prev = g1ptr;
    g2ptr->next = g1r;
    g1ptr->prev = g2l;
    if ( g1r != NULL )
      g1r->prev = g2ptr;
    if ( g2l != NULL )
      g2l->next = g1ptr;
  
  }  
  else {  /* g1ptr is not adjacent to g2ptr */
    /*  Reset g1ptr*/
	  g1ptr->prev = g2l;
	  if ( g2l != NULL ) 
	    g2l->next = g1ptr;
	  g1ptr->next = g2r;
	  if ( g2r != NULL ) 
	    g2r->prev = g1ptr;

	/* Reset g2ptr*/
	  g2ptr->prev = g1l;
	  if ( g1l != NULL ) 
	    g1l->next = g2ptr;
	  g2ptr->next = g1r;
	  if ( g1r != NULL ) 
	    g1r->prev = g2ptr;
  }  
  


  return(1);
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Genome.c
------------------------------------------------------------------ */


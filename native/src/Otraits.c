/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Otraits.c) is part of QTL Cartographer
         
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
/*  Functions to process and analyze categorical variables. These are called
    other traits.    */

/*
  Go through all the other traits and determine how many different values
  each one has (sort of like determining the number of alleles for each
  of the Other or categorical traits.  
  
  I suppose it is better to think of these as Categorical traits, and this
  function determines the number of categories for each trait.
  
  After this is done, each individual will have an integer value for the value
  of each of their categorical traits.
*/
void process_otraits(individual *iptr,int n,markermap *themap)
{
  int ii,tr,k;
  otnode *first,*optr;
  k=0;
  if ( themap->otypes == NULL ) 
    themap->otypes = ivector(1, themap->otraits);
  for ( ii = 1 ; ii <= n ; ii++ )
    iptr[ii].oyt = ivector(1,themap->otraits);

  for ( tr = 1 ; tr <= themap->otraits ; tr++ ) {
    first = alloc_otnode(NULL,NULL);
    strcpy( first->name, ".");
    for ( ii = 1 ; ii <= n ; ii++ )
      iptr[ii].oyt[tr] = search_otnodes(&iptr[ii],tr,first);
    for ( optr = first ; optr->next != NULL ; optr = optr->next ) k+=1;
    themap->otypes[tr] = optr->which;
    dealloc_otnode(first);
  }
}


/* 
  This is used when determining how many different values the Other trait
  can have.  It searches through the linked list of types already seen.  If
  it finds one, it returns which value it is (sort of like which allele).
  If this individual has an Other trait not yet seen, a new node is created
  and the number of 'alleles' increases by one.
*/
int search_otnodes(individual *iptr,int tr,otnode *optr)
{
  int found;
  otnode *ptr;
  found = 1;
  ptr = optr; 
  do {
    if ( !strcmp( ptr->name, iptr->oy[tr] ) )
      return(ptr->which);
    if ( ptr->next != NULL )
      ptr = ptr->next;
    else
      found = 0;
  }  while ( found == 1 );
  ptr = alloc_otnode(ptr,NULL);
  ptr->prev->next = ptr;
  ptr->which = ptr->prev->which+1;
  strcpy( ptr->name, iptr->oy[tr] );
  return(ptr->which);
}


/* 
  Deallocate memory for a node that defines an Other Trait.
*/
void dealloc_otnode(otnode *optr)
{
  otnode *ptr;
  ptr = optr ;
  while ( ptr->next != NULL ) {
    ptr = ptr->next;
    free_cvector(ptr->prev->name,0,MAXNAME);
    if ( debugging > 2 ) {
        sprintf(gwarn,"In dealloc_otnode(), deallocated 1 otnode  at %p\n",ptr->prev);
        MemoryAccount(gwarn);
    }
    free((char *) ptr->prev);
  }
  free_cvector(ptr->name,0,MAXNAME);
    if ( debugging > 2 ) {
        sprintf(gwarn,"In dealloc_otnode(), deallocated 1 otnode  at %p\n",ptr);
        MemoryAccount(gwarn);
    }
  free((char *) ptr);
}


/* 
  Allocate memory for a node that defines an Other Trait.
*/
otnode *alloc_otnode(otnode *prev,otnode *next)
{
  otnode *ptr;
#if defined(MACWARRIOR) || defined(WINWARRIOR)
  ptr = (otnode *) malloc( (size_t) sizeof(otnode));
#else
  ptr = (otnode *) malloc( (unsigned) sizeof(otnode));
#endif
 if ( debugging > 2 ) {
        sprintf(gwarn,"In alloc_otnode(), allocated 1 otnode  at %p\n",ptr);
        MemoryAccount(gwarn);
 }

  ptr->name = cvector(0,MAXNAME);
  ptr->which = 0;
  ptr->prev = prev;
  ptr->next = next;
  return(ptr);
}


/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Otraits.c
------------------------------------------------------------------ */


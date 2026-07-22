/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Genome.h) is part of QTL Cartographer
         
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



genome *firstone(genome *gptr);
void show_genome(genome *thegenome,  char *outfile );
genome *elim_gnode(genome *gptr);

genome *del_genome_node(genome *gptr);
genome *create_genome(markermap *themap);
void clear_genome(genome *gptr);
aqtl *qtlvector(markermap *themap);
void place_qtls(params *theparams, genome *first);
genome *alloc_genome(void);
genome *mv_genome_front(genome *gptr);
void cleanse_genome(genome *gptr);
void free_qtlvector(aqtl *qtlptr, markermap *themap);

void CleanQTL(aqtl *qtlptr,int kk, markermap *themap);
genome *aqtl_genomenode( genome *gptr, aqtl *qtlptr);

void free_genomevector(genome **m, int nrl, int nrh );
genome **vectorizegenome(genome *gptr) ;
genome **genomevector(int nrl, int nrh );
int genomelength(genome *gptr);
void  addmarkernames(markermap *themap,genome *thegenome);
void mv_genome_node(genome *gptr, genome *rptr);
int swap_gnodes( genome *g1ptr, genome *g2ptr);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Genome.h
------------------------------------------------------------------ */


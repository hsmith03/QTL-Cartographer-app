/* ------------------------------------------------------ XCutXCodeXSkip
     This file (Mdatain.h) is part of QTL Cartographer
         
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
  These subroutines allow one to read in a map of markers in the format 
  that is specified in the output subroutines of Rmap.
*/
void GetTheMap(params *theparams, char *infile);
int read_marker_params(int *mm,int  *ll, FPN *ssigl, FPN *ss, FPN *ssigs, FPN *brdrs, char *inputf);
int read_marker_numbers(markermap *themap, char *inputf);
int what_chrom(char *lastline);
int det_if_mm(char *infile);
void read_markers(markermap *themap,char *inputf);
void deallocate_markermap(markermap *themap);
void print_map(markermap *themap, char *outfile);
void print_map_inp(markermap *themap,char *outfile);
void print_plabqtl_map(markermap *themap,FILE *outf);
void do_map_stats(markermap *themap);
void print_marker_names(markermap *themap, FILE *fptr);
void get_markernames(markermap *themap, char *inputf);
void get_traitnames(params *theparams,  markermap *themap);
markermap *trans_stdm(char *infile);
markermap *allocate_markermap(void);
markermap *trans_maps(char *infile);
markermap *initialize_markermap(int kk, int mm, int ll, FPN ssigl, FPN ss, FPN ssigs, FPN bbrdrs, char *infile);
markermap *get_markermap(int knum,char *infile);
markermap *trans_mm(int marks, char *infile, char *outfile);

void    create_bogus_otraits(params *theparams, markermap *themap);
void    create_bogus_traits(params *theparams, markermap *themap);
void    create_bogus_markers(params *theparams, markermap *themap);
markermap *mapfromgenome(genome *thegenome);
markermap *bogus_markermap(int chrom, int mark, int maxl, int traits, FPN dist);
markermap *rcrossbogusmap(char *iinfile);
markermap *crossbogusmap(char *iinfile);
void print_genome_map(genome *thegenome,params *theparams, char *outfile,char *progname,char *chptr);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file Mdatain.h
------------------------------------------------------------------ */


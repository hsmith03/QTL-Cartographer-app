/* ------------------------------------------------------ XCutXCodeXSkip
     This file (EQfunc.c) is part of QTL Cartographer
         
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


/*  Functions called by Eqtl and Preplot */

#include "Main.h"


/*

   Some new code to create a linked list of where the various results in
   the Zmapqtl output file are.  
   
   
   zresult is  a structure that tells which trait, chromosome and model.  It
   indicates where the data are.  
   
   zresult_node creates, destroys and initializes zresult nodes.
   zresult_list goes through the Zmapqtl output and collects the positions of the output.
   zresult_rank ranks the nodes and sorts them by traits, chromosomes or models.  
   
   To do:  
   
   Preplot
     1.  Get and rank results by trait, then chromosome then model
     2.  Go down chain.  If a new chromosome...
            a.  Start a new graph
            b.  if there is a qtl file, process it
            c.  if there is an lr file, process it
            d.  if there is an eqt file, process it
            e.  put in a significance level.
         Do the node.
     
     
  Eqtl  Should be an option.  will forgo all else and do this stuff.
  
    1. Rank by model, trait, chromosome
    2. Eliminate all nodes that don't have the correct model
    3. Go down chain, count number of peaks for each trait.  Put number in themap->knum
    4. Allocate qtlvector
    5. Go down nodes.  At peaks, record in qtlvector
    6. collect results for lrfile?
    
    Output format flag:
    
    1. Rqtl.out
    2. bootstrap  
    3. Hypertext  html
    4. latex2e    tex
 

struct ZmapqtlResults {
	int chrom;
	int trait;
	int model;
	FPN window;
	int nbp;
	long rank;
	long start_offset;
	long end_offset;
	struct ZmapqtlResults prev;
	struct ZmapqtlResults next;
}
*/

void calc_recomb_dist(aqtl *theqtls,markermap *themap,int numqtls)
{
  FPN rdist,ldist;
  int i,k;
  genome *gptr,*first;
  first = create_genome(themap);
  k=0;
  
  for ( i = 1 ; i <= numqtls ; i++ ) {
    for ( gptr = first ; gptr != NULL && (gptr->chrom != theqtls[i].chrm || gptr->markr != theqtls[i].mrk) ; gptr = gptr->next ) k+=1;
    ldist = theqtls[i].c1/(FPN) 100.0 - gptr->pos;
    rdist = gptr->dist - ldist;
    theqtls[i].c1 = mapfunc(ldist, -1);
    theqtls[i].c2 = mapfunc(rdist, -1);
  }
  clear_genome(first);
}
/*

Estimate the QTLs

Out              Fx        Bi,Ri
Col  estimate  Quant       Quant.
 1              c           c
 2              m           m
 3              position    position
 4   1          H0:H3       H0:H1
 5   2          H1:H3       R2(0:1)
 6   3          H2:H3       TR2(0:1)
 7   4          H1:a        H1:a 
 8   5          H3:a        S
 9   6          H2:d       
10   7          H3:d      
11   8          H0:H1      
12   9          H0:H2
13  10          R2(0:3)      
14  11          R2(1:3)      
15  12          R2(2:3)      
16  13          TR2(0:3)    
17  14          TR2(1:3)     
18  15          TR2(2:3)     
19              S1
20              S2
21              S3


*/
int znode_QTLs_estimate(zresult *znode,params *theparams,aqtl *theqtls,markermap *themap,int which)
{
  FILE *inputf;
  FPN lrvalue,newlrvalue;
  FPN odominance,oadditive,oleft_r,oright_r,or2,otr2,os;
  int ochrom,omark,assignem;   
  int error,ch,numqtl,start, lrcol,acol,dcol,updown,chrom,r2col,tr2col,scol;
  numqtl = which;
  updown = 1;
  lrvalue = newlrvalue = (FPN) 0.0;  
  assign_zfilecols(theparams,&lrcol,&acol,&dcol,&r2col,&tr2col,&scol) ;

  inputf = fileopen(theparams->zfile, "r");
  if (inputf == NULL)
    return(numqtl);
  error = fseek(inputf,znode->start_offset,SEEK_SET);
  if ( error != 0 ) {
    printf("\nProblem in znode_QTLs");
    start = 1;
    updown = 0;
  }
  else 
    start = 0;
  assignem = 0;
  while (start == 0) {
    ch =  get_next_line(gbuffer, MAXLINE, inputf);
    if (ch != EOF ) {
      if ( gbuffer[0] != '-' ) {
        get_field(1, gname, gbuffer);
        chrom = is_pinteger(gname);
        if ( chrom != znode->chrom ) {
          start = 1;
          if (updown == 1 && lrvalue > theparams->siglevel )   /*This is for LRs that max at the telomere.*/
            assignem = 1;
        }
        else {
          get_field(lrcol,gname,gbuffer);
          newlrvalue = (FPN) atof(gname);
          if ( updown == 1 && newlrvalue < lrvalue ) {
            updown = -1;
            if ( lrvalue > theparams->siglevel ) 
              assignem = 1;
          }
          else if ( updown == -1 && newlrvalue > lrvalue )
            updown = 1;
          lrvalue = newlrvalue;      
        }
      }
      else if (  gbuffer[1] == 'e'  ) {
        if ( updown == 1 && lrvalue > theparams->siglevel ) 
          assignem = 1;
        start = 1;
      }
      if (assignem == 1  ) {  
          theqtls[numqtl].chrm =  ochrom;
          theqtls[numqtl].mrk =  omark;
          theqtls[numqtl].c1 =  (FPN) 100.0 * oleft_r;
          if (theparams->lodflag == 1 )
            theqtls[numqtl].c2 =  lrtolod(oright_r);
          else
            theqtls[numqtl].c2 =  oright_r;
          if (theqtls[numqtl].c2 > theparams->maxlr)
            theparams->maxlr = theqtls[numqtl].c2;
          theqtls[numqtl].a =  oadditive;
          if (  theparams->maxeffect < (FPN) fabs(theqtls[numqtl].a)  )
            theparams->maxeffect = (FPN) fabs(theqtls[numqtl].a) ;
          if ( theparams->ngt == 3 && theparams->ihypo == 20 )
            theqtls[numqtl].a = (FPN) 0.0;
          theqtls[numqtl].d =  odominance;
          if ( theparams->ngt == 3 && theparams->ihypo == 10 )
            theqtls[numqtl].d = (FPN) 0.0;
          theqtls[numqtl].map = themap;
          theqtls[numqtl].r2 =  or2;
          theqtls[numqtl].tr2 =  otr2;
          theqtls[numqtl].s =  os;
          numqtl = numqtl+1;
          assignem = 0;
      }
      if ( gbuffer[0] != '-' ) {      
      get_field(1, gname, gbuffer);
      ochrom = atoi(gname);
      get_field(2, gname, gbuffer);
      omark = atoi(gname);
      get_field(3, gname, gbuffer);
      oleft_r = (FPN) atof(gname);
      get_field(lrcol, gname, gbuffer);
      oright_r = (FPN) atof(gname);
      get_field(acol, gname, gbuffer);
      oadditive = (FPN) atof(gname);
      get_field(dcol, gname, gbuffer);
      odominance = (FPN) atof(gname);      
      get_field(r2col, gname, gbuffer);
      or2 = (FPN) atof(gname);      
      get_field(tr2col, gname, gbuffer);
      otr2 = (FPN) atof(gname);      
      get_field(scol, gname, gbuffer);
      os = (FPN) atof(gname);     
      } 
    }
    else
      start = 1;
  }
  fileclose(theparams->zfile, inputf);
  return(numqtl);  
}

/*this routine knows which columns in the Zmapqtl.out file are which.


   3                2
1  c                c
2  m                m
3  position         position
4  H0:H3            H1:H0
5  H1:H3            R2(0:1)
6  H2:H3            TR2(0:1)
7  H1:a             a        
8  H3:a             S1
9  H2:d           
10 H3:d          
11 H0:H1          
12 H0:H2         
13 R2(0:1)        
14 R2(0:2)        
15 R2(0:3)         
16 TR2(0:1)       
17 TR2(0:2)       
18 TR2(0:3)         
19 S1             
20 S2             
21 S3 



*/
void assign_zfilecols(params *theparams, int *lrcol,int *acol,int *dcol, int *r2col, int *tr2col, int *scol) 
{
/* default for BC or Ri lines ... */  
    *lrcol = 4;	/* LR H1 v H0 */
    *acol = 7;	/* est(a) under H1 */
    *dcol = 9;	/* will be 0.0... */
    *r2col = 5;  /*r2 for H1:H0*/
    *tr2col = 6;  /*tr2 for H1:H0*/
    *scol = 8;  /*S for H1 */
/* 
     If the cross produced three marker genotypes, there are esimates for
     dominance as well.  We need to take this into account, as well as
     use the proper hypothesis tests. 
     
     The default is ihypo == 1 or 30 which means H3 v H0.  Only change if 
     directed to do so.   
*/
  if ( theparams->ngt == 3 ) {
    *r2col  = 15;  /* 3:0 */
    *tr2col = 18;  /* 3:0 */
    *scol   = 21;  /* 3 */
    *acol = 8;	/* est(a) under H3 */
    *dcol = 10;	/* est(d) under H3 */
    switch (theparams->ihypo) {
      case 1: 
      case 30:   theparams->ihypo = 30;
      case 34: 
      default:
        *lrcol = 4;	/* LR H3 v H0 */
        break;
      case 2:
      case 31:  theparams->ihypo = 31;
        *lrcol = 5;	/* LR H3 v H1 */
        break;
      case 3:
      case 32:  theparams->ihypo = 32;
        *lrcol = 6;	/* LR H3 v H2 */        
        break;
      case 10:
        *lrcol = 11;	/* LR H1 v H0 */        
        *r2col  = 13;
        *tr2col = 16;
        *scol   = 19;
        *acol = 7;	   /* est(a) under H1 */
        *dcol = 9;	   /* est(d) under H1  should be zero */
        break;
      case 20:
        *lrcol = 12;	/* LR H2 v H0 */        
        *r2col  = 14;
        *tr2col = 17;
        *scol   = 20;
        *acol = 7;	  /* est(a) under H1 */
        *dcol = 9;	  /* est(d) under H1  should be zero */
        break;
        
    }
  }
}



/*
  Count the number of QTL on each node.  Return an integer.
*/
int znode_QTLs(zresult *znode,params *theparams)
{
  FILE *inputf;
  FPN lrvalue,newlrvalue;
  int error,ch,numqtl,start, lrcol,acol,dcol,r2col,tr2col,scol,updown,chrom;
  numqtl = 0;
  updown = 1;
  lrvalue = newlrvalue = (FPN) 0.0;
  assign_zfilecols(theparams,&lrcol,&acol,&dcol,&r2col,&tr2col,&scol) ;
  inputf = fileopen(theparams->zfile, "r");
  if (inputf == NULL)
    return(0);
  error = fseek(inputf,znode->start_offset,SEEK_SET);
  if (   error != 0 ) {
    printf("\nProblem in znode_QTLs");
    start = 1;
    updown = 0;
  }
  else
    start = 0;
  while (start == 0) {
    ch =  get_next_line(gbuffer, MAXLINE, inputf);
    if (ch != EOF && gbuffer[0] != '-' ) {
      get_field(1, gname, gbuffer);
      chrom = is_pinteger(gname);
      if ( chrom != znode->chrom ) {
        start = 1;
        if (updown == 1 && lrvalue > theparams->siglevel )
          numqtl = numqtl+1;
      }
      else {
        get_field(lrcol,gname,gbuffer);
        newlrvalue = (FPN) atof(gname);
        if ( updown == 1 && newlrvalue < lrvalue ) {
          updown = -1;
          if ( lrvalue > theparams->siglevel )
            numqtl = numqtl+1;
        }
        else if ( updown == -1 && newlrvalue > lrvalue )
          updown = 1;
        lrvalue = newlrvalue;      
      }
    }
    else {
      start = 1;
      if (updown == 1 && lrvalue > theparams->siglevel )
          numqtl = numqtl+1;
    }
  }
  fileclose(theparams->zfile, inputf);
  return(numqtl);
  
}


/*
  if the node doesn't match the model, eliminate it.
*/
zresult *zresult_elim_nonmodel(zresult *first,params *theparams)
{
  zresult *ptr,*tptr,*pptr;
  int k;
  k=0;
  ptr = first;
  while ( ptr != NULL )      
    if ( ptr->model != theparams->Model ) {
      tptr = ptr->next;
      if ( tptr == NULL )
        pptr = ptr->prev;
      ptr = zresult_node(0,1,1,1,0L,0L,ptr,ptr->prev,ptr->next);
      ptr = tptr;
    }
    else {
      if ( ptr->next == NULL )
        pptr = ptr;
      ptr = ptr->next;
    }
  for ( tptr = pptr ; tptr->prev != NULL ; tptr = tptr->prev )  k+=1;

  return(tptr);
}




/*
	rearrange the nodes...flag indicates how (chrom, model, trait, etc)
*/
zresult *zresult_rank(zresult *first,int flag)
{
  zresult *tfirst,*tptr,*tptr2,*newzresult, *biggest;
  int numnodes,k;
  long cc,mm,tt,mil,thou;
  k=0;
  cc = mm = tt = 1L;
  mil = 1000000L;
  thou = 1000L;
  if ( flag > 300 ) {
    tt = mil;
    if ( flag > 312 ) 
      mm = thou;
    else 
      cc = thou;
  }
  else if ( flag > 200 ) {
    mm = mil;
    if ( flag > 213 ) 
      tt = thou;
    else 
      cc = thou;  
  }
  else if ( flag > 100 ) {
    cc = mil;
    if ( flag > 123 ) 
      tt = thou;
    else 
      mm = thou;    
  }
  numnodes = 0;
  for ( tptr = first ; tptr != NULL ; tptr = tptr->next ) {
    numnodes = numnodes + 1;
    tptr->rank = cc * (long) tptr->chrom + mm * (long) tptr->model + tt * (long) tptr->trait;
  }
  tfirst = first;
  
  while ( tfirst != NULL ) {
    biggest = tfirst;
    for ( tptr = tfirst ; tptr != NULL ; tptr = tptr->next ) {
      if ( tptr->rank > biggest->rank )
        biggest = tptr;
    }    
    if ( biggest != tfirst )
      tfirst = znode_switch(tfirst,biggest);
    tfirst = tfirst->next;
  }
  for ( tptr = first ; tptr->prev != NULL ; tptr=tptr->prev ) k+=1;
  newzresult = tptr;
  for ( tptr=newzresult ; tptr != NULL ; tptr=tptr->prev ) { /*reverse direction of list*/
    tptr2 = tptr->next;
    tptr->next = tptr->prev;
    tptr->prev = tptr2;
    tfirst = tptr;
  }
  return(tfirst);
}
  
  
  
/*switch a pair of result nodes*/
zresult *znode_switch(zresult *tfirst,zresult *biggest)
{
  zresult *tptr,zptr;
  tptr = &zptr;
  if ( biggest != tfirst->next && biggest != tfirst->prev ) { 
	  if ( tfirst->next != NULL )
	    tfirst->next->prev = biggest;
	  if ( tfirst->prev != NULL )
	    tfirst->prev->next = biggest;
	  
	  if ( biggest->next != NULL )
	    biggest->next->prev = tfirst;
	  if ( biggest->prev != NULL )
	    biggest->prev->next = tfirst;
	  
	  tptr->prev = tfirst->prev;
	  tptr->next = tfirst->next;
	  tfirst->prev = biggest->prev;
	  tfirst->next = biggest->next;
	
	  biggest->prev = tptr->prev;
	  biggest->next = tptr->next;  
  }
  else if ( biggest == tfirst->next ) {

	  if ( tfirst->prev != NULL )
	    tfirst->prev->next = biggest;
	  
	  if ( biggest->next != NULL )
	    biggest->next->prev = tfirst;

	  
	  tptr->prev = tfirst->prev;
	  tptr->next = tfirst;
	  tfirst->prev = biggest;
	  tfirst->next = biggest->next;
	
	  biggest->prev = tptr->prev;
	  biggest->next = tptr->next;    
  }
  else if ( biggest == tfirst->prev ) {

	  if ( biggest->prev != NULL )
	    biggest->prev->next = tfirst;
	  
	  if ( tfirst->next != NULL )
	    tfirst->next->prev = biggest;

	  
	  tptr->prev = biggest->prev;
	  tptr->next = biggest;
	  biggest->prev = tfirst;
	  biggest->next = tfirst->next;
	
	  tfirst->prev = tptr->prev;
	  tfirst->next = tptr->next;    
  }
  return(biggest);
}

/*
        -1    reinitialize
  flag = 0    deallocate
         1    allocate
*/
zresult *zresult_node(int flag,int c,int t,int m,long start,long end,zresult *celle,zresult *prev,zresult *next)
{
  zresult *new_node;
  if ( flag == 0 ) {  /* deallocate the node, connect prev to next and vice versa */
      if ( next != NULL )
        next->prev = prev;
      if ( prev != NULL )
        prev->next = next;    
    if ( debugging > 2 ) {
        sprintf(gwarn,"In zresult_node(), deallocated 1 znode at %p\n", celle);
        MemoryAccount(gwarn);
    }
    free((char *) celle);
    new_node = NULL;
  }
  else {
    if ( flag == 1 ) {  /* allocate new node */
#if defined(MACWARRIOR) || defined(WINWARRIOR)
      new_node = (zresult *) malloc( (size_t) sizeof(zresult));
#else
      new_node = (zresult *) malloc( (unsigned) sizeof(zresult));
#endif
      if ( debugging > 2 ) {
        sprintf(gwarn,"In zresult_node(), allocated 1 znodes at %p\n",new_node);
        MemoryAccount(gwarn);
      }
      if ( next != NULL )
        next->prev = new_node;
      if ( prev != NULL )
        prev->next = new_node;    
    }
    else
      new_node = celle;
    new_node->chrom = c;
    new_node->trait = t;
    new_node->model = m;
    new_node->start_offset = start;
    new_node->end_offset = end;
    new_node->prev = prev;
    new_node->next = next;
    new_node->window = (FPN) 0.0;
    new_node->nbp = 0;
    new_node->rank = 0;
  }
  return(new_node);
}

/*Print the node.  Used in debugging*/
void print_zresult(zresult *new_node)
{
    printf("\n\nChromosome %d", new_node->chrom) ;
    printf("\nTrait  %d",new_node->trait) ;
    printf("\nModel  %d",new_node->model) ;
    printf("\nStart  %ld",new_node->start_offset) ;
    printf("\nEnd    %ld",new_node->end_offset) ;
    printf("\nWindow %f",new_node->window) ;
    printf("\nNBP    %d",new_node->nbp) ;
    printf("\nRank   %d",new_node->rank) ;
}

/*  This gets rid of the list */
zresult *zresult_list_abolish(zresult *first)
{
  zresult *ptr,*tptr;
  ptr = first ; 
  while ( ptr != NULL ) {
    tptr = ptr->next ;
    ptr = zresult_node(0,1,1,1,0L,0L,ptr,ptr->prev,ptr->next);
    ptr = tptr;
  }
  return(NULL);

}

/*

   This will open up the results file of Zmapqtl.  It will go through it
   and produce a node for each set of results with a unique Model, trait
   and chromosome.  It will then sort the nodes based on ranker.  
   
            123  chromosome, Model, trait
            132  chromosome, trait, Model
            213  Model, chromosome, trait
   ranker = 231  Model, trait, chromosome
            312  trait, chromosome, Model
            321  trait, Model, chromosome
*/

zresult *zmapqtl_list(params *theparams,int ranker)
{
  FILE *inputf;
  zresult *first, *celle,*prev,*next;
  int k,c,t,m,nbp,flag,ch,chrom;
  FPN window;
  long start_o;
  first = celle = prev = next = NULL;
  inputf = fileopen(theparams->zfile, "r");
  if (inputf == NULL)
    return(NULL);
  chrom = flag = 0;
  start_o = ftell(inputf);
  k = 0;
  while ( (ch=get_next_token(gname,MAXNAME,inputf)) != EOF ) {
    switch (flag) {
      case 0:
        if ( !strcmp(gname,"-window") )
          flag = -1;
        else if ( !strcmp(gname,"-background") )
          flag = -2;
        else if ( !strcmp(gname,"-Model") )
          flag = -3;
        else if ( !strcmp(gname,"-trait") )
          flag = -4;
        else if ( !strcmp(gname,"-s") )
          flag = 1;
        break;
      case -1:
        window = (FPN) atof(gname);
        flag = 0;
        break;
      case -2:
        nbp = atoi(gname);
        flag = 0;
        break;
      case -3:
        m = atoi(gname);
        flag = 0;
        break;
      case -4:
        t = atoi(gname);
        flag = 0;
        break;
      default :
        if ( !strcmp(gname,"-e") ) 
          chrom = flag = 0;
        else {
          if ( flag == 1 ) {
            c = is_pinteger(gname);
             
            if ( c != -1 && c != chrom ) {
              chrom = c;
              if ( first == NULL ) 
                prev = first = zresult_node(1,c,t,m,start_o,0L,NULL,NULL,NULL);
              else
                prev = celle = zresult_node(1,c,t,m,start_o,0L,NULL,prev,NULL); 
              prev->window = window;
              prev->nbp = nbp;             
            }    
            if ( c == -1 )
              chrom = flag = 0;        
          }
          if ( flag == 12 ) {
            while ( (ch=fgetc(inputf))   != '\n') k+=1 ;
            flag = 1;
          }
          else
            flag = flag+1;          
        }      
        break;
    }
    start_o = ftell(inputf);
  }  
  fileclose(theparams->zfile, inputf);
  
  first = zresult_rank(first,ranker);
  return(first);
}

/*
  Calculate the mean and sample standard deviation from the 
  bootstrap or jackknife.
*/
int process_bootfile(params *theparams,char *infile,char *outfile,char *progname,char *chptr)    
{
  FILE *outf,*inputf;
  int model,trait,start,ch,n,chrom,mark,ft;
  FPN lr,lr2,add,add2,dom,dom2,pos,idn,idnd;
  
  inputf = fileopen(infile,"r");
	  model = trait = -1;
	  start = 0;
	  do {
	    ch = get_next_token(gname, MAXNAME, inputf);
	    if ( ch == EOF ) {
	      fileclose(infile, inputf);  
	      return(-2);
	    }
	    else if ( !strcmp(gname,"-filetype") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      ft  = file_to_int(gname);
	      if ( !(ft==63||ft==64) ) {
	        fileclose(infile, inputf);  
	        return(-1);
	      }	      
	    }
	    else if ( !strcmp(gname,"-Model") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      model = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-trait") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      trait = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-boots") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      n = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-start") ) {
	      if ( model == theparams->Model  && trait == theparams->whichtrait )
	        start = 1;
	    }
	  } while ( start != 1 && ch != EOF );  /*We should be at the right place... */
  idn = (FPN) 1.0 / (FPN) n ;
  idnd = (FPN) 1.0 / (FPN) (n-1);
  print_head(progname,outfile,chptr,0,ft+10,theparams); 
  write_zheader2(outfile, theparams, 4,n);


  outf = fileopen(outfile, "a");  
 

  do {    
	ch = get_next_token(gname, MAXNAME, inputf);/*c*/
	if ( !strcmp(gname,"-e") )  
	   start = 0;
	else {   
	  chrom = atoi(gname);
	  ch = get_next_token(gname, MAXNAME, inputf);/*m*/
	  mark = atoi(gname);
	  ch = get_next_token(gname, MAXNAME, inputf);/*p*/
	  pos = (FPN) atof(gname);
	  ch = get_next_token(gname, MAXNAME, inputf);/*lr*/
	  lr = (FPN) atof(gname)  ;
	  ch = get_next_token(gname, MAXNAME, inputf);/*lr2*/
	  lr2 = (FPN) atof(gname) ;
	  ch = get_next_token(gname, MAXNAME, inputf);/*a*/
	  add = (FPN) atof(gname) ;
	  ch = get_next_token(gname, MAXNAME, inputf);/*a2*/
	  add2 = (FPN) atof(gname) ;
	  ch = get_next_token(gname, MAXNAME, inputf);/*d*/
	  dom = (FPN) atof(gname)  ;
	  ch = get_next_token(gname, MAXNAME, inputf);/*d2*/
	  dom2 = (FPN) atof(gname);
	  lr2 = (lr2 - idn*lr*lr)*idnd;
	  lr2 = (FPN) sqrt(lr2);
	  add2 = (add2 - idn*add*add)*idnd;
	  add2 = (FPN) sqrt(add2);
	  dom2 = (dom2 - idn*dom*dom)*idnd;
	  dom2 = (FPN) sqrt(dom2);
	  lr = lr*idn;
	  add = add*idn;
	  dom = dom*idn;	 
	  if ( theparams->lodflag == 1 ) {
	    lr = lrtolod(lr);
	    lr2 = lrtolod(lr2);
	  } 
	  fprintf(outf, "\n%2d %2d %7.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f ",  chrom, mark, pos,lr,lr2,add,add2,dom,dom2);
    }
  } while ( start != 0 && ch != EOF );  
  fileclose(infile, inputf);
  fileclose(outfile, outf);
  return(n);
}


/*
  Use the results of the permutation test to calculate the 
  Experimentwise significance level.
*/
int process_permEfile(params *theparams,char *efile,char *cfile )    
{
  FILE *outf,*inputf;
  int ii,jj,model,trait,start,ch,n,nchroms;
  FPN **maxlr,l1,l05,l025,l01;
  
  inputf = fileopen(cfile,"r");
	  model = trait = -1;
	  start = 0;
	  do {
	    ch = get_next_token(gname, MAXNAME, inputf);
	    if ( ch == EOF ) {
	      fileclose(cfile, inputf);  
	      return(-2);
	    }
	    else if ( !strcmp(gname,"-Model") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      model = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-trait") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      trait = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-rwd") ) {
	      theparams->rwd = 1;
	    }
	    else if ( !strcmp(gname,"-perm") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      n = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-start") ) {
	      if ( model == theparams->Model  && trait == theparams->whichtrait )
	        start = 1;
	    }
	  } while ( start != 1 && ch != EOF );  /*We should be at the right place... */
  fileclose(cfile, inputf);


  inputf = fileopen(efile,"r");
	  model = trait = -1;
	  start = 0;
	  do {
	    ch = get_next_token(gname, MAXNAME, inputf);
	    if ( ch == EOF ) {
	      fileclose(efile, inputf);  
	      return(-2);
	    }
	    else if ( !strcmp(gname,"-Model") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      model = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-trait") ) {
	      ch = get_next_token(gname, MAXNAME, inputf);
	      trait = atoi(gname);
	    }
	    else if ( !strcmp(gname,"-start") ) {
	      if ( model == theparams->Model  && trait == theparams->whichtrait )
	        start = 1;
	    }
	  } while ( start != 1 && ch != EOF );  /*We should be at the right place... */
  if ( theparams->rwd == 1 )
    nchroms = theparams->chrom;
  else
    nchroms = 0;
  maxlr = dmatrix(0,theparams->chrom,1,n);
  for ( ii = 1 ; ii <= n ; ii++ )  { 
	ch = get_next_token(gname, MAXNAME, inputf);/* row */
	for ( jj = 0; jj<= nchroms; jj++ ) {
	  ch = get_next_token(gname, MAXNAME, inputf);/*maxlr*/
      maxlr[jj][ii] = (FPN) atof(gname);
    }
  }
  fileclose(efile, inputf);
  outf = fileopen(theparams->error, "a");

	for ( ii = 0; ii <= nchroms; ii++ ) {
	  sort(n, maxlr[ii]);
	  if (outf != NULL) {
	    l1 = maxlr[ii][ (int) (0.9 * (FPN) n) + 1];
	    l05 = maxlr[ii][ (int) (0.95 * (FPN) n) + 1];
	    l025 = maxlr[ii][ (int) (0.975 * (FPN) n) + 1];
	    l01 = maxlr[ii][ (int) (0.99 * (FPN) n) + 1];
	    if ( theparams->lodflag == 1 ) {
	      l1 = lrtolod(l1);
	      l05 = lrtolod(l05);
	      l025 = lrtolod(l025);
	      l01 = lrtolod(l01);    
	    }
	    fprintf(outf, "\n\n-start\n Performed %d permutations of the phenotypes and genotypes",n);
	    fprintf(outf, "\nHere are the Experimentwise significance levels for different sizes alpha for chromosome %d",ii);
	    fprintf(outf, "\nPermutation significance level for alpha = 0.1   : %7.4f", l1);
	    fprintf(outf, "\nPermutation significance level for alpha = 0.05  : %7.4f", l05);
	    fprintf(outf, "\nPermutation significance level for alpha = 0.025 : %7.4f", l025);
	    fprintf(outf, "\nPermutation significance level for alpha = 0.01  : %7.4f", l01);
	    fprintf(outf, "\n-end of shuffling results\n\n");
	  }
	  if ( ii == 0 ) {
	    theparams->siglevel = maxlr[ii][(int) ((1.0-theparams->size) * (FPN) n) + 1];
/*	    if ( theparams->lodflag == 1 )
	      theparams->siglevel = lrtolod(theparams->siglevel);*/
	  }
	}  
  if (outf != NULL)
    fileclose(theparams->error, outf);
  free_dmatrix(maxlr,0,theparams->chrom,1,n);
  return(n);
}

/*  Now for some code to handle the JZmapqtl output.   First, we need to 
    determine if JZmapqtl was run.   There seems little point to running it
    if there is only one trait, so this section will be entered into only if
    the number of traits is greater than 1.  
*/  
int process_jzfiles(params *theparams,char *mainfile,   genome *first)
{
  FPN window;
  int ii,nbp,model,traits,ihypo,ch,*whtraits,go_on,i;
  int nqtls,pfact;
  jzresult *thejzresults,*lptr;
  long offset;
  FILE *fptr;
  thejzresults = NULL;
  pfact = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    pfact = theparams->tcross;
  if (pfact < 3 )
    pfact = 1;
  else /*Need to know how many parameters. */
    pfact = 2;

  whtraits = NULL;
  fptr = fileopen(mainfile, "r");
  go_on = 1;
  while ( go_on == 1) {
    ch=get_next_token(gname,MAXNAME,fptr);
    if ( !strcmp(gname,"-window") ) {
      ch=get_next_token(gname,MAXNAME,fptr);
      window = (FPN) atof(gname);
    }
    else if (!strcmp(gname,"-background") ) {
      ch=get_next_token(gname,MAXNAME,fptr);
      nbp =   atoi(gname);
    }
    else if (!strcmp(gname,"-Model") ) {
      ch=get_next_token(gname,MAXNAME,fptr);
      model =   atoi(gname);
    }
    else if (!strcmp(gname,"-traits") ) {
      ch=get_next_token(gname,MAXNAME,fptr);
      traits =   atoi(gname);
      whtraits = ivector(1,traits);
    }
    else if (!strcmp(gname,"-ihypo") ) {
      ch=get_next_token(gname,MAXNAME,fptr);
      ihypo =   atoi(gname);
    }
    else if (!strcmp(gname,"-Names") ) {
      for ( i = 1; i<= traits; i++ ) {
        ch=get_next_token(gname,MAXNAME,fptr);
        whtraits[i] =   atoi(gname);
        ch=get_next_token(gname,MAXNAME,fptr);
      }
    }
    else if ( !strcmp(gname,"-s") && model == theparams->Model && ihypo == theparams->ihypo ) 
      go_on = 0;
    if ( ch == EOF )
      return(-1);
  }
  offset = ftell(fptr);
  nqtls = 0;
  thejzresults = jz_qtlnum(theparams,fptr,&nqtls,traits);
  fileclose(mainfile, fptr);
  if ( theparams->verbosity == 1 )
    printf("\nModel: %d, Window: %f, NBP: %d\nNumber of QTL in joint analysis %d \n\n",model,window,nbp, nqtls);
  if ( whtraits != NULL ) {
    for ( ii = 1; ii<= traits ; ii++ ) 
      get_jzestimates(theparams,mainfile,thejzresults,whtraits,ii);
  }
  if ( theparams->lodflag == 1 ) 
    for (lptr=thejzresults; lptr != NULL ; lptr = lptr->next ) {
	        lptr->lr10 = lrtolod(lptr->lr10);
	        lptr->lr20 = lrtolod(lptr->lr20);
	        lptr->lr30 = lrtolod(lptr->lr30);
	        lptr->lr31 = lrtolod(lptr->lr31);
	        lptr->lr32 = lrtolod(lptr->lr32);
	        lptr->lrGxE = lrtolod(lptr->lrGxE);  
    }
  
  if ( theparams->verbosity == 1 ) 
    print_jzestimates(theparams,traits,whtraits,thejzresults,stdout);

  fptr = fileopen(theparams->eqtl, "a");
  print_jzestimates(theparams,traits,whtraits,thejzresults,fptr);
  fileclose(theparams->eqtl, fptr);

  write_jzmapqtl_ranks(thejzresults, theparams, first, theparams->srfile);
  
  if ( whtraits != NULL )
    free_ivector(whtraits,1,traits);
  return(0);
}


void    print_jzestimates(params *theparams,int traits,int *whtraits,jzresult *thejzresults,FILE *fptr)
{
  int jj,ii,pfact,linelength,k;
  jzresult *backward,*forward;
  k=0;
  if ( thejzresults == NULL )
    return;
  for (forward = thejzresults; forward->next != NULL ; forward = forward->next ) k+=1;/*start with last node.*/


  fprintf(fptr,"\n\n          JZmapqtl analysis summary for Model %d, Hypothesis %d.",theparams->Model,theparams->ihypo);


  if ( theparams->ihypo < 15 )
    pfact = 1;
  else 
    pfact = 2;
  linelength = 24+9*traits*pfact;
  if ( theparams->ihypo == 14 )
    linelength = linelength + 20;
  if ( theparams->ihypo == 34 )
    linelength = linelength + 30;
  putline(fptr,'=',linelength);
  if (theparams->lodflag == 1 )
    fprintf(fptr,"\nChrom   Pos.     JointLOD ");
  else
    fprintf(fptr,"\nChrom   Pos.     JointLR  ");
  for ( jj=1; jj<= traits ; jj++ ) {
      if (pfact == 1 )
        fprintf(fptr,"  a1(%d)  ",whtraits[jj]);
      else if (pfact == 2)
        fprintf(fptr,"  a3(%d)  ",whtraits[jj]);
      if (pfact == 2   )
        fprintf(fptr,"  d3(%d)  ",whtraits[jj]);
  }  
  
  if ( theparams->ihypo == 14) {
    if ( theparams->lodflag == 1)
      fprintf(fptr,"  LOD(GxE)      a*  "); 
    else
      fprintf(fptr,"   LR(GxE)      a*  "); 
  } 
  
  if ( theparams->ihypo == 34 ) {
    if ( theparams->lodflag == 1)
      fprintf(fptr,"  LOD(GxE)      a*      d*  ");  
    else
      fprintf(fptr,"   LR(GxE)      a*      d*  ");  
  }
   
  putline(fptr,'-',linelength);

  for (  backward = forward; backward != NULL ; backward = backward->prev ) {
    fprintf(fptr,"\n    %2d %7.4f",backward->chrom,backward->position);
    if ( theparams->ihypo < 15 ) {
      fprintf(fptr," %10.4f",backward->lr10);
      for ( ii = 1; ii <= traits; ii++ )
        fprintf(fptr," %7.4f ",backward->gparams[1][whtraits[ii]]);
      if ( theparams->ihypo == 14)
        fprintf(fptr," %10.4f %7.4f",backward->lrGxE,backward->gparams[1][0]);
    }
    else {
      fprintf(fptr," %10.4f",backward->lr30);
      for ( ii = 1; ii <= traits; ii++ )
        fprintf(fptr," %7.4f  %7.4f ",backward->gparams[2][whtraits[ii]],backward->gparams[3][whtraits[ii]]);
      if ( theparams->ihypo == 34)
        fprintf(fptr," %10.4f %7.4f %7.4f",backward->lrGxE,backward->gparams[2][0],backward->gparams[3][0]);
    
    
    }
  
  
  }
  putline(fptr,'=',linelength);
  fprintf(fptr,"\n");
}
/*

Format of the inidividual files:  
-ihypo          14      Hypothesis test
# Chromosome   Marker   MarkerName  Position    col5 col6 .....
-s
     1          1            ve001  0.0001000   ..............  

Chromosome   Marker   MarkerName  Position  THEN
------------------------------------------------------
------------------------------------------------------
Hypo        5        6       7       8     9      10
------------------------------------------------------
10       LR(1:0)     a
14       LR(1:0)  LR(GxE)    a                
30       LR(3:0)    a3       d3
31       LR(3:0)  LR(3:1)  LR(1:0)   a3    d3      a1
32       LR(3:0)  LR(3:1)  LR(2:0)   a3    d3      d2
34       LR(3:0)  LR(GxE)    a3      d3
------------------------------------------------------
------------------------------------------------------

*/ 
void  get_jzestimates(params *theparams,char *mainfile,jzresult *thejzresults,int *whtraits,int wt )
{
  int ihypo,model,j,ch,go_on,thischrom,ii,k;
  FPN pos;
  long lpos,ljpos;
  jzresult *lptr;
  FILE *fptr;
  k=0;
  ii = whtraits[wt];
  for ( j=1;j<=MAXNAME;j++ )
    mainfile[j] = '\0';
  sprintf(mainfile,"%s.z%d",theparams->stem,whtraits[wt]);
  fptr = fileopen(mainfile, "r");
  go_on = 1;
  while ( go_on == 1) {
    ch=get_next_token(gname,MAXNAME,fptr);
    if (!strcmp(gname,"-Model") ) {
      ch=get_next_token(gname,MAXNAME,fptr);
      model = atoi(gname);
    } 
    else if (!strcmp(gname,"-ihypo") ) {
      ch=get_next_token(gname,MAXNAME,fptr);
      ihypo = atoi(gname);
    }
    else if ( !strcmp(gname,"-s") && model == theparams->Model && ihypo == theparams->ihypo) 
      go_on = 0;
    if ( ch == EOF )
      return;    
  }
  for ( lptr = thejzresults; lptr->next != NULL ; lptr = lptr->next ) k+=1;/*go to last node*/
  go_on = 1;
  while ( go_on == 1) {
    ch=get_next_line(gbuffer,(int) MAXLINE,fptr);

    if ( ch == EOF )
      return;    
    if ( !strcmp(gname,"-e") || lptr == NULL ) 
      go_on = 0;
    else {
      get_field(1, gname, gbuffer);
      thischrom = (int) atoi(gname);
      get_field(4, gname, gbuffer);
      pos = (FPN) atof(gname);
      lpos = (long) (1000000.0 *  pos);
      ljpos = (long) (1000000.0 * lptr->position) ;
      if ( thischrom ==lptr->chrom && lpos == ljpos ) {
        get_field(5, gname, gbuffer);
        if ( theparams->ihypo > 29 ) {
          if ( theparams->ihypo == 30 ) {
            get_field(6, gname, gbuffer);
            lptr->gparams[2][ii] = (FPN) atof(gname);         
            get_field(7, gname, gbuffer);
            lptr->gparams[3][ii] = (FPN) atof(gname);         
          
          }
          else if  ( theparams->ihypo == 31 ) {
            get_field(8, gname, gbuffer);
            lptr->gparams[2][ii] = (FPN) atof(gname);         
            get_field(9, gname, gbuffer);
            lptr->gparams[3][ii] = (FPN) atof(gname);         
            get_field(10, gname, gbuffer);
            lptr->gparams[1][ii] = (FPN) atof(gname);         
          }
          else if  ( theparams->ihypo == 32 ) {
            get_field(8, gname, gbuffer);
            lptr->gparams[2][ii] = (FPN) atof(gname);         
            get_field(9, gname, gbuffer);
            lptr->gparams[3][ii] = (FPN) atof(gname);         
            get_field(10, gname, gbuffer);
            lptr->gparams[4][ii] = (FPN) atof(gname);         
          }
          else if  ( theparams->ihypo == 34 ) {
            get_field(7, gname, gbuffer);
            lptr->gparams[2][ii] = (FPN) atof(gname);         
            get_field(8, gname, gbuffer);
            lptr->gparams[3][ii] = (FPN) atof(gname);         
          }
        }
        else {
          get_field(6, gname, gbuffer);
          if ( theparams->ihypo == 10 )
            lptr->gparams[1][ii] = (FPN) atof(gname);
          else {
            get_field(7, gname, gbuffer);
            lptr->gparams[1][ii] = (FPN) atof(gname);         
          }
        }
         
      
        lptr = lptr->prev;
      }
    }       
  }
  
  fileclose(mainfile, fptr);

}

/*

         1 if rise
pup,up = 0 if telomere
        -1 if fall


Format of the JZmapqtl joint mapping file (*.z0)

-ihypo          14      Hypothesis test
# Chromosome   Marker   MarkerName  Position    col5 col6 .....
-s
     1          1            ve001  0.0001000   ..............  

c m MN pos then
------------------------------------------------------
------------------------------------------------------
Hypo        5        6       7       8      
------------------------------------------------------
10       LR(1:0)
14       LR(1:0)  LR(GxE)    a*
30       LR(3:0)                          JOINT (*.z0)
31       LR(3:0)  LR(3:1)  LR(1:0)
32       LR(3:0)  LR(3:2)  LR(2:0)
34       LR(3:0)  LR(GxE)    a*      d*
------------------------------------------------------
------------------------------------------------------
a*, d* are average parameter values for GxE mapping
*/ 
jzresult *jz_qtlnum(params *theparams,FILE *fptr,int *nqtls,int traits)
{
  int cb,cn,ca,ch,go_on,doit;
  FPN lrb,lrn,lra,lr30,lr31,lr32,lr10,lr20,lrGxE,a1,a3,d3,d2,pos;
  jzresult *thejzresults;
  thejzresults = NULL;
  go_on = 1;
  cb = -3;
  cn = -2;
  ca = -1;
  lr30=lr31=lr32=lr10=lr20=lrGxE=a1=a3=d3=d2=lrb=lrn=lra=(FPN) 0.0;
  
  do {
    ch=get_next_line(gbuffer,(int) MAXLINE,fptr);
    cb = cn;
    cn = ca;
    lrb = lrn;
    lrn = lra;
    if ( !strcmp(gbuffer,"-e") || ch==EOF) {
      go_on = 0;
      ca = -1;
      lra = (FPN) 0.0;
    }
    else {
        get_field(1, gname, gbuffer);
        ca = atoi(gname);
        get_field(5, gname, gbuffer);
        lra = (FPN) atof(gname);
    }      
	    if ( lrn > theparams->siglevel ) {
		  doit = 0;
		  if ( cb == cn && cn == ca ) { /* Int. on chrom. */
		    if ( lrn > lrb && lrn > lra  )
		      doit = 1;
		  }
		  else if ( cn == ca && cb != cn ) { /*  Left telomere*/
		    if ( lrn > lra )
		      doit = 1;
		  }
		  else if ( cn == cb && cn != ca ) {/*  Right telomere telomere*/
		    if ( lrn > lrb )
		      doit = 1;	  
		  }
	      if ( doit == 1 ) {
	        *nqtls +=1;   
	        thejzresults = jzresult_node(1,thejzresults, traits);
	        thejzresults->chrom = cn; 
	        thejzresults->position = pos; 
	        thejzresults->lr10 = lr10;
	        thejzresults->lr20 = lr20;
	        thejzresults->lr30 = lr30;
	        thejzresults->lr31 = lr31;
	        thejzresults->lr32 = lr32;
	        thejzresults->lrGxE = lrGxE;
	        thejzresults->gparams[1][0] = a1;
	        thejzresults->gparams[2][0] = a3;
	        thejzresults->gparams[3][0] = d3;
	        thejzresults->gparams[4][0] = d2;
	      }
	    } 
	     
    if ( !strcmp(gbuffer,"-e") || ch==EOF) {
      go_on = 0;
      ca = -1;
      lra = (FPN) 0.0;
    }
    else {
        get_field(4, gname, gbuffer);
        pos = (FPN) atof(gname);
        if ( theparams->ihypo == 10   ) 
          lr10 = lra;
        else if ( theparams->ihypo == 14 ) {
          lr10 = lra;
          get_field(6, gname, gbuffer);
          lrGxE = (FPN) atof(gname);
          get_field(7, gname, gbuffer);
          a1 = (FPN) atof(gname);
        }
        else if ( theparams->ihypo == 30 )  
          lr30 = lra;
        else if ( theparams->ihypo == 31 ) {
          lr30 = lra;
          get_field(6, gname, gbuffer);
          lr31 = (FPN) atof(gname);
          get_field(7, gname, gbuffer);
          lr10 = (FPN) atof(gname);
        }
        else if ( theparams->ihypo == 32 ) {
          lr30 = lra;
          get_field(6, gname, gbuffer);
          lr32 = (FPN) atof(gname);
          get_field(7, gname, gbuffer);
          lr20 = (FPN) atof(gname);
        }
        else if ( theparams->ihypo == 34 ) {
          lr30 = lra;
          get_field(6, gname, gbuffer);
          lrGxE = (FPN) atof(gname);
          get_field(7, gname, gbuffer);
          a3 = (FPN) atof(gname);
          get_field(8, gname, gbuffer);
          d3 = (FPN) atof(gname);
        }
    }
  }    while (go_on > 0);



  return(thejzresults);
}

/*

  flag = 0    deallocate
         1    allocate
*/
jzresult *jzresult_node(int flag, jzresult *thisnode,int traits)
{
  jzresult *new_node;

  if ( flag == 1 ) {  /* allocate new node */
#if defined(MACWARRIOR) || defined(WINWARRIOR)
      new_node = (jzresult *) malloc( (size_t) sizeof(jzresult));
#else
      new_node = (jzresult *) malloc( (unsigned) sizeof(jzresult));
#endif
      if ( debugging > 2 ) {
        sprintf(gwarn,"In jzresult_node(), allocated 1 jznodes at %p\n",new_node);
        MemoryAccount(gwarn);
      }
      new_node->chrom = 0;
      new_node->traits = traits;
      new_node->position = (FPN) 0.0;
      new_node->lr10 = (FPN) 0.0;
      new_node->lr20 = (FPN) 0.0;
      new_node->lr30 = (FPN) 0.0;
      new_node->lr31 = (FPN) 0.0;
      new_node->lr32 = (FPN) 0.0;
      new_node->lrGxE = (FPN) 0.0;
      new_node->gparams = dmatrix(1,4,0,traits);
      new_node->prev = NULL;
      new_node->next = thisnode;
      if ( thisnode != NULL )
        thisnode->prev = new_node;
  }
  else { /* deallocate all nodes */
    new_node = thisnode;
    while ( new_node != NULL ) {    
      if ( debugging > 2 ) {
        sprintf(gwarn,"In jzresult_node(), deallocated 1 jznode at %p\n", new_node);
        MemoryAccount(gwarn);
      }
      free_dmatrix(new_node->gparams,1,4,0,traits);
      if ( new_node->next != NULL ) {
        new_node = new_node->next;
        free((char *) new_node->prev);
      }
      else {
        free((char *) new_node);
        new_node = NULL;
      }
    }
  }

  return( new_node);

}

/*
  Write header for the Zmapqtl output files.  
           0    Zmapqtl.out
  eorc =   1    ZpermC.out
           2    ZpermE.out
           3    Zboot.out
           
           
   boots is only used in the eorc==3 case.
*/         
void write_zheader2(char *outfile,params *theparams, int eorc,int boots)
{
  int cross;
  FILE *outf;
  cross = theparams->cross;
  if ( theparams->tcross == 1 || theparams->tcross == 2 )
    cross = theparams->tcross;
  outf = fileopen(outfile, "a");
  if (outf == NULL)
    outf = fileopen(outfile, "w");
  if (outf != NULL) {
    fprintf(outf, "\n#The position is from the left telomere on the chromosome");
    fprintf(outf, "\n-window     %6.2f      Window size for models 5 and 6",theparams->window);
    fprintf(outf, "\n-background %6d      Background parameters in model 6",theparams->nbp);
    fprintf(outf, "\n-Model      %6d      Model number",theparams->Model);
    fprintf(outf, "\n-trait      %6d      Trait analyzed", theparams->whichtrait);
    fprintf(outf, "\n-cross      %6s      Cross", theparams->thecross);
    if ( eorc == 0 ) {
      fprintf(outf, "\n#  Test Site   *");
      fprintf(outf, " Like. Ratio Test Statistics   *     Additive       *     Dominance        * Misc. HT");
      fprintf(outf, "\n c  m position");
      if ( cross == 1 || cross == 2 || cross == 5 )
        fprintf(outf, "     H0:H1      .....      .....       H1:a       ....       ....       ....      .....      .....\n-s");
      else 
        fprintf(outf, "     H0:H3      H1:H3      H2:H3       H1:a       H3:a       H2:d       H3:d      H0:H1      H0:H2\n-s");
    }   
    else if ( eorc == 1 ) {
      fprintf(outf, "\n#Position is from the left telomere on the chromosome");
      fprintf(outf, "\n#Row Chrom Mark Position  LR(Samp)    P-Val     Count     ");
    }
    else if ( eorc == 2)  
      fprintf(outf, "\n#  Rep.        Max           -start  ");
    else if (eorc == 4) {
      fprintf(outf, "\n#Position is from the left telomere on the chromosome");
      fprintf(outf, "\n#                           LR                 Add                   Dom");
      fprintf(outf, "\n#C  M  Position     Mean      StDev       Mean      StDev       Mean      StDev -boots %d  -start",boots);
    }
    fileclose(outfile, outf);
  }
}

/*
Write the the results in SRmapqtl format
*/
void write_zmapqtl_ranks(params *theparams, aqtl *qtlptr,char *srfile)
{
  FILE *outf;
  int jj,kk,k,ii,marker,wo,printit,prevmarker,prevchrom;
  
  
    /*  Write the markers to the output file....*/
    outf = fileopen(srfile, "a");
    fprintf(outf,"\n#The markers below are not ranked, but they are above the threshold given to Eqtl.\n#");
    for (kk = 1; kk <= theparams->themap->traits; kk++) {
	  fprintf(outf,"\n-Zmapqtl model %d peaks converted to nearest markers for -trait %d",theparams->Model,kk);
      if ( theparams->themap->tnames != NULL )
        fprintf(outf, " named %s \n",theparams->themap->tnames[kk]);
      else 
        fprintf(outf, "\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf, "\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf,"\nChromosome Marker   WhichQTL     C1      C2\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf,"              -start");
      wo = 0;
      k = theparams->themap->knum[kk];
      prevmarker = prevchrom = -1;
      for (ii = 1; ii <= k; ii++) {
        printit = 1;
	    wo = wo + 1;
	    if ( qtlptr[wo].c1 > qtlptr[wo].c2 )
	      marker = qtlptr[wo].mrk + 1;
	    else
	      marker = qtlptr[wo].mrk;
	    if ( marker == prevmarker && qtlptr[wo].chrm == prevchrom)
	        printit = 0;
	    prevmarker = marker;
	    prevchrom =  qtlptr[wo].chrm;
	    if ( printit == 1 )
	      fprintf(outf, "\n   %5d  %5d  %5d   %12.4f  %12.4f  ",  qtlptr[wo].chrm, marker, ii, qtlptr[wo].c1, qtlptr[wo].c2);
      }
      fprintf(outf,"                -end\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf,"\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf, "\n");
    }

    fileclose(srfile, outf);
}


/*
Write the the results in SRmapqtl format
*/
void write_jzmapqtl_ranks(jzresult *thejzresults, params *theparams,genome *first,char *srfile)
{
  FILE *outf;
  int jj,kk,k;
  jzresult *ljz;
  genome *lgen;
  
      for ( lgen = first; lgen != NULL ; lgen = lgen->next ) {
        lgen->whichqtl = 0; 
        lgen->pxo = (FPN) 0.0;
      }
      k = 0;
      for ( ljz = thejzresults; ljz != NULL ; ljz = ljz->next ) {
        lgen = NearMark(first, ljz->chrom, ljz->position);
        if ( lgen != NULL ) {
          lgen->whichqtl += 1;
          if ( theparams->ihypo > 29 && ljz->lr30 > lgen->pxo )
            lgen->pxo = ljz->lr30;
          else if ( ljz->lr10 > lgen->pxo )
            lgen->pxo = ljz->lr10;
          
          
        }   
      }
      for ( lgen = first; lgen != NULL ; lgen = lgen->next )
        if ( lgen->whichqtl > 0 ) 
          k +=1;
  
    /*  Write the markers to the output file....*/
    outf = fileopen(srfile, "a");
    fprintf(outf,"\n#The markers below are not ranked, but they are above the threshold given to Eqtl.\n#");
    fprintf(outf,"\n#Since these are based on a joint likelihood, the same set of markers are listed for each trait.\n#");
    for (kk = 1; kk <= theparams->traits; kk++) {
	  fprintf(outf,"\n-JZmapqtl model %d peaks converted to nearest markers for -trait %d\n",theparams->Model,kk);

      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf, "\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf,"\nChromosome Marker   WhichQTL  LargestJointLR  Num Peaks\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf,"              -start");
      for ( lgen=first; lgen != NULL ; lgen = lgen->next )
        if ( lgen->whichqtl > 0 ) 
          fprintf(outf, "\n   %5d  %5d      1   %12.4f            %d  ",  lgen->chrom,lgen->markr, lgen->pxo, lgen->whichqtl);

      fprintf(outf,"                -end\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf,"\n");
      for ( jj = 1 ; jj<= 55 ; jj++ )
        fprintf(outf,"-");
      fprintf(outf, "\n");
    }

    fileclose(srfile, outf);
}


/*
  For a QTL position, return the nearest marker.   
*/
genome *NearMark(genome *first, int chrom, FPN cMposition) {
  genome *gptr;
  FPN position;
  position =   cMposition;
  
  
  for (gptr = first; gptr != NULL  ; gptr = gptr->next ) {
    if ( gptr->chrom == chrom ) {
      if ( gptr->next != NULL  ) {
        if ( gptr->next->chrom == chrom && gptr->next->pos > position ) { /*  still on same chromosome.  */
          if ( gptr->next->pos - position > position - gptr->pos ) {
            if ( gptr->markr == 0 )
              return(gptr->next);
            else
              return(gptr);
          }
          else
            return(gptr->next); 
        }
        else if ( gptr->next->chrom == chrom && gptr->next->pos < position )
          continue;
        else /*  at the end of the chromosome */
          return(gptr);      
      }
      else /* we are on the correct chromosome, but at the end of the genome  */
        return(gptr);    
    }    
  }
  return(NULL);
}

/*
  Allocate space for the results and return a pointer.   This one will be the first
  in line.  
*/
aresults *AllocAnalysisResults(aresults *first, int rows, FPN window, int nbp, int trait,int filetype, char *filename, int model,int ihypo, long offset) {
  aresults *new_node;

#if defined(MACWARRIOR) || defined(WINWARRIOR)
      new_node = (aresults *) malloc( (size_t) sizeof(aresults));
#else
      new_node = (aresults *) malloc( (unsigned) sizeof(aresults));
#endif

  new_node->filename = cvector(0,MAXLINE);
  strcpy(new_node->filename, filename);
  new_node->rows = rows;
  new_node->cols = 19;
  new_node->window = window;
  new_node->nbp = nbp;
  new_node->trait = trait;
  new_node->ihypo = ihypo;
  new_node->model = model;
  new_node->offset = offset;
  new_node->filetype =  filetype;
  if ( filetype == 65   ) {
    new_node->cols = 2;  /* this will change depending on ihypo.   10, 14, 30, 31, 32, 34 */
    if ( ihypo == 14 || ihypo == 31 || ihypo == 32)
      new_node->cols = 4;     
    else if ( ihypo == 34 ) 
      new_node->cols = 5;     
    
  }
  else if ( filetype == 68   )
    new_node->cols = 13;
  new_node->cm = imatrix(1,rows,0,2);
  new_node->results = dmatrix(1,rows,1,new_node->cols);
  new_node->headers = cmatrix(1,new_node->cols,0,(int) MAXNAME);
  new_node->prev = NULL;
  new_node->next = first;
  if ( first !=NULL )
    first->prev = new_node;
  return(new_node);
}
void UnAllocAnalysisResults(aresults *first) {
  aresults *cnode;
  cnode = first;
  while ( cnode != NULL ) {
    if ( cnode->filename != NULL )
      free_cvector(cnode->filename, 0, MAXLINE);
    if (cnode->results != NULL )
      free_dmatrix(cnode->results,1, cnode->rows, 1, cnode->cols );
    if ( cnode->cm != NULL )
      free_imatrix(cnode->cm, 1, cnode->rows, 0, 2 );
    if ( cnode->headers != NULL )
      free_cmatrix(cnode->headers,1,cnode->cols,0,(int) MAXNAME);
    if ( cnode->next != NULL ) {
      cnode = cnode->next;
      free((char *) cnode->prev);
    }
    else {  /* This is the terminus.  */
      free((char *) cnode);
      cnode = NULL;
    }
  } 
}

aresults *AnalysisFirstPass(char *filename) {

  aresults *cnode;
  int rowcount, nbp, trait,model,ihypo, filetype;
  FPN window;
  int ch;
  long offset;
  
  FILE *input;
  input = fileopen(filename, "r");
  cnode = NULL;
  nbp=trait=model=ihypo=0;
  do {
    ch = get_next_token( gbuffer, (int) MAXLINE, input );
    if ( !strcmp(gbuffer,"-window") ) {
      ch=get_next_token(gname,MAXNAME,input);
      window = (FPN) atof(gname);
    }
    else if (!strcmp(gbuffer,"-background") ) {
      ch=get_next_token(gname,MAXNAME,input);
      nbp =   atoi(gname);
    }
    else if (!strcmp(gbuffer,"-Model") ) {
      ch=get_next_token(gname,MAXNAME,input);
      model =   atoi(gname);
    }
    else if (!strncmp(gbuffer,"-trait",6) ) {
      ch=get_next_token(gname,MAXNAME,input);
      trait =   atoi(gname);
    }
    else if (!strcmp(gbuffer,"-ihypo") ) {
      ch=get_next_token(gname,MAXNAME,input);
      ihypo =   atoi(gname);
    }
    else if (!strcmp(gbuffer,"-cross") ) {
      ch=get_next_token(gname,MAXNAME,input);
      ihypo =   atoi(gname);
    }
    else if (!strcmp(gbuffer,"-filetype") ) {
      ch=get_next_token(gname,MAXNAME,input);
      filetype = file_to_int(gname);
    }
    else if (!strcmp(gbuffer,"-s") ) {
      offset = ftell(input);
      rowcount = 0;
      do {
        ch = get_next_line(gbuffer,(int) MAXLINE, input);
        if ( ch == EOF || ( gbuffer[0] == '-' && gbuffer[1] == 'e' ) ) {
          cnode = AllocAnalysisResults(cnode, rowcount, window, nbp, trait,filetype, filename, model, ihypo, offset );        
          rowcount = -1;
        }
        else
          rowcount +=1;      
      } while (ch != EOF  && rowcount >= 0 );
    }
  } while (ch != EOF);
  fileclose(filename,input);
  return(cnode);
}
/*
  Collect all the numbers into the structures that have been created.  
   
              Zmapqtl                         JZmapqtl                       Pleiotopy
Column    bc        sf           10    14     30     31      32   34       
  1        ---------------------  position in Morgans ------------------------------
  2      H0:H1     H0:H3       H0:H1  H0:H1  H0:H3  H0:H3  H0:H3  H0:H3        H11:H0
  3      R2(0:1)   H1:H3              GxE           H1:H3  H2:H3  GxE          H11:H10
  4     TR2(0:1)   H2:H3               a            H0:H1  H0:H2   a           H11:H01
  5      H1:a      H1:a                                            d           l11
  6      S1        H3:a                                                        H11:a1
  7       .        H2:d                                                        H11:d1
  8       .        H3:d                                                        H11:a2
  9       .        H0:H1                                                       H11:d2
 10       .        H0:H2                                                       H10:a1
 11       .        R2(0:1)                                                     H10:d1
 12       .        R2(0:2)                                                     H01:a1
 13       .        R2(0:3)                                                     H01:d1
 14       .       TR2(0:1)                                                      
 15       .       TR2(0:2)                                                       
 16       .       TR2(0:3)                                                       
 17       .        S1                                                       
 18       .        S2                                                       
 19       .        S3                                                       
-------------------------------------------------------------------------------------------  
        19        19             2      4        2       4     4     5           13
-------------------------------------------------------------------------------------------  
 ft     60        60             65    65       65       65    65   65           68 
  
  
*/
void AnalysisSecondPass(aresults *first) {
  aresults *cnode;
  FILE *input;
  int i,j,ch;
  
  for ( cnode=first; cnode!=NULL; cnode=cnode->next ) {
    input = fileopen(cnode->filename, "r");
    fseek(input,cnode->offset, 0L);

    for ( i=1; i<=cnode->rows; i++ ) {
      ch=get_next_token(gname,MAXNAME,input);    
      cnode->cm[i][1] = atoi(gname);
      ch=get_next_token(gname,MAXNAME,input);    
      cnode->cm[i][2] = atoi(gname);
      if ( cnode->filetype == 68 ||  cnode->filetype == 65)
        ch=get_next_token(gname,MAXNAME,input);
      for ( j=1; j<=cnode->cols; j++ ) {
        ch=get_next_token(gname,MAXNAME,input);    
        cnode->results[i][j] = (FPN) atof(gname);
      
      }
    }    
    
    
    fileclose(cnode->filename,input);
  }
}
/*
              Zmapqtl                         JZmapqtl                       Pleiotopy
Column    bc        sf           10    14     30     31      32   34       
  1        ---------------------  position in Morgans ------------------------------
  2      H0:H1     H0:H3       H0:H1  H0:H1  H0:H3  H0:H3  H0:H3  H0:H3        H11:H0
  3      R2(0:1)   H1:H3              GxE           H1:H3  H2:H3  GxE          H11:H10
  4     TR2(0:1)   H2:H3               a            H0:H1  H0:H2   a           H11:H01
  5      H1:a      H1:a                                            d           l11
  6      S1        H3:a                                                        H11:a1
  7       .        H2:d                                                        H11:d1
  8       .        H3:d                                                        H11:a2
  9       .        H0:H1                                                       H11:d2
 10       .        H0:H2                                                       H10:a1
 11       .        R2(0:1)                                                     H10:d1
 12       .        R2(0:2)                                                     H01:a1
 13       .        R2(0:3)                                                     H01:d1
 14       .       TR2(0:1)                                                      
 15       .       TR2(0:2)                                                       
 16       .       TR2(0:3)                                                       
 17       .        S1                                                       
 18       .        S2                                                       
 19       .        S3                                                       
-------------------------------------------------------------------------------------------  
        19        19             2      4        2       4     4     5           13
-------------------------------------------------------------------------------------------  
 ft     60        60             65    65       65       65    65   65           68 
  

theparams->ngt = 2, 3 for (BC, RI) v (Fx)
*/
void SetAnalysisHeaders(aresults *cnode,params *theparams) {

  int i;
  
  for (i=1; i<=cnode->cols; i++ )
    strcpy( cnode->headers[i], "N/A");
  strcpy(cnode->headers[1], "Postions");
  if ( cnode->filetype == 60) {
    if (theparams->ngt == 2 ) {
      strcpy(cnode->headers[2], "l(H1/H0)");
      strcpy(cnode->headers[3], "r2(H1/H0)");
      strcpy(cnode->headers[4], "tr2(H1/H0)");
      strcpy(cnode->headers[5], "est(a)");
      strcpy(cnode->headers[6], "S1");
    }
    else {
      strcpy(cnode->headers[2], "l(H3/H0)");
      strcpy(cnode->headers[3], "l(H3/H1)");
      strcpy(cnode->headers[4], "l(H3/H2)");
      strcpy(cnode->headers[5], "H1est(a)");
      strcpy(cnode->headers[6], "H3est(a)");
      strcpy(cnode->headers[7], "H2est(d)");
      strcpy(cnode->headers[8], "H3est(d)");
      strcpy(cnode->headers[9], "l(H1/H0)");
      strcpy(cnode->headers[10], "l(H2/H0)");
      strcpy(cnode->headers[11], "r2(H1/H0)");
      strcpy(cnode->headers[12], "r2(H2/H0)");
      strcpy(cnode->headers[13], "r2(H3/H0)");
      strcpy(cnode->headers[14], "tr2(H1/H0)");
      strcpy(cnode->headers[15], "tr2(H2/H0)");
      strcpy(cnode->headers[16], "tr2(H3/H0)");
      strcpy(cnode->headers[17], "S1");
      strcpy(cnode->headers[18], "S2");
      strcpy(cnode->headers[19], "S3");
    }
  }
  else if ( cnode->filetype == 65 ) {
    if ( theparams->ngt == 2 )
      strcpy(cnode->headers[2], "l(H1/H0)");
    else
      strcpy(cnode->headers[2], "l(H3/H0)");
    if ( cnode->ihypo == 14 || cnode->ihypo == 34 ) {
      strcpy(cnode->headers[3], "l(GxE/H1)");
      strcpy(cnode->headers[4], "GxEest(a)");
      if ( cnode->ihypo == 34 )
         strcpy(cnode->headers[4], "GxEest(d)");
    
    } 
    else if (  cnode->ihypo == 31 ) {
      strcpy(cnode->headers[3], "l(H3/H1)");
      strcpy(cnode->headers[4], "l(H1/H0)");
    
    }   
    else if (  cnode->ihypo == 32 ) {
      strcpy(cnode->headers[3], "l(H3/H2)");
      strcpy(cnode->headers[4], "l(H2/H0)");
    
    }   
  
  } 
  else if ( cnode->filetype == 68 ) {
      strcpy(cnode->headers[2], "l(H11/H0)");
      strcpy(cnode->headers[3], "l(H11/H10)");
      strcpy(cnode->headers[4], "l(H11/H01)");
      strcpy(cnode->headers[5], "l11");
      strcpy(cnode->headers[6], "H11est(a1)");
      strcpy(cnode->headers[7], "H11est(d1)");
      strcpy(cnode->headers[8], "H11est(a2)");
      strcpy(cnode->headers[9], "H11est(d2)");
      strcpy(cnode->headers[10], "H10est(a)");
      strcpy(cnode->headers[11], "H10est(d)");
      strcpy(cnode->headers[12], "H01est(a)");
      strcpy(cnode->headers[13], "H01est(d)");
  
  
  }

}
void ShowAnalysisResults(FILE *out, aresults *cnode, params *theparams) {
  int i,j;
  
  fprintf(out, "\n#    Here are the results from file %s at position %ld with ", cnode->filename,cnode->offset);
  write_file_type(cnode->filetype, out);
  fprintf(out, "\n-window %f", cnode->window);
  fprintf(out, "\n-nbp         %d", cnode->nbp);
  fprintf(out, "\n-model       %d", cnode->model);
  fprintf(out, "\n-ihypo       %d", cnode->ihypo);
  fprintf(out, "\n-trait       %d", cnode->trait);
  fprintf(out, "\n-Cross       %s", theparams->thecross);
  fprintf(out, "\n# Ch Mk");
  for (i=1;i<=cnode->cols; i++) 
    fprintf(out, " %12s", cnode->headers[i]);
  fprintf(out, "\n-s");
  for (i=1; i<=cnode->rows; i++ ) {
    fprintf(out,"\n%3d %3d", cnode->cm[i][1], cnode->cm[i][2]);
    for (j=1; j<=cnode->cols; j++ )
      fprintf(out, " %12.6f", cnode->results[i][j] );
  
  
  }
  
  fprintf(out,"\n-e\n\n");
}
/*
 attach second onto first.   
*/
void AppendResults(aresults *first, aresults *second) {

  aresults *cnode;
  if ( first != NULL && second != NULL ) {
    for ( cnode = first; cnode->next != NULL; cnode=cnode->next);
  
    cnode->next = second;
    second->prev = cnode;
  }
}

/*
   Zip down a column and 

*/

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file EQfunc.c
------------------------------------------------------------------ */


/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MRfunc.c) is part of QTL Cartographer
         
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



/*  Multiple regression subroutines for use with MultiRegress.*/

#include "Main.h"



void  AppendAtoB(char *filea,char *fileb) {
  int ch;
  FILE *outputf,*inputf;
  /*  copy over results*/
  outputf = fileopen(fileb, "a");
  inputf  = fileopen(filea, "r");
  while ( (ch=fgetc(inputf))  != EOF )
    fputc(ch,outputf);
  fileclose(fileb, outputf);
  fileclose(filea, inputf);
}


/*
   Show the parameters used in this run.   Can write to stdout or file. 
   Mainly for header in output file.
*/
void ShowMRParams(params *theparams, char *outfile, linpakws *lnpk,char **otnames ) {
  int i;
  FILE *outf;
   if ( outfile[0] == '-' )
     outf = stdout;
   else 
     outf = fileopen(outfile, "a");

   fprintf(outf,"\n#     Here are the parameters used in this run\n#");
   fprintf(outf,"\n-inputfile  %15s  #  Data file",theparams->mrinput);
   fprintf(outf,"\n-whichtrait %10d       #  Trait to analyze  ",  theparams->whichtrait );
   fprintf(outf,"\n-cattraits  %10d       #  Include [1] or exclude [0] categorical traits",  theparams->boots );
   if ( theparams->boots == 1  )  {
     for (i=1; i<=lnpk->k; i++ )   
       fprintf(outf,"\n\t-factor %d named %s has %d levels",i,otnames[i],lnpk->bp[i][0]);
     fprintf(outf,"\n-catparams %10d       #  Total number of parameters used by categorical traits",  lnpk->ldx );
   }   
   fprintf(outf,"\n-size       %10.5f       #  Type I error rate (size) of tests",  theparams->size );
   fprintf(outf,"\n-maxsites   %10d       #  Max. #  parameters allowed in model",  theparams->srupper );
   fprintf(outf,"\n-window     %10.5f       #  Blocked window around sites in the model",  theparams->window  );
   fprintf(outf,"\n-hypothesis %10d       #  Dominance inclusion code:  10 => no,  30 => yes",  theparams->ihypo );
   fprintf(outf,"\n#\n#" );
   if ( outfile[0] != '-' )
     fileclose(outfile, outf);

}


/*

#   12345789    -filetype qtls.inp             # Documentation at end
-Units       cM   
-named      yes   

-start qtls 3
Trait_1 4
    1   9.1  0.75  0.0
    1  89.1  0.5   0.0
    3  68.4  0.22  0.0
    4  43.2  0.95  0.0
Trait_2 2
    2  93.4  0.42  0.0
    4  33.2  0.90  0.0
Trait_3 1
    1  33.4  0.84  0.2
-stop qtls



*/
void ShowQTLSinp(genome *first,params *theparams,char *outfile, int i, char *tnames,int where) {
   genome *aptr,*dptr;
   int cntr,k;
   FILE *outf;
   if ( outfile[0] == '-' )
     outf = stdout;
   else 
     outf = fileopen(outfile, "a");
  k=0;
  if ( where == 1 ) {
	fprintf(outf,"\n-Units       cM          for Morgans  ");  
	fprintf(outf,"\n-named      yes ");  
	fprintf(outf,"\n-start qtls %d", theparams->traits); 
  }
  
  else if ( where == 2 ) {
    if ( i == 0 ) 
      fprintf(outf,"\n%s 0", tnames  );
    else {
      cntr = 0;
      for ( aptr=first; aptr != NULL ; aptr = aptr->next )  
        if ( aptr->whichqtl > 0 ) 
          cntr +=1;
      fprintf(outf,"\n%s %d", tnames,cntr );
      for ( aptr=first; aptr != NULL ; aptr = aptr->next )  
        if ( aptr->whichqtl > 0 )  {
          if ( theparams->ihypo == 30 )
            for ( dptr=first; dptr !=NULL && dptr->whichqtl != -aptr->whichqtl; dptr = dptr->next ) k+=1;
          else
            dptr = NULL;
          if ( dptr == NULL )         
            fprintf(outf,"\n%4d  %8.4f  %8.4f  0.0", aptr->chrom,  100.0 * aptr->pos, aptr->mxo);
          else
            fprintf(outf,"\n%4d  %8.4f  %8.4f %8.4f ", aptr->chrom,  100.0 *   aptr->pos,  aptr->mxo,  dptr->mxo);

        }
    }
  }  
  else if ( where == 3 ) {
	fprintf(outf,"\n-stop qtls\n-end\n");
  }
  
      
  if ( outfile[0] != '-' )
     fileclose(outfile, outf);

}

/*
    Print the results contained in the  genome structure.  This can go to a file or
    standard output.
*/
void ShowMultiRegress(genome *first,int *pstatus, params *theparams,char *outfile,int i,char *tnames) {
   genome *aptr,*dptr;
   FILE *outf;
   int k;
   k=0;
   if ( outfile[0] == '-' )
     outf = stdout;
   else 
     outf = fileopen(outfile, "a");
     
   if ( first != NULL ) {
      fprintf(outf,"\n  Here is the model for -trait %d -named %s...",i,tnames);
      putline(outf, '-', 60);
      putline(outf, '-', 60);
      fprintf(outf,"\nQTL  Chrom Marker  Position  P-value  Additive Dominance");
      putline(outf, '-', 60);
      for ( aptr=first; aptr != NULL ; aptr = aptr->next )  
        if ( aptr->whichqtl > 0 )  {
          if ( theparams->ihypo == 30 )
            for ( dptr=first; dptr !=NULL && dptr->whichqtl != -aptr->whichqtl; dptr = dptr->next ) k+=1;
          else
            dptr = NULL;
          if ( dptr == NULL )         
            fprintf(outf,"\n%4d  %4d %4d   %8.4f  %8.4f %8.4f       na",pstatus[aptr->whichqtl], aptr->chrom, aptr->markr, aptr->pos, aptr->pxo,aptr->mxo);
          else
            fprintf(outf,"\n%4d  %4d %4d   %8.4f  %8.4f %8.4f %8.4f ",pstatus[aptr->whichqtl], aptr->chrom, aptr->markr, aptr->pos, aptr->pxo,aptr->mxo,dptr->mxo);

        }

   
      putline(outf, '-', 60);
    
      fprintf(outf,"\n");


   }


  if ( outfile[0] != '-' )
     fileclose(outfile, outf);
}




/*
   Estimate Effects
*/
void EstimateEffects(genome *first,linpakws *lnpk, params *theparams, int trait) {
  genome *tptr;
  int rows,j,k,error;
  FPN SSe1;
  
  for ( k=1; k<=theparams->nn; k++ )
    lnpk->xx[1][k] = (FPN)1.0;   /*  This is a row for the means. */
  
  rows = j = SetUpMatrix(lnpk,theparams);
  for ( tptr= first; tptr!=NULL ; tptr = tptr->next ) {
    j +=1;
    scopy(theparams->nn,tptr->values,1,lnpk->xx[j],1);
  }
  error = sqrst(lnpk->xx,theparams->nn,theparams->nn,j,lnpk->y[trait],(FPN)0.0,lnpk->bb[trait],lnpk->rsd[trait],&k,lnpk->jpvt[trait],lnpk->qraux);
  SSe1 = sdot(theparams->nn,lnpk->rsd[trait],1,lnpk->rsd[trait],1);
  j=rows;
  for ( tptr= first; tptr!=NULL ; tptr = tptr->next ) {
    j +=1;
    tptr->mxo = lnpk->bb[trait][j];
  }

}

/*
    Do stepwise multiple regression of trait on the positions.   
*/
genome *multiregress(char *infile,  params *theparams, linpakws *lnpk,long *positions,long *dpositions, int trait) {
  int i,k, go_on, *pstatus,ii,whichi,*chromosomes,thestep,hmsteps,hardbound; 
  FPN pvalue, max,min,ysum,ysum2,syy,*locations;
  genome *ttgptr,*tgptr,*gptr,*first, *aptr,*dptr;
  k=0;
  /*  lnpk->pcnts[i]  =  { -1 => blocked,   0 => try,  >0 => already in model  */
  /*  lnpk->samplesize[i]  =  chromosome of position i */
  /*  lnpk->lratio[i]  =  location from left telomere of position i*/
  
  pstatus = lnpk->pcnts;
  chromosomes = lnpk->samplesize;
  locations = lnpk->lratio;
  go_on = 1;
  first = NULL;
  for ( i=1; i<= (int) positions[0] ; i++ )
        pstatus[i] = 0;   /* Start with all at zero */
  
  ysum = ysum2 = (FPN)0.0;
  for ( ii = 1 ; ii <= theparams->nn ; ii++ )  {
      ysum = ysum + lnpk->y[trait][ii];
      ysum2 = ysum2 + lnpk->y[trait][ii] * lnpk->y[trait][ii];
  }
  syy = ysum2 - ysum*ysum/((FPN) theparams->nn-1) ;

/*

Do the regression ...   first try to find a new site
                        if new one found, test all sites
                        if no new one found, end process
                        if max num already found, quit
*/
  hmsteps = theparams->srupper/theparams->ngt - lnpk->ldx - 1;
  hardbound = 2 * (int) ( sqrt( (FPN) theparams->nn) );
  if ( hardbound < hmsteps )
    hmsteps = hardbound; 
  while ( go_on > 0 && go_on <= hmsteps ) {  
    update_pstatus(pstatus,chromosomes,locations,positions,first,theparams);
    aptr = alloc_genome();  /* Allocate a genome node, allocate space for data and name */
    aptr->values = dvector(1,theparams->nn);
    aptr->n = theparams->nn;
    aptr->markername = cvector(0,MAXNAME);
    aptr->next = first;  /* this node points to the existing chain */
    if ( first != NULL )
       first->prev = aptr; /*  the existing chain points back to it */
    first = aptr;        /*  it is now the beginning of the chain.  */
    min = (FPN)1.0;
    i = 0;
    thestep = 0;
    if ( theparams->ihypo == 30 ) {    
	    dptr = alloc_genome();  /* Allocate a genome node, allocate space for data and name */
	    dptr->values = dvector(1,theparams->nn);
	    dptr->n = theparams->nn;
	    dptr->markername = cvector(0,MAXNAME);
	    dptr->next = first;  /* this node points to the existing chain */
	    if ( first != NULL )
	       first->prev = dptr; /*  the existing chain points back to it */
	    first = dptr;        /*  it is now the beginning of the chain.  */    
    }
    else
      dptr = aptr;
    for ( i=1; i<= (int) positions[0]; i++ ) 
      if ( pstatus[i] == 0 ) {  /*  get it and test it.  */
	  	    if ( theparams->verbosity == 1 ) 
	  	      thestep = Rotator(thestep);
        GetSiteExpecteds(infile,theparams->nn,aptr,positions,i);
        if ( theparams->ihypo == 30 )
          GetSiteExpecteds(infile,theparams->nn,dptr,dpositions,i);
        pvalue = TestThisModel( aptr, dptr, first, lnpk, theparams,trait,syy);
        if ( pvalue > 0.0 && pvalue < min ) {
          min = pvalue;
          whichi = i;
          
        }
      }
    if ( min < theparams->size ) { /* We can add another site to the chain */  
      aptr->pxo = min;
      aptr->whichqtl = whichi;
      if ( whichi < 0 ) 
        mypause();
      pstatus[whichi] = go_on;
      GetSiteExpecteds(infile,theparams->nn,aptr,positions,whichi);
      if ( theparams->ihypo == 30 ) {
        dptr->pxo = min;
        dptr->whichqtl = -aptr->whichqtl;
        GetSiteExpecteds(infile,theparams->nn,dptr,dpositions,whichi);
      }
      if ( theparams->verbosity == 1 )
        printf("   Add QTL %d (lfm %s on chrom %d, mark %d, pos %7.3f cM)\n",go_on,aptr->markername,aptr->chrom,aptr->markr, 100.0*aptr->pos);
      go_on += 1;
 /*  Check all additive nodes for significance...  */ 
      do {   
        max = (FPN)0.0;
        for ( gptr=first; gptr!=NULL ; gptr=gptr->next ) 
          if ( gptr != aptr  && gptr->whichqtl > 0 ) {            
            for (tgptr=first; tgptr != NULL && tgptr->whichqtl != -gptr->whichqtl; tgptr = tgptr->next ) k +=1;;
            if ( tgptr == NULL )  /*  if there is a dominance effect, find it and test it jointly.  */
              tgptr = gptr;
                     
            gptr->mxo = TestThisModel( gptr, tgptr, first, lnpk, theparams,trait,syy);
            if ( gptr->mxo > max ) {
              ttgptr = gptr;
              max = gptr->mxo;        
            }
          }
        if ( max > theparams->size ) {
          if ( theparams->verbosity == 1 )          
            printf("    --> Delete QTL %d (lfm %s on chrom %d, mark %d, pos %7.3f cM)\n", pstatus[ttgptr->whichqtl],ttgptr->markername,ttgptr->chrom,ttgptr->markr, 100.0*ttgptr->pos);
          for (tgptr=first; tgptr != NULL && tgptr->whichqtl != -ttgptr->whichqtl; tgptr = tgptr->next ) k -=1;
          first = del_genome_node(ttgptr);  /*  if there is a dominance effect, find it and pull it too. */
          if ( tgptr != NULL )
            first = del_genome_node(tgptr);   
              
        }
      } while ( max > theparams->size ) ;  
    }
    else {  /* get rid of the new node(s) */
      go_on = 0;
      first = del_genome_node(aptr);
      if ( theparams->ihypo == 30 )  
        first = del_genome_node(dptr);
      if ( theparams->verbosity == 1 )
        printf(" ------>  Finished searching\n");      
    }  
  }

  return(first); /*  This will return a list of sites...each with a rank and a pvalue*/
}

/*
     Get chromosome and locations from teleomers for all positions

*/
void  GetChromLocales(int *chromosomes,FPN *locations,long *positions,char *infile) {
  FILE *fptr;
  int ch,i;
  
  fptr = fileopen(infile, "r");
    
  for ( i=1; i<= (int) positions[0] ; i++ ) {
    fseek(fptr,positions[i],SEEK_SET);
    
    ch = 1;
    while ( ch != EOF ) {
	    ch = get_next_token(gname, MAXNAME, fptr);
	    if  (!strncmp("-chromosome",gname,11)) {
	      ch = get_next_token(gname, MAXNAME, fptr);
	      chromosomes[i]  =   atoi( gname );
	    }
	    else if (!strncmp("-position",gname,9)) {
	      ch = get_next_token(gname, MAXNAME, fptr);
	      locations[i] = (FPN)100.0 *  (FPN)  atof( gname );
	      ch = EOF;
	    }
    }    
  }
  fileclose(infile, fptr);
}


/*
   Status of each position.   
                    -1 => blocked  (within the window of a position in the model)
   pstatus[i] =      0 => try this in the model
                     >0 => already in model   
   This routine will ratchet down the values in pstatus...That is, check all 0 values.  If in
   the window around any additive QTL, then make it -1.
*/
void     update_pstatus(int *pstatus,int *chromosomes,FPN *locations,long *positions,genome *first,params *theparams) {
  int i;
  genome *gptr;

  if ( first != NULL )                 /*  If there is a model, then modify the status of each site. */
    for ( i=1; i<=positions[0]; i++ ) /*  Look at all sites. */
      if ( pstatus[i] <= 0 )         /*  Check a site if it is not already in the model */
        for ( gptr = first; gptr != NULL ; gptr = gptr->next ) /* look at all nodes */
          if ( gptr->whichqtl > 0 && chromosomes[i] == chromosomes[gptr->whichqtl] )  /* If this node is additive and on the same chromosome,*/
            if ( fabs(locations[gptr->whichqtl]-locations[i]) <= theparams->window ) /* and if it is within the window around the QTL,      */
              pstatus[i] = -1;                                                      /* bar it from being considered in future analyses.    */


}

/*
   Calculate a p-value for the pair of models:  
   
     H0:  all nodes in first structure minus aptr and dptr    (if it exists)
     H1:  all nodes in first structure
*/
FPN TestThisModel(genome *aptr, genome *dptr,genome *first,linpakws *lnpk, params *theparams, int trait,FPN syy) {
  genome *tptr;
  int j,k,error;
  FPN SSe1,SSr1,SSe0,SSr0,fstat0,v1,v2,pvalue;
  

  
  j=SetUpMatrix(lnpk,theparams);
  
  for ( tptr= first; tptr!=NULL ; tptr = tptr->next ) {
    j +=1;
    scopy(theparams->nn,tptr->values,1,lnpk->xx[j],1);
  }
  /*
    ShowDVector(lnpk->y[trait], 1, theparams->nn, 7);
  for ( i=1; i<= j ; i++ )
    ShowDVector(lnpk->xx[i], 1, theparams->nn, 6);
    */
  error = sqrst(lnpk->xx,theparams->nn,theparams->nn,j,lnpk->y[trait],(FPN)0.0,lnpk->bb[trait],lnpk->rsd[trait],&k,lnpk->jpvt[trait],lnpk->qraux);
  SSe1 = sdot(theparams->nn,lnpk->rsd[trait],1,lnpk->rsd[trait],1);
  SSr1 = syy - SSe1;
    
  j = SetUpMatrix(lnpk,theparams);
  for ( tptr= first; tptr!=NULL ; tptr = tptr->next ) 
    if ( tptr != aptr && tptr != dptr) {
      j +=1;
      scopy(theparams->nn,tptr->values,1,lnpk->xx[j],1);
    }
  
  error = sqrst(lnpk->xx,theparams->nn,theparams->nn,j,lnpk->y[trait],(FPN)0.0,lnpk->bb[trait],lnpk->rsd[trait],&k,lnpk->jpvt[trait],lnpk->qraux);
  SSe0 = sdot(theparams->nn,lnpk->rsd[trait],1,lnpk->rsd[trait],1);
  SSr0 = syy - SSe0;
  
  fstat0 = (SSr1-SSr0) * (FPN) (theparams->nn-j-2) / SSe1 ;
  v2 = (FPN) (theparams->nn-j-2) ;
  v1 = (FPN)1.0;
  /*  fstat0 has an F distribution with v1 = 1 and v2 = n-j-2 dof 
  if ( fstat0 < 0.0 ) 
    pvalue = 1.0;
  else*/
    pvalue = betai( (FPN)0.5* v2, (FPN)0.5*v1 ,  v2 / ( v2 + v1*fstat0 ) );
  return(pvalue);
}

/*
  Set up the matrix for analysis.   The first row is the mean, and the next rows will
  contain indicator variables for the categorical traits.   
*/
int SetUpMatrix(linpakws *lnpk, params *theparams) {
  int i,j,k,l;
  j=1;
  for ( i=1; i<=theparams->nn; i++ )
    lnpk->xx[1][i] = (FPN)1.0;   /*  This is a row for the means. */
  if ( theparams->boots == 1 ) /* Use categorical traits as well.  */
    for ( k=1; k<=lnpk->k; k++ ) 
      for ( l=1; l<lnpk->bp[k][0]; l++ ) {  /* each factor has bp[k][0]-1 levels*/
        j +=1;
        for ( i=1; i<=theparams->nn; i++ )
          if ( lnpk->bp[k][i] == l )
            lnpk->xx[j][i] = 1;
          else
            lnpk->xx[j][i] = 0;      
      }
  return(j);
}


/*  Mainly for debugging*/
void ShowSiteExpecteds(genome *gptr) {
  printf("\n Got  Expected values for marker %d on chromosome %d named %s, postion %f",gptr->markr,gptr->chrom, gptr->markername, gptr->pos);
/*  for (i=1; i<=gptr->n; i++ )
     printf("\n %f",gptr->values[i]);*/

}


/* Get the expected genotype values at the indicated site.
*/
void GetSiteExpecteds(char *infile,int n, genome *gptr,long *positions, int whichsite) {
  FILE *fptr;
  int ch,i;
 
  fptr = fileopen(infile, "r");
  fseek(fptr,positions[whichsite],SEEK_SET);
  ch = 1;
  while ( ch != EOF ) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if (!strncmp("-chromosome",gname,11)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      gptr->chrom = atoi(gname);
    }
    else if (!strncmp("-name",gname,5)) {
      ch = get_next_token(gptr->markername, MAXNAME, fptr);
    }
    else if (!strncmp("-marker",gname,7)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      gptr->markr = atoi(gname);
    }
    else if (!strncmp("-position",gname,9)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      gptr->pos = (FPN) atof(gname);
    }
    else if (!strncmp("-values",gname,7)) {        
      /*ch = get_next_token(names[whichtrait], MAXNAME, fptr);*/
      for ( i=1; i<=n; i++ ) {
        ch = get_next_token(gname, MAXNAME, fptr);
        gptr->values[i] = (FPN) atof(gname);          
      }
    }
    else if (!strncmp("-Site",gname,5))
      ch = EOF;
    else if (!strncmp("-end",gname,4))
      ch = EOF;
  }
  fileclose(infile, fptr);


}

/*  Get the parameters at the beginning of the file.   Stop when the
    trait data have been reached.   
*/
FPN  **GetTraitData(char *infile, int ntraits, int n,char **names) {
  FILE *fptr;
  int ch,i,whichtrait;
  FPN **traits;
  
  fptr = fileopen(infile, "r");
  if ( fptr == NULL )
    return(NULL);
  traits = dmatrix(0,ntraits,1,n);
  ch = 1;
  while ( ch != EOF ) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if (!strncmp("-Trait",gname,6)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      whichtrait = atoi(gname);
      ch = get_next_token(names[whichtrait], MAXNAME, fptr);
      for ( i=1; i<=n; i++ ) {
        ch = get_next_token(gname, MAXNAME, fptr);
        traits[whichtrait][i] = (FPN) atof(gname);          
      }
    }
    else if (!strncmp("-Site",gname,5))
      ch = EOF;
  }
  fileclose(infile, fptr);
  return(traits);
}

/*  
   Get the categorical trait values.   
*/
int  **GetOTraitData(char *infile, int ntraits, int n,char **names) {
  FILE *fptr;
  int ch,i,whichtrait,max;
  int **traits;
  
  fptr = fileopen(infile, "r");
  if ( fptr == NULL )
    return(NULL);
  traits = imatrix(1,ntraits,0,n);
  ch = 1;
  while ( ch != EOF ) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if (!strncmp("-Otrait",gname,7)) {
      max = 0;
      ch = get_next_token(gname, MAXNAME, fptr);
      whichtrait = atoi(gname);
      ch = get_next_token(names[whichtrait], MAXNAME, fptr);
      for ( i=1; i<=n; i++ ) {
        ch = get_next_token(gname, MAXNAME, fptr);
        traits[whichtrait][i] =   atoi(gname);  
        if (  traits[whichtrait][i] > max )
          max = traits[whichtrait][i];       
      }
      traits[whichtrait][0] = max;
    }
    else if (!strncmp("-Site",gname,5))
      ch = EOF;
  }
  fileclose(infile, fptr);
  return(traits);
}


/*  Get the parameters at the beginning of the file.   Stop when the
    trait data have been reached.   
*/
void GetMRParams(char *infile, FPN *walk, int *otraits, int *traits, int *n) {
  FILE *fptr;
  int ch;
  
  fptr = fileopen(infile, "r");
  ch = 1;
  while ( ch != EOF ) {
    ch = get_next_token(gname, MAXNAME, fptr);
    if  (!strncmp("-walk",gname,5)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      *walk = (FPN) atof( gname );
    }
    else if (!strncmp("-otraits",gname,8)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      *otraits =   atoi( gname );
    }
    else if (!strncmp("-traits",gname,7)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      *traits =   atoi( gname );
    }
    else if (!strncmp("-n",gname,2)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      *n =   atoi( gname );
    }
    else if (!strncmp("-Trait",gname,6)) 
      ch = EOF;
  }
  fileclose(infile, fptr);

}

/*  Create an index vector of where the expected values
are in the file.   In this way, we won't need to 
load all the data into memory.   */
long *GetPositions(char *infile,int ad) {
  FILE *fptr;
  long  fileposition, *positions;
  int whpos,ch;
  
  fptr = fileopen(infile, "r");
  if ( fptr == NULL )
    return(NULL);
  ch = 1;
  while ( ch != EOF ) {
    ch = get_next_token(gname, MAXNAME, fptr);
    
    if  (!strncmp("-positions",gname,10)) {
      ch = get_next_token(gname, MAXNAME, fptr);
      whpos =  atoi( gname );
      positions = lvector(0,whpos);
      positions[0] =  (long)  whpos;
    }
    else if (!strncmp("-Site",gname,5)) {
      fileposition = ftell(fptr);
      ch = get_next_token(gname, MAXNAME, fptr);
      whpos = atoi( gname );
    }
    else if (!strncmp("-parameter",gname,10) ) {    
      ch = get_next_token(gname, MAXNAME, fptr);
      if ( ( ad == 1 && !strncmp("additive",gname,8)  ) || ( ad == 2 && !strncmp("dominance",gname,8) ) ) {
	      if ( positions == NULL ) {
	        fileclose(infile, fptr);
	        return(NULL);
	      }
	      else if ( whpos <= positions[0] ) 
	        positions[whpos] = fileposition; 
	      else 
	        ch = EOF ;
      }
    }
  }
  fileclose(infile, fptr);
  return( positions );
}

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MRfunc.c
------------------------------------------------------------------ */


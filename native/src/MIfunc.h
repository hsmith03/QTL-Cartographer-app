/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MIfunc.h) is part of QTL Cartographer
         
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


void WhichMultilocusGT(int loci, long j, int *genotype, int ngt);
long LongMultilocusGT(int loci, int *genotype, int ngt);
void copy_traits(params *theparams,markermap *themap,individual *individs,linpakws *lnpk);

void UpdateGTindex(mimparam *mimptr);
mimparam *ConstructMIMparams(params *theparams,markermap *themap,aqtl *theqtls);
mimparam *insert_mimparam(mimparam *inptr, int abdtype, FPN initvalue, int qtl1, int qtl2);
mimparam *delete_mimparam(mimparam *inptr);
void insert_epistatic(mimparam *tmimptr, mimparam *rmimptr, mimparam *smimptr, int abdtype );
void PushMIMnode(mimparam *lptr);
void PullMIMnode(mimparam *lptr);
mimparam *WhichMIMnode(mimparam *mimptr, int qtl1,int qtl2, int abd,int indic);
void QTLpostions(params *theparams,markermap *themap,mimparam *mimptr,genome *agptr);

void init_design(mimparam *gptr, int *genotype, int crosstype);
void mimplace_qtls(params *theparams,mimparam *mimptr,genome *first);
int HowManyQTL(mimparam *mimptr,int abd);
mimparam *TerminalMIMnode(mimparam *mimptr);
int CheckForEpistatics(mimparam *mimptr, mimparam *aptr);


void write_mimheader(char *outfile, params *theparams, char *onamae, char *chptr, markermap *themap, int oc);
void do_mimanalysis(params *theparams, markermap *themap, aqtl *theqtls, individual *individs,  genome *agptr, linpakws *lnpk);

FPN EvalModel(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);
void CalculatePriors(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);
void EstimateParameters(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);

void FindMoreQTL(FILE *outf,params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk,int nparams );
int ScanGenomeQTL(FILE *outf,params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);
void SearchForEpistasis(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk,int nparams,FPN initlogln);
int BackElimEpistasis(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk,int nparams,FPN initlogln);
int SearchForDominance(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);
FPN CheckEpistasis(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,mimparam *rmimptr, mimparam *smimptr,mimparam *lmimptr,genome *agptr,linpakws *lnpk,int epitype);
void CheckThisModel(FILE *outf,params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) ;
void ConfirmCurrentModel(FILE *outf,params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk) ;
int RefineEstimates(params *theparams,mimparam *mimptr, linpakws *lnpk);
void  JitterPositions(params *theparams,markermap *themap ,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk );
FPN  DoAJitter(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,mimparam *lptr,genome *agptr,linpakws *lnpk);
void  TestAdjacentIntervals(params *theparams,markermap *themap ,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);
int isclearint(params *theparams, mimparam *mimptr,mimparam *lptr,genome *agptr) ;

void PrintModelHeader(FILE *outf);
void PrintFinalModel(FILE *outf, params *theparams,markermap *themap,mimparam *mimptr,char *outstring,int nn);
void PrintMIMnode(FILE *outf,mimparam *gptr );
void PrintRqtlout(FILE *outf, params *theparams,markermap *themap,mimparam *mimptr,int position);

void CalculateResiduals(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);
void CalculatePriors(params *theparams,markermap *themap,individual *individs,mimparam *mimptr,genome *agptr,linpakws *lnpk);
void CalcPriorMatrix(params *theparams,genome *gptr,individual *individs,int jj, linpakws *lnpk, mimparam *mptr);
void ShowPriorMatrix( params *theparams, linpakws *lnpk, mimparam *mptr);
void ShowSigGenotypes(  params *theparams, linpakws *lnpk);

void   CalcMeanVar(params *theparams,mimparam *mimptr,linpakws *lnpk);
void   CalcGTvalues(params *theparams,mimparam *mimptr,linpakws *lnpk);
void   ShowGTvalues(FILE *outf,params *theparams,linpakws *lnpk,individual *individs);
FPN CalcVarCovar(params *theparams,mimparam *mimptr,linpakws *lnpk);
void   ShowVarCovar(FILE *outf, params *theparams,markermap *themap,mimparam *mimptr,linpakws *lnpk,FPN genvar);

FPN CalcIC( FPN logln, int nparams, int nn,int whoseic);
void   MStepE(params *theparams,mimparam *mimptr,linpakws *lnpk);
void   EStepPi(params *theparams,mimparam *mimptr,linpakws *lnpk) ;
int    ShouldWeContinue(mimparam *mimptr,int jj);

/*  PROTOTYPES    from Otto and Jones*/
int OttoJones(params *theparams,FILE *fileptr, FPN *fParentalDiff, char *scUseEstThress,mimparam *mimptr );
FPN esttotalqtl(FPN fParentalDiff, FPN fMeanEffect, FPN dMinEffectUsed);
FPN estimatethreshold(FPN fNumQtlFound, FPN fMeanEffect, FPN fMinEffect);
FPN estimatelowerbound(FPN fNumQtlFound, FPN fMeanEffect, FPN dMinEffectUsed, FPN fParentalDiff, FPN fTotNumQTl);
FPN estimateupperbound(FPN fNumQtlFound, FPN fMeanEffect, FPN dMinEffectUsed, FPN fParentalDiff);
FPN meaneffectundetected(FPN fMeanEffect, FPN fTau);
FPN segvarundetected(FPN fTau);
FPN calctau(FPN fMeanEffect, FPN dMinEffectUsed);

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MIfunc.h
------------------------------------------------------------------ */


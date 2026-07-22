/* ------------------------------------------------------ XCutXCodeXSkip
     This file (MissMark.h) is part of QTL Cartographer
         
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



void a_Mv(FPN *a, FPN (*M)[3], FPN *v);
void AinR(FPN (*rr)[3], FPN (*aa)[3]);
void AdotB(FPN (*rr)[3], FPN (*aa)[3], FPN (*bb)[3]);

void cond_prob(params *theparams, FPN (*Muv)[3], FPN *quk, FPN *pRuk, FPN *pLuk);
void calc_cond_prob(params *theparams, genome *gptr, individual *individs, int kk, FPN *pAA, FPN *pAa, FPN *paa);
void markov_rec(params *theparams, FPN (*mm)[3], FPN rec, int rl);
void assign_q(FPN *ql, params *theparams);

FPN expect_mark(params *theparams, individual *individs, int kk, genome *gptr);
void calc_priors(params *theparams, FPN *pp1, FPN *pp2, int nn, individual *individs, genome *gptr, FPN thetaL, FPN thetaR);
void selfed_f_tpm(FPN (*sftpm)[3], FPN r, int t);
void bselfed_f_tpm(FPN (*sftpm)[3], FPN r, int t, int bc);

void add_virtual_node(genome *tgptr, genome *lgptr, genome *rgptr, FPN thetaL, FPN thetaR);
void del_virtual_node(genome *tgptr);

 

/* ------------------------------------------------------- XCutXCodeXSkip
             End of file MissMark.h
------------------------------------------------------------------ */


(* ::Package:: *)

(*****************************
* dmformfactor.m
* Package written by Nikhil Anand, A. Liam Fitzpatrick, and Wick Haxton.
* 
* For documentation, see http://arxiv.org/abs/arXiv:13xx.xxxx
* 
* *)

BeginPackage["dmformfactor`"]

(* ::Text:: *)
(*Wick and Celia's Library *)


(******************
* Section for Wick and Cecilia's Operators
*)

(***Triangular inequalities and 3,6, and 9-J Symbols***)
minus=-1;
JNorm[j_]:=2j+1;
QNorm[j_]:=Sqrt[2 j+1];
Triangular[x_,y_,z_]:=Abs[x-y]<=z<=x+y
RangeM[j_,m_]:=-Abs[j]<=m<=Abs[j]

(*Conditions for non-zero 6-J Symbol*)
SixJConditionTriad[j1_,j2_,j3_]:=IntegerQ[j1+j2+j3]&&Triangular[j1,j2,j3]
SixJCondition[j1_,j2_,j3_,J1_,J2_,J3_]:=SixJConditionTriad[j1,j2,j3]&&SixJConditionTriad[j1,J2,J3]&&SixJConditionTriad[J1,j2,J3]&&SixJConditionTriad[J1,J2,j3]
SixJ[{j1_,j2_,j3_},{J1_,J2_,J3_}]:=Which[SixJCondition[j1,j2,j3,J1,J2,J3],SixJSymbol[{j1,j2,j3},{J1,J2,J3}],True,0]

(*Likewise for 3-J*)
ThreeJCondition[{j1_,m1_},{j2_,m2_},{j3_,m3_}]:=Triangular[j1,j2,j3]&&RangeM[j1,m1]&&RangeM[j2,m2]&&RangeM[j3,m3]&&m1+m2+m3==0
ThreeJ[{j1_,m1_},{j2_,m2_},{j3_,m3_}]:=Which[ThreeJCondition[{j1,m1},{j2,m2},{j3,m3}],ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}],True,0]

(*And the 9-J*)
cutsum=50;
SummandNineJ[j1_,j2_,j12_,j3_,j4_,j34_,j13_,j24_,j_,g_]:=minus^(2*g) JNorm[g] SixJ[{j1,j2,j12},{j34,j,g}]*SixJ[{j3,j4,j34},{j2,g,j24}]*SixJ[{j13,j24,j},{g,j1,j3}]
NineJSymbol[{j1_,j2_,j12_},{j3_,j4_,j34_},{j13_,j24_,j_}]:=Sum[SummandNineJ[j1,j2,j12,j3,j4,j34,j13,j24,j,g],{g,0,cutsum+1/2,1/2}]

(**Auxillary functions**)
BesselFactor1[y_,{np_,lp_},{n_,l_},lcap_]:=(2^lcap)/(JNorm[lcap]!!) y^(lcap/2) Exp[-y] Sqrt[(np-1)! (n-1)!]
BesselFactor2[y_,{np_,lp_},{n_,l_},lcap_]:=Sqrt[Gamma[np+lp+1/2]*Gamma[n+l+1/2]]
Summand1[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=minus^(m+mp)/(m! mp! (n-1-m)! (np-1-mp)!)
Summand2[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=Gamma[(l+lp+lcap+2*m+2*mp+3)/2]/(Gamma[l+m+3/2]*Gamma[lp+mp+3/2])
Summand3[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=Hypergeometric1F1[(lcap-lp-l-2*mp-2*m)/2,lcap+3/2,y]
BesselFactor3[y_,{np_,lp_},{n_,l_},lcap_]:=Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*Summand2[y,{np,lp,mp},{n,l,m},lcap]*Summand3[y,{np,lp,mp},{n,l,m},lcap]),{m,0,n-1},{mp,0,np-1}]
BesselElement[y_,{np_,lp_},{n_,l_},lcap_]:=BesselFactor1[y,{np,lp},{n,l},lcap]*BesselFactor2[y,{np,lp},{n,l},lcap]*BesselFactor3[y,{np,lp},{n,l},lcap]
BesselFactor1A[y_,{np_,lp_},{n_,l_},lcap_]:=(2^(lcap-1))/(JNorm[lcap]!!) y^((lcap-1)/2) Exp[-y] Sqrt[(np-1)! (n-1)!]
BesselFactor2[y_,{np_,lp_},{n_,l_},lcap_]:=Sqrt[Gamma[np+lp+1/2]*Gamma[n+l+1/2]]

(**)
Summand1[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=(-1)^(m+mp)/(m! mp! (n-1-m)! (np-1-mp)!)
Summand2A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=Gamma[(l+lp+lcap+2*m+2*mp+2)/2]/(Gamma[l+m+3/2]*Gamma[lp+mp+3/2])
Summand3A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=-(l+lp+lcap+2*m+2*mp+2)/2*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m-1)/2,lcap+3/2,y]+2*m*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m+1)/2,lcap+3/2,y]
BesselFactor3A[y_,{np_,lp_},{n_,l_},lcap_]:=Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*Summand2A[y,{np,lp,mp},{n,l,m},lcap]*Summand3A[y,{np,lp,mp},{n,l,m},lcap]),{m,0,n-1},{mp,0,np-1}]

(**)
BesselElementMinus[y_,{np_,lp_},{n_,l_},lcap_]:=BesselFactor1A[y,{np,lp},{n,l},lcap]*BesselFactor2[y,{np,lp},{n,l},lcap]*BesselFactor3A[y,{np,lp},{n,l},lcap]

(**)
Summand4A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=-(l+lp+lcap+2*m+2*mp+2)/2*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m-1)/2,lcap+3/2,y]+(2*l+2*m+1)*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m+1)/2,lcap+3/2,y]
BesselFactor4A[y_,{np_,lp_},{n_,l_},lcap_]:=Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*Summand2A[y,{np,lp,mp},{n,l,m},lcap]*Summand4A[y,{np,lp,mp},{n,l,m},lcap]),{m,0,n-1},{mp,0,np-1}]
BesselElementPlus[y_,{np_,lp_},{n_,l_},lcap_]:=BesselFactor1A[y,{np,lp},{n,l},lcap]*BesselFactor2[y,{np,lp},{n,l},lcap]*BesselFactor4A[y,{np,lp},{n,l},lcap]
Lnumber[NPrincipal_,j_]:=Which[EvenQ[NPrincipal-(j+1/2)],(j+1/2),True,j-1/2]
Nodal[NPrincipal_,j_]:=(NPrincipal-Lnumber[NPrincipal,j])/2+1


(***Parity and Parity Conservation***)
ParityState[NPrincipal_,j_]:=minus^(Lnumber[NPrincipal,j])
ParityNormal[jcap_]:=minus^(jcap)
ParityConsNormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityState[NPrincipal,j]*ParityState[NPrincipalp,jp]*ParityNormal[jcap]
PhysicalConditionsNormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityConsNormal[NPrincipal,j,NPrincipalp,jp,jcap]==1&&Nodal[NPrincipal,j]>0&&Nodal[NPrincipalp,jp]>0&&Abs[j-jp]<=jcap<=j+jp
ParityAbnormal[jcap_]:=minus^(jcap+1)
ParityConsAbnormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityState[NPrincipal,j]*ParityState[NPrincipalp,jp]*ParityAbnormal[jcap]
PhysicalConditionsAbnormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityConsAbnormal[NPrincipal,j,NPrincipalp,jp,jcap]==1&&Nodal[NPrincipal,j]>0&&Nodal[NPrincipalp,jp]>0&&Abs[j-jp]<=jcap<=j+jp


(***Normal Parity Operators***)

(**MJ**)
mjelement[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=minus^(1/2+j+jcap) Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap])/(4 Pi)]*ThreeJ[{lp,0},{jcap,0},{l,0}]*SixJ[{lp,jp,1/2},{j,l,jcap}]*BesselElement[y,{np,lp},{n,l},jcap]
MJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=mjelement[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]
MJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],MJE[y,{ncapp,jp},{ncap,j},jcap],0]

(**SigmaJ**)
MJLSigma[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=minus^lp Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap]*JNorm[lcap])/(4 Pi)]*Sqrt[6]*ThreeJ[{lp,0},{lcap,0},{l,0}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}]*BesselElement[y,{np,lp},{n,l},lcap]
SigmaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=MJLSigma[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap,jcap]
SigmaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],SigmaJE[y,{ncapp,jp},{ncap,j},jcap],0]

(**PhiPPJ**)
PhiPPoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=minus^(lp+1)*6 *QNorm[lp]*QNorm[jp]*QNorm[j]/Sqrt[4 Pi]
PhiPPsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=QNorm[l+1]*Sqrt[l+1]*ThreeJ[{lp,0},{jcap+1,0},{l+1,0}]*BesselElementMinus[y,{np,lp},{n,l},jcap+1]*Sum[minus^(jcap-lcap+1)*JNorm[lcap]*SixJ[{jcap+1,1,lcap},{1,jcap,1}]*SixJ[{jcap+1,1,lcap},{l,lp,l+1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap,jcap+1}]
PhiPPsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=If[l==0,0,QNorm[l-1]* Sqrt[l]*ThreeJ[{lp,0},{jcap+1,0},{l-1,0}]*BesselElementPlus[y,{np,lp},{n,l},jcap+1]*Sum[minus^(jcap-lcap)*JNorm[lcap]*SixJ[{jcap+1,1,lcap},{1,jcap,1}]*SixJ[{jcap+1,1,lcap},{l,lp,l-1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap,jcap+1}]]
PhiPPsummand3[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=If[jcap==0,0,QNorm[l+1]*Sqrt[l+1]*ThreeJ[{lp,0},{jcap-1,0},{l+1,0}]*BesselElementMinus[y,{np,lp},{n,l},jcap-1]*Sum[minus^(jcap-lcap+1)*JNorm[lcap]*SixJ[{jcap-1,1,lcap},{1,jcap,1}]*SixJ[{jcap-1,1,lcap},{l,lp,l+1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap-1,jcap}]]
PhiPPsummand4[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=If[jcap==0,0,If[l==0,0,QNorm[l-1]* Sqrt[l]*ThreeJ[{lp,0},{jcap-1,0},{l-1,0}]*BesselElementPlus[y,{np,lp},{n,l},jcap-1]*Sum[minus^(jcap-lcap)*JNorm[lcap]*SixJ[{jcap-1,1,lcap},{1,jcap,1}]*SixJ[{jcap-1,1,lcap},{l,lp,l-1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap-1,jcap}]]]
PhiPPJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=PhiPPoverall[y,{np,lp,jp},{n,l,j},jcap]*(QNorm[jcap+1]*Sqrt[jcap+1]*(PhiPPsummand1[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand2[y,{np,lp,jp},{n,l,j},jcap])+If[jcap==0,0,QNorm[jcap-1]*Sqrt[jcap]*(PhiPPsummand3[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand4[y,{np,lp,jp},{n,l,j},jcap])])
PhiPPJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiPPJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]
PhiPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],PhiPPJE[y,{ncapp,jp},{ncap,j},jcap],0]

(**PhiTPJ**)
PhiPJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=PhiPPoverall[y,{np,lp,jp},{n,l,j},jcap]*(-QNorm[jcap+1]*Sqrt[jcap]*(PhiPPsummand1[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand2[y,{np,lp,jp},{n,l,j},jcap])+If[jcap==0,0,QNorm[jcap-1]*Sqrt[jcap+1]*(PhiPPsummand3[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand4[y,{np,lp,jp},{n,l,j},jcap])])
PhiPJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiPJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]
PhiPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],PhiPJE[y,{ncapp,jp},{ncap,j},jcap],0]
PhiTPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiPJ[y,{ncapp,jp},{ncap,j},jcap]+SigmaJ[y,{ncapp,jp},{ncap,j},jcap]/2


(***Operators of Abnormal Parity***)

(**DeltaJ**)
MJLDivQoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=minus^(lcap+j+1/2)*QNorm[lp]*QNorm[jp]*QNorm[j]*QNorm[jcap]*QNorm[lcap]*SixJ[{lp,jp,1/2},{j,l,jcap}]/Sqrt[4 Pi]
MJLDivQsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=-Sqrt[l+1] QNorm[l+1] SixJ[{lcap,1,jcap},{l,lp,l+1}] ThreeJ[{lp,0},{lcap,0},{l+1,0}] BesselElementMinus[y,{np,lp},{n,l},lcap]
MJLDivQsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=Sqrt[l] QNorm[l-1] SixJ[{lcap,1,jcap},{l,lp,l-1}] ThreeJ[{lp,0},{lcap,0},{l-1,0}] BesselElementPlus[y,{np,lp},{n,l},lcap]
MJLDivQ[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=MJLDivQoverall[y,{np,lp,jp},{n,l,j},jcap,lcap]*(MJLDivQsummand1[y,{np,lp,jp},{n,l,j},jcap,lcap]+MJLDivQsummand2[y,{np,lp,jp},{n,l,j},jcap,lcap])

DeltaPrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=-Sqrt[jcap]/QNorm[jcap]*MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+Sqrt[jcap+1]/QNorm[jcap]*MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap-1]
Deltaop[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap]
DeltaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=Deltaop[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]
DeltaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],DeltaJE[y,{ncapp,jp},{ncap,j},jcap],0]

(**SigmaPJ**)
SigmaPrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=-Sqrt[jcap]/QNorm[jcap]*MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+Sqrt[jcap+1]/QNorm[jcap]*MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap-1]
SigmaPrimeJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaPrime[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]
SigmaPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaPrimeJE[y,{ncapp,jp},{ncap,j},jcap],0]

(**SigmaPPJ**)
SigmaSecond[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=(Sqrt[jcap+1] MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+Sqrt[jcap] MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap-1])/QNorm[jcap]
SigmaSecondJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaSecond[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]
SigmaPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaSecondJE[y,{ncapp,jp},{ncap,j},jcap],0]




















(* ::Text:: *)
(*Section for the Responses*)


(* ::Input:: *)
(*(* "Response" is a function of the nucleus and the operator whose matrix elements we want to calculate.  It is convenient to keep the size "DMmatrixLength" of the density matrix and the isospin "T" as additional arguments, though these could conceivably be merged with "DMmatrix" into a single argument. *)
(**)
(*"Response" needs to be evaluated the first time with "DMmatrixLength" and "T" as inputs in addition to "DMmatrix" and "Operator", but from then on it can be evaluated with just "DMmatrix" and "Operator". *)*)


Jbase[0]=j0;Jbase[1]=j1;Jbase[2]=j2;Jbase[3]=j3;Jbase[4]=j4;
Jbase[5]=j5;Jbase[6]=j6;Jbase[7]=j7;Jbase[8]=j8;Jbase[9]=j9;
Jbase[10]=j10;Jbase[11]=j11;Jbase[12]=j12;Jbase[13]=j13;Jbase[14]=j14;
Jmax=14;


(*Initially, the J and T for the isotope set to 0. This is changed
once the density matrix is read in.*)
JIso=0;
TIso=0;



(*******************
* Section for basic nuclear form factors
*)
(*Basic form factors for non-interference and interference between operators*)

FF[DMmatrix_,Operator_,J_,T_]:=Block[{Response},

If[Operator==MJ&&UseHelm==True,Return[FHelm={{ZZ^2,ZZ(M-ZZ)},{ZZ(M-ZZ),(M-ZZ)^2}}*(3SphericalBesselJ[1,2y^(1/2) r0Helm/b]/(2y^(1/2) r0Helm/b))^2*Exp[-4y s^2/b^2]];];

Response[iDMmatrix_,iOperator_,iT_]:=Response[iDMmatrix,iOperator,iT]=Simplify[Block[{},
Sum[iDMmatrix[i][[1]](If[iDMmatrix[i][[3]]==0,Sqrt[2]/Sqrt[2iT+1] ((an+ap)/2),-(Sqrt[6iT]/Sqrt[(2iT+1)(iT+1)])((ap-an)/2)])*E^y (iOperator[y,{iDMmatrix[i][[4]],iDMmatrix[i][[5]]/2},{iDMmatrix[i][[6]],iDMmatrix[i][[7]]/2},iDMmatrix[i][[2]]])Jbase[iDMmatrix[i][[2]]],{i,DMmatrixLength}]]];

Simplify[(4\[Pi])/(2J+1) Sum[1/2 Table[D[(Coefficient[Response[DMmatrix,Operator,T],Jbase[n]]^2),{ap,an}[[ii]],{ap,an}[[jj]]],{ii,2},{jj,2}],{n,0,Jmax}]]
];

FF[DMmatrix_,Operator1_,Operator2_,J_,T_]:=Block[{Response},

Response[iDMmatrix_,iOperator_,iT_]:=Response[iDMmatrix,iOperator,iT]=Simplify[Block[{},
Sum[iDMmatrix[i][[1]](If[iDMmatrix[i][[3]]==0,Sqrt[2]/Sqrt[2iT+1] ((an+ap)/2),-(Sqrt[6iT]/Sqrt[(2iT+1)(iT+1)])((ap-an)/2)])*E^y (iOperator[y,{iDMmatrix[i][[4]],iDMmatrix[i][[5]]/2},{iDMmatrix[i][[6]],iDMmatrix[i][[7]]/2},iDMmatrix[i][[2]]])Jbase[iDMmatrix[i][[2]]],{i,DMmatrixLength}]]];


Simplify[(4\[Pi])/(2J+1) Sum[Table[D[(Coefficient[Response[DMmatrix,Operator1,T]/.{ap->ap1,an->an1},Jbase[n]])(Coefficient[Response[DMmatrix,Operator2,T]/.{ap->ap2,an->an2},Jbase[n]]),{ap1,an1}[[ii]],{ap2,an2}[[jj]]],{ii,2},{jj,2}],{n,0,Jmax}]]
];



ResponseNuclear[yyy_,jjj_,tau1_,tau2_]:=Block[{Operator1, Operator2},
Operator1=0;
Operator2=0;
If[jjj<1||jjj>8||(!(jjj\[Element]Integers)), Print["Error: Operator index must be an integer between 1 to 8, corresponding to the eight response functions in Eqn. 39 of documentation."];
Return[];];
If[!((tau1===0||tau1===1)&&(tau2===0||tau2===1)), Print["Error: isospin indices must be either 0 or 1."];
Return[];];

If[jjj==1,Operator1=MJ;];
If[jjj==2,Operator1=SigmaPPJ;];
If[jjj==3,Operator1=SigmaPJ;];
If[jjj==4,Operator1=PhiPPJ;];
If[jjj==5,Operator1=PhiTPJ;];
If[jjj==6,Operator1=DeltaJ;];
If[jjj==7,Operator1=MJ; Operator2=PhiPPJ;];
If[jjj==8,Operator1=SigmaPJ; Operator2=DeltaJ;];

If[1<=jjj<=6,
If[tau1==0&&tau2==0,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(FF[DensityMatrix,Operator1,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,JIso,TIso][[1,1]]+FF[DensityMatrix,Operator1,JIso,TIso][[1,2]]+FF[DensityMatrix,Operator1,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];
If[tau1==1&&tau2==1,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(FF[DensityMatrix,Operator1,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,JIso,TIso][[1,1]]-FF[DensityMatrix,Operator1,JIso,TIso][[1,2]]-FF[DensityMatrix,Operator1,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];
If[tau1==0&&tau2==1,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(-FF[DensityMatrix,Operator1,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,JIso,TIso][[1,1]]-FF[DensityMatrix,Operator1,JIso,TIso][[1,2]]+FF[DensityMatrix,Operator1,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];
If[tau1==1&&tau2==0,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(-FF[DensityMatrix,Operator1,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,JIso,TIso][[1,1]]+FF[DensityMatrix,Operator1,JIso,TIso][[1,2]]-FF[DensityMatrix,Operator1,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];];

If[7<=jjj<=8,
If[tau1==0&&tau2==0,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,1]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,2]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];
If[tau1==1&&tau2==1,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,1]]-FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,2]]-FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];
If[tau1==0&&tau2==1,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(-FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,1]]-FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,2]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];
If[tau1==1&&tau2==0,Return[E^(-2y)*(2JIso+1)/(4Pi)*1/4(-FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,2]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,1]]+FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[1,2]]-FF[DensityMatrix,Operator1,Operator2,JIso,TIso][[2,1]])/.y->yyy/.FormalReplace//Simplify//MyChop]
];];

Return[];
];

(*When the user wants to turn on coefficients, they do it using the SetCoeffsNonrel function*)

UsingRelCoeffs=False;

ZeroCoeffs[]:=Block[{iii},
For[iii=1,iii<=12,iii++,
cpvector[[iii]]=0;
cnvector[[iii]]=0;
];
For[iii=1,iii<=20,iii++,
dpvector[[iii]]=0;
dnvector[[iii]]=0;
];
ellpnucl={{0}}~Join~Table[{0,0,0},{i,2,4}];
ellnnucl={{0}}~Join~Table[{0,0,0},{i,2,4}];
];


SetCoeffsNonrel[Op_,coeffdimless_,nucleon_]:=Block[{c0,c1,coeff},

If[UsingRelCoeffs,

Print[Style["WARNING: You are switching from using relativistic operators to non-relativistic operators.",32,Orange]];
Print[Style["All previous settings for coefficients will be discarded.",32,Orange]];
ZeroCoeffs[];
UsingRelCoeffs=False;
];

(*Convert from the old normalization of c coefficients to the new one *)
If[Op==1, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
If[Op==2, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
If[Op==3, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
If[Op==4, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
If[Op==5, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
If[Op==6, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN^2;];
If[Op==7, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
If[Op==8, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
If[Op==9, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
If[Op==10, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
If[Op==11, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
If[Op==12, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless;];
If[Op==13, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
If[Op==14, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN;];
If[Op==15, coeff=(4mN*mchiFORMAL/mV^2)coeffdimless/mN^2;];

If[nucleon!="n"&&nucleon!="p"&&nucleon!=0&&nucleon!=1,
Print["3rd argument must be one of \"p\",\"n\",0, or 1"];Return[];
];
If[!(Integer[Op]&&1<=Op&&Op<=15),
Print["1st argument must be an integer between 1 and 15"];Return[];
];
If[nucleon=="p",cpvector[[Op]]=coeff; Return[];];
If[nucleon=="n",cnvector[[Op]]=coeff; Return[];];
If[nucleon==0,
(*Store the old value of isoscalar,isovector coefficient*)
c0=cpvector[[Op]]+cnvector[[Op]];c1=cpvector[[Op]]-cnvector[[Op]];
(*Set the new value for isoscalar*)
c0=coeff;
(*Convert back to neutron, proton basis*)
cpvector[[Op]]=(c0+c1)/2//Expand;
cnvector[[Op]]=(c0-c1)/2//Expand;
Return[];];
If[nucleon==1,
(*Store the old value of isoscalar,isovector coefficient*)
c0=cpvector[[Op]]+cnvector[[Op]];c1=cpvector[[Op]]-cnvector[[Op]];
(*Set the new value for isoscalar*)
c1=coeff;
(*Convert back to neutron, proton basis*)
cpvector[[Op]]=(c0+c1)/2//Expand;
cnvector[[Op]]=(c0-c1)/2//Expand;
Return[];];
];

(*Same procedure, except if the user prefers the relativistic coefficients*)
SetCoeffsRel[Op_,coeffReldimless_,nucleon_]:=Block[{d0,d1,coeffRel},

(*If[(Op==12)||(Op==16)||(Op==18)||(Op==19),
 Print[Style["Error: Relativistic operators 12,16,18, and 19 have not been implemented! Ignoring function call.",32,Red]];
Return[];
];*)

If[!UsingRelCoeffs,

Print[Style["WARNING: You are switching from using non-relativistic operators to relativistic operators.",32,Orange]];
Print[Style["All previous settings for coefficients will be discarded and DM spin will be set to 1/2.",32,Orange]];
ZeroCoeffs[];
SetJChi[1/2];
UsingRelCoeffs=True;
];

If[Op==1, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==2, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==3, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==4, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==5, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==6, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==7, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==8, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==9, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==10, coeffRel=coeffReldimless/mMFORMAL^4;];
If[Op==11, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==12, coeffRel=coeffReldimless/mMFORMAL^4;];
If[Op==13, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==14, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==15, coeffRel=coeffReldimless/mMFORMAL^2;];
If[Op==16, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==17, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==18, coeffRel=coeffReldimless/mMFORMAL^4;];
If[Op==19, coeffRel=coeffReldimless/mMFORMAL^3;];
If[Op==20, coeffRel=coeffReldimless/mMFORMAL^4;];

If[nucleon!="n"&&nucleon!="p"&&nucleon!=0&&nucleon!=1,
Print["3rd argument must be one of \"p\",\"n\",0, or 1"];Return[];
];
If[!(Integer[Op]&&1<=Op&&Op<=20),
Print["1st argument must be an integer between 1 and 20"];Return[];
];
If[nucleon=="p",dpvector[[Op]]=coeffRel; ConvertRel[]; Return[];];
If[nucleon=="n",dnvector[[Op]]=coeffRel; ConvertRel[]; Return[];];
If[nucleon==0,
d0=dpvector[[Op]]+dnvector[[Op]];
d1=dpvector[[Op]]-dnvector[[Op]];
d0=coeffRel;
dpvector[[Op]]=(d0+d1)/2//Expand;
dnvector[[Op]]=(d0-d1)/2//Expand;
ConvertRel[];
Return[];];
If[nucleon==1,
d0=dpvector[[Op]]+dnvector[[Op]];d1=dpvector[[Op]]-dnvector[[Op]];
d1=coeffRel;
dpvector[[Op]]=(d0+d1)/2//Expand;
dnvector[[Op]]=(d0-d1)/2//Expand;
ConvertRel[];
Return[];];
];

ConvertRel[]:=Block[{cpvectorTEMP,cnvectorTEMP},
cnvectorTEMP=Table[0,{i,15}];
cpvectorTEMP=Table[0,{i,15}];
(*Now determine the non-relativistic coefficients*)
cpvectorTEMP[[1]]=(4 mchiFORMAL mN)(dpvector[[1]]+dpvector[[5]]+q^2/(2mN)dpvector[[6]]-q^2/(2mchiFORMAL)dpvector[[9]]);
cnvectorTEMP[[1]]=(4 mchiFORMAL mN)(dnvector[[1]]+dnvector[[5]]+q^2/(2mN)dnvector[[6]]-q^2/(2mchiFORMAL)dnvector[[9]]);
cpvectorTEMP[[2]]=0;
cnvectorTEMP[[2]]=0;
cpvectorTEMP[[3]]=(1/mN)(4 mchiFORMAL mN)(-2mN dpvector[[6]]);
cnvectorTEMP[[3]]=(1/mN)(4 mchiFORMAL mN)(-2mN dnvector[[6]]);
cpvectorTEMP[[4]]=(4 mchiFORMAL mN)(2q^2/mchiFORMAL dpvector[[6]]-2q^2/mN dpvector[[9]]+4q^2 dpvector[[10]]-4dpvector[[15]]);
cnvectorTEMP[[4]]=(4 mchiFORMAL mN)(2q^2/mchiFORMAL dnvector[[6]]-2q^2/mN dnvector[[9]]+4q^2 dnvector[[10]]-4dnvector[[15]]);
cpvectorTEMP[[5]]=(1/mN)(4 mchiFORMAL mN)(2mN dpvector[[9]]);
cnvectorTEMP[[5]]=(1/mN)(4 mchiFORMAL mN)(2mN dnvector[[9]]);
cpvectorTEMP[[6]]=(1/mN^2)(4 mchiFORMAL mN)(-mN/mchiFORMAL dpvector[[4]]-2mN^2/mchiFORMAL dpvector[[6]]+2mN dpvector[[9]]-4mN^2 dpvector[[10]]+4mN^2dpvector[[20]]);
cnvectorTEMP[[6]]=(1/mN^2)(4 mchiFORMAL mN)(-mN/mchiFORMAL dnvector[[4]]-2mN^2/mchiFORMAL dnvector[[6]]+2mN dnvector[[9]]-4mN^2 dnvector[[10]]+4mN^2dnvector[[20]]);
cpvectorTEMP[[7]]=(4 mchiFORMAL mN)(-2dpvector[[7]]);
cnvectorTEMP[[7]]=(4 mchiFORMAL mN)(-2dnvector[[7]]);
cpvectorTEMP[[8]]=(4 mchiFORMAL mN)(2dpvector[[13]]);
cnvectorTEMP[[8]]=(4 mchiFORMAL mN)(2dnvector[[13]]);
cpvectorTEMP[[9]]=(1/mN)(4 mchiFORMAL mN)(2mN/mchiFORMAL dpvector[[7]]+4mN dpvector[[11]]+2dpvector[[13]]-4mN dpvector[[14]]);
cnvectorTEMP[[9]]=(1/mN)(4 mchiFORMAL mN)(2mN/mchiFORMAL dnvector[[7]]+4mN dnvector[[11]]+2dnvector[[13]]-4mN dnvector[[14]]);
cpvectorTEMP[[10]]=(1/mN)(4 mchiFORMAL mN)(dpvector[[2]]+2mN dpvector[[8]]-mN/mchiFORMAL q^2 dpvector[[12]]);
cnvectorTEMP[[10]]=(1/mN)(4 mchiFORMAL mN)(dnvector[[2]]+2mN dnvector[[8]]-mN/mchiFORMAL q^2 dnvector[[12]]);
cpvectorTEMP[[11]]=(1/mN)(4 mchiFORMAL mN)(-mN/mchiFORMAL dpvector[[3]]+2mN dpvector[[17]]+q^2 dpvector[[18]]);
cnvectorTEMP[[11]]=(1/mN)(4 mchiFORMAL mN)(-mN/mchiFORMAL dnvector[[3]]+2mN dnvector[[17]]+q^2 dnvector[[18]]);
cpvectorTEMP[[12]]=(4 mchiFORMAL mN)(-4q^2 dpvector[[12]]);
cnvectorTEMP[[12]]=(4 mchiFORMAL mN)(-4q^2 dnvector[[12]]);
cpvectorTEMP[[13]]=(1/mN)(4 mchiFORMAL mN)(4mN dpvector[[16]]);
cnvectorTEMP[[13]]=(1/mN)(4 mchiFORMAL mN)(4mN dnvector[[16]]);
cpvectorTEMP[[14]]=(1/mN)(4 mchiFORMAL mN)(-4mN dpvector[[19]]);
cnvectorTEMP[[14]]=(1/mN)(4 mchiFORMAL mN)(-4mN dnvector[[19]]);
cpvectorTEMP[[15]]=(1/mN^2)(4 mchiFORMAL mN)(-4mN^2 dpvector[[12]]+4mN^2 dpvector[[18]]);
cnvectorTEMP[[15]]=(1/mN^2)(4 mchiFORMAL mN)(-4mN^2 dnvector[[12]]+4mN^2 dnvector[[18]]);
(*c's used internally are in the old normalization so there is no need to change*)
cpvector=cpvectorTEMP;
cnvector=cnvectorTEMP;
Return[];];


(*Define constants*)

mN:=0.938272 GeV
(*mn and mp are close enough together that we just use a single variable mN for them*)
mV:=246.2GeV;

(*For a dark matter particle with general spin, Cl refers to an effective quadratic casimir*)
Cl[jin_]:=(4jin(jin+1))/3
(*M refers to isotope mass number (# of nucleons), 
mchi is the dark matter particle mass and \[Mu]T is the reduced mass
of the target*)
\[Mu]T[mchi_, M_]:= mchi M mN/(mchi+M mN)


(*
	O. Gorton Notes
	This is equation (40)
	The ResponseCoeff equations are related to equation (38)
*)
FFfinal[DMmatrix_,J_,T_]:=Block[{iiFF,jjFF,ResponseCoeff,iiResponse},

cvector[1]=cpvector;
cvector[2]=cnvector;

(*These response coefficients differ in their definition from those in the paper 
because they have old c's instead of new c's, resulting in an overall factor of 
1/(4mN mchi)^2.  This will show up as an additional factor that we multiply
FF by before printing out StrucFunction, TransitionProbability, and 
ResponseNucl *)

ResponseCoeff[MJ]=Table[
	1/4 Cl[jchi](
		(
			cvector[iiFF][[5]]cvector[jjFF][[5]] q^2
			+ cvector[iiFF][[8]]cvector[jjFF][[8]]
		)
		( v^2-q^2/(4(\[Mu]T[mchi,M])^2) )
		+ cvector[iiFF][[11]]cvector[jjFF][[11]]q^2
	)
	+ (
		cvector[iiFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2))
		+ cvector[iiFF][[1]]
	)
	(
		cvector[jjFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2))
		+ cvector[jjFF][[1]]
	),
	{iiFF,2},{jjFF,2}
];

ResponseCoeff[PhiPPJ]=Table[
	q^2/(16mN^2) Cl[jchi](
		cvector[iiFF][[12]]
		-cvector[iiFF][[15]]q^2
	)(
		cvector[jjFF][[12]]
		-cvector[jjFF][[15]]q^2
	)
	+ q^4/(4mN^2)cvector[iiFF][[3]]cvector[jjFF][[3]],
	{iiFF,2},{jjFF,2}
];

ResponseCoeff[MJ,PhiPPJ]=Table[
	q^2/(4mN)Cl[jchi]cvector[iiFF][[11]](
		cvector[jjFF][[12]]
		-cvector[jjFF][[15]]q^2
	) + q^2/(mN)cvector[jjFF][[3]](
		cvector[iiFF][[1]]
		+cvector[iiFF][[2]](
			v^2-q^2/(4(\[Mu]T[mchi,M])^2)
		)
	),
	{iiFF,2},{jjFF,2}
];

ResponseCoeff[PhiTPJ]=Table[
	q^2/(16mN^2)Cl[jchi](
		cvector[iiFF][[13]]cvector[jjFF][[13]]q^2
		+cvector[iiFF][[12]]cvector[jjFF][[12]]
	),
	{iiFF,2},{jjFF,2}
];
			
ResponseCoeff[SigmaPPJ]=Table[
		1/16Cl[jchi](
			cvector[iiFF][[6]]cvector[jjFF][[6]]q^4
			+ (
				cvector[iiFF][[13]]cvector[jjFF][[13]]q^2
				+cvector[iiFF][[12]]cvector[jjFF][[12]]
			)(v^2-q^2/(4(\[Mu]T[mchi,M])^2))
			+ 2cvector[iiFF][[4]]cvector[jjFF][[6]]q^2
			+ cvector[iiFF][[4]]cvector[jjFF][[4]]
		)
		+ 1/4 cvector[iiFF][[10]]cvector[jjFF][[10]]q^2,
		{iiFF,2},{jjFF,2}
];

ResponseCoeff[SigmaPJ]=Table[
	1/32 Cl[jchi](
		2cvector[iiFF][[9]]cvector[jjFF][[9]]q^2
		+(
			cvector[iiFF][[15]]cvector[jjFF][[15]]q^4
			+cvector[iiFF][[14]]cvector[jjFF][[14]]q^2
			-2cvector[iiFF][[12]]cvector[jjFF][[15]]q^2
			+cvector[iiFF][[12]]cvector[jjFF][[12]]
		)(v^2-q^2/(4(\[Mu]T[mchi,M])^2))
		+ 2cvector[iiFF][[4]]cvector[jjFF][[4]]
	)
	+1/8 (
		cvector[iiFF][[3]]cvector[jjFF][[3]]q^2+cvector[iiFF][[7]]cvector[jjFF][[7]]
	)(v^2-q^2/(4(\[Mu]T[mchi,M])^2)),
	{iiFF,2},{jjFF,2}
];

ResponseCoeff[DeltaJ]=Table[
	q^2/(4mN^2) Cl[jchi](
		cvector[iiFF][[5]]cvector[jjFF][[5]]q^2
		+cvector[iiFF][[8]]cvector[jjFF][[8]])
		+2 q^2/mN^2 cvector[iiFF][[2]]cvector[jjFF][[2]](
			v^2-q^2/(4(\[Mu]T[mchi,M])^2)
		),
	{iiFF,2},{jjFF,2}
];

ResponseCoeff[SigmaPJ,DeltaJ]=Table[
	q^2/(4mN)Cl[jchi](
		cvector[iiFF][[4]]cvector[jjFF][[5]]
		-cvector[jjFF][[8]]cvector[iiFF][[9]]
	) - q^2/(mN)cvector[jjFF][[2]]cvector[iiFF][[3]](
		v^2
		-q^2/(4(\[Mu]T[mchi,M])^2)
	),
	{iiFF,2},{jjFF,2}
];

Sum[
	ResponseCoeff[MJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,J,T][[iiFF,jjFF]]
	+ ResponseCoeff[SigmaPPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPPJ,J,T][[iiFF,jjFF]]
	+ ResponseCoeff[SigmaPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,J,T][[iiFF,jjFF]]
	+ ResponseCoeff[DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,DeltaJ,J,T][[iiFF,jjFF]]
	+ ResponseCoeff[PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiPPJ,J,T][[iiFF,jjFF]]
	+ ResponseCoeff[PhiTPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiTPJ,J,T][[iiFF,jjFF]]
	+ ResponseCoeff[MJ,PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,PhiPPJ,J,T][[iiFF,jjFF]]
	+ ResponseCoeff[SigmaPJ,DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,DeltaJ,J,T][[iiFF,jjFF]],
	{iiFF,2},{jjFF,2}
]

]



(**Define another nuclear form factor for cross-check purposes(?), via Eqns. 36,38,39 of Wick's version**)

FFfinalnucl2[DMmatrix_,J_,T_]:=Block[{iiFF,jjFF,DMResponseCoeff,iiResponse},

cvector[1]=cpvector;
cvector[2]=cnvector;
(*These response coefficients differ in their definition from those in the paper 
because they have c's instead of a's, resulting in an overall factor of 
1/(4mN mchi)^2.  This will show up as an additional factor that we multiply
FF by before printing out StrucFunction, TransitionProbability, and 
ResponseNucl *)


DMResponseCoeff[MJ] = Table[ 
	(cvector[iiFF][[1]]+v^2 cvector[iiFF][[2]])
	(cvector[jjFF][[1]]+v^2 cvector[jjFF][[2]])
		+ Cl[jchi]/4 ( 
			q^2/mN^2 v^2 cvector[iiFF][[5]]cvector[jjFF][[5]] 
			+ v^2 cvector[iiFF][[8]]cvector[jjFF][[8]] 
			+ q^2/mN^2 cvector[iiFF][[11]]cvector[jjFF][[11]]
		),
	{iiFF,2},{jjFF,2}];

DMResponseCoeff[PhiPPJ] = Table[ 
	q^2/(4mN^2)cvector[iiFF][[3]] cvector[jjFF][[3]] 
	+ Cl[jchi]/16 cvector[iiFF][[12]]cvector[jjFF][[12]],
	{iiFF,2},{jjFF,2}
];

DMResponseCoeff[MJ,PhiPPJ] = Table[
	cvector[iiFF][[3]](
		cvector[jjFF][[1]]
		+ v^2 cvector[jjFF][[2]])
		+ Cl[jchi]/4 cvector[iiFF][[12]]cvector[jjFF][[11]],
	{iiFF,2},{jjFF,2}
];

(*Matches (38)*)
DMResponseCoeff[PhiTPJ] = Table[
	1/16 Cl[jchi](
		cvector[iiFF][[12]] cvector[jjFF][[12]] 
		+ q^2/mN^2 cvector[iiFF][[13]]cvector[jjFF][[13]]
	),
	{iiFF,2},{jjFF,2}
];

DMResponseCoeff[SigmaPPJ] = Table[
	q^2/(4mN^2) cvector[iiFF][[10]]cvector[jjFF][[10]] 
	+ 1/16 Cl[jchi](
		cvector[iiFF][[4]]cvector[jjFF][[4]]
		+ q^2/mN^2(
			cvector[iiFF][[4]]cvector[jjFF][[6]]
			+ cvector[iiFF][[6]]cvector[jjFF][[4]]
		) 
		+ q^4/mN^4 cvector[iiFF][[6]]cvector[jjFF][[6]] 
		+ v^2 cvector[iiFF][[12]]cvector[jjFF][[12]] 
		+ q^2/mN^2 v^2cvector[iiFF][[13]]cvector[jjFF][[13]]
	),
	{iiFF,2},{jjFF,2}
];

DMResponseCoeff[SigmaPJ] = Table[
	1/8(
		q^2/mN^2 v^2 cvector[iiFF][[3]]cvector[jjFF][[3]] 
		+ v^2 cvector[iiFF][[7]]cvector[jjFF][[7]]
	) 
	+ 1/16 Cl[jchi](
		cvector[iiFF][[4]]cvector[jjFF][[4]] 
		+ q^2/mN^2 cvector[iiFF][[9]]cvector[jjFF][[9]] 
		+ v^2/2 cvector[iiFF][[12]]cvector[jjFF][[12]]
		+ q^2/(2mN^2)cvector[iiFF][[14]]cvector[jjFF][[14]]
	),
	{iiFF,2},{jjFF,2}
];

DMResponseCoeff[DeltaJ] = Table[
	2v^2 cvector[iiFF][[2]]cvector[jjFF][[2]] 
	+ Cl[jchi]/4 (
		q^2/mN^2 cvector[iiFF][[5]]cvector[jjFF][[5]] 
		+ cvector[iiFF][[8]]cvector[jjFF][[8]]
	),
	{iiFF,2},{jjFF,2}
];

DMResponseCoeff[SigmaPJ, DeltaJ] = Table[
	-v^2 cvector[iiFF][[2]]cvector[jjFF][[3]] 
	+ Cl[jchi]/4 (
		cvector[iiFF][[5]]cvector[jjFF][[4]]
		-cvector[iiFF][[8]]cvector[jjFF][[9]]
	),
	{iiFF,2},{jjFF,2}
];

Sum[	
	DMResponseCoeff[MJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,J,T][[iiFF,jjFF]]
	+ DMResponseCoeff[SigmaPPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPPJ,J,T][[iiFF,jjFF]]
	+ DMResponseCoeff[SigmaPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,J,T][[iiFF,jjFF]]
	+ DMResponseCoeff[DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,DeltaJ,J,T][[iiFF,jjFF]]
	+ DMResponseCoeff[PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiPPJ,J,T][[iiFF,jjFF]]
	+ DMResponseCoeff[PhiTPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiTPJ,J,T][[iiFF,jjFF]]
	+ DMResponseCoeff[MJ,PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,PhiPPJ,J,T][[iiFF,jjFF]]
	+ DMResponseCoeff[SigmaPJ,DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,DeltaJ,J,T][[iiFF,jjFF]],
	{iiFF,2},{jjFF,2}
]

]



(* ::Text:: *)
(*Velocity Distributions*)


Fv0MB[\[Eta]_,xmin_]:=(Erf[xmin+\[Eta]]-Erf[xmin-\[Eta]])/(2 \[Eta]);
FvsqMB[\[Eta]_,xmin_]:=((E^-(-xmin+\[Eta])^2 (xmin+\[Eta])-E^-(xmin+\[Eta])^2 (xmin-\[Eta]))/(2 Sqrt[\[Pi]] \[Eta])+((1+2 \[Eta]^2) (Erf[xmin+\[Eta]]-Erf[xmin-\[Eta]]))/(4 \[Eta]));
Fv0MBcutoff[\[Eta]_,xmin_,xesc_]:=(-4 \[Eta] (-3-3 xesc^2+3 xmin^2+\[Eta]^2)+3 E^xesc^2 Sqrt[\[Pi]] (Erf[xmin-\[Eta]]-Erf[xmin+\[Eta]]))/(2 \[Eta] (6 xesc+4 xesc^3-3 E^xesc^2 Sqrt[\[Pi]] Erf[xesc])) UnitStep[xesc-xmin-\[Eta]]+(2 (xesc-xmin+\[Eta]) (3+(2 xesc+xmin-\[Eta]) (xesc-xmin+\[Eta]))+3 E^xesc^2 Sqrt[\[Pi]] (-Erf[xesc]+Erf[xmin-\[Eta]]))/(2 \[Eta] (6 xesc+4 xesc^3-3 E^xesc^2 Sqrt[\[Pi]] Erf[xesc])) UnitStep[xesc-xmin+\[Eta],-xesc+xmin+\[Eta]];
FvsqMBcutoff[\[Eta]_,xmin_,xesc_]:=(-(E^(-2 (xmin^2+\[Eta]^2)) (-30 E^(xesc^2+(xmin-\[Eta])^2) (xmin-\[Eta])+30 E^(xesc^2+(xmin+\[Eta])^2) (xmin+\[Eta])+4 E^(2 (xmin^2+\[Eta]^2)) \[Eta] (-15 (2+2 xesc^2+xesc^4-xmin^4)-10 (1+xesc^2) \[Eta]^2+\[Eta]^4)+15 E^(xesc^2+2 (xmin^2+\[Eta]^2)) Sqrt[\[Pi]] (1+2 \[Eta]^2) (-Erf[xmin-\[Eta]]+Erf[xmin+\[Eta]])))/(20 \[Eta] (6 xesc+4 xesc^3-3 E^xesc^2 Sqrt[\[Pi]] Erf[xesc])))UnitStep[xesc-xmin-\[Eta]]+((E^(-xmin^2-\[Eta]^2) (2 E^(xmin^2+\[Eta]^2) (15 xesc+10 xesc^3+4 xesc^5-10 xmin^3-10 xesc^2 xmin^3+6 xmin^5+15 (2+2 xesc^2+xesc^4-xmin^4) \[Eta]+10 (3 xesc+2 xesc^3+xmin^3) \[Eta]^2+10 (1+xesc^2) \[Eta]^3-\[Eta]^5)+15 E^xesc^2 (-2 E^(2 xmin \[Eta]) (xmin+\[Eta])+E^(xmin^2+\[Eta]^2) Sqrt[\[Pi]] (1+2 \[Eta]^2) (-Erf[xesc]+Erf[xmin-\[Eta]]))))/(20 \[Eta] (6 xesc+4 xesc^3-3 E^xesc^2 Sqrt[\[Pi]] Erf[xesc])))UnitStep[xesc-xmin+\[Eta],-xesc+xmin+\[Eta]];


(* ::Text:: *)
(*Here is the function for the user*)


(*the optional argument for formfactor if the user wants to output the form
factor to a TeX source, myFormFactor.tex. This file, by default, will appear
in the Users/username directory*)
StrucFunction[vv_,qqdimless_,ifOutFile_:0]:=Block[{ii, myFF},


PrintLag[];


Print["Your structure function is"];
(* The factor of 1/(4mN mchi)^2 is necessary to convert the c coefficients to a coefficients *)
myFF=E^(-2y) FFfinal[DensityMatrix,JIso,TIso]/.FormalReplace;
If[ifOutFile!=1,Return[1/(4mN mchi)^2*myFF/.y->(q bHO/2)^2/.q->qqdimless*GeV/.v->vv/.b->bHO//Expand//MyChop//Simplify]];

If[ifOutFile==1,
Print["Writing to TeX file..."];
TeXForm[1/(4mN mchi)^2*myFF/.y->(q bHO/2)^2/.q->qqdimless*GeV/.v->vv/.b->bHO//Expand//MyChop//Simplify]>>"myFormFactor.tex";
];];

TransitionProbability[vv_,qqdimless_,IfRel_:False]:=Block[{ii, myFF,ANonrelToRel},

ANonrelToRel=1;
If[IfRel,ANonrelToRel=(4mchi*M*mN)^2];

PrintLag[];


(* The factor of 1/(4mN mchi)^2 is necessary to convert the c coefficients to a coefficients *)
Print["Your transition probability is"];
myFF=E^(-2y) FFfinal[DensityMatrix,JIso,TIso]/.FormalReplace;
Return[ANonrelToRel/(4mN mchi)^2*myFF/.y->(q bHO/2)^2/.q->qqdimless*GeV/.v->vv/.b->bHO//Expand//MyChop//Simplify];
];


DiffCrossSection[ERkeV_,vv_]:=Block[{FFTemp,ER,bb,qq},
bb=bHO;
ER=ERkeV*10^(-6)*GeV;
qq=Sqrt[2M*(mN/GeV)*ERkeV]*10^(-3)GeV;

PrintLag[];


FFTemp=E^(-2y) FFfinal[DensityMatrix,JIso,TIso]/.q->qq/.y->((qq bb)/2)^2/.b->bb/.v->vv;


Print["Your cross-section is"];
Return[M/(32\[Pi] vv^2 mchi^2 mN) FFTemp(*/.vmin->qq/(2\[Mu]T[mchi,M])*)]/.FormalReplace//MyChop
];


ApproxTotalCrossSection[vv_]:=Block[{FFTempNoExp,FFTempNoExp2,ER,qqmax,bb,qq},
bb=bHO;
qqmax=2\[Mu]T[mchi,M]*vv;

PrintLag[];


FFTempNoExp=FFfinal[DensityMatrix,JIso,TIso]/.q->qq/.y->((qq bb)/2)^2/.b->bb/.v->vv;
(*Change integration variables using dER=(q/mT)dq and integrate*)
FFTempNoExp2=Collect[(qq/(M*mN))*FFTempNoExp,qq]/.{qq^{a_Integer} -> Superscript[qqmax, a+1]/(a+1),qq->qqmax/2} (*/.qq^(a_)->(qqmax^(a+1))/(a+1)*);

Print["Here is your small-nucleus-radius approximate cross-section.
This approximation is valid only at v \[LessLess] 1/(2\!\(\*SubscriptBox[\(\[Mu]\), \(T\)]\)b)."];
Return[M/(32\[Pi] vv^2 mchi^2 mN) FFTempNoExp2]/.FormalReplace//MyChop//Expand
];


EventRate[NT_,\[Rho]chi_,qqGeV_,vve_,vv0_]:=EventRate[NT,\[Rho]chi,qqGeV,vve,vv0,12*vv0];
EventRate[NT_,\[Rho]chi_,qqGeV_,vve_,vv0_,vvesc_]:=EventRateInelastic[NT,\[Rho]chi,0,qqGeV,vve,vv0,vvesc];


EventRateInelastic[NT_,\[Rho]chi_,\[Delta]\[Delta]_,qqGeV_,vve_,vv0_]:=EventRateInelastic[NT,\[Rho]chi,\[Delta]\[Delta],qqGeV,vve,vv0,12vv0];
EventRateInelastic[NT_,\[Rho]chi_,\[Delta]\[Delta]_,qqGeV_,vve_,vv0_,vvesc_]:=Block[{FFTemp,ERTemp,bb,qq},

bb=bHO;
qq=qqGeV*GeV;
PrintLag[];



FFTemp=E^(-2y) FFfinal[DensityMatrix,JIso,TIso]/.q->qq/.y->((qq bb)/2)^2/.b->bb;

(*FFTemp=FormFactor[((qq bb)/2)^2,bb,v];*)
If[(HALO!="MB")&&(HALO!="MBcutoff"),Print["Warning: Halo option not recognized.  Setting to Maxwell-Boltzmann."]; HALO="MB";];
If[HALO=="MB",ERTemp=Coefficient[FFTemp,v,0] Fv0MB[vve/vv0,vmin/vv0]/vv0+Coefficient[FFTemp,v,2]FvsqMB[vve/vv0,vmin/vv0]*vv0];
If[HALO=="MBcutoff",ERTemp=Coefficient[FFTemp,v,0] Fv0MBcutoff[vve/vv0,vmin/vv0,vvesc/vv0]/vv0+Coefficient[FFTemp,v,2]FvsqMBcutoff[vve/vv0,vmin/vv0,vvesc/vv0]*vv0];
Print["Your event rate is"];
Return[NT \[Rho]chi M/(32\[Pi] mchi^3 mN) ERTemp/.vmin->qq/(2\[Mu]T[mchi,M])+\[Delta]\[Delta]/qq]/.FormalReplace
];




ReadInDMFile[DensityMatrixFile_]:=Block[{},
inputstring=Import[DensityMatrixFile,"TSV"];
DensityMatrixVector={};
For[i=1,i<=Length[inputstring],i++,
JT=StringCases[inputstring[[i,1]],RegularExpression[".+ONE.+JO.+\\s(\\d+)\\s.+(\\d+).*"]->{"$1","$2"}]//ToExpression;
dens=StringCases[inputstring[[i,1]],RegularExpression[".+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(.+\\d+.\\d+).+"]->{"$1","$2","$3","$4","$5"}]//ToExpression;
If[{}!=JT,oJT=JT[[1]]/2;];
If[{}!=dens,oDens={dens[[1,5]],oJT[[1]],oJT[[2]],dens[[1,1]],dens[[1,2]],dens[[1,3]],dens[[1,4]]};DensityMatrixVector=DensityMatrixVector~Join~{oDens}];
];

For[i=1,i<=Length[DensityMatrixVector],i++,
DensityMatrix[i]=DensityMatrixVector[[i]];
];

SetDMmatrixLength[Length[DensityMatrixVector]];
(*DensityMatrix=sdF;*)
Return[];];


End[]

EndPackage[]

Print[Style["Welcome to DMFormFactor version 1.1.",32,Blue]];
Print["Functions are SetCoeffsNonrel, SetCoeffsRel, SetCoeffsNucl, ZeroCoeffs, SetJChi, SetMchi, SetIsotope, SetHALO, SetHelm, TransitionProbability, ResponseNucl, DiffCrossSection, ApproxTotalCrossSection, and EventRate."];




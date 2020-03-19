(* ::Package:: *)

GeV;
Femtometer=5.0677/GeV;
KilometerPerSecond=3.3356*10^-6;
KilogramDay=7.3634*10^55

(*SetCoeffsNonrel::usage="SetCoeffsNonrel[Op,coeff,nucleon] with Op=1,...,12, 
  coeff=EFT coefficient, and nucleon=n,p,0, or 1."*)

SetCoeffsNucl::usage="SetCoeffsNucl[var,ell,nucleon] with var=1,...,9, 
  ell=nuclear coefficient, and nucleon=n,p,0, or 1."

SetCoeffsRel::usage="SetCoeffsRel[Op,coeff,nucleon] with Op=1,...,20,
  coeff=EFT relativistic coefficient, and nucleon=n,p,0, or 1."

ZeroCoeffs::usage="ZeroCoeffs[] resets all operator coefficients to zero."

Begin["`Private`"]

SetDMmatrixLength[x_]:=Block[{},
DMmatrixLength=x;
];
(*The user can set the dark matter spin via SetJChi[..]*)
jchi=0;
SetJChi[jx_]:=Block[{},
jchi=jx;
];
(*And likewise for the dark matter mass*)
mchi=m\[Chi];
SetMChi[mx_]:=Block[{},
mchi=mx;
];
(*And likewise for the fiducial interaction scale*)
mM=mV;
SetMM[mm_]:=Block[{},
mM=mm;
];
(*And also the numbers needed to compute the target reduced mass*)
M=A;
SetM[Num_]:=Block[{},
M=Num;
];
ZZ=Z;
SetZ[Num_]:=Block[{},
ZZ=Num;
];
bHO=HarmonicOscillatorParameter;
SetbHO[Num_]:=Block[{},
bHO=Num;
];
(*The variables mchiFORMAL and mMFORMAL are internal to the program and are 
protected so that they can never be set to any value.  This way, if the user 
changes them externally, all formulae will be updated appropriately*)
SetAttributes[mchiFORMAL,Protected];
SetAttributes[mMFORMAL,Protected];
FormalReplace:={mchiFORMAL->mchi,mMFORMAL->mM};
(*The Helm form factor for the density depends on the following parameters *)
s=Femtometer;
rHelm:=1.2 M^(1/3) Femtometer;
r0Helm:=Sqrt[(rHelm^2-5s^2)GeV^2]/GeV;
UseHelm=False;
SetHelm[bool_]:=Block[{},
UseHelm=bool;
];
MyChop[t_]:=Chop[t,10.^-99.];

HALO="MB";
SetHALO[halomodel_]:=Block[{},
HALO=halomodel;
];


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



CoeffDimNonrel[1]=0;CoeffDimNonrel[2]=0;CoeffDimNonrel[3]=-1;CoeffDimNonrel[4]=0;CoeffDimNonrel[5]=-1;CoeffDimNonrel[6]=-2;CoeffDimNonrel[7]=0;CoeffDimNonrel[8]=0;CoeffDimNonrel[9]=-1;CoeffDimNonrel[10]=-1;CoeffDimNonrel[11]=-1;CoeffDimNonrel[12]=0;CoeffDimNonrel[13]=-1;CoeffDimNonrel[14]=-1;CoeffDimNonrel[15]=-2;
CoeffDimRel[1]=-2;CoeffDimRel[2]=-2;CoeffDimRel[3]=-2;CoeffDimRel[4]=-4;CoeffDimRel[5]=-2;CoeffDimRel[6]=-4;CoeffDimRel[7]=-3;CoeffDimRel[8]=-4;CoeffDimRel[9]=-4;CoeffDimRel[10]=-4;CoeffDimRel[11]=-3;CoeffDimRel[12]=-4;CoeffDimRel[13]=-3;CoeffDimRel[14]=-3;CoeffDimRel[15]=-2;CoeffDimRel[16]=-3;CoeffDimRel[17]=-4;CoeffDimRel[18]=-4;CoeffDimRel[19]=-3;CoeffDimRel[20]=-4;


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

PrintLag[]:=Block[{},
OldToNewNorm=1/(4mN*mchi);


Print["Your Lagrangian is"];
If[!UsingRelCoeffs,
Print[" \!\(\*SubscriptBox[\(L\), \(prot\)]\)=",Sum[OldToNewNorm*cpvector[[ii]]Op[ii]/.FormalReplace,{ii,1,15}]];
Print[" \!\(\*SubscriptBox[\(L\), \(neut\)]\)=",Sum[OldToNewNorm*cnvector[[ii]]Op[ii]/.FormalReplace,{ii,1,15}]];
If[!(cpvector[[2]]===0&&cnvector[[2]]===0),
   Print[Style["Warning: Implementation of \!\(\*SubscriptBox[\(O\), \(2\)]\) discards \!\(\*SuperscriptBox[SubscriptBox[\(v\), \(N\)], \(2\)]\) piece.  ",32,Orange]];
];
,
Print[" \!\(\*SubscriptBox[\(L\), \(prot\)]\)=",Sum[dpvector[[ii]]OpRel[ii]/.FormalReplace,{ii,1,20}]];
Print[" \!\(\*SubscriptBox[\(L\), \(neut\)]\)=",Sum[dnvector[[ii]]OpRel[ii]/.FormalReplace,{ii,1,20}]];;
];

]

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



FFfinal[DMmatrix_,J_,T_]:=Block[{iiFF,jjFF,ResponseCoeff,iiResponse},

cvector[1]=cpvector;
cvector[2]=cnvector;
(*These response coefficients differ in their definition from those in the paper 
because they have old c's instead of new c's, resulting in an overall factor of 
1/(4mN mchi)^2.  This will show up as an additional factor that we multiply
FF by before printing out StrucFunction, TransitionProbability, and 
ResponseNucl *)
ResponseCoeff[MJ]=Table[1/4 Cl[jchi]((cvector[iiFF][[5]]cvector[jjFF][[5]] q^2+cvector[iiFF][[8]]cvector[jjFF][[8]])(v^2-q^2/(4(\[Mu]T[mchi,M])^2))+cvector[iiFF][[11]]cvector[jjFF][[11]]q^2)+(cvector[iiFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2))+cvector[iiFF][[1]])(cvector[jjFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2))+cvector[jjFF][[1]]),{iiFF,2},{jjFF,2}];
ResponseCoeff[SigmaPPJ]=Table[1/16 Cl[jchi](cvector[iiFF][[6]]cvector[jjFF][[6]]q^4+(cvector[iiFF][[13]]cvector[jjFF][[13]]q^2+cvector[iiFF][[12]]cvector[jjFF][[12]])(v^2-q^2/(4(\[Mu]T[mchi,M])^2))+2cvector[iiFF][[4]]cvector[jjFF][[6]]q^2+cvector[iiFF][[4]]cvector[jjFF][[4]])+1/4 cvector[iiFF][[10]]cvector[jjFF][[10]]q^2,{iiFF,2},{jjFF,2}];
ResponseCoeff[SigmaPJ]=Table[1/32 Cl[jchi](2cvector[iiFF][[9]]cvector[jjFF][[9]]q^2+(cvector[iiFF][[15]]cvector[jjFF][[15]]q^4+cvector[iiFF][[14]]cvector[jjFF][[14]]q^2-2cvector[iiFF][[12]]cvector[jjFF][[15]]q^2+cvector[iiFF][[12]]cvector[jjFF][[12]])(v^2-q^2/(4(\[Mu]T[mchi,M])^2))+2cvector[iiFF][[4]]cvector[jjFF][[4]])+1/8 (cvector[iiFF][[3]]cvector[jjFF][[3]]q^2+cvector[iiFF][[7]]cvector[jjFF][[7]])(v^2-q^2/(4(\[Mu]T[mchi,M])^2)),{iiFF,2},{jjFF,2}];
ResponseCoeff[DeltaJ]=Table[q^2/(4mN^2) Cl[jchi](cvector[iiFF][[5]]cvector[jjFF][[5]]q^2+cvector[iiFF][[8]]cvector[jjFF][[8]])+2 q^2/mN^2 cvector[iiFF][[2]]cvector[jjFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2)),{iiFF,2},{jjFF,2}];
ResponseCoeff[PhiPPJ]=Table[q^2/(16mN^2) Cl[jchi](cvector[iiFF][[12]]-cvector[iiFF][[15]]q^2)(cvector[jjFF][[12]]-cvector[jjFF][[15]]q^2)+q^4/(4mN^2)cvector[iiFF][[3]]cvector[jjFF][[3]],{iiFF,2},{jjFF,2}];
ResponseCoeff[PhiTPJ]=Table[q^2/(16mN^2)Cl[jchi](cvector[iiFF][[13]]cvector[jjFF][[13]]q^2+cvector[iiFF][[12]]cvector[jjFF][[12]]),{iiFF,2},{jjFF,2}];
ResponseCoeff[MJ,PhiPPJ]=Table[q^2/(4mN)Cl[jchi]cvector[iiFF][[11]](cvector[jjFF][[12]]-cvector[jjFF][[15]]q^2)+q^2/(mN)cvector[jjFF][[3]](cvector[iiFF][[1]]+cvector[iiFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2))),{iiFF,2},{jjFF,2}];
ResponseCoeff[SigmaPJ,DeltaJ]=Table[q^2/(4mN)Cl[jchi](cvector[iiFF][[4]]cvector[jjFF][[5]]-cvector[jjFF][[8]]cvector[iiFF][[9]])-q^2/(mN)cvector[jjFF][[2]]cvector[iiFF][[3]](v^2-q^2/(4(\[Mu]T[mchi,M])^2)),{iiFF,2},{jjFF,2}];
Sum[ResponseCoeff[MJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,J,T][[iiFF,jjFF]]+ResponseCoeff[SigmaPPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[SigmaPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,DeltaJ,J,T][[iiFF,jjFF]]+ResponseCoeff[PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiPPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[PhiTPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiTPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[MJ,PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,PhiPPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[SigmaPJ,DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,DeltaJ,J,T][[iiFF,jjFF]],{iiFF,2},{jjFF,2}]
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


(*SetIsotope sets the density matrix to the one in filename. If filename
is specified as "default", then it will load
one of the default density matrices built-in to the program. Which one it
loads depends on M and Z - the available isotopes are presented in
Fitizpatrick et al. arXiv:1203.3542v2*)

SetIsotope[Ziso_,Mxx_, bHOdimlessxx_,filename_]:=Block[{},

JIso=IsoJLookupTable[Ziso,Mxx];
TIso=Mxx/2-Ziso//Abs;
SetM[Mxx];
SetZ[Ziso];
If[bHOdimlessxx==="default"||bHOdimlessxx==="Default"||bHOdimlessxx==="DEFAULT",
SetbHO[(41.467/(45Mxx^(-1/3)-25Mxx^(-2/3)))^(1/2)Femtometer];,
SetbHO[bHOdimlessxx Femtometer];
]


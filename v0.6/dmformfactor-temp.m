(* ::Package:: *)

(*****************************
* dmformfactor.m
* Package written by Nikhil Anand, A. Liam Fitzpatrick, and Wick Haxton.
* 
* For documentation, see http://arxiv.org/abs/arXiv:13xx.xxxx
* 
* *)

BeginPackage["dmformfactor`"]


(*Provide GeV as a unit for the user to check dimensions.*)

GeV;
Femtometer=5.0677/GeV;
KilometerPerSecond=3.3356*10^-6;
KilogramDay=7.3634*10^55
(*********)



SetCoeffsNonrel::usage="SetCoeffsNonrel[Op,coeff,nucleon] with Op=1,...,12, 
  coeff=EFT coefficient, and nucleon=n,p,0, or 1."

SetCoeffsNucl::usage="SetCoeffsNucl[var,ell,nucleon] with var=1,...,9, 
  ell=nuclear coefficient, and nucleon=n,p,0, or 1."

SetCoeffsRel::usage="SetCoeffsRel[Op,coeff,nucleon] with Op=1,...,20,
  coeff=EFT relativistic coefficient, and nucleon=n,p,0, or 1."

ZeroCoeffs::usage="ZeroCoeffs[] resets all operator coefficients to zero."

CheckCoeffDims::usage="CheckCoeffDims[]."

SetDMmatrixLength::usage="SetDMMatrixLength[DMMatrixLength]."

SetJChi::usage="SetJChi[Jchi], where Jchi is the dark matter spin."

SetMChi::usage="SetMChi[Mchi], where Mchi is the dark matter particle mass."

SetMM::usage="SetMM[MM], where MM is a fiducial scale for the relativistic interactions."

SetM::usage="SetM[Num], where Num is the isotope mass number."

SetHALO::usage="SetHALO[HALO], where HALO=\"MB\" or \"MBcutoff\". "

SetHelm::usage="SetHelm[UseHelm], where UseHelm=True or False."

(*To come? Setb, which sets the harmonic oscillator parameter*)

StrucFunction::usage="StrucFunction[y,b,v] outputs the final structure function."

TransitionProbability::usage="TransitionProbabilty[v,q] outputs the final transition probability."

ResponseNuclear::usage="ResponseNuclear[y,i,tau1,tau2] outputs the nuclear response functions, specified 0<=i<=8 (which correspond to the eight operators and interferences in Eqn. 39) and the isospins tau1 and tau2."

FormFactor2::usage="FormFactor2[y,b,v] outputs the form factors."

DiffCrossSection::usage="DiffCrossSection[ER,v] outputs the differential cross-section."

ApproxTotalCrossSection::usage="ApproxTotalCrossSection[v] outputs the total cross-section
integrate over all recoil energies in an approximation that the nucleus harmonic 
oscillator parameter b is small, so that the Exp[-2y] term in the response functions 
can be neglected."

EventRate::usage="EventRate[NT,rhoDM,q,ve,v0] outputs the event rate."

(*EventRateInelastic::usage="EventRateInelastic[NT,rhoDM,delta,q,b,ve,v0] outputs the event rate for ...."*)

(*ReadInDMFile::usage="ReadInDMFile[DensityMatrix] reads in the density matrix from file DensityMatrix."*)

SetIsotope::usage="SetIsotope[Z,M,b, filename] reads the density matrix for
the isotope stored in filename. Z is the charge, M is the mass number, and b is
the harmonic oscillator parameter. If b is set to default (in quotes), it will
be set to \!\(\*SqrtBox[\(41.467/\((45 \*SuperscriptBox[\(M\), \(\(-1\)/3\)] - 25 \*SuperscriptBox[\(M\), \(\(-2\)/3\)])\)\)]\)femtometers. If filename is set
to default (in quotes), SetIsotope will load the default density matrix."




Begin["`Private`"]


(**************
* Section for Setting the EFT coefficients and model parameters
*)

(*By default, make all the coefficients vanish*)
cnvector=Table[0,{i,15}];
cpvector=Table[0,{i,15}];
(*and likewise for the relativisitc operator coefficients*)
dnvector=Table[0,{i,20}];
dpvector=Table[0,{i,20}];
(*and likewise for the nuclear Hamiltonian coefficients*)
ellpnucl={{0}}~Join~Table[{0,0,0},{i,2,4}];
ellnnucl={{0}}~Join~Table[{0,0,0},{i,2,4}];

(*Function to keep track of the length of the density matrix. *)

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



(* ::Text:: *)
(*Halo options*)


(*By default, use Maxwell-Boltzmann*)
HALO="MB";
SetHALO[halomodel_]:=Block[{},
HALO=halomodel;
];


(* ::Text:: *)
(*Default Density Matrices*)


(*Library of default density matrices*)

(* rho, Jt, Tt , nbra, 2jbra, nket, 2jket *)
sHe4[1]={sHe00x0101=2.0000000,0,0,0,1,0,1};
sHe4length=1;

sO16[1]={sO1600x0101=2.0000000,0,0,0,1,0,1};
sO16[2]={sO1600x1111=2.0000000,0,0,1,1,1,1};
sO16[3]={sO1600x1313=2.8284300,0,0,1,3,1,3};
sO16length=3;

sdF[1]={sdF00x0101=4.0000000,0,0,0,1,0,1};
sdF[2]={sdF00x1111=4.0000000,0,0,1,1,1,1};
sdF[3]={sdF00x1313=5.6568542,0,0,1,3,1,3};
sdF[4]={sdF00x2121=1.2252593,0,0,2,1,2,1};
sdF[5]={sdF00x2323=0.2036611,0,0,2,3,2,3};
sdF[6]={sdF00x2525=0.8583583,0,0,2,5,2,5};
sdF[7]={sdF01x2121=0.3698483,0,1,2,1,2,1};
sdF[8]={sdF01x2323=0.0479437,0,1,2,3,2,3};
sdF[9]={sdF01x2525=0.3246722,0,1,2,5,2,5};
sdF[10]={sdF10x2121=0.4451426,1,0,2,1,2,1};
sdF[11]={sdF10x2321=-0.0119775,1,0,2,3,2,1};
sdF[12]={sdF10x2123=0.0119775,1,0,2,1,2,3};
sdF[13]={sdF10x2323=-0.0542883,1,0,2,3,2,3};
sdF[14]={sdF10x2523=-0.1217257,1,0,2,5,2,3};
sdF[15]={sdF10x2325=0.1217257,1,0,2,3,2,5};
sdF[16]={sdF10x2525=0.1228063,1,0,2,5,2,5};
sdF[17]={sdF11x2121=-0.4078034,1,1,2,1,2,1};
sdF[18]={sdF11x2321=-0.0127852,1,1,2,3,2,1};
sdF[19]={sdF11x2123=0.0127852,1,1,2,1,2,3};
sdF[20]={sdF11x2323=0.0120967,1,1,2,3,2,3};
sdF[21]={sdF11x2523=0.1054748,1,1,2,5,2,3};
sdF[22]={sdF11x2325=-0.1054748,1,1,2,3,2,5};
sdF[23]={sdF11x2525=-0.2411054,1,1,2,5,2,5};
sdFlength=23;

sdSi28[1]={sdSi2800x0101=2.,0,0,0,1,0,1};
sdSi28[2]={sdSi2800x1111=2.,0,0,1,1,1,1};
sdSi28[3]={sdSi2800x1313=2.82843,0,0,1,3,1,3};
sdSi28[4]={sdSi2800x2121=0.7038300,0,0,2,1,2,1};
sdSi28[5]={sdSi2800x2323=0.4757098,0,0,2,3,2,3};
sdSi28[6]={sdSi2800x2525=2.6691689,0,0,2,5,2,5};
sdSi28length=6;

sdSi29[1]={sdSi2800x0101=4.,0,0,0,1,0,1};
sdSi29[2]={sdSi2800x1111=4.,0,0,1,1,1,1};
sdSi29[3]={sdSi2800x1313=5.65685,0,0,1,3,1,3};
sdSi29[4]={sdSi2900x2121=1.6791266,0,0,2,1,2,1};
sdSi29[5]={sdSi2900x2323=0.8589358,0,0,2,3,2,3};
sdSi29[6]={sdSi2900x2525=5.8344745,0,0,2,5,2,5};
sdSi29[7]={sdSi2901x2121=0.5065078,0,1,2,1,2,1};
sdSi29[8]={sdSi2901x2323=0.1242974,0,1,2,3,2,3};
sdSi29[9]={sdSi2901x2525=0.1833635,0,1,2,5,2,5};
sdSi29[10]={sdSi2910x2121=0.5720814,1,0,2,1,2,1};
sdSi29[11]={sdSi2910x2321=0.0249074,1,0,2,3,2,1};
sdSi29[12]={sdSi2910x2123=-0.0249074,1,0,2,1,2,3};
sdSi29[13]={sdSi2910x2323=-0.0803281,1,0,2,3,2,3};
sdSi29[14]={sdSi2910x2523=0.1966218,1,0,2,5,2,3};
sdSi29[15]={sdSi2910x2325=-0.1966218,1,0,2,3,2,5};
sdSi29[16]={sdSi2910x2525=0.1152554,1,0,2,5,2,5};
sdSi29[17]={sdSi2911x2121=0.5528034,1,1,2,1,2,1};
sdSi29[18]={sdSi2911x2321=0.0142132,1,1,2,3,2,1};
sdSi29[19]={sdSi2911x2123=-0.0142132,1,1,2,1,2,3};
sdSi29[20]={sdSi2911x2323=-0.0875390,1,1,2,3,2,3};
sdSi29[21]={sdSi2911x2523=0.1856207,1,1,2,5,2,3};
sdSi29[22]={sdSi2911x2325=-0.1856207,1,1,2,3,2,5};
sdSi29[23]={sdSi2911x2525=0.1110496,1,1,2,5,2,5};
sdSi29length=23;

sdSi30[1]={sdSi2800x0101=3.4641,0,0,0,1,0,1};
sdSi30[2]={sdSi2800x1111=3.4641,0,0,1,1,1,1};
sdSi30[3]={sdSi2800x1313=4.89898,0,0,1,3,1,3};
sdSi30[4]={sdSi3000x2121=1.6386966,0,0,2,1,2,1};
sdSi30[5]={sdSi3000x2323=1.0012031,0,0,2,3,2,3};
sdSi30[6]={sdSi3000x2525=5.2361869,0,0,2,5,2,5};
sdSi30[7]={sdSi3001x2121=0.5920034,0,1,2,1,2,1};
sdSi30[8]={sdSi3001x2323=0.3893739,0,1,2,3,2,3};
sdSi30[9]={sdSi3001x2525=0.1567129,0,1,2,5,2,5};
sdSi30length=9;

sdgI[1]={sdgI00x0101=22.9782505,0,0,0,1,0,1};
sdgI[2]={sdgI00x1111=22.9782505,0,0,1,1,1,1};
sdgI[3]={sdgI00x1313=32.4961536,0,0,1,3,1,3};
sdgI[4]={sdgI00x2121=22.9782505,0,0,2,1,2,1};
sdgI[5]={sdgI00x2323=32.4961536,0,0,2,3,2,3};
sdgI[6]={sdgI00x2525=39.7994974,0,0,2,5,2,5};
sdgI[7]={sdgI00x3131=22.9782505,0,0,3,1,3,1};
sdgI[8]={sdgI00x3333=32.4961536,0,0,3,3,3,3};
sdgI[9]={sdgI00x3535=39.7994974,0,0,3,5,3,5};
sdgI[10]={sdgI00x3737=45.9565011,0,0,3,7,3,7};
sdgI[11]={sdgI00x4141=12.2143316,0,0,4,1,4,1};
sdgI[12]={sdgI00x4343=17.2413539,0,0,4,3,4,3};
sdgI[13]={sdgI00x4545=27.6852741,0,0,4,5,4,5};
sdgI[14]={sdgI00x4747=23.7869625,0,0,4,7,4,7};
sdgI[15]={sdgI00x4949=51.3809303,0,0,4,9,4,9};
sdgI[16]={sdgI00x511511=9.3807897,0,0,5,11,5,11};
sdgI[17]={sdgI01x4141=6.5035815,0,1,4,1,4,1};
sdgI[18]={sdgI01x4343=9.2169838,0,1,4,3,4,3};
sdgI[19]={sdgI01x4545=7.3193388,0,1,4,5,4,5};
sdgI[20]={sdgI01x4747=13.3948955,0,1,4,7,4,7};
sdgI[21]={sdgI01x511511=5.6680417,0,1,5,11,5,11};
sdgI[22]={sdgI10x4141=0.0237469,1,0,4,1,4,1};
sdgI[23]={sdgI10x4341=-0.0894454,1,0,4,3,4,1};
sdgI[24]={sdgI10x4143=0.0894454,1,0,4,1,4,3};
sdgI[25]={sdgI10x4343=0.1124782,1,0,4,3,4,3};
sdgI[26]={sdgI10x4543=0.1683758,1,0,4,5,4,3};
sdgI[27]={sdgI10x4345=-0.1683758,1,0,4,3,4,5};
sdgI[28]={sdgI10x4545=2.2032630,1,0,4,5,4,5};
sdgI[29]={sdgI10x4745=0.1716681,1,0,4,7,4,5};
sdgI[30]={sdgI10x4547=-0.1716681,1,0,4,5,4,7};
sdgI[31]={sdgI10x4747=0.0601998,1,0,4,7,4,7};
sdgI[32]={sdgI10x511511=0.3344252,1,0,5,11,5,11};
sdgI[33]={sdgI11x4141=-0.0143481,1,1,4,1,4,1};
sdgI[34]={sdgI11x4341=0.0540445,1,1,4,3,4,1};
sdgI[35]={sdgI11x4143=-0.0540445,1,1,4,1,4,3};
sdgI[36]={sdgI11x4343=-0.0679603,1,1,4,3,4,3};
sdgI[37]={sdgI11x4543=-0.1017357,1,1,4,5,4,3};
sdgI[38]={sdgI11x4345=0.1017357,1,1,4,3,4,5};
sdgI[39]={sdgI11x4545=-1.3312444,1,1,4,5,4,5};
sdgI[40]={sdgI11x4745=-0.1037250,1,1,4,7,4,5};
sdgI[41]={sdgI11x4547=0.1037250,1,1,4,5,4,7};
sdgI[42]={sdgI11x4747=-0.0363690,1,1,4,7,4,7};
sdgI[43]={sdgI11x511511=0.2020650,1,1,5,11,5,11};
sdgI[44]={sdgI20x4341=-0.0367672,2,0,4,3,4,1};
sdgI[45]={sdgI20x4541=0.0331361,2,0,4,5,4,1};
sdgI[46]={sdgI20x4143=0.0367672,2,0,4,1,4,3};
sdgI[47]={sdgI20x4343=0.0016208,2,0,4,3,4,3};
sdgI[48]={sdgI20x4543=0.3246385,2,0,4,5,4,3};
sdgI[49]={sdgI20x4743=0.0094186,2,0,4,7,4,3};
sdgI[50]={sdgI20x4145=0.0331361,2,0,4,1,4,5};
sdgI[51]={sdgI20x4345=-0.3246385,2,0,4,3,4,5};
sdgI[52]={sdgI20x4545=0.4757186,2,0,4,5,4,5};
sdgI[53]={sdgI20x4745=0.0257897,2,0,4,7,4,5};
sdgI[54]={sdgI20x4347=0.0094186,2,0,4,3,4,7};
sdgI[55]={sdgI20x4547=-0.0257897,2,0,4,5,4,7};
sdgI[56]={sdgI20x4747=0.0655241,2,0,4,7,4,7};
sdgI[57]={sdgI20x511511=0.0840819,2,0,5,11,5,11};
sdgI[58]={sdgI21x4341=0.0222154,2,1,4,3,4,1};
sdgI[59]={sdgI21x4541=-0.0200214,2,1,4,5,4,1};
sdgI[60]={sdgI21x4143=-0.0222154,2,1,4,1,4,3};
sdgI[61]={sdgI21x4343=-0.0009804,2,1,4,3,4,3};
sdgI[62]={sdgI21x4543=-0.1961524,2,1,4,5,4,3};
sdgI[63]={sdgI21x4743=-0.0056909,2,1,4,7,4,3};
sdgI[64]={sdgI21x4145=-0.0200214,2,1,4,1,4,5};
sdgI[65]={sdgI21x4345=0.1961524,2,1,4,3,4,5};
sdgI[66]={sdgI21x4545=-0.2874359,2,1,4,5,4,5};
sdgI[67]={sdgI21x4745=-0.0155826,2,1,4,7,4,5};
sdgI[68]={sdgI21x4347=-0.0056909,2,1,4,3,4,7};
sdgI[69]={sdgI21x4547=0.0155826,2,1,4,5,4,7};
sdgI[70]={sdgI21x4747=-0.0395901,2,1,4,7,4,7};
sdgI[71]={sdgI21x511511=0.0508037,2,1,5,11,5,11};
sdgI[72]={sdgI30x4541=-0.0619387,3,0,4,5,4,1};
sdgI[73]={sdgI30x4741=-0.0134526,3,0,4,7,4,1};
sdgI[74]={sdgI30x4343=-0.1709088,3,0,4,3,4,3};
sdgI[75]={sdgI30x4543=0.0448811,3,0,4,5,4,3};
sdgI[76]={sdgI30x4743=0.0183765,3,0,4,7,4,3};
sdgI[77]={sdgI30x4145=-0.0619387,3,0,4,1,4,5};
sdgI[78]={sdgI30x4345=-0.0448811,3,0,4,3,4,5};
sdgI[79]={sdgI30x4545=2.5543096,3,0,4,5,4,5};
sdgI[80]={sdgI30x4745=-0.1067677,3,0,4,7,4,5};
sdgI[81]={sdgI30x4147=0.0134526,3,0,4,1,4,7};
sdgI[82]={sdgI30x4347=0.0183765,3,0,4,3,4,7};
sdgI[83]={sdgI30x4547=0.1067677,3,0,4,5,4,7};
sdgI[84]={sdgI30x4747=-0.0290251,3,0,4,7,4,7};
sdgI[85]={sdgI30x511511=0.1036858,3,0,5,11,5,11};
sdgI[86]={sdgI31x4541=0.0374245,3,1,4,5,4,1};
sdgI[87]={sdgI31x4741=0.0081283,3,1,4,7,4,1};
sdgI[88]={sdgI31x4343=0.1032660,3,1,4,3,4,3};
sdgI[89]={sdgI31x4543=-0.0271179,3,1,4,5,4,3};
sdgI[90]={sdgI31x4743=-0.0111034,3,1,4,7,4,3};
sdgI[91]={sdgI31x4145=0.0374245,3,1,4,1,4,5};
sdgI[92]={sdgI31x4345=0.0271179,3,1,4,3,4,5};
sdgI[93]={sdgI31x4545=-1.5433611,3,1,4,5,4,5};
sdgI[94]={sdgI31x4745=0.0645109,3,1,4,7,4,5};
sdgI[95]={sdgI31x4147=-0.0081283,3,1,4,1,4,7};
sdgI[96]={sdgI31x4347=-0.0111034,3,1,4,3,4,7};
sdgI[97]={sdgI31x4547=-0.0645109,3,1,4,5,4,7};
sdgI[98]={sdgI31x4747=0.0175367,3,1,4,7,4,7};
sdgI[99]={sdgI31x511511=0.0626488,3,1,5,11,5,11};
sdgI[100]={sdgI40x4741=0.1385490,4,0,4,7,4,1};
sdgI[101]={sdgI40x4543=-0.4515199,4,0,4,5,4,3};
sdgI[102]={sdgI40x4743=-0.0719338,4,0,4,7,4,3};
sdgI[103]={sdgI40x4345=0.4515199,4,0,4,3,4,5};
sdgI[104]={sdgI40x4545=0.4072144,4,0,4,5,4,5};
sdgI[105]={sdgI40x4745=0.2055684,4,0,4,7,4,5};
sdgI[106]={sdgI40x4147=-0.1385490,4,0,4,1,4,7};
sdgI[107]={sdgI40x4347=-0.0719338,4,0,4,3,4,7};
sdgI[108]={sdgI40x4547=-0.2055684,4,0,4,5,4,7};
sdgI[109]={sdgI40x4747=-0.0768913,4,0,4,7,4,7};
sdgI[110]={sdgI40x511511=-0.7408876,4,0,5,11,5,11};
sdgI[111]={sdgI41x4741=-0.0837138,4,1,4,7,4,1};
sdgI[112]={sdgI41x4543=0.2728165,4,1,4,5,4,3};
sdgI[113]={sdgI41x4743=0.0434637,4,1,4,7,4,3};
sdgI[114]={sdgI41x4345=-0.2728165,4,1,4,3,4,5};
sdgI[115]={sdgI41x4545=-0.2460467,4,1,4,5,4,5};
sdgI[116]={sdgI41x4745=-0.1242081,4,1,4,7,4,5};
sdgI[117]={sdgI41x4147=0.0837138,4,1,4,1,4,7};
sdgI[118]={sdgI41x4347=0.0434637,4,1,4,3,4,7};
sdgI[119]={sdgI41x4547=0.1242081,4,1,4,5,4,7};
sdgI[120]={sdgI41x4747=0.0464610,4,1,4,7,4,7};
sdgI[121]={sdgI41x511511=-0.4476574,4,1,5,11,5,11};
sdgI[122]={sdgI50x4743=-0.0157470,5,0,4,7,4,3};
sdgI[123]={sdgI50x4545=2.1534283,5,0,4,5,4,5};
sdgI[124]={sdgI50x4745=-0.0039983,5,0,4,7,4,5};
sdgI[125]={sdgI50x4347=-0.0157470,5,0,4,3,4,7};
sdgI[126]={sdgI50x4547=0.0039983,5,0,4,5,4,7};
sdgI[127]={sdgI50x4747=-0.0466446,5,0,4,7,4,7};
sdgI[128]={sdgI50x511511=0.1435908,5,0,5,11,5,11};
sdgI[129]={sdgI51x4743=0.0095146,5,1,4,7,4,3};
sdgI[130]={sdgI51x4545=-1.3011404,5,1,4,5,4,5};
sdgI[131]={sdgI51x4745=0.0024158,5,1,4,7,4,5};
sdgI[132]={sdgI51x4347=0.0095146,5,1,4,3,4,7};
sdgI[133]={sdgI51x4547=-0.0024158,5,1,4,5,4,7};
sdgI[134]={sdgI51x4747=0.0281839,5,1,4,7,4,7};
sdgI[135]={sdgI51x511511=0.0867601,5,1,5,11,5,11};
sdgIlength=135;

pfGe70[1]={pfGe00x0101=5.291502,0,0,0,1,0,1};
pfGe70[2]={pfGe00x1111=5.291502,0,0,1,1,1,1};
pfGe70[3]={pfGe00x1313=7.483314,0,0,1,3,1,3};
pfGe70[4]={pfGe00x2121=5.291502,0,0,2,1,2,1};
pfGe70[5]={pfGe00x2323=7.483314,0,0,2,3,2,3};
pfGe70[6]={pfGe00x2525=9.165151,0,0,2,5,2,5};
pfGe70[7]={pfGe00x3131=2.5701432,0,0,3,1,3,1};
pfGe70[8]={pfGe00x3333=4.4035658,0,0,3,3,3,3};
pfGe70[9]={pfGe00x3535=4.5331907,0,0,3,5,3,5};
pfGe70[10]={pfGe00x3737=10.58300,0,0,3,7,3,7};
pfGe70[11]={pfGe00x4949=0.8274453,0,0,4,9,4,9};
pfGe70[12]={pfGe01x3131=0.3360468,0,1,3,1,3,1};
pfGe70[13]={pfGe01x3333=0.9629148,0,1,3,3,3,3};
pfGe70[14]={pfGe01x3535=1.4160823,0,1,3,5,3,5};
pfGe70[15]={pfGe01x4949=0.5058569,0,1,4,9,4,9};
pfGe70length=15;

pfGe72[1]={pfGe00x0101=6.000000,0,0,0,1,0,1};
pfGe72[2]={pfGe00x1111=6.000000,0,0,1,1,1,1};
pfGe72[3]={pfGe00x1313=8.485281,0,0,1,3,1,3};
pfGe72[4]={pfGe00x2121=6.000000,0,0,2,1,2,1};
pfGe72[5]={pfGe00x2323=8.485281,0,0,2,3,2,3};
pfGe72[6]={pfGe00x2525=10.39230,0,0,2,5,2,5};
pfGe72[7]={pfGe00x3131=3.5592558,0,0,3,1,3,1};
pfGe72[8]={pfGe00x3333=7.2327270,0,0,3,3,3,3};
pfGe72[9]={pfGe00x3535=5.4209780,0,0,3,5,3,5};
pfGe72[10]={pfGe00x3737=12.00000,0,0,3,7,3,7};
pfGe72[11]={pfGe00x4949=0.3664625,0,0,4,9,4,9};
pfGe72[12]={pfGe01x3131=1.4485508,0,1,3,1,3,1};
pfGe72[13]={pfGe01x3333=0.6166944,0,1,3,3,3,3};
pfGe72[14]={pfGe01x3535=2.9687550,0,1,3,5,3,5};
pfGe72[15]={pfGe01x4949=0.1259596,0,1,4,9,4,9};
pfGe72length=15;

pfGe73[1]={pfGe00x0101=20.00000,0,0,0,1,0,1};
pfGe73[2]={pfGe00x1111=20.00000,0,0,1,1,1,1};
pfGe73[3]={pfGe00x1313=28.28427,0,0,1,3,1,3};
pfGe73[4]={pfGe00x2121=20.00000,0,0,2,1,2,1};
pfGe73[5]={pfGe00x2323=28.28427,0,0,2,3,2,3};
pfGe73[6]={pfGe00x2525=34.64101,0,0,2,5,2,5};
pfGe73[7]={pfGe00x3131=12.4173870,0,0,3,1,3,1};
pfGe73[8]={pfGe00x3333=21.5749149,0,0,3,3,3,3};
pfGe73[9]={pfGe00x3535=19.4369812,0,0,3,5,3,5};
pfGe73[10]={pfGe00x3737=40.00000,0,0,3,7,3,7};
pfGe73[11]={pfGe00x4949=3.7456708,0,0,4,9,4,9};
pfGe73[12]={pfGe01x3131=4.1355395,0,1,3,1,3,1};
pfGe73[13]={pfGe01x3333=3.5267536,0,1,3,3,3,3};
pfGe73[14]={pfGe01x3535=8.6396446,0,1,3,5,3,5};
pfGe73[15]={pfGe01x4949=2.0663068,0,1,4,9,4,9};
pfGe73[16]={pfGe10x3131=0.0224806,1,0,3,1,3,1};
pfGe73[17]={pfGe10x3331=0.0169551,1,0,3,3,3,1};
pfGe73[18]={pfGe10x3133=-0.0169551,1,0,3,1,3,3};
pfGe73[19]={pfGe10x3333=0.1571694,1,0,3,3,3,3};
pfGe73[20]={pfGe10x3533=0.0053792,1,0,3,5,3,3};
pfGe73[21]={pfGe10x3335=-0.0053792,1,0,3,3,3,5};
pfGe73[22]={pfGe10x3535=0.1165488,1,0,3,5,3,5};
pfGe73[23]={pfGe10x4949=2.1414728,1,0,4,9,4,9};
pfGe73[24]={pfGe11x3131=-0.0163453,1,1,3,1,3,1};
pfGe73[25]={pfGe11x3331=-0.0088582,1,1,3,3,3,1};
pfGe73[26]={pfGe11x3133=0.0088582,1,1,3,1,3,3};
pfGe73[27]={pfGe11x3333=-0.0709766,1,1,3,3,3,3};
pfGe73[28]={pfGe11x3533=-0.0109537,1,1,3,5,3,3};
pfGe73[29]={pfGe11x3335=0.0109537,1,1,3,3,3,5};
pfGe73[30]={pfGe11x3535=-0.0443756,1,1,3,5,3,5};
pfGe73[31]={pfGe11x4949=1.3447058,1,1,4,9,4,9};
pfGe73[32]={pfGe20x3331=0.9937153,2,0,3,3,3,1};
pfGe73[33]={pfGe20x3531=0.4791123,2,0,3,5,3,1};
pfGe73[34]={pfGe20x3133=-0.9937153,2,0,3,1,3,3};
pfGe73[35]={pfGe20x3333=0.7912950,2,0,3,3,3,3};
pfGe73[36]={pfGe20x3533=-0.3084388,2,0,3,5,3,3};
pfGe73[37]={pfGe20x3135=0.4791123,2,0,3,1,3,5};
pfGe73[38]={pfGe20x3335=0.3084388,2,0,3,3,3,5};
pfGe73[39]={pfGe20x3535=0.9502520,2,0,3,5,3,5};
pfGe73[40]={pfGe20x4949=1.7906727,2,0,4,9,4,9};
pfGe73[41]={pfGe21x3331=-0.5251137,2,1,3,3,3,1};
pfGe73[42]={pfGe21x3531=-0.1589314,2,1,3,5,3,1};
pfGe73[43]={pfGe21x3133=0.5251137,2,1,3,1,3,3};
pfGe73[44]={pfGe21x3333=-0.4086682,2,1,3,3,3,3};
pfGe73[45]={pfGe21x3533=0.1220642,2,1,3,5,3,3};
pfGe73[46]={pfGe21x3135=-0.1589314,2,1,3,1,3,5};
pfGe73[47]={pfGe21x3335=-0.1220642,2,1,3,3,3,5};
pfGe73[48]={pfGe21x3535=-0.4476262,2,1,3,5,3,5};
pfGe73[49]={pfGe21x4949=1.1154441,2,1,4,9,4,9};
pfGe73[50]={pfGe30x3531=0.0106724,3,0,3,5,3,1};
pfGe73[51]={pfGe30x3333=-0.0092243,3,0,3,3,3,3};
pfGe73[52]={pfGe30x3533=-0.0063204,3,0,3,5,3,3};
pfGe73[53]={pfGe30x3135=0.0106724,3,0,3,1,3,5};
pfGe73[54]={pfGe30x3335=0.0063204,3,0,3,3,3,5};
pfGe73[55]={pfGe30x3535=0.0843433,3,0,3,5,3,5};
pfGe73[56]={pfGe30x4949=1.8186721,3,0,4,9,4,9};
pfGe73[57]={pfGe31x3531=-0.0093953,3,1,3,5,3,1};
pfGe73[58]={pfGe31x3333=0.0144374,3,1,3,3,3,3};
pfGe73[59]={pfGe31x3533=0.0036504,3,1,3,5,3,3};
pfGe73[60]={pfGe31x3135=-0.0093953,3,1,3,1,3,5};
pfGe73[61]={pfGe31x3335=-0.0036504,3,1,3,3,3,5};
pfGe73[62]={pfGe31x3535=-0.0575719,3,1,3,5,3,5};
pfGe73[63]={pfGe31x4949=1.1527181,3,1,4,9,4,9};
pfGe73[64]={pfGe40x3533=-0.3054407,4,0,3,5,3,3};
pfGe73[65]={pfGe40x3335=0.3054407,4,0,3,3,3,5};
pfGe73[66]={pfGe40x3535=0.1956187,4,0,3,5,3,5};
pfGe73[67]={pfGe40x4949=1.4346717,4,0,4,9,4,9};
pfGe73[68]={pfGe41x3533=0.1603219,4,1,3,5,3,3};
pfGe73[69]={pfGe41x3335=-0.1603219,4,1,3,3,3,5};
pfGe73[70]={pfGe41x3535=-0.1069728,4,1,3,5,3,5};
pfGe73[71]={pfGe41x4949=0.9149673,4,1,4,9,4,9};
pfGe73[72]={pfGe50x3535=0.0239712,5,0,3,5,3,5};
pfGe73[73]={pfGe50x4949=1.5819857,5,0,4,9,4,9};
pfGe73[74]={pfGe51x3535=-0.0229403,5,1,3,5,3,5};
pfGe73[75]={pfGe51x4949=1.0073045,5,1,4,9,4,9};
pfGe73[76]={pfGe60x4949=1.3630590,6,0,4,9,4,9};
pfGe73[77]={pfGe61x4949=0.8746149,6,1,4,9,4,9};
pfGe73[78]={pfGe70x4949=1.5237149,7,0,4,9,4,9};
pfGe73[79]={pfGe71x4949=0.9716610,7,1,4,9,4,9};
pfGe73[80]={pfGe80x4949=1.3834725,8,0,4,9,4,9};
pfGe73[81]={pfGe81x4949=0.8862418,8,1,4,9,4,9};
pfGe73[82]={pfGe90x4949=1.9384257,9,0,4,9,4,9};
pfGe73[83]={pfGe91x4949=1.2358058,9,1,4,9,4,9};
pfGe73length=83;

pfGe74[1]={pfGe00x0101=6.633249,0,0,0,1,0,1};
pfGe74[2]={pfGe00x1111=6.633249,0,0,1,1,1,1};
pfGe74[3]={pfGe00x1313=9.380831,0,0,1,3,1,3};
pfGe74[4]={pfGe00x2121=6.633249,0,0,2,1,2,1};
pfGe74[5]={pfGe00x2323=9.380831,0,0,2,3,2,3};
pfGe74[6]={pfGe00x2525=11.48912,0,0,2,5,2,5};
pfGe74[7]={pfGe00x3131=4.1213665,0,0,3,1,3,1};
pfGe74[8]={pfGe00x3333=6.2867717,0,0,3,3,3,3};
pfGe74[9]={pfGe00x3535=6.9126963,0,0,3,5,3,5};
pfGe74[10]={pfGe00x3737=13.26649,0,0,3,7,3,7};
pfGe74[11]={pfGe00x4949=2.1649847,0,0,4,9,4,9};
pfGe74[12]={pfGe01x3131=1.1887984,0,1,3,1,3,1};
pfGe74[13]={pfGe01x3333=1.6090327,0,1,3,3,3,3};
pfGe74[14]={pfGe01x3535=2.3665826,0,1,3,5,3,5};
pfGe74[15]={pfGe01x4949=1.3026152,0,1,4,9,4,9};
pfGe74length=15;

pfGe76[1]={pfGe00x0101=7.211102,0,0,0,1,0,1};
pfGe76[2]={pfGe00x1111=7.211102,0,0,1,1,1,1};
pfGe76[3]={pfGe00x1313=10.19803,0,0,1,3,1,3};
pfGe76[4]={pfGe00x2121=7.211102,0,0,2,1,2,1};
pfGe76[5]={pfGe00x2323=10.19803,0,0,2,3,2,3};
pfGe76[6]={pfGe00x2525=12.48999,0,0,2,5,2,5};
pfGe76[7]={pfGe00x3131=4.4121060,0,0,3,1,3,1};
pfGe76[8]={pfGe00x3333=6.4540329,0,0,3,3,3,3};
pfGe76[9]={pfGe00x3535=7.8256316,0,0,3,5,3,5};
pfGe76[10]={pfGe00x3737=14.42220,0,0,3,7,3,7};
pfGe76[11]={pfGe00x4949=3.9898245,0,0,4,9,4,9};
pfGe76[12]={pfGe01x3131=1.2103980,0,1,3,1,3,1};
pfGe76[13]={pfGe01x3333=1.8588700,0,1,3,3,3,3};
pfGe76[14]={pfGe01x3535=2.4187805,0,1,3,5,3,5};
pfGe76[15]={pfGe01x4949=2.4333946,0,1,4,9,4,9};
pfGe76length=15;

sdNa[1]={sdNa00x0101=5.6568542,0,0,0,1,0,1};
sdNa[2]={sdNa00x1111=5.6568542,0,0,1,1,1,1};
sdNa[3]={sdNa00x1313=8.0000000,0,0,1,3,1,3};
sdNa[4]={sdNa00x2121=1.2899591,0,0,2,1,2,1};
sdNa[5]={sdNa00x2323=0.7714725,0,0,2,3,2,3};
sdNa[6]={sdNa00x2525=4.3408126,0,0,2,5,2,5};
sdNa[7]={sdNa01x2121=-0.0364236,0,1,2,1,2,1};
sdNa[8]={sdNa01x2323=0.2040327,0,1,2,3,2,3};
sdNa[9]={sdNa01x2525=0.6709343,0,1,2,5,2,5};
sdNa[10]={sdNa10x2121=-0.1075534,1,0,2,1,2,1};
sdNa[11]={sdNa10x2321=-0.0171936,1,0,2,3,2,1};
sdNa[12]={sdNa10x2123=0.0171936,1,0,2,1,2,3};
sdNa[13]={sdNa10x2323=0.1014424,1,0,2,3,2,3};
sdNa[14]={sdNa10x2523=-0.0577659,1,0,2,5,2,3};
sdNa[15]={sdNa10x2325=0.0577659,1,0,2,3,2,5};
sdNa[16]={sdNa10x2525=0.4984805,1,0,2,5,2,5};
sdNa[17]={sdNa11x2121=0.0969449,1,1,2,1,2,1};
sdNa[18]={sdNa11x2321=0.0724723,1,1,2,3,2,1};
sdNa[19]={sdNa11x2123=-0.0724723,1,1,2,1,2,3};
sdNa[20]={sdNa11x2323=-0.0308415,1,1,2,3,2,3};
sdNa[21]={sdNa11x2523=0.0994210,1,1,2,5,2,3};
sdNa[22]={sdNa11x2325=-0.0994210,1,1,2,3,2,5};
sdNa[23]={sdNa11x2525=-0.2917088,1,1,2,5,2,5};
sdNa[24]={sdNa20x2321=0.1636153,2,0,2,3,2,1};
sdNa[25]={sdNa20x2521=-0.5865685,2,0,2,5,2,1};
sdNa[26]={sdNa20x2123=-0.1636153,2,0,2,1,2,3};
sdNa[27]={sdNa20x2323=-0.0314851,2,0,2,3,2,3};
sdNa[28]={sdNa20x2523=-0.2891304,2,0,2,5,2,3};
sdNa[29]={sdNa20x2125=-0.5865685,2,0,2,1,2,5};
sdNa[30]={sdNa20x2325=0.2891304,2,0,2,3,2,5};
sdNa[31]={sdNa20x2525=-0.5494576,2,0,2,5,2,5};
sdNa[32]={sdNa21x2321=0.0462096,2,1,2,3,2,1};
sdNa[33]={sdNa21x2521=0.0988956,2,1,2,5,2,1};
sdNa[34]={sdNa21x2123=-0.0462096,2,1,2,1,2,3};
sdNa[35]={sdNa21x2323=0.0149520,2,1,2,3,2,3};
sdNa[36]={sdNa21x2523=-0.0731602,2,1,2,5,2,3};
sdNa[37]={sdNa21x2125=0.0988956,2,1,2,1,2,5};
sdNa[38]={sdNa21x2325=0.0731602,2,1,2,3,2,5};
sdNa[39]={sdNa21x2525=-0.1166088,2,1,2,5,2,5};
sdNa[40]={sdNa30x2521=0.0429231,3,0,2,5,2,1};
sdNa[41]={sdNa30x2323=0.0594685,3,0,2,3,2,3};
sdNa[42]={sdNa30x2523=0.0558610,3,0,2,5,2,3};
sdNa[43]={sdNa30x2125=0.0429231,3,0,2,1,2,5};
sdNa[44]={sdNa30x2325=-0.0558610,3,0,2,3,2,5};
sdNa[45]={sdNa30x2525=-0.5268011,3,0,2,5,2,5};
sdNa[46]={sdNa31x2521=-0.0570678,3,1,2,5,2,1};
sdNa[47]={sdNa31x2323=-0.0853394,3,1,2,3,2,3};
sdNa[48]={sdNa31x2523=-0.0815626,3,1,2,5,2,3};
sdNa[49]={sdNa31x2125=-0.0570678,3,1,2,1,2,5};
sdNa[50]={sdNa31x2325=0.0815626,3,1,2,3,2,5};
sdNa[51]={sdNa31x2525=0.5382747,3,1,2,5,2,5};
sdNalength=51;

sdgXe128[1]={sdgXe00x0101=9.1651,0,0,0,1,0,1};
sdgXe128[2]={sdgXe00x1111=9.1651,0,0,1,1,1,1};
sdgXe128[3]={sdgXe00x1313=12.961,0,0,1,3,1,3};
sdgXe128[4]={sdgXe00x2121=9.1651,0,0,2,1,2,1};
sdgXe128[5]={sdgXe00x2323=12.961,0,0,2,3,2,3};
sdgXe128[6]={sdgXe00x2525=15.874,0,0,2,5,2,5};
sdgXe128[7]={sdgXe00x3131=9.1651,0,0,3,1,3,1};
sdgXe128[8]={sdgXe00x3333=12.961,0,0,3,3,3,3};
sdgXe128[9]={sdgXe00x3535=15.874,0,0,3,5,3,5};
sdgXe128[10]={sdgXe00x3737=18.330,0,0,3,7,3,7};
sdgXe128[11]={sdgXe00x4141=5.0163569,0,0,4,1,4,1};
sdgXe128[12]={sdgXe00x4343=7.1678524,0,0,4,3,4,3};
sdgXe128[13]={sdgXe00x4545=11.5199689,0,0,4,5,4,5};
sdgXe128[14]={sdgXe00x4747=9.9389543,0,0,4,7,4,7};
sdgXe128[15]={sdgXe00x4949=20.493,0,0,4,9,4,9};
sdgXe128[16]={sdgXe00x511511=3.7415199,0,0,5,11,5,11};
sdgXe128[17]={sdgXe01x4141=2.5114359,0,1,4,1,4,1};
sdgXe128[18]={sdgXe01x4343=3.5071084,0,1,4,3,4,3};
sdgXe128[19]={sdgXe01x4545=2.6354470,0,1,4,5,4,5};
sdgXe128[20]={sdgXe01x4747=5.0796470,0,1,4,7,4,7};
sdgXe128[21]={sdgXe01x511511=2.2656025,0,1,5,11,5,11};
sdgXe128length=21;

sdgXe129[1]={sdgXe00x0101=13.26649,0,0,0,1,0,1};
sdgXe129[2]={sdgXe00x1111=13.26649,0,0,1,1,1,1};
sdgXe129[3]={sdgXe00x1313=18.76166,0,0,1,3,1,3};
sdgXe129[4]={sdgXe00x2121=13.26649,0,0,2,1,2,1};
sdgXe129[5]={sdgXe00x2323=18.76166,0,0,2,3,2,3};
sdgXe129[6]={sdgXe00x2525=22.97825,0,0,2,5,2,5};
sdgXe129[7]={sdgXe00x3131=13.26649,0,0,3,1,3,1};
sdgXe129[8]={sdgXe00x3333=18.76166,0,0,3,3,3,3};
sdgXe129[9]={sdgXe00x3535=22.97825,0,0,3,5,3,5};
sdgXe129[10]={sdgXe00x3737=26.53299,0,0,3,7,3,7};
sdgXe129[11]={sdgXe00x4141=5.7196121,0,0,4,1,4,1};
sdgXe129[12]={sdgXe00x4343=7.6784286,0,0,4,3,4,3};
sdgXe129[13]={sdgXe00x4545=15.2459268,0,0,4,5,4,5};
sdgXe129[14]={sdgXe00x4747=16.6382923,0,0,4,7,4,7};
sdgXe129[15]={sdgXe00x4949=29.66479,0,0,4,9,4,9};
sdgXe129[16]={sdgXe00x511511=8.1228113,0,0,5,11,5,11};
sdgXe129[17]={sdgXe01x4141=3.4558932,0,1,4,1,4,1};
sdgXe129[18]={sdgXe01x4343=4.6394448,0,1,4,3,4,3};
sdgXe129[19]={sdgXe01x4545=4.6675825,0,1,4,5,4,5};
sdgXe129[20]={sdgXe01x4747=5.9734473,0,1,4,7,4,7};
sdgXe129[21]={sdgXe01x511511=4.9079494,0,1,5,11,5,11};
sdgXe129[22]={sdgXe10x4141=0.9042038,1,0,4,1,4,1};
sdgXe129[23]={sdgXe10x4341=-0.0024811,1,0,4,3,4,1};
sdgXe129[24]={sdgXe10x4143=0.0024811,1,0,4,1,4,3};
sdgXe129[25]={sdgXe10x4343=-0.7546895,1,0,4,3,4,3};
sdgXe129[26]={sdgXe10x4545=0.1387504,1,0,4,5,4,5};
sdgXe129[27]={sdgXe10x4745=-0.0121087,1,0,4,7,4,5};
sdgXe129[28]={sdgXe10x4547=0.0121087,1,0,4,5,4,7};
sdgXe129[29]={sdgXe10x4747=0.1139526,1,0,4,7,4,7};
sdgXe129[30]={sdgXe10x511511=0.1732518,1,0,5,11,5,11};
sdgXe129[31]={sdgXe11x4141=0.5463362,1,1,4,1,4,1};
sdgXe129[32]={sdgXe11x4341=-0.0014991,1,1,4,3,4,1};
sdgXe129[33]={sdgXe11x4143=0.0014991,1,1,4,1,4,3};
sdgXe129[34]={sdgXe11x4343=-0.4559973,1,1,4,3,4,3};
sdgXe129[35]={sdgXe11x4545=-0.0838353,1,1,4,5,4,5};
sdgXe129[36]={sdgXe11x4745=0.0073162,1,1,4,7,4,5};
sdgXe129[37]={sdgXe11x4547=-0.0073162,1,1,4,5,4,7};
sdgXe129[38]={sdgXe11x4747=-0.0688515,1,1,4,7,4,7};
sdgXe129[39]={sdgXe11x511511=0.1046819,1,1,5,11,5,11};
sdgXe129length=39;

sdgXe130[1]={sdgXe00x0101=9.5916,0,0,0,1,0,1};
sdgXe130[2]={sdgXe00x1111=9.5916,0,0,1,1,1,1};
sdgXe130[3]={sdgXe00x1313=13.564,0,0,1,3,1,3};
sdgXe130[4]={sdgXe00x2121=9.5916,0,0,2,1,2,1};
sdgXe130[5]={sdgXe00x2323=13.564,0,0,2,3,2,3};
sdgXe130[6]={sdgXe00x2525=16.613,0,0,2,5,2,5};
sdgXe130[7]={sdgXe00x3131=9.5916,0,0,3,1,3,1};
sdgXe130[8]={sdgXe00x3333=13.564,0,0,3,3,3,3};
sdgXe130[9]={sdgXe00x3535=16.613,0,0,3,5,3,5};
sdgXe130[10]={sdgXe00x3737=19.183,0,0,3,7,3,7};
sdgXe130[11]={sdgXe00x4141=5.2442660,0,0,4,1,4,1};
sdgXe130[12]={sdgXe00x4343=7.4567446,0,0,4,3,4,3};
sdgXe130[13]={sdgXe00x4545=11.3239822,0,0,4,5,4,5};
sdgXe130[14]={sdgXe00x4747=11.0671272,0,0,4,7,4,7};
sdgXe130[15]={sdgXe00x4949=21.447,0,0,4,9,4,9};
sdgXe130[16]={sdgXe00x511511=5.8730974,0,0,5,11,5,11};
sdgXe130[17]={sdgXe01x4141=2.6201291,0,1,4,1,4,1};
sdgXe130[18]={sdgXe01x4343=3.6811604,0,1,4,3,4,3};
sdgXe130[19]={sdgXe01x4545=3.1870365,0,1,4,5,4,5};
sdgXe130[20]={sdgXe01x4747=4.8913521,0,1,4,7,4,7};
sdgXe130[21]={sdgXe01x511511=3.5416116,0,1,5,11,5,11};
sdgXe130length=21;

sdgXe131[1]={sdgXe00x0101=19.59591,0,0,0,1,0,1};
sdgXe131[2]={sdgXe00x1111=19.59591,0,0,1,1,1,1};
sdgXe131[3]={sdgXe00x1313=27.71281,0,0,1,3,1,3};
sdgXe131[4]={sdgXe00x2121=19.59591,0,0,2,1,2,1};
sdgXe131[5]={sdgXe00x2323=27.71281,0,0,2,3,2,3};
sdgXe131[6]={sdgXe00x2525=33.94112,0,0,2,5,2,5};
sdgXe131[7]={sdgXe00x3131=19.59591,0,0,3,1,3,1};
sdgXe131[8]={sdgXe00x3333=27.71281,0,0,3,3,3,3};
sdgXe131[9]={sdgXe00x3535=33.94112,0,0,3,5,3,5};
sdgXe131[10]={sdgXe00x3737=39.19183,0,0,3,7,3,7};
sdgXe131[11]={sdgXe00x4141=9.0950002,0,0,4,1,4,1};
sdgXe131[12]={sdgXe00x4343=10.8862533,0,0,4,3,4,3};
sdgXe131[13]={sdgXe00x4545=21.2939682,0,0,4,5,4,5};
sdgXe131[14]={sdgXe00x4747=25.6417293,0,0,4,7,4,7};
sdgXe131[15]={sdgXe00x4949=43.81780,0,0,4,9,4,9};
sdgXe131[16]={sdgXe00x511511=15.9975309,0,0,5,11,5,11};
sdgXe131[17]={sdgXe01x4141=5.4745469,0,1,4,1,4,1};
sdgXe131[18]={sdgXe01x4343=6.5527534,0,1,4,3,4,3};
sdgXe131[19]={sdgXe01x4545=7.6082420,0,1,4,5,4,5};
sdgXe131[20]={sdgXe01x4747=8.1510658,0,1,4,7,4,7};
sdgXe131[21]={sdgXe01x511511=9.6293792,0,1,5,11,5,11};
sdgXe131[22]={sdgXe10x4141=-0.3117988,1,0,4,1,4,1};
sdgXe131[23]={sdgXe10x4341=0.6189256,1,0,4,3,4,1};
sdgXe131[24]={sdgXe10x4143=-0.6189256,1,0,4,1,4,3};
sdgXe131[25]={sdgXe10x4343=2.1510624,1,0,4,3,4,3};
sdgXe131[26]={sdgXe10x4545=0.1025763,1,0,4,5,4,5};
sdgXe131[27]={sdgXe10x4745=-0.0276369,1,0,4,7,4,5};
sdgXe131[28]={sdgXe10x4547=0.0276369,1,0,4,5,4,7};
sdgXe131[29]={sdgXe10x4747=0.1557883,1,0,4,7,4,7};
sdgXe131[30]={sdgXe10x511511=0.1432732,1,0,5,11,5,11};
sdgXe131[31]={sdgXe11x4141=-0.1876811,1,1,4,1,4,1};
sdgXe131[32]={sdgXe11x4341=0.3725494,1,1,4,3,4,1};
sdgXe131[33]={sdgXe11x4143=-0.3725494,1,1,4,1,4,3};
sdgXe131[34]={sdgXe11x4343=1.2947871,1,1,4,3,4,3};
sdgXe131[35]={sdgXe11x4545=-0.0617405,1,1,4,5,4,5};
sdgXe131[36]={sdgXe11x4745=0.0166354,1,1,4,7,4,5};
sdgXe131[37]={sdgXe11x4547=-0.0166354,1,1,4,5,4,7};
sdgXe131[38]={sdgXe11x4747=-0.0937705,1,1,4,7,4,7};
sdgXe131[39]={sdgXe11x511511=0.0862400,1,1,5,11,5,11};
sdgXe131[40]={sdgXe20x4341=0.2741695,2,0,4,3,4,1};
sdgXe131[41]={sdgXe20x4143=-0.2741695,2,0,4,1,4,3};
sdgXe131[42]={sdgXe20x4343=-1.3365351,2,0,4,3,4,3};
sdgXe131[43]={sdgXe20x4545=-1.1115330,2,0,4,5,4,5};
sdgXe131[44]={sdgXe20x4745=0.2959931,2,0,4,7,4,5};
sdgXe131[45]={sdgXe20x4547=-0.2959931,2,0,4,5,4,7};
sdgXe131[46]={sdgXe20x4747=-1.4990736,2,0,4,7,4,7};
sdgXe131[47]={sdgXe20x511511=-1.7205836,2,0,5,11,5,11};
sdgXe131[48]={sdgXe21x4341=0.1650306,2,1,4,3,4,1};
sdgXe131[49]={sdgXe21x4143=-0.1650306,2,1,4,1,4,3};
sdgXe131[50]={sdgXe21x4343=-0.8044994,2,1,4,3,4,3};
sdgXe131[51]={sdgXe21x4545=0.6690664,2,1,4,5,4,5};
sdgXe131[52]={sdgXe21x4745=-0.1781669,2,1,4,7,4,5};
sdgXe131[53]={sdgXe21x4547=0.1781669,2,1,4,5,4,7};
sdgXe131[54]={sdgXe21x4747=0.9023374,2,1,4,7,4,7};
sdgXe131[55]={sdgXe21x511511=-1.0356696,2,1,5,11,5,11};
sdgXe131[56]={sdgXe30x4343=2.1419019,3,0,4,3,4,3};
sdgXe131[57]={sdgXe30x4545=0.0558688,3,0,4,5,4,5};
sdgXe131[58]={sdgXe30x4745=-0.0105682,3,0,4,7,4,5};
sdgXe131[59]={sdgXe30x4547=0.0105682,3,0,4,5,4,7};
sdgXe131[60]={sdgXe30x4747=0.0474560,3,0,4,7,4,7};
sdgXe131[61]={sdgXe30x511511=0.0672938,3,0,5,11,5,11};
sdgXe131[62]={sdgXe31x4343=1.2892731,3,1,4,3,4,3};
sdgXe131[63]={sdgXe31x4545=-0.0336297,3,1,4,5,4,5};
sdgXe131[64]={sdgXe31x4745=0.0063613,3,1,4,7,4,5};
sdgXe131[65]={sdgXe31x4547=-0.0063613,3,1,4,5,4,7};
sdgXe131[66]={sdgXe31x4747=-0.0285656,3,1,4,7,4,7};
sdgXe131[67]={sdgXe31x511511=0.0405060,3,1,5,11,5,11};
sdgXe131length=67;

sdgXe132[1]={sdgXe00x0101=10.00000,0,0,0,1,0,1};
sdgXe132[2]={sdgXe00x1111=10.00000,0,0,1,1,1,1};
sdgXe132[3]={sdgXe00x1313=14.14213,0,0,1,3,1,3};
sdgXe132[4]={sdgXe00x2121=10.00000,0,0,2,1,2,1};
sdgXe132[5]={sdgXe00x2323=14.14213,0,0,2,3,2,3};
sdgXe132[6]={sdgXe00x2525=17.32050,0,0,2,5,2,5};
sdgXe132[7]={sdgXe00x3131=10.00000,0,0,3,1,3,1};
sdgXe132[8]={sdgXe00x3333=14.14213,0,0,3,3,3,3};
sdgXe132[9]={sdgXe00x3535=17.32050,0,0,3,5,3,5};
sdgXe132[10]={sdgXe00x3737=20.00000,0,0,3,7,3,7};
sdgXe132[11]={sdgXe00x4141=5.4107909,0,0,4,1,4,1};
sdgXe132[12]={sdgXe00x4343=7.6701021,0,0,4,3,4,3};
sdgXe132[13]={sdgXe00x4545=11.3384866,0,0,4,5,4,5};
sdgXe132[14]={sdgXe00x4747=12.0487699,0,0,4,7,4,7};
sdgXe132[15]={sdgXe00x4949=22.36068,0,0,4,9,4,9};
sdgXe132[16]={sdgXe00x511511=8.1643047,0,0,5,11,5,11};
sdgXe132[17]={sdgXe01x4141=2.7571082,0,1,4,1,4,1};
sdgXe132[18]={sdgXe01x4343=3.8882699,0,1,4,3,4,3};
sdgXe132[19]={sdgXe01x4545=3.5936000,0,1,4,5,4,5};
sdgXe132[20]={sdgXe01x4747=4.7767691,0,1,4,7,4,7};
sdgXe132[21]={sdgXe01x511511=4.9061365,0,1,5,11,5,11};
sdgXe132length=21;

sdgXe134[1]={sdgXe00x0101=10.39230,0,0,0,1,0,1};
sdgXe134[2]={sdgXe00x1111=10.39230,0,0,1,1,1,1};
sdgXe134[3]={sdgXe00x1313=14.69693,0,0,1,3,1,3};
sdgXe134[4]={sdgXe00x2121=10.39230,0,0,2,1,2,1};
sdgXe134[5]={sdgXe00x2323=14.69693,0,0,2,3,2,3};
sdgXe134[6]={sdgXe00x2525=18.00000,0,0,2,5,2,5};
sdgXe134[7]={sdgXe00x3131=10.39230,0,0,3,1,3,1};
sdgXe134[8]={sdgXe00x3333=14.69693,0,0,3,3,3,3};
sdgXe134[9]={sdgXe00x3535=18.00000,0,0,3,5,3,5};
sdgXe134[10]={sdgXe00x3737=20.78461,0,0,3,7,3,7};
sdgXe134[11]={sdgXe00x4141=4.4653115,0,0,4,1,4,1};
sdgXe134[12]={sdgXe00x4343=6.2070241,0,0,4,3,4,3};
sdgXe134[13]={sdgXe00x4545=9.8459224,0,0,4,5,4,5};
sdgXe134[14]={sdgXe00x4747=14.1426544,0,0,4,7,4,7};
sdgXe134[15]={sdgXe00x4949=23.23790,0,0,4,9,4,9};
sdgXe134[16]={sdgXe00x511511=12.1316423,0,0,5,11,5,11};
sdgXe134[17]={sdgXe01x4141=2.5399527,0,1,4,1,4,1};
sdgXe134[18]={sdgXe01x4343=3.4498467,0,1,4,3,4,3};
sdgXe134[19]={sdgXe01x4545=4.6284947,0,1,4,5,4,5};
sdgXe134[20]={sdgXe01x4747=3.8375825,0,1,4,7,4,7};
sdgXe134[21]={sdgXe01x511511=7.0794744,0,1,5,11,5,11};
sdgXe134length=21;

sdgXe136[1]={sdgXe00x0101=10.77033,0,0,0,1,0,1};
sdgXe136[2]={sdgXe00x1111=10.77033,0,0,1,1,1,1};
sdgXe136[3]={sdgXe00x1313=15.23154,0,0,1,3,1,3};
sdgXe136[4]={sdgXe00x2121=10.77033,0,0,2,1,2,1};
sdgXe136[5]={sdgXe00x2323=15.23154,0,0,2,3,2,3};
sdgXe136[6]={sdgXe00x2525=18.65475,0,0,2,5,2,5};
sdgXe136[7]={sdgXe00x3131=10.77033,0,0,3,1,3,1};
sdgXe136[8]={sdgXe00x3333=15.23154,0,0,3,3,3,3};
sdgXe136[9]={sdgXe00x3535=18.65475,0,0,3,5,3,5};
sdgXe136[10]={sdgXe00x3737=21.54065,0,0,3,7,3,7};
sdgXe136[11]={sdgXe00x4141=5.4804468,0,0,4,1,4,1};
sdgXe136[12]={sdgXe00x4343=7.8045373,0,0,4,3,4,3};
sdgXe136[13]={sdgXe00x4545=10.4772539,0,0,4,5,4,5};
sdgXe136[14]={sdgXe00x4747=14.7770242,0,0,4,7,4,7};
sdgXe136[15]={sdgXe00x4949=24.08318,0,0,4,9,4,9};
sdgXe136[16]={sdgXe00x511511=13.3553810,0,0,5,11,5,11};
sdgXe136[17]={sdgXe01x4141=3.1612961,0,1,4,1,4,1};
sdgXe136[18]={sdgXe01x4343=4.4384675,0,1,4,3,4,3};
sdgXe136[19]={sdgXe01x4545=4.8869710,0,1,4,5,4,5};
sdgXe136[20]={sdgXe01x4747=4.0420193,0,1,4,7,4,7};
sdgXe136[21]={sdgXe01x511511=7.7847480,0,1,5,11,5,11};
sdgXe136length=21;


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


IsotopicSpins={{1,1,1/2},{1,2,1},{1,3,1/2},{2,3,1/2},{2,4,0},{3,6,1},{3,7,3/2},{4,9,3/2},{5,10,3},{5,11,3/2},{6,12,0},{6,13,1/2},{6,14,3},{7,14,1},{7,15,1/2},{8,16,0},{8,17,5/2},{8,18,0},{9,19,1/2},{10,20,0},{10,21,3/2},{10,22,0},{11,22,3},{11,23,3/2},{12,24,0},{12,25,5/2},{12,26,0},{13,27,5/2},{14,28,0},{14,29,1/2},{14,30,0},{15,31,1/2},{16,32,0},{16,33,3/2},{16,34,0},{16,36,0},{17,35,3/2},{17,36,2},{17,37,3/2},{18,36,0},{18,38,0},{18,39,7/2},{18,40,0},{19,39,3/2},{19,40,4},{19,41,3/2},{20,40,0},{20,41,7/2},{20,42,0},{20,43,7/2},{20,44,0},{20,46,0},{20,48,0},{21,45,7/2},{22,46,0},{22,47,5/2},{22,48,0},{22,49,7/2},{22,50,0},{23,50,6},{23,51,7/2},{24,50,0},{24,52,0},{24,53,3/2},{24,54,0},{25,53,7/2},{25,55,5/2},{26,54,0},{26,56,0},{26,57,1/2},{26,58,0},{27,59,7/2},{27,60,5},{28,58,0},{28,60,0},{28,61,3/2},{28,62,0},{28,64,0},{29,63,3/2},{29,65,3/2},{30,64,0},{30,66,0},{30,67,5/2},{30,68,0},{30,70,0},{31,69,3/2},{31,71,3/2},{32,70,0},{32,72,0},{32,73,9/2},{32,74,0},{32,76,0},{33,75,3/2},{34,74,0},{34,76,0},{34,77,1/2},{34,78,0},{34,79,7/2},{34,80,0},{34,82,0},{35,79,3/2},{35,81,3/2},{36,78,0},{36,80,0},{36,82,0},{36,83,9/2},{36,84,0},{36,85,9/2},{36,86,0},{37,85,5/2},{37,87,3/2},{38,84,0},{38,86,0},{38,87,9/2},{38,88,0},{39,89,1/2},{40,90,0},{40,91,5/2},{40,92,0},{40,94,0},{40,96,0},{41,93,9/2},{42,92,0},{42,94,0},{42,95,5/2},{42,96,0},{42,97,5/2},{42,98,0},{42,100,0},{43,99,9/2},{44,96,0},{44,98,0},{44,99,5/2},{44,100,0},{44,101,5/2},{44,102,0},{44,104,0},{45,102,6},{45,103,1/2},{46,102,0},{46,104,0},{46,105,5/2},{46,106,0},{46,108,0},{46,110,0},{47,107,1/2},{47,109,1/2},{48,106,0},{48,108,0},{48,110,0},{48,111,1/2},{48,112,0},{48,113,1/2},{48,114,0},{48,116,0},{49,113,9/2},{49,115,9/2},{50,112,0},{50,114,0},{50,115,1/2},{50,116,0},{50,117,1/2},{50,118,0},{50,119,1/2},{50,120,0},{50,122,0},{50,124,0},{51,121,5/2},{51,123,7/2},{51,125,7/2},{52,120,0},{52,122,0},{52,123,1/2},{52,124,0},{52,125,1/2},{52,126,0},{52,128,0},{52,130,0},{53,127,5/2},{53,129,7/2},{54,124,0},{54,126,0},{54,128,0},{54,129,1/2},{54,130,0},{54,131,3/2},{54,132,0},{54,134,0},{54,136,0},{55,133,7/2},{55,134,4},{55,135,7/2},{55,137,7/2},{56,130,0},{56,132,0},{56,133,1/2},{56,134,0},{56,135,3/2},{56,136,0},{56,137,3/2},{56,138,0},{57,137,7/2},{57,138,5},{57,139,7/2},{58,136,0},{58,138,0},{58,140,0},{58,142,0},{59,141,5/2},{60,142,0},{60,143,7/2},{60,144,0},{60,145,7/2},{60,146,0},{60,148,0},{60,150,0},{61,147,7/2},{62,144,0},{62,147,7/2},{62,148,0},{62,149,7/2},{62,150,0},{62,151,5/2},{62,152,0},{62,154,0},{63,151,5/2},{63,152,3},{63,153,5/2},{63,154,3},{63,155,5/2},{64,152,0},{64,154,0},{64,155,3/2},{64,156,0},{64,157,3/2},{64,158,0},{64,160,0},{65,157,3/2},{65,159,3/2},{65,160,3},{66,156,0},{66,158,0},{66,160,0},{66,161,5/2},{66,162,0},{66,163,5/2},{66,164,0},{67,165,7/2},{68,162,0},{68,164,0},{68,166,0},{68,167,7/2},{68,168,0},{68,170,0},{69,169,1/2},{69,171,1/2},{70,168,0},{70,170,0},{70,171,1/2},{70,172,0},{70,173,5/2},{70,174,0},{70,176,0},{71,173,7/2},{71,174,1},{71,175,7/2},{71,176,7},{72,174,0},{72,176,0},{72,177,7/2},{72,178,0},{72,179,9/2},{72,180,0},{73,180,0},{73,181,7/2},{74,180,0},{74,182,0},{74,183,1/2},{74,184,0},{74,186,0},{75,185,5/2},{75,187,5/2},{76,184,0},{76,186,0},{76,187,1/2},{76,188,0},{76,189,3/2},{76,190,0},{76,192,0},{77,191,3/2},{77,193,3/2},{78,190,0},{78,192,0},{78,194,0},{78,195,1/2},{78,196,0},{78,198,0},{79,197,3/2},{80,196,0},{80,198,0},{80,199,1/2},{80,200,0},{80,201,3/2},{80,202,0},{80,204,0},{81,203,1/2},{81,204,2},{81,205,1/2},{82,204,0},{82,206,0},{82,207,1/2},{82,208,0},{83,207,9/2},{83,209,9/2},{84,209,1/2},{89,227,3/2},{90,229,5/2},{90,232,0},{92,234,0},{92,235,7/2},{92,238,0}};
For[i=1,i<=Length[IsotopicSpins],i++,
IsoJLookupTable[IsotopicSpins[[i,1]],IsotopicSpins[[i,2]]]=IsotopicSpins[[i,3]];
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


(*******
*  Here are the names of the operators - 
 this is useful for printing out the effective Lagrangian
*********)

Op[1]="1";
Op[2]="(\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\)\!\(\*SuperscriptBox[\()\), \(2\)]\)";
Op[3]="\[ImaginaryI] \!\(\*SubscriptBox[\(S\), \(N\)]\)\[CenterDot](q\[Times]\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\))";
Op[4]="\!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\)\[CenterDot]\!\(\*SubscriptBox[\(S\), \(N\)]\)";
Op[5]="\[ImaginaryI] \!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\)\[CenterDot](q\[Times]\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\))";
Op[6]="(\!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\)\[CenterDot]q)(\!\(\*SubscriptBox[\(S\), \(N\)]\)\[CenterDot]q)";
Op[7]="\!\(\*SubscriptBox[\(S\), \(N\)]\)\[CenterDot]\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\)";
Op[8]="\!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\)\[CenterDot]\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\)";
Op[9]="\[ImaginaryI] \!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\).(\!\(\*SubscriptBox[\(S\), \(N\)]\)\[Times]q)";
Op[10]="\[ImaginaryI] \!\(\*SubscriptBox[\(S\), \(N\)]\)\[CenterDot]q";
Op[11]="\[ImaginaryI] \!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\)\[CenterDot]q";
Op[12]="\!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\)\[CenterDot](\!\(\*SubscriptBox[\(S\), \(N\)]\)\[Times]\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\))";
Op[13]="(\[ImaginaryI]q\[Times]\!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\))\[CenterDot](\!\(\*SubscriptBox[\(S\), \(N\)]\)\[CenterDot]\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\))";
Op[14]="(\[ImaginaryI] q\[CenterDot]\!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\))(\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\)\[CenterDot]\!\(\*SubscriptBox[\(S\), \(N\)]\))";
Op[15]="(\[ImaginaryI] q\[CenterDot]\!\(\*SubscriptBox[\(S\), \(\[Chi]\)]\))(\!\(\*SuperscriptBox[\(v\), \(\[Perpendicular]\)]\).(\[ImaginaryI]q\[Times]\!\(\*SubscriptBox[\(S\), \(N\)]\)))";
OpRel[1]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)N";
OpRel[2]="\[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[3]="\[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)N";
OpRel[4]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[5]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N";
OpRel[6]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)N";
OpRel[7]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[8]="\[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[9]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N";
OpRel[10]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)N";
OpRel[11]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[12]="\[ImaginaryI] \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[13]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N";
OpRel[14]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)N";
OpRel[15]="\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[16]="\[ImaginaryI] \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[17]="\[ImaginaryI] \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N";
OpRel[18]="\[ImaginaryI] \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)N";
OpRel[19]="\[ImaginaryI] \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
OpRel[20]="\[ImaginaryI] \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(q\), \(\[Nu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi] \!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*SuperscriptBox[\(q\), \(\[Alpha]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N";
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

If[Op==1,Print["Relativistic Operator 1 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)N"];]
If[Op==2,Print["Relativistic Operator 2 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==3,Print["Relativistic Operator 3 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)N"];]
If[Op==4,Print["Relativistic Operator 4 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==5,Print["Relativistic Operator 5 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N"];]
If[Op==6,Print["Relativistic Operator 6 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)N"];]
If[Op==7,Print["Relativistic Operator 7 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==8,Print["Relativistic Operator 8 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==9,Print["Relativistic Operator 9 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N"];]
If[Op==10,Print["Relativistic Operator 10 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)N"];]
If[Op==11,Print["Relativistic Operator 11 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==12,Print["Relativistic Operator 12 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==13,Print["Relativistic Operator 13 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N"];]
If[Op==14,Print["Relativistic Operator 14 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)N"];]
If[Op==15,Print["Relativistic Operator 15 is \!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==16,Print["Relativistic Operator 16 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==17,Print["Relativistic Operator 17 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)N"];]
If[Op==18,Print["Relativistic Operator 18 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)N"];]
If[Op==19,Print["Relativistic Operator 19 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]
If[Op==20,Print["Relativistic Operator 20 is \[ImaginaryI]\!\(\*OverscriptBox[\(\[Chi]\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\)\!\(\*FractionBox[SubscriptBox[\(q\), \(\[Nu]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)\[Chi]\!\(\*OverscriptBox[\(N\), \(_\)]\)\!\(\*SubscriptBox[\(\[ImaginaryI]\[Sigma]\), \(\[Mu]\[Alpha]\)]\)\!\(\*FractionBox[SuperscriptBox[\(q\), \(\[Alpha]\)], SubscriptBox[\(m\), \(M\)]]\)\!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)N"];]

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


(* ::Text:: *)
(*Section for the Form Factors*)


(* Underlying functions
*)

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

(* DEBUG)
Print["mj coef"]
Print[ResponseCoeff[MJ]]
*)
ResponseCoeff[SigmaPPJ]=Table[1/16 Cl[jchi](cvector[iiFF][[6]]cvector[jjFF][[6]]q^4+(cvector[iiFF][[13]]cvector[jjFF][[13]]q^2+cvector[iiFF][[12]]cvector[jjFF][[12]])(v^2-q^2/(4(\[Mu]T[mchi,M])^2))+2cvector[iiFF][[4]]cvector[jjFF][[6]]q^2+cvector[iiFF][[4]]cvector[jjFF][[4]])+1/4 cvector[iiFF][[10]]cvector[jjFF][[10]]q^2,{iiFF,2},{jjFF,2}];
ResponseCoeff[SigmaPJ]=Table[1/32 Cl[jchi](2cvector[iiFF][[9]]cvector[jjFF][[9]]q^2+(cvector[iiFF][[15]]cvector[jjFF][[15]]q^4+cvector[iiFF][[14]]cvector[jjFF][[14]]q^2-2cvector[iiFF][[12]]cvector[jjFF][[15]]q^2+cvector[iiFF][[12]]cvector[jjFF][[12]])(v^2-q^2/(4(\[Mu]T[mchi,M])^2))+2cvector[iiFF][[4]]cvector[jjFF][[4]])+1/8 (cvector[iiFF][[3]]cvector[jjFF][[3]]q^2+cvector[iiFF][[7]]cvector[jjFF][[7]])(v^2-q^2/(4(\[Mu]T[mchi,M])^2)),{iiFF,2},{jjFF,2}];
ResponseCoeff[DeltaJ]=Table[q^2/(4mN^2) Cl[jchi](cvector[iiFF][[5]]cvector[jjFF][[5]]q^2+cvector[iiFF][[8]]cvector[jjFF][[8]])+2 q^2/mN^2 cvector[iiFF][[2]]cvector[jjFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2)),{iiFF,2},{jjFF,2}];
ResponseCoeff[PhiPPJ]=Table[q^2/(16mN^2) Cl[jchi](cvector[iiFF][[12]]-cvector[iiFF][[15]]q^2)(cvector[jjFF][[12]]-cvector[jjFF][[15]]q^2)+q^4/(4mN^2)cvector[iiFF][[3]]cvector[jjFF][[3]],{iiFF,2},{jjFF,2}];
ResponseCoeff[PhiTPJ]=Table[q^2/(16mN^2)Cl[jchi](cvector[iiFF][[13]]cvector[jjFF][[13]]q^2+cvector[iiFF][[12]]cvector[jjFF][[12]]),{iiFF,2},{jjFF,2}];
ResponseCoeff[MJ,PhiPPJ]=Table[q^2/(4mN)Cl[jchi]cvector[iiFF][[11]](cvector[jjFF][[12]]-cvector[jjFF][[15]]q^2)+q^2/(mN)cvector[jjFF][[3]](cvector[iiFF][[1]]+cvector[iiFF][[2]](v^2-q^2/(4(\[Mu]T[mchi,M])^2))),{iiFF,2},{jjFF,2}];
ResponseCoeff[SigmaPJ,DeltaJ]=Table[q^2/(4mN)Cl[jchi](cvector[iiFF][[4]]cvector[jjFF][[5]]-cvector[jjFF][[8]]cvector[iiFF][[9]])-q^2/(mN)cvector[jjFF][[2]]cvector[iiFF][[3]](v^2-q^2/(4(\[Mu]T[mchi,M])^2)),{iiFF,2},{jjFF,2}];
Sum[ResponseCoeff[MJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,J,T][[iiFF,jjFF]]+ResponseCoeff[SigmaPPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[SigmaPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,DeltaJ,J,T][[iiFF,jjFF]]+ResponseCoeff[PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiPPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[PhiTPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiTPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[MJ,PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,PhiPPJ,J,T][[iiFF,jjFF]]+ResponseCoeff[SigmaPJ,DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,DeltaJ,J,T][[iiFF,jjFF]],{iiFF,2},{jjFF,2}]
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
DMResponseCoeff[MJ] = Table[(cvector[iiFF][[1]]+v^2 cvector[iiFF][[2]])(cvector[jjFF][[1]]+v^2 cvector[jjFF][[2]])+Cl[jchi]/4 (q^2/mN^2 v^2 cvector[iiFF][[5]]cvector[jjFF][[5]] + v^2 cvector[iiFF][[8]]cvector[jjFF][[8]] + q^2/mN^2 cvector[iiFF][[11]]cvector[jjFF][[11]]),{iiFF,2},{jjFF,2}];
DMResponseCoeff[PhiPPJ] = Table[q^2/(4mN^2)cvector[iiFF][[3]] cvector[jjFF][[3]] + Cl[jchi]/16 cvector[iiFF][[12]]cvector[jjFF][[12]],{iiFF,2},{jjFF,2}];
DMResponseCoeff[MJ,PhiPPJ] = Table[cvector[iiFF][[3]](cvector[jjFF][[1]]+v^2 cvector[jjFF][[2]])+Cl[jchi]/4 cvector[iiFF][[12]]cvector[jjFF][[11]],{iiFF,2},{jjFF,2}];
DMResponseCoeff[PhiTPJ] = Table[1/16 Cl[jchi](cvector[iiFF][[12]] cvector[jjFF][[12]] + q^2/mN^2 cvector[iiFF][[13]]cvector[jjFF][[13]]),{iiFF,2},{jjFF,2}];
DMResponseCoeff[SigmaPPJ] = Table[q^2/(4mN^2) cvector[iiFF][[10]]cvector[jjFF][[10]] + 1/16 Cl[jchi](cvector[iiFF][[4]]cvector[jjFF][[4]]+q^2/mN^2(cvector[iiFF][[4]]cvector[jjFF][[6]]+cvector[iiFF][[6]]cvector[jjFF][[4]]) + q^4/mN^4 cvector[iiFF][[6]]cvector[jjFF][[6]] + v^2 cvector[iiFF][[12]]cvector[jjFF][[12]] + q^2/mN^2 v^2cvector[iiFF][[13]]cvector[jjFF][[13]]),{iiFF,2},{jjFF,2}];
DMResponseCoeff[SigmaPJ] = Table[1/8(q^2/mN^2 v^2 cvector[iiFF][[3]]cvector[jjFF][[3]] + v^2 cvector[iiFF][[7]]cvector[jjFF][[7]]) + 1/16 Cl[jchi](cvector[iiFF][[4]]cvector[jjFF][[4]] + q^2/mN^2 cvector[iiFF][[9]]cvector[jjFF][[9]] + v^2/2 cvector[iiFF][[12]]cvector[jjFF][[12]]+q^2/(2mN^2)cvector[iiFF][[14]]cvector[jjFF][[14]]),{iiFF,2},{jjFF,2}];
DMResponseCoeff[DeltaJ] = Table[2v^2 cvector[iiFF][[2]]cvector[jjFF][[2]] + Cl[jchi]/4 (q^2/mN^2 cvector[iiFF][[5]]cvector[jjFF][[5]] + cvector[iiFF][[8]]cvector[jjFF][[8]]),{iiFF,2},{jjFF,2}];
DMResponseCoeff[SigmaPJ, DeltaJ] = Table[-v^2 cvector[iiFF][[2]]cvector[jjFF][[3]] + Cl[jchi]/4 (cvector[iiFF][[5]]cvector[jjFF][[4]]-cvector[iiFF][[8]]cvector[jjFF][[9]]),{iiFF,2},{jjFF,2}];
Sum[DMResponseCoeff[MJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,J,T][[iiFF,jjFF]]+DMResponseCoeff[SigmaPPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPPJ,J,T][[iiFF,jjFF]]+DMResponseCoeff[SigmaPJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,J,T][[iiFF,jjFF]]+DMResponseCoeff[DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,DeltaJ,J,T][[iiFF,jjFF]]+DMResponseCoeff[PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiPPJ,J,T][[iiFF,jjFF]]+DMResponseCoeff[PhiTPJ][[iiFF,jjFF]]*FF[DMmatrix,PhiTPJ,J,T][[iiFF,jjFF]]+DMResponseCoeff[MJ,PhiPPJ][[iiFF,jjFF]]*FF[DMmatrix,MJ,PhiPPJ,J,T][[iiFF,jjFF]]+DMResponseCoeff[SigmaPJ,DeltaJ][[iiFF,jjFF]]*FF[DMmatrix,SigmaPJ,DeltaJ,J,T][[iiFF,jjFF]],{iiFF,2},{jjFF,2}]
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





(***************
* Section for Density Matrix I/O 
*)

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


If[filename=="default"||filename=="Default"||filename=="DEFAULT",
(*set it to the default isotope, in this case sdF*)
Print["Getting default matrix..."];

If[Ziso==2&&Mxx==4,
Print["Setting isotope to helium-4 (default.)"];
DensityMatrix=sHe4;
SetDMmatrixLength[sHe4length];
Return[];];

If[Ziso==8&&Mxx==16,
Print["Setting isotope to oxygen-16 (default.)"];
DensityMatrix=sO16;
SetDMmatrixLength[sO16length];
Return[];];

If[Ziso==9&&Mxx==19,
Print["Setting isotope to fluorine-19 (default.)"];
DensityMatrix=sdF;
SetDMmatrixLength[sdFlength];
Return[];];

(*another test density matrix*)
If[Ziso==53&&Mxx==127,
Print["Setting isotope to iodine-127."];
DensityMatrix=sdgI;
SetDMmatrixLength[sdgIlength];
Return[];];

If[Ziso==32&&Mxx==70,
Print["Setting isotope to germanium-70."];
DensityMatrix=pfGe70;
SetDMmatrixLength[pfGe70length];
Return[];];

If[Ziso==32&&Mxx==72,
Print["Setting isotope to germanium-72."];
DensityMatrix=pfGe72;
SetDMmatrixLength[pfGe72length];
Return[];];

If[Ziso==32&&Mxx==73,
Print["Setting isotope to germanium-73."];
DensityMatrix=pfGe73;
SetDMmatrixLength[pfGe73length];
Return[];];

If[Ziso==32&&Mxx==74,
Print["Setting isotope to germanium-74."];
DensityMatrix=pfGe74;
SetDMmatrixLength[pfGe74length];
Return[];];

If[Ziso==32&&Mxx==76,
Print["Setting isotope to germanium-76."];
DensityMatrix=pfGe76;
SetDMmatrixLength[pfGe76length];
Return[];];

If[Ziso==11,
Print["Setting isotope to sodium-23."];
DensityMatrix=sdNa;
SetDMmatrixLength[sdNalength];
Return[];];

If[Ziso==14&&Mxx==28,
Print["Setting isotope to silicon-28."];
DensityMatrix=sdSi28;
SetDMmatrixLength[sdSi28length];
Return[];];

If[Ziso==14&&Mxx==29,
Print["Setting isotope to silicon-29."];
DensityMatrix=sdSi29;
SetDMmatrixLength[sdSi29length];
Return[];];

If[Ziso==14&&Mxx==30,
Print["Setting isotope to silicon-30."];
DensityMatrix=sdSi30;
SetDMmatrixLength[sdSi30length];
Return[];];

If[Ziso==54&&Mxx==128,
Print["Setting isotope to xenon-128."];
DensityMatrix=sdgXe128;
SetDMmatrixLength[sdgXe128length];
Return[];];

If[Ziso==54&&Mxx==129,
Print["Setting isotope to xenon-129."];
DensityMatrix=sdgXe129;
SetDMmatrixLength[sdgXe129length];
Return[];];

If[Ziso==54&&Mxx==130,
Print["Setting isotope to xenon-130."];
DensityMatrix=sdgXe130;
SetDMmatrixLength[sdgXe130length];
Return[];];

If[Ziso==54&&Mxx==131,
Print["Setting isotope to xenon-131."];
DensityMatrix=sdgXe131;
SetDMmatrixLength[sdgXe131length];
Return[];];

If[Ziso==54&&Mxx==132,
Print["Setting isotope to xenon-132."];
DensityMatrix=sdgXe132;
SetDMmatrixLength[sdgXe132length];
Return[];];

If[Ziso==54&&Mxx==134,
Print["Setting isotope to xenon-134."];
DensityMatrix=sdgXe134;
SetDMmatrixLength[sdgXe134length];
Return[];];

If[Ziso==54&&Mxx==136,
Print["Setting isotope to xenon-136."];
DensityMatrix=sdgXe136;
SetDMmatrixLength[sdgXe136length];
Return[];];

Print["Cannot find default matrix which matches input M and Z!"];
Return[];
,
(*If not default, read specified file*)
SetM[Mxx];
ReadInDMFile[filename];

Return[];];
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




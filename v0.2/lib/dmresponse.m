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
		+cvector[iiFF][[8]]cvector[jjFF][[8]]
	)
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

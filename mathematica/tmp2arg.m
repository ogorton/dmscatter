FF[DMmatrix_,Operator1_,Operator2_,J_,T_]:=Block[{Response},
	Response[iDMmatrix_,iOperator_,iT_]:=Response[iDMmatrix,iOperator,iT]
	= Simplify[
		Block[
			{},
			Sum[
				iDMmatrix[i][[1]]
				(
					If[
						iDMmatrix[i][[3]]==0,
						Sqrt[2]/Sqrt[2iT+1] ((an+ap)/2),
						-(Sqrt[6iT]/Sqrt[(2iT+1)(iT+1)])((ap-an)/2)
					]
				)
				* E^y (
					iOperator[y,{iDMmatrix[i][[4]],iDMmatrix[i][[5]]/2},{iDMmatrix[i][[6]],iDMmatrix[i][[7]]/2},iDMmatrix[i][[2]]]
				)
				Jbase[iDMmatrix[i][[2]]],
		
				{i,DMmatrixLength}
			]
		]
	];

	Simplify[
		(4\[Pi])/(2J+1) 
		Sum[
			Table[
				D[
					(Coefficient[Response[DMmatrix,Operator1,T]/.{ap->ap1,an->an1},Jbase[n]]) (Coefficient[Response[DMmatrix,Operator2,T]/.{ap->ap2,an->an2},Jbase[n]]),
					{ap1,an1}[[ii]],
					{ap2,an2}[[jj]]
				],
	
				{ii,2},
				{jj,2}
			],
	
			{n,0,Jmax}
		]
	]
];

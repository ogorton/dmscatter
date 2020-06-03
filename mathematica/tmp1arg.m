FF[DMmatrix_,Operator_,J_,T_]:=Block[{Response},
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
			1/2 Table[
				D[
					(Coefficient[Response[DMmatrix,Operator,T],Jbase[n]]^2),
					{ap,an}[[ii]],
					{ap,an}[[jj]]
				],
	
				{ii,2},
				{jj,2}
			],
			
			{n,0,Jmax}
		]
	]
];

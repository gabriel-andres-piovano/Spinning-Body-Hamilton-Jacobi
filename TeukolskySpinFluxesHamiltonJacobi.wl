(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["TeukolskySpinFluxes`",
  {
  "SpinWeightedSpheroidalHarmonics`",
  "SpinningOrbit`"
  }
];


(* ::Text:: *)
(*If you make use of this package, please acknowledge *)
(*- "Piovano, Pantelidou, Mac Uilliam and Witzany" arXiv:2410.05769 (https://arxiv.org/abs/2410.05769v1 ) *)
(*- "Piovano, Brito, Maselli, Pani" PRD 104, 124019 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.124019 )*)
(* - "Skoup\[YAcute], Lukes-Gerakopoulos, Drummond, Hughes" PRD 108, 044041 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.044041 )*)
(**)
(*IMPORTANT NOTE: the following package is a modified version of the code used in Skoup\[YAcute] et al, PRD 108, 044041 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.044041 ) kindly provided by Viktor Skoup\[YAcute]. The Mathematica functions in the sections "Functions to compute linearized angular Teukolsky equations" and "Functions to compute linearized radial Teukolsky equations" were developed in Piovano et al PRD 104, 124019 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.124019 ).*)
(*The main functions*)


(* ::ItemNumbered:: *)
(*TeukolskySpinModeCorrectionDH*)


(* ::ItemNumbered:: *)
(*TeukolskySpinModeCorrectionFC*)


(* ::Text:: *)
(*uses the solutions of the spin-corrections to the orbits computed with the Hamilton-Jacobi formalism. See arXiv:2410.05769 (https://arxiv.org/abs/2410.05769v1 ). The source term and the algorithm for computing the amplitudes and fluxes are mostly unchanged from the code of PRD 108, 044041, with the following exceptions:*)
(*- this version of code returns the fully linearize spin-corrections to the amplitudes and fluxes. It does NOT use finite difference methods (see Eq. 66 of PRD 108, 044041)*)
(*- we performed some tweaks to the original code to improve its performance (for instance, we exploit the vectorization of functions in Mathematica to speed up some calculations).*)


(* ::Text:: *)
(*The main functions*)


(* ::ItemNumbered:: *)
(*TeukolskySpinModeCorrectionOrtPlus*)


(* ::ItemNumbered:: *)
(*TeukolskySpinModeCorrectionOrtMinus*)


(* ::Text:: *)
(*employs the solutions of the spin-corrections to the orbits computed with the Hamilton-Jacobi formalism. The source term is modified to include the orthogonal component of the secondary spin, whereas the algorithm for computing the amplitudes are mostly unchanged (minus minor tweaks)  from the code of PRD 108, 044041. *)


(* ::Text:: *)
(*PREREQUISITEs: this package requires the following packages:*)
(*- SpinOrbitsHamiltonJacobi*)
(*- mod_SpinWeightedSpheroidalHarmonics2 (which can be found in the folder of the same name)*)
(**)
(*mod_SpinWeightedSpheroidalHarmonics2 was developed in PRD 104, 124019, and it is a modified version of the SpinWeightedSpheroidalHarmonics package, version 0.3.0 (https://zenodo.org/records/8091168 ) of the  BHPToolkit (https://bhptoolkit.org/ ).  See the notebook Amplitudes_and_fluxes_corrections _v3 .nb on how to load the package mod_SpinWeightedSpheroidalHarmonics2. To load it properly, make sure that mod_SpinWeightedSpheroidalHarmonics2 has the same version number or higher compare to the vanilla SpinWeightedSpheroidalHarmonics package. The version can be checked in PacletInfo.wl contain inside the folder "mod_SpinWeightedSpheroidalHarmonics2". To update the version of the mod_SpinWeightedSpheroidalHarmonics2 , open PacletInfo.wl , increase the number in "Version", then run the command PacletDataRebuild[]*)


TeukolskySpinModeCorrectionDH::usage = "TeukolskySpinModeCorrectionDH[l,m,n,k,orbitCorrection] calculates linear correction to the fluxes and amplitudes in the fixed turning points (on average) parametrization";


TeukolskySpinModeCorrectionFC::usage = "TeukolskySpinModeCorrectionFC[l,m,n,k,orbitCorrection] calculates linear correction to the fluxes and amplitudes in the fixed constants of motion parametrization";


TeukolskySpinModeCorrectionOrtPlus::usage = "TeukolskySpinModeCorrectionOrtPlus[l,m,n,k,orbitCorrection] calculates the corrections to amplitudes due to the orthogonal component of the secondary spin for the +1 precession mode";


TeukolskySpinModeCorrectionOrtMinus::usage = "TeukolskySpinModeCorrectionOrtMinus[l,m,n,k,orbitCorrection] calculates the corrections to amplitudes due to the orthogonal component of the secondary spin for the -1 precession mode";


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Radial and polar geodesic trajectories*)


rgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,jSN},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi] ellK wr,krg];
	
	(r3g(r1g-r2g)jSN^2-r2g(r1g-r3g))/((r1g-r2g)jSN^2-(r1g-r3g))
]


zgfun[wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,z1g,z2g,kzg,ellK},
	EEg=const[[1]];
	Lzg=const[[2]];
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	z1g JacobiSN[2/\[Pi] ellK(wz+\[Pi]/2),kzg]
]


(* ::Subsection::Closed:: *)
(*Radial and polar geodesic velocities*)


dzgd\[Lambda]fun[wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,zg,z1g,z2g,kzg,ellK,Yzg},
	EEg=const[[1]];
	Lzg=const[[2]];
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	Yzg=Sqrt[-a^2 (1-EEg^2)zg^2+z2g^2];
	
	z1g*JacobiCN[2/\[Pi] ellK(wz+\[Pi]/2),kzg]Yzg
]


(* ::Subsection::Closed:: *)
(*Spin precession phase *)


\[Psi]rfun[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,krg,ellK,xaux},
 {EEg,Lzg,KKg}=const;
 r1g=p/(1-e);
 r2g=p/(1+e);
 r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
 r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));

 krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
 ellK=EllipticK[krg];
 xaux=JacobiAmplitude[wr/\[Pi] ellK,krg];

 (2(r2g-r3g)((KKg-a^2)EEg+a Lzg))/((KKg+r2g^2)(KKg+r3g^2)Sqrt[(1-EEg^2)(r1g-r3g)(r2g-r4g)]) (Im[(r2g+I Sqrt[KKg])(r3g+I Sqrt[KKg])EllipticPi[(r1g-r2g)/(r1g-r3g)(r3g-I Sqrt[KKg])/(r2g-I Sqrt[KKg]),xaux,krg]]-wr/(\[Pi])Im[(r2g+I Sqrt[KKg])(r3g+I Sqrt[KKg])EllipticPi[(r1g-r2g)/(r1g-r3g)(r3g-I Sqrt[KKg])/(r2g-I Sqrt[KKg]),krg]])
]


\[Psi]zfun[wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,z2g,kzg,ellK,xaux},
 {EEg,Lzg,KKg}=const;
 r1g=p/(1-e);
 r2g=p/(1+e);
 z1g=Sqrt[1-xg^2];

 z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
 kzg=a^2(1-EEg^2)z1g^2/z2g^2;
 ellK=EllipticK[kzg];
 xaux=JacobiAmplitude[(2/\[Pi])*ellK(wz+\[Pi]/2),kzg];

 -(((a^2-KKg)EEg-a*Lzg)/(Sqrt[KKg]z2g))(EllipticPi[(a^2z1g^2)/KKg,xaux,kzg]-2/\[Pi] ellK(wz+\[Pi]/2)EllipticPi[a^2z1g^2/KKg,kzg])
]


(* ::Section::Closed:: *)
(*Functions to compute linearized angular Teukolsky equations*)


(* ::Subsection::Closed:: *)
(*SWSH eigenvalues and eigenfunctions*)


angpar[s_,l_,m_,a_,\[Omega]0_,\[Omega]1_]:=Module[{\[Lambda]0,\[Lambda]1,\[Gamma]0,\[Gamma]1,swshY0,swsh0all,swsh1all},
	 \[Gamma]0 = a \[Omega]0;
	 \[Gamma]1 = a \[Omega]1;
	If[a==0,
	 \[Lambda]0=l (1+l)-s (1+s);
	 swshY0 = Evaluate[SpinWeightedSphericalHarmonicY[s,l,m,#,0]]&;
	 swsh0all = Function[{\[Theta]},{swshY0[\[Theta]],swshY0'[\[Theta]],swshY0''[\[Theta]]}];
	 {\[Lambda]1,swsh1all} = {0,{0,0,0}&};
	 ,
	 {\[Lambda]0,\[Lambda]1,swsh0all,swsh1all}= angparKerr[s,l,m,\[Gamma]0,\[Gamma]1]//Quiet;
	];
	 {\[Lambda]0,\[Lambda]1,swsh0all,swsh1all}
]


angparLeading[s_,l_,m_,a_,\[Omega]0_]:=Module[{\[Lambda]0,\[Gamma]0,swshY0,swsh0all,swshcoeff,swshS},
	 \[Gamma]0 = a \[Omega]0;
	If[a==0,
	 \[Lambda]0=l (1+l)-s (1+s);
	 swshY0 = Evaluate[SpinWeightedSphericalHarmonicY[s,l,m,#,0]]&;
	 swsh0all = Function[{\[Theta]},{swshY0[\[Theta]],swshY0'[\[Theta]],swshY0''[\[Theta]]}];
	 ,
	 \[Lambda]0 = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0]//Quiet;
	 swshcoeff = SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0];
	 swshS = swshcoeff[#,0]&;
	 swsh0all = Function[{\[Theta]},{swshS[\[Theta]]//Quiet,swshS'[\[Theta]]//Quiet,swshS''[\[Theta]]//Quiet}];
	];
	 {\[Lambda]0,swsh0all}
]


angparKerr[s_,l_,m_,\[Gamma]0_,\[Gamma]1_]:=Module[{swshcoeff,\[Lambda]0,\[Lambda]1,swsh0all,swsh1all,swshS},
	 \[Lambda]0 = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0]//Quiet;
	 swshcoeff = SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0];
	 swshS = swshcoeff[#,0]&;
	 swsh0all = Function[{\[Theta]},{swshS[\[Theta]]//Quiet,swshS'[\[Theta]]//Quiet,swshS''[\[Theta]]//Quiet}];
	
	{\[Lambda]1,swsh1all} = angpar1corr[s,l,m,\[Gamma]0,\[Gamma]1,\[Lambda]0,swshcoeff[[5,2;;]]];
	{\[Lambda]0,\[Lambda]1,swsh0all,swsh1all}
]


(* ::Subsubsection::Closed:: *)
(*Kernel linear corrections SWSH and their eigenvalues*)


angpar1corr[s_,l_,m_,\[Gamma]0_,\[Gamma]1_,\[Lambda]0_,coeff_]:=Module[{km, kp, ks,an1, nmin,nmax,\[Alpha]n,an0Tab,an1Tab,an0square,an1square,an1an0prod,\[Beta]n0,\[Beta]n1,\[Gamma]n0,\[Gamma]n1,\[Lambda]1,swsh1all,normterm0,norm0,\[ScriptCapitalA],\[ScriptCapitalB],coeff1order,sign,sum0,sum1,nprec1,\[Theta]},
	 km = Abs[m-s]/2;
	 kp = Abs[m+s]/2;
	 ks = 2(kp+km);
	
	 an0Tab = coeff[[1]];
	 nmin = coeff[[2]];
	 nmax = coeff[[3]];
	
	 an0square = Table[Sum[an0Tab[[j]]an0Tab[[i+2-nmin-j]], {j, 1, i-nmin+1}],{i,nmin,nmax}];
	 normterm0[i_]:= 2^i Pochhammer[2+ks+i,-2kp-1] Hypergeometric1F1[1+2 km+i,2+ks+i,4 \[Gamma]0]an0square[[i+1]];
	 sum0 = Sum[normterm0[i],{i,nmin,nmax}];
	
	 (* Normalisation such that \[Integral]\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubsuperscriptBox[\(S\), \(lm\), \(*\)]\)\)(\[Theta],\[Phi];\[Gamma])\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubscriptBox[\(S\), \(l'm'\)]\)\)(\[Theta],\[Phi];\[Gamma])d\[CapitalOmega] = Subscript[\[Delta], ll']Subscript[\[Delta], mm'] *)
	 norm0=Sqrt[2\[Pi]](2^(1+ks)E^(-2 \[Gamma]0) Gamma[1+2 kp] sum0)^(1/2);
	
	 (* Overall sign such that the \[Gamma]\[Rule]0 limit is continuous *)
	 sign = If[(OddQ[l]&&EvenQ[m+s])||(OddQ[l]&&OddQ[m+s]&&m>=s)||(EvenQ[l]&&OddQ[m-s]&&m<=s), -1, 1];
	
	 coeff1order[n_]:= 2\[Pi](Gamma[1+2 km+n]/ Gamma[3+ks+n])2^(2+ks+n)E^(-2\[Gamma]0) Gamma[1+2kp] ((1+2kp) (2+ks+n+2 s) Hypergeometric1F1[1+2 km+n,3+ks+n,4 \[Gamma]0]-(2+ks+n)(1+2 kp-m+s)Hypergeometric1F1[1+2 km+n,2+ks+n,4 \[Gamma]0]);
	
	 (* First order correction to the eigenvalues: \[Lambda]1 = \[Gamma]1(Alm\[Gamma]1 + 2(\[Gamma]0 - m)) *)
	 \[Lambda]1 = -(1/norm0)^2Sum[coeff1order[i]an0square[[i+1]], {i,nmin, nmax}];
	
	(* The interpolation method introduces an exta \[Gamma]1 that is not needed in the computation of the eigenfunctions. This factor is reintroduced at the end *)
	 If[Precision[\[Lambda]1]<Max[Precision[\[Gamma]1],16],
	  \[Lambda]1 = (1/\[Gamma]1)angeigenvalue1corrint[s,l,m,\[Gamma]0,\[Gamma]1];
	 ];
	
	 (* Coefficients of first order expansion in \[Sigma] \[Gamma]1*);
	 \[Alpha]n[n_]:= -2(n+1)(n+2 km+1);
	 \[Beta]n0[n_]:=  n(n-1)+2n(km+kp+1-2\[Gamma]0)-(2\[Gamma]0(2km+s+1)-(km+kp)(km+kp+1))-(s(s+1)+\[Lambda]0+2m \[Gamma]0);
	 \[Beta]n1[n_]:= -2(1+2 km+m+2 n+s)-\[Lambda]1;
	 \[Gamma]n0[n_]:= 2\[Gamma]0(n+kp+km+s);
	 \[Gamma]n1[n_]:= (ks+2 (n+s));
	
	 an1[nmin] = 0;
	 an1[nmin+1] = -(\[Beta]n1[nmin]/\[Alpha]n[nmin])an0Tab[[nmin+1]];
	 an1[n_]:= an1[n] =(-an1[n-1] \[Beta]n0[n-1]-an0Tab[[n]] \[Beta]n1[n-1]-an1[n-2] \[Gamma]n0[n-1]-an0Tab[[n-1]] \[Gamma]n1[n-1])/\[Alpha]n[n-1];     
	 an1Tab = Table[an1[i],{i,nmin,nmax}];
	
	(*This algorithm is taken from "Convergence acceleration for the numerical solution of second-order linear recurrence relations"
	   by P. Levrie and A. Bultheel *)
	 Block[{rs1,rs2,fs},
	 (*It is possible to reuse the rs[n] of the zeroth order correction *)
	  rs1[n_]:= rs1[n]= -(\[Gamma]n0[n+1]/\[Alpha]n[n+1])/((\[Beta]n0[n+1]/\[Alpha]n[n+1])+rs1[n+1]);
	  rs1[nmax] = 0;
	  Table[rs1[j],{j,nmax,1,-1}];
	  rs2[n_]:= rs2[n]= -((an0Tab[[n+2]] \[Beta]n1[n+1]+an0Tab[[n+1]] \[Gamma]n1[n+1])/\[Alpha]n[n+1]+rs2[n+1])/((\[Beta]n0[n+1]/\[Alpha]n[n+1])+rs1[n+1]);
	  rs2[nmax] = 0;
	  Table[rs2[j],{j,nmax,1,-1}];
	  (*Number of high precision elements obtained by solving the three-terms recurrence relation directly. These terms are the seeds for the backward algorithm*)
	  fs[n_] := fs[n] = rs1[n-1]fs[n-1]+rs2[n-1];
	  nprec1 = Length[Position[an1Tab,_?(Precision[#]>(Precision[\[Gamma]0]-5)&)]];
	  Table[fs[i] = an1[i],{i,nmin,nprec1}];
	  an1Tab = Table[fs[j],{j,nmin,nmax}];
	 ];
	
	 an1an0prod = Table[Sum[an1Tab[[j]]an0Tab[[i+2-nmin-j]], {j, 1, i-nmin+1}],{i,nmin,nmax}];
	
	 \[ScriptCapitalA][i_]:= 2^i Pochhammer[2+i+ks,-1-2 kp](4 (1+i+2 km)/(2+i+ks))Hypergeometric1F1[2+i+2 km,3+i+ks,4 \[Gamma]0];
	 \[ScriptCapitalB][i_]:= 2^i Pochhammer[2+i+ks,-1-2 kp]Hypergeometric1F1[1+i+2 km,2+i+ks,4 \[Gamma]0];
	 
	 sum1 = Sum[\[ScriptCapitalB][i]2an1an0prod[[i+1]]+\[ScriptCapitalA][i] an0square[[i+1]],{i,nmin,nmax}];
	
	swsh1all=Function[{\[Theta]},Module[{u,oneplusu,oneminusu,\[Theta]var1,swsh1,dswsh1,ddswsh1},
	 u = Cos[\[Theta]var1];
	 oneplusu = 2Cos[\[Theta]var1/2]^2;
	 oneminusu = 2Sin[\[Theta]var1/2]^2;
	
	 (*Compensated summation does not improve the precision of the corrections to the eigenfunctions*)
	swsh1 = \[Gamma]1(sign/norm0)E^(\[Gamma]0 u)If[km==0,1,oneplusu^km] If[kp==0,1,oneminusu^kp] ((-sum1/(2sum0)+oneplusu)an0Tab+an1Tab) . oneplusu^Range[0,nmax]/.\[Theta]var1->\[Theta];
	dswsh1 = \[Gamma]1(sign/norm0)D[E^(\[Gamma]0 u)If[km==0,1,oneplusu^km] If[kp==0,1,oneminusu^kp]((-sum1/(2sum0)+oneplusu) an0Tab+an1Tab) . oneplusu^Range[0,nmax],\[Theta]var1]/.\[Theta]var1->\[Theta];ddswsh1=\[Gamma]1(sign/norm0)  D[E^(\[Gamma]0 u)If[km==0,1,oneplusu^km] If[kp==0,1,oneminusu^kp]((-sum1/(2sum0)+oneplusu)an0Tab+an1Tab) . oneplusu^Range[0,nmax],{\[Theta]var1,2}]/.\[Theta]var1->\[Theta];
	(*(*Quite demanding, it can be adjusted*)
	If[(Precision[swsh1]<Min[Precision[\[Gamma]1],precODE])||(Precision[dswsh1]<Min[Precision[\[Gamma]1],precODE])||(swsh1===Indeterminate)||(dswsh1===Indeterminate)||(ddswsh1===Indeterminate),
	  {swsh1,dswsh1,ddswsh1} = angeigenfunction1corrint[s,l,m,\[Gamma]0,\[Gamma]1,\[Theta]]
	 ];*)
	
	(*Quite demanding, it can be adjusted*)
	If[(Precision[swsh1]<20)||(Precision[dswsh1]<20)||(swsh1===Indeterminate)||(dswsh1===Indeterminate)||(ddswsh1===Indeterminate),
	  {swsh1,dswsh1,ddswsh1} = angeigenfunction1corrint[s,l,m,\[Gamma]0,\[Gamma]1,\[Theta]]
	 ];
	{swsh1,dswsh1,ddswsh1}]
	];
	
	 Remove[\[Alpha]n,\[Beta]n0,\[Beta]n1,\[Gamma]n0,\[Gamma]n1];
	{\[Gamma]1 \[Lambda]1,swsh1all}
]


(* ::Subsubsection::Closed:: *)
(*Corrections to eigenvalues and eigenfunctions through numerical derivative*)


angeigenvalue1corrint[s_,l_,m_,\[Gamma]0_,\[Gamma]1_]:=Module[{qratio,\[Lambda]p,\[Lambda]m},
  qratio = 1 10^-6;
  \[Lambda]p = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0+qratio \[Gamma]1];
  \[Lambda]m = SpinWeightedSpheroidalEigenvalue[s,l,m,\[Gamma]0-qratio \[Gamma]1];
  (\[Lambda]p-\[Lambda]m)/(2qratio)
]


angeigenfunction1corrint[s_,l_,m_,\[Gamma]0_,\[Gamma]1_,\[Theta]_]:=Module[{qratio,\[Lambda]p,\[Lambda]m,swshcoeffp,swshcoeffm,swshSp,swshSm,swshp,swshm,dswshp,dswshm,ddswshp,ddswshm,swsh1,dswsh1,ddswsh1},
   qratio = 1 10^-6;
   swshcoeffp = SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0+qratio \[Gamma]1];
   swshcoeffm = SpinWeightedSpheroidalHarmonicS[s,l,m,\[Gamma]0-qratio \[Gamma]1];
   swshSp = swshcoeffp[#,0]&;
   swshSm = swshcoeffm[#,0]&;

   swshp = swshSp[\[Theta]]//Quiet;
   dswshp = swshSp'[\[Theta]]//Quiet;
   ddswshp = swshSp''[\[Theta]]//Quiet;
   swshm = swshSm[\[Theta]]//Quiet;
   dswshm = swshSm'[\[Theta]]//Quiet;
   ddswshm = swshSm''[\[Theta]]//Quiet;

  {(swshp-swshm)/(2qratio),(dswshp-dswshm)/(2qratio),(ddswshp-ddswshm)/(2qratio)}
]


(* ::Section::Closed:: *)
(*Functions to compute linearized radial Teukolsky equations*)


(* ::Subsection::Closed:: *)
(*Coefficients Teukolsky equation in hyperboloidal slicing - linear in the secondary spin*)


TeukolskyHSCoeff1spin := 
 Module[{M=1,r,m ,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1,\[CapitalDelta],s,f,H,Gtilda0,Gtilda1,Utilda0,Utilda1,p0,p1,q0,q1, coefficients0,coefficients1},
	  \[CapitalDelta] = r^2-2M r+a^2;
	  f = \[CapitalDelta]/(r^2+a^2); (*dr/drstar*)
	  Gtilda0 = a^2 \[CapitalDelta]+(r^2+a^2)(r s (r-M)-I r((r^2+a^2)\[Omega]0 H+m a));
	  Gtilda1 = -I H r (r^2+a^2)^2 \[Omega]1;
	  Utilda0 = \[CapitalDelta] (2 a^2-2 M r (1+s)-r^2 \[Lambda]0)-2 a (1+H) m r^2 (r^2+a^2) \[Omega]0+2 I r^2 s (-(1+H) M (-a^2+r^2)+(1-H) r \[CapitalDelta]) \[Omega]0-(-1+H^2) r^2 (r^2+a^2)^2 \[Omega]0^2-2 I a r \[CapitalDelta] (m+a H \[Omega]0);
	  Utilda1 = r (-r \[CapitalDelta] \[Lambda]1-2 a (1+H) m r (r^2+a^2) \[Omega]1-2 I a^2 H \[CapitalDelta] \[Omega]1+2 I r s ((1+H) M (a-r) (a+r)-(-1+H) r \[CapitalDelta]) \[Omega]1-2 (-1+H^2) r (r^2+a^2)^2 \[Omega]0 \[Omega]1);
	 
	  p0 = ((r^2+a^2)/\[CapitalDelta])D[f,r]-(1/(\[CapitalDelta](r^2+a^2)))(1/r)2Gtilda0;
	  q0 = Utilda0/(r^2 \[CapitalDelta]^2);
	 
	  p1 = -(1/(\[CapitalDelta](r^2+a^2)))(1/r)2Gtilda1;
	  q1 = Utilda1/(r^2 \[CapitalDelta]^2);
	 
	  coefficients0 = {Function@@{p0}/.Thread[{r,H,s,m,a,\[Omega]0,\[Lambda]0}->Array[Slot,7]], Function@@{q0}/.Thread[{r,H,s,m,a,\[Omega]0,\[Lambda]0}->Array[Slot,7]]};
	  coefficients1 = {Function@@{p1}/.Thread[{r,H,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1}->Array[Slot,9]], Function@@{q1}/.Thread[{r,H,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1}->Array[Slot,9]]};
	
	  Remove[r,H,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1];
	  {coefficients0[[1]],coefficients1[[1]],coefficients0[[2]],coefficients1[[2]]}
]


(* ::Subsection::Closed:: *)
(*Boundary conditions*)


(* ::Subsubsection::Closed:: *)
(*BCs at the horizon*)


bchorHS0spin[workingprecision_,s_,m_,a_,\[Omega]0_,\[Lambda]0_]:=Module[{M=1,deltarp,An0,rp,rm,rin,err,errold,i,chor0,p,q,phor0,qhor0,dphor0,dqhor0,\[Psi]hor0,an0,cInHor0},
	rp=M+Sqrt[M^2-a^2];
	rm=M-Sqrt[M^2-a^2];
	deltarp=(rp-rm)/50;
	rin=rp+deltarp; (*Closest r to the horizon *)
	chor0=2I ( 2M rp)/(rp-rm) (\[Omega]0-(a m)/(2 M rp))+ s;
	
	dqhor0[n_]:=Which[
		n==0,
		0,
		n==1,
		(2 I a m+2 M (-1+s)-2 I a^2 \[Omega]0+rp (2+\[Lambda]0-4 I rp s \[Omega]0))/((rm-rp) rp),
		n>1,
		2 (-1+n) (-rp)^(-n)+(rm-rp)^(-n)(2+\[Lambda]0-4 I rm s \[Omega]0)+1/rp 2 n (rm-rp)^-n (M (-1+s)+I a (m-a \[Omega]0)) Hypergeometric2F1[1,1-n,2,rm/rp]
	];
		
	dphor0[n_]:=Which[
		n==0,
		1-chor0,
		n==1,
		 1/(rp (rm-rp)^2) (-2rm^2+a^2(3+2s+4I*rp*\[Omega]0)+I*rp(-2a*m+2a^2\[Omega]0+I(rp+2M*s+2I*rp^2*\[Omega]0))),
		n>1,
		2 (-rp)^(-n)-(rm-rp)^(-n)+(rm-rp)^(-1-n)(rm (+2s)+2I*rm^2\[Omega]0+2I(-a*m+I*M*s+a^2\[Omega]0))
	];
		
	an0[0]=1;
	An0[n_]:= -(1/(n(n-chor0)))Sum[(j*dphor0[n-j]+dqhor0[n-j])an0[j],{j,0,n-1}];
	
	err=1;
	i=1;
	{p,q}=TeukolskyHSCoeff1spin[[{1,3}]];
	phor0=p[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
	qhor0=q[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
	
	While[err > 10^(-workingprecision),
	(*Apparently,it works better when comparing to the MST method by adding more terms instead of decreasing the starting radius. I am not sure why.*)
		If[Mod[i,30]==0,
			deltarp=deltarp/2;
			rin=rp+deltarp;
			phor0=p[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
			qhor0=q[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
		];
	    an0[i]=An0[i];
		\[Psi]hor0=Evaluate[1+Sum[an0[k](#-rp)^k,{k,i}]]&;
		err=Abs[\[Psi]hor0''[rin]+phor0 \[Psi]hor0'[rin]+qhor0 \[Psi]hor0[rin]];
	   
	   If[errold<=err&&errold> 10^(-workingprecision)&&i>5,
			Break[];
			,
			errold=err;
		];
		i++;
		
		If[i > 100, Break[]]  (*Safeguard to avoid ruwaway computation*)  
	];
	
	cInHor0=Table[an0[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
	
	Remove[chor0,M,dqhor0,dphor0,an0];
	{cInHor0,rin}
]


(* ::Subsubsection::Closed:: *)
(*BCs  at infinity*)


bcinfHS0spin[workingprecision_,s_,m_,a_,\[Omega]0_,\[Lambda]0_]:=Module[{M=1,rp,rm,rout,err,i,p,q,pinf0,qinf0,dpinf0,dqinf0,\[Psi]inf0,bn0,Bn0,cOutinf0},
	rp=M+Sqrt[M^2-a^2];
	rm=M-Sqrt[M^2-a^2];
	rout=10\[Pi](1/Abs[\[Omega]0]+Abs[\[Omega]0]/(1+Abs[\[Omega]0]));
	
	dqinf0[n_]:=Which[
		n==0,
		0,
		n==1,
		0,
		n==2,
		-(4 a m \[Omega]0+4 I M s \[Omega]0+\[Lambda]0),
		n>2,
		1/(rm-rp)^3 ((rm-rp)^2(2rm^(-2+n)rp-2rm rp^(-2+n)+2I*a*m(-rm^(-2+n)+rp^(-2+n))+2M(-rm^(-2+n)+rp^(-2+n))(1+s)+(-rm^(-1+n)+rp^(-1+n))\[Lambda]0)+2I (rm-rp)^2(-rm^(-1+n) rp+rm rp^(-1+n))\[Omega]0+4(a*m(-rp^n(2M+n rm-n rp)+rm^n(2M-n rm+n rp)+a^2(rp^(-2+n)(rm-n rm+(-3+n)rp)+rm^(-2+n)(-(-3+n)rm+(-1+n)rp)))+I*M(-rp^n(2M+n rm-n rp)+rm^n(2M-n rm+n rp)+a^2(rp^(-2+n)((-1+n) rm-(-3+n)rp)+rm^(-2+n)((-3+n) rm+rp-n rp)))s)\[Omega]0)
	];
		
	dpinf0[n_]:=Which[
		n==0,
		2I*\[Omega]0,
		n==1,
		-2s+2I*2M*\[Omega]0 ,
		n>1,
		rm^(-1+n)+rp^(-1+n)+(2 rm^(-1+n)((M-rm)s+I (a m+(a^2+rm^2)\[Omega]0)))/(rm-rp)-(2 rp^(-1+n)((M-rp)s+I (a*m+(a^2+rp^2) \[Omega]0)))/(rm-rp)
	];
		
	err=1;
	i=1;
	bn0[0]=1;
	Bn0[n_]:=(n-1)/(2I*\[Omega]0 ) bn0[n-1]+1/(2I*\[Omega]0*n) Sum[(dqinf0[j+1]-(n-j)dpinf0[j])bn0[n-j],{j,1,n}];
	
	{p,q}=TeukolskyHSCoeff1spin[[{1,3}]];
	pinf0=p[rout,1,s,m,a,\[Omega]0,\[Lambda]0];
	qinf0=q[rout,1,s,m,a,\[Omega]0,\[Lambda]0];
	
	While[err > 10^(-workingprecision),
		(*Apparently, it works better when comparing to the MST method by adding more terms instead of increasing the starting radius. Probably because the longest is the integration interval, the greater are the errors.*)
		(*If[Mod[i,50]\[Equal]0,
		rout=2*rout;
		pinf=p[rout,1,s,a,\[Omega],m,\[Lambda]];
		qinf=q[rout,1,s,a,\[Omega],m,\[Lambda]];
		];
		*)
	     bn0[i]=Bn0[i];
		\[Psi]inf0=Evaluate[1+Sum[bn0[k](#)^-k,{k,i}]]&;
		err=Abs[\[Psi]inf0''[rout]+pinf0 \[Psi]inf0'[rout]+qinf0 \[Psi]inf0[rout]];
	    i++;
		If[i > 100, Break[]]  (*Asymptotic expansions are not convergent*)    
	];
	cOutinf0=Table[bn0[k],{k,0,i-1}]; (*Coefficients for outgoing waves at \[Infinity] (+)*)
		
	Remove[M,dqinf0,dpinf0,bn0];
	
	{cOutinf0,rout}
]


(* ::Subsubsection::Closed:: *)
(*Linearized BCs in the spin at the horizon*)


bchorHS1spin[workingprecision_,s_,m_,a_,\[Omega]0_,\[Omega]1_,\[Lambda]0_,\[Lambda]1_]:=Module[{M=1,deltarp,An0,rp,rm,rin,err,i,chor0,chor1,p,q,phor0,qhor0,dphor0,dphor1,dqhor0,dqhor1,\[Psi]hor0,an0,an1,cInHor0,cInHor1},
	rp=M+Sqrt[M^2-a^2];
	rm=M-Sqrt[M^2-a^2];
	deltarp=(rp-rm)/50;
	rin=rp+deltarp; (*Closest r to the horizon *)
	chor0=2I ( 2M rp)/(rp-rm) (\[Omega]0-(a m)/(2 M rp))+ s;
	chor1=(4 I M rp \[Omega]1)/(rp-rm);
	
	dqhor0[n_]:=Which[
		n==0,
		0,
		n==1,
		(2 I a m+2 M (-1+s)-2 I a^2 \[Omega]0+rp (2+\[Lambda]0-4 I rp s \[Omega]0))/((rm-rp) rp),
		n>1,
		2 (-1+n) (-rp)^(-n)+(rm-rp)^(-n)(2+\[Lambda]0-4 I rm s \[Omega]0)+1/rp 2 n (rm-rp)^-n (M (-1+s)+I a (m-a \[Omega]0)) Hypergeometric2F1[1,1-n,2,rm/rp]
	];
	
	dqhor1[n_]:=Which[
		n==0,
		0,
		n==1,
		(-2 I a^2 \[Omega]1+rp (\[Lambda]1-4 I rp s \[Omega]1))/((rm-rp) rp),
		n>1,
		1/rp (rm-rp)^(-n)(rp*\[Lambda]1-4I*a^2*s*\[Omega]1-2I*a^2*n*\[Omega]1 Hypergeometric2F1[1,1-n,2,rm/rp])
	];
	
	dphor0[n_]:=Which[
		n==0,
		1-chor0,
		n==1,
		 1/(rp (rm-rp)^2) (-2rm^2+a^2(3+2s+4I*rp*\[Omega]0)+I*rp(-2a*m+2a^2\[Omega]0+I(rp+2M*s+2I*rp^2*\[Omega]0))),
		n>1,
		2 (-rp)^(-n)-(rm-rp)^(-n)+(rm-rp)^(-1-n)(rm (+2s)+2I*rm^2\[Omega]0+2I(-a*m+I*M*s+a^2\[Omega]0))
	];
	
	dphor1[n_]:=Which[
		n==0,
		-chor1,
		n==1,
		 -((2 I (-3a^2+rp^2)\[Omega]1)/(rm-rp)^2),
		n>1,
		2 I (a^2+rm^2) (rm-rp)^(-1-n)\[Omega]1
	];
	
	an0[0]=1;
	An0[n_]:= -(1/(n(n-chor0)))Sum[(j*dphor0[n-j]+dqhor0[n-j])an0[j],{j,0,n-1}];
	
	err=1;
	i=1;
	{p,q}=TeukolskyHSCoeff1spin[[{1,3}]] ;
	phor0=p[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
	qhor0=q[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
	
	While[err > 10^(-workingprecision),
	(*Apparently,it works better when comparing to the MST method by adding more terms instead of decreasing the starting radius. I am not sure why.*)
		If[Mod[i,30]==0,
			deltarp=deltarp/2;
			rin=rp+deltarp;
			phor0=p[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
			qhor0=q[rin,-1,s,m,a,\[Omega]0,\[Lambda]0];
		];
	        an0[i]=An0[i];
		\[Psi]hor0=Evaluate[1+Sum[an0[k](#-rp)^k,{k,i}]]&;
		err=Abs[\[Psi]hor0''[rin]+phor0 \[Psi]hor0'[rin]+qhor0 \[Psi]hor0[rin]];
	    i++;
	];
	
	an1[0]=0;
	an1[n_]:=an1[n]=-(1/(n(n-chor0)))Sum[(chor1 /(n-chor0) an0[j]+ an1[j]) (j dphor0[n-j]+dqhor0[n-j])+ (j dphor1[-j+n]+dqhor1[-j+n])an0[j],{j,0,n-1}];
	
	cInHor0=Table[an0[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
	cInHor1=Table[an1[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
	
	Remove[chor0,chor1,M,dqhor0,dphor0,dqhor1,dphor1,an0];
	{cInHor0,cInHor1,rin}
]


(* ::Subsubsection::Closed:: *)
(*Linearized BCs in the spin at infinity*)


bcinfHS1spin[workingprecision_,s_,m_,a_,\[Omega]0_,\[Omega]1_,\[Lambda]0_,\[Lambda]1_]:=Module[{M=1,rp,rm,rout,err,i,p,q,pinf0,qinf0,dpinf0,dpinf1,dqinf0,dqinf1,\[Psi]inf0,bn0,bn1,Bn0,cOutinf0,cOutinf1},
	rp=M+Sqrt[M^2-a^2];
	rm=M-Sqrt[M^2-a^2];
	rout=10\[Pi](1/Abs[\[Omega]0]+Abs[\[Omega]0]/(1+Abs[\[Omega]0]));
	
	dqinf0[n_]:=Which[
		n==0,
		0,
		n==1,
		0,
		n==2,
		-(4 a m \[Omega]0+4 I M s \[Omega]0+\[Lambda]0),
		n>2,
		1/(rm-rp)^3 ((rm-rp)^2(2rm^(-2+n)rp-2rm rp^(-2+n)+2I*a*m(-rm^(-2+n)+rp^(-2+n))+2M(-rm^(-2+n)+rp^(-2+n))(1+s)+(-rm^(-1+n)+rp^(-1+n))\[Lambda]0)+2I (rm-rp)^2(-rm^(-1+n) rp+rm rp^(-1+n))\[Omega]0+4(a*m(-rp^n(2M+n rm-n rp)+rm^n(2M-n rm+n rp)+a^2(rp^(-2+n)(rm-n rm+(-3+n)rp)+rm^(-2+n)(-(-3+n)rm+(-1+n)rp)))+I*M(-rp^n(2M+n rm-n rp)+rm^n(2M-n rm+n rp)+a^2(rp^(-2+n)((-1+n) rm-(-3+n)rp)+rm^(-2+n)((-3+n) rm+rp-n rp)))s)\[Omega]0)
	];
	
	dqinf1[n_]:=Which[
		n==0,
		0,
		n==1,
		0,
		n==2,
		-\[Lambda]1-4(a*m+I*M*s)\[Omega]1,
		n>2,
		1/(rm-rp)^3 ((rm-rp)^2(-rm^(-1+n)+rp^(-1+n))\[Lambda]1+2I (rm-rp)^2(-rm^(-1+n)rp+rm rp^(-1+n))\[Omega]1+4a*m(rp^n(-2M+n(-rm+rp))+rm^n(2 M+n (-rm+rp))+a^2(rp^(-2+n)(rm-n rm+(-3+n)rp)+rm^(-2+n)(-(-3+n)rm+(-1+n)rp)))\[Omega]1+4I*M(rp^n(-2M+n (-rm+rp))+rm^n(2M+n (-rm+rp))+a^2(rp^(-2+n)((-1+n) rm-(-3+n) rp)+rm^(-2+n)((-3+n)rm+rp-n*rp)))s*\[Omega]1)
	];
	
	dpinf0[n_]:=Which[
		n==0,
		2I*\[Omega]0,
		n==1,
		-2s+2I*2M*\[Omega]0 ,
		n>1,
		rm^(-1+n)+rp^(-1+n)+(2 rm^(-1+n)((M-rm)s+I (a m+(a^2+rm^2)\[Omega]0)))/(rm-rp)-(2 rp^(-1+n)((M-rp)s+I (a*m+(a^2+rp^2) \[Omega]0)))/(rm-rp)
	];
	
	dpinf1[n_]:=Which[
		n==0,
		2I*\[Omega]1,
		n==1,
		4I*M*\[Omega]1,
		n>1,
		(4I*M(rm^n -rp^n)\[Omega]1)/(rm-rp)
	];
	
	err=1;
	i=1;
	bn0[0]=1;
	Bn0[n_]:=(n-1)/(2I*\[Omega]0 ) bn0[n-1]+1/(2I*\[Omega]0*n) Sum[(dqinf0[j+1]-(n-j)dpinf0[j])bn0[n-j],{j,1,n}];
	
	{p,q}=TeukolskyHSCoeff1spin[[{1,3}]];
	pinf0=p[rout,1,s,m,a,\[Omega]0,\[Lambda]0];
	qinf0=q[rout,1,s,m,a,\[Omega]0,\[Lambda]0];
	
	While[err > 10^(-workingprecision),
		(*Apparently, it works better when comparing to the MST method by adding more terms instead of increasing the starting radius. Probably because the longest is the integration interval, the greater are the errors.*)
		(*If[Mod[i,50]\[Equal]0,
		rout=2*rout;
		pinf=p[rout,1,s,a,\[Omega],m,\[Lambda]];
		qinf=q[rout,1,s,a,\[Omega],m,\[Lambda]];
		];
		*)
	     bn0[i]=Bn0[i];
		\[Psi]inf0=Evaluate[1+Sum[bn0[k](#)^-k,{k,i}]]&;
		err=Abs[\[Psi]inf0''[rout]+pinf0 \[Psi]inf0'[rout]+qinf0 \[Psi]inf0[rout]];
	    i++;
		If[i > 100, Break[]]  (*Asymptotic expansions are not convergent*)    
	];
	cOutinf0=Table[bn0[k],{k,0,i-1}]; (*Coefficients for outgoing waves at \[Infinity] (+)*)
	
	bn1[0]=0;
	
	bn1[n_]:=bn1[n]= (n-1) /(2I \[Omega]0) (bn1[n-1]-\[Omega]1/\[Omega]0 bn0[n-1])+1/(2I*\[Omega]0*n) Sum[(bn1[n-j]-\[Omega]1/\[Omega]0 bn0[n-j]) (dqinf0[j+1]-(n-j) dpinf0[j])+(dqinf1[j+1]-(n-j) dpinf1[j])bn0[n-j],{j,1,n}];
	
	cOutinf1=Table[bn1[k],{k,0,i-1}];
	
	Remove[M,dqinf0,dqinf1,dpinf0,dpinf1,bn0];
	
	{cOutinf0,cOutinf1,rout}
]


(* ::Subsection::Closed:: *)
(*Teukolsky solver in hyperboloidal slicing coordinates - linear in the secondary spin and specialized to eccentric orbits*)


(* ::Text:: *)
(*Remember to specify "InterpolationOrder->All"! It is needed when you are not evaluating your solution at one the extreme of the integrating interval, but in one of the intermediate points.*)


(* ::Text:: *)
(*TeukolskyHS2 is specialized for eccentric orbits, and it solve the In (Up) ODE up to (down to) the apostraton (periastron). This code provide the solution to the HST version of the radial Teukolsky equation. *)


erre[ecc_,semilatus_,\[Chi]_]:= semilatus/(1 + ecc Cos[\[Chi]]);


(* ::Subsubsection::Closed:: *)
(*Leading order term only*)


TeukolskySolverHS0spinecc[s_,l_,m_,a_,\[Omega]0_,\[Lambda]0_,ecc_,semilatus_]:=Module[{M=1,prec,r,rin,rout,rp,rm,p0,q0,cInH0,cOutinf0,\[Psi]hor0,\[Psi]inf0,eqhor,eqinf,X0,Y0,\[Psi]in0,d\[Psi]in0,\[Psi]up0,d\[Psi]up0,nmaxhor,nmaxinf,r0in,r0up},
   rp = M+Sqrt[M^2-a^2];
   rm = M-Sqrt[M^2-a^2];
 
   prec = Max[Precision[\[Omega]0]-8,MachinePrecision];
   
   {p0,q0}=TeukolskyHSCoeff1spin[[{1,3}]];
   
   r0in=erre[ecc,semilatus,-\[Pi]];
   r0up=erre[ecc,semilatus,0];
 
  {cInH0,rin} = bchorHS0spin[prec,s,m,a,\[Omega]0,\[Lambda]0];
  {cOutinf0,rout} = bcinfHS0spin[prec,s,m,a,\[Omega]0,\[Lambda]0];
  nmaxhor = Length[cInH0];
  nmaxinf = Length[cOutinf0];
  \[Psi]hor0 = Evaluate[Sum[cInH0[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;
  \[Psi]inf0[r_]:= Sum[cOutinf0[[i]]r^(-i+1),{i,nmaxinf}];
 
  eqhor={
   X0'[r] == Y0[r],
   Y0'[r] == -p0[r,-1,s,m,a,\[Omega]0,\[Lambda]0]Y0[r]-q0[r,-1,s,m,a,\[Omega]0,\[Lambda]0] X0[r],
   X0[rin] == \[Psi]hor0[rin],Y0[rin] == \[Psi]hor0'[rin]
  };
  
  eqinf={
   X0'[r] == Y0[r],
   Y0'[r] == -p0[r,1,s,m,a,\[Omega]0,\[Lambda]0]Y0[r]-q0[r,1,s,m,a,\[Omega]0,\[Lambda]0] X0[r],
   X0[rout] == \[Psi]inf0[rout],Y0[rout] == \[Psi]inf0'[rout]
 };  
  {\[Psi]in0,d\[Psi]in0} = {X0,Y0}/.First@NDSolve[eqhor,{X0,Y0},{r,rin,r0in},Method->"StiffnessSwitching",WorkingPrecision->prec,InterpolationOrder->All];
  {\[Psi]up0,d\[Psi]up0} = {X0,Y0}/.First@NDSolve[eqinf,{X0,Y0},{r,r0up,rout},Method->"StiffnessSwitching",WorkingPrecision->prec,InterpolationOrder->All];
  
   (* Remove local variables not garbage collected*)
  Remove[X0,Y0,r];
  ClearSystemCache[];

{{\[Psi]in0,d\[Psi]in0},{\[Psi]up0,d\[Psi]up0}}
]


TeukolskyRadialHSLeading[r_,s_,l_,m_,a_,\[Omega]0_,\[Psi]_]:=Module[{M=1,rp,rm,rs,\[Psi]in0,d\[Psi]in0,\[Psi]up0,d\[Psi]up0,\[CapitalDelta],dfacexp0,Rin0,dRin0,Rup0,dRup0,resfac},
	{{\[Psi]in0,d\[Psi]in0},{\[Psi]up0,d\[Psi]up0}}=\[Psi];
	rp = M+Sqrt[M^2-a^2];
	rm = M-Sqrt[M^2-a^2];
	rs=((2M rp)/(rp-rm) Log[(r-rp)/(2M)]-(2M rm)/(rp-rm) Log[(r-rm)/(2M)]+r);

	 \[CapitalDelta] = r^2+a^2-2 M r;
	
	dfacexp0 = -(1/r +2 s (r-M)/\[CapitalDelta])+(I/\[CapitalDelta])(# (r^2+a^2)\[Omega]0+a m)&;
	
	Rin0 = rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])] \[Psi]in0[r];
	dRin0 = rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])](\[Psi]in0[r] dfacexp0[-1]+ d\[Psi]in0[r]);
	
	Rup0 = \[Psi]up0[r];
	dRup0 = \[Psi]up0[r] dfacexp0[1]+d\[Psi]up0[r];
	
	resfac[H_]:=r^(-1) \[CapitalDelta]^(-s)Exp[H I \[Omega]0 rs]Exp[I m a/(rp-rm) Log[(r-rp)/(r-rm)]];
	
	{resfac[-1]{Rin0,dRin0},resfac[1]{Rup0,dRup0}}
]


(* ::Text:: *)
(*Function which returns radial function and its correction*)


RLeading[s_,l_,m_,a_,\[Omega]_,\[Lambda]_,e_,p_]:=Module[{\[Psi]},
	\[Psi]=TeukolskySolverHS0spinecc[s,l,m,a,\[Omega],\[Lambda],e,p];
	Function[{r},
		Module[{R},
				R=TeukolskyRadialHSLeading[r,s,l,m,a,\[Omega],\[Psi]];
				Association[
							"R0"-><|
							"In"-><|"R"->R[[1,1]],"dR"->R[[1,2]]|>,
							"Up"-><|"R"->R[[2,1]],"dR"->R[[2,2]]|>
								|>
						]
				]
	]
]


(* ::Subsubsection::Closed:: *)
(*Leading order and spin correction terms*)


TeukolskySolverHS1spinecc[s_,l_,m_,a_,\[Omega]0_,\[Omega]1_,\[Lambda]0_,\[Lambda]1_,ecc_,semilatus_]:=Module[{M=1,prec,r,rin,rout,rp,rm,p0,p1,q0,q1,cInH0,cInH1,cOutinf0,cOutinf1,\[Psi]hor0,\[Psi]hor1,\[Psi]inf0,\[Psi]inf1,eqhor,eqinf,X0,X1,Y0,Y1,\[Psi]in0,d\[Psi]in0,\[Psi]in1,d\[Psi]in1,\[Psi]up0,d\[Psi]up0,\[Psi]up1,d\[Psi]up1,nmaxhor,nmaxinf,r0in,r0up},
   rp = M+Sqrt[M^2-a^2];
   rm = M-Sqrt[M^2-a^2];
 
   prec = Max[Precision[\[Omega]0]-8,MachinePrecision];
   
   {p0,p1,q0,q1}=TeukolskyHSCoeff1spin;
   r0in=erre[ecc,semilatus,-\[Pi]];
   r0up=erre[ecc,semilatus,0];
 
  {cInH0,cInH1,rin} = bchorHS1spin[prec,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1];
  {cOutinf0,cOutinf1,rout} = bcinfHS1spin[prec,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1];
  nmaxhor = Length[cInH0];
  nmaxinf = Length[cOutinf0];
  \[Psi]hor0 = Evaluate[Sum[cInH0[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;
  \[Psi]hor1 = Evaluate[Sum[cInH1[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;
  \[Psi]inf0[r_]:= Sum[cOutinf0[[i]]r^(-i+1),{i,nmaxinf}];
  \[Psi]inf1[r_]:= Sum[cOutinf1[[i]]r^(-i+1),{i,nmaxinf}];
 
  eqhor={
   X0'[r] == Y0[r],
   Y0'[r] == -p0[r,-1,s,m,a,\[Omega]0,\[Lambda]0]Y0[r]-q0[r,-1,s,m,a,\[Omega]0,\[Lambda]0] X0[r],
   X1'[r] == Y1[r],
   Y1'[r] == -(p0[r,-1,s,m,a,\[Omega]0,\[Lambda]0]Y1[r]+q0[r,-1,s,m,a,\[Omega]0,\[Lambda]0] X1[r])-(p1[r,-1,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1]Y0[r]+q1[r,-1,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1] X0[r]),
   X0[rin] == \[Psi]hor0[rin],Y0[rin] == \[Psi]hor0'[rin],X1[rin] == \[Psi]hor1[rin],Y1[rin] == \[Psi]hor1'[rin]
  };
  
  eqinf={
   X0'[r] == Y0[r],
   Y0'[r] == -p0[r,1,s,m,a,\[Omega]0,\[Lambda]0]Y0[r]-q0[r,1,s,m,a,\[Omega]0,\[Lambda]0] X0[r],
   X1'[r] == Y1[r],
   Y1'[r] == -(p0[r,1,s,m,a,\[Omega]0,\[Lambda]0] Y1[r]+q0[r,1,s,m,a,\[Omega]0,\[Lambda]0] X1[r])-(p1[r,1,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1] Y0[r]+q1[r,1,s,m,a,\[Omega]0,\[Omega]1,\[Lambda]0,\[Lambda]1] X0[r]),
   X0[rout] == \[Psi]inf0[rout],Y0[rout] == \[Psi]inf0'[rout],X1[rout] == \[Psi]inf1[rout],Y1[rout] == \[Psi]inf1'[rout]
 };  
  {\[Psi]in0,d\[Psi]in0,\[Psi]in1,d\[Psi]in1} = {X0,Y0,X1,Y1}/.First@NDSolve[eqhor,{X0,Y0,X1,Y1},{r,rin,r0in},Method->"StiffnessSwitching",WorkingPrecision->prec,InterpolationOrder->All];
  {\[Psi]up0,d\[Psi]up0,\[Psi]up1,d\[Psi]up1} = {X0,Y0,X1,Y1}/.First@NDSolve[eqinf,{X0,Y0,X1,Y1},{r,r0up,rout},Method->"StiffnessSwitching",WorkingPrecision->prec,InterpolationOrder->All];
  
   (* Remove local variables not garbage collected*)
  Remove[X0,Y0,X1,Y1,r];
  ClearSystemCache[];

{{\[Psi]in0,d\[Psi]in0,\[Psi]in1,d\[Psi]in1},{\[Psi]up0,d\[Psi]up0,\[Psi]up1,d\[Psi]up1}}
]


TeukolskyRadialHS[r_,s_,l_,m_,a_,\[Omega]0_,\[Omega]1_,\[Psi]_]:=Module[{M=1,rp,rm,rs,\[Psi]in0,d\[Psi]in0,\[Psi]in1,d\[Psi]in1,\[Psi]up0,d\[Psi]up0,\[Psi]up1,d\[Psi]up1,\[CapitalDelta],dfacexp0,dfacexp1,Rin0,Rin1,dRin0,dRin1,Rup0,Rup1,dRup0,dRup1,resfac},
	{{\[Psi]in0,d\[Psi]in0,\[Psi]in1,d\[Psi]in1},{\[Psi]up0,d\[Psi]up0,\[Psi]up1,d\[Psi]up1}}=\[Psi];
	rp = M+Sqrt[M^2-a^2];
	rm = M-Sqrt[M^2-a^2];
	rs=((2M rp)/(rp-rm) Log[(r-rp)/(2M)]-(2M rm)/(rp-rm) Log[(r-rm)/(2M)]+r);

	 \[CapitalDelta] = r^2+a^2-2 M r;
	
	dfacexp0 = -(1/r +2 s (r-M)/\[CapitalDelta])+(I/\[CapitalDelta])(# (r^2+a^2)\[Omega]0+a m)&;
	
	dfacexp1[H_]:= I H \[Omega]1((r^2+a^2)/\[CapitalDelta]+rs dfacexp0[H]);
	
	Rin0 = rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])] \[Psi]in0[r];
	Rin1 = rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])](\[Psi]in1[r]-I rs \[Omega]1  \[Psi]in0[r]);
	
	dRin0 = rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])](\[Psi]in0[r] dfacexp0[-1]+ d\[Psi]in0[r]);
	dRin1 = rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])](\[Psi]in1[r] dfacexp0[-1]+ \[Psi]in0[r] dfacexp1[-1]+d\[Psi]in1[r]-I rs \[Omega]1 d\[Psi]in0[r]);
	
	Rup0 = \[Psi]up0[r];
	Rup1 = \[Psi]up1[r]+I rs \[Omega]1 Rup0;
	
	dRup0 = \[Psi]up0[r] dfacexp0[1]+d\[Psi]up0[r];
	dRup1 = \[Psi]up1[r] dfacexp0[1]+\[Psi]up0[r] dfacexp1[1]+d\[Psi]up1[r] +I rs \[Omega]1 d\[Psi]up0[r];
	
	resfac[H_]:=r^(-1) \[CapitalDelta]^(-s)Exp[H I \[Omega]0 rs]Exp[I m a/(rp-rm) Log[(r-rp)/(r-rm)]];
	
	{resfac[-1]{Rin0,dRin0,Rin1,dRin1},resfac[1]{Rup0,dRup0,Rup1,dRup1}}
]


(* ::Text:: *)
(*Function which returns radial function and its correction*)


RCorrection[s_,l_,m_,a_,\[Omega]_,\[Omega]S_,\[Lambda]_,\[Lambda]S_,e_,p_]:=Module[{\[Psi]},
	\[Psi]=TeukolskySolverHS1spinecc[s,l,m,a,\[Omega],\[Omega]S,\[Lambda],\[Lambda]S,e,p];
	Function[{r},
		Module[{R},
				R=TeukolskyRadialHS[r,s,l,m,a,\[Omega],\[Omega]S,\[Psi]];
				Association[
							"R0"-><|
							"In"-><|"R"->R[[1,1]],"dR"->R[[1,2]]|>,
							"Up"-><|"R"->R[[2,1]],"dR"->R[[2,2]]|>
								|>,
							"R1"-><|
							"In"-><|"R"->R[[1,3]],"dR"->R[[1,4]]|>,
							"Up"-><|"R"->R[[2,3]],"dR"->R[[2,4]]|>
							|>
						]
				]
	]
]


(* ::Section::Closed:: *)
(*Spinning fluxes*)


(* ::Subsection::Closed:: *)
(*Functions f_ab^(i)*)


fabi[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,a_]:={-2*\[Zeta]^2/\[CapitalDelta]^2*(L1L2S-2I*a/\[Zeta]*sin\[Theta]*L2S),
                                                       4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*((I*K/\[CapitalDelta]+1/\[Zeta]+1/\[Zeta]bar)*L2S-a*sin\[Theta]*K/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S),
                                                       4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*(L2S+I*a*sin\[Theta]*(1/\[Zeta]bar-1/\[Zeta])*S),
                                                       \[Zeta]^2/\[Zeta]bar^2*(I*(dK*\[CapitalDelta]-K*d\[CapitalDelta])-2/\[Zeta]*I*K*\[CapitalDelta]+K^2)/\[CapitalDelta]^2*S,
                                                       -2*\[Zeta]^2/\[Zeta]bar^2*(1/\[Zeta]+I*K/\[CapitalDelta])*S,
                                                       -\[Zeta]^2/\[Zeta]bar^2*S};
                                                       
dfabidS[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,K1_,dK1_,S_,L2S_,L1L2S_,S1_,L2S1_,L1L2S1_,a_]:={
      -2*\[Zeta]^2/\[CapitalDelta]^2*(L1L2S1-2I*a/\[Zeta]*sin\[Theta]*L2S1),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*((I*K/\[CapitalDelta]+1/\[Zeta]+1/\[Zeta]bar)*L2S1+(I*K1/\[CapitalDelta])*L2S-a*sin\[Theta]*K/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S1-a*sin\[Theta]*K1/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*(L2S1+I*a*sin\[Theta]*(1/\[Zeta]bar-1/\[Zeta])*S1),
      \[Zeta]^2/\[Zeta]bar^2*((I*(dK*\[CapitalDelta]-K*d\[CapitalDelta])-2/\[Zeta]*I*K*\[CapitalDelta]+K^2)/\[CapitalDelta]^2*S1+(I*(dK1*\[CapitalDelta]-K1*d\[CapitalDelta])-2/\[Zeta]*I*K1*\[CapitalDelta]+2*K*K1)/\[CapitalDelta]^2*S),
      -2*\[Zeta]^2/\[Zeta]bar^2*((1/\[Zeta]+I*K/\[CapitalDelta])*S1+(I*K1/\[CapitalDelta])*S),
      -\[Zeta]^2/\[Zeta]bar^2*S1};

dfabidr[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,a_,\[Omega]_]:={
      1/\[CapitalDelta]^3 (\[CapitalDelta] (4 I a L2S sin\[Theta]-4 L1L2S \[Zeta])+4 \[Zeta] (-2 I a L2S sin\[Theta]+L1L2S \[Zeta]) d\[CapitalDelta]),
      1/(\[Zeta]bar^3 \[CapitalDelta]^3) 2 Sqrt[2] (a S (\[Zeta]-\[Zeta]bar) sin\[Theta] (2 d\[CapitalDelta] \[Zeta] \[Zeta]bar K-(dK \[Zeta] \[Zeta]bar+(-2 \[Zeta]+\[Zeta]bar) K) \[CapitalDelta])+L2S (-I \[Zeta] \[Zeta]bar K (2 d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-2 \[Zeta]bar) \[CapitalDelta])+\[CapitalDelta] (\[Zeta] \[Zeta]bar (I dK \[Zeta] \[Zeta]bar-d\[CapitalDelta] (\[Zeta]+\[Zeta]bar))-(\[Zeta]-\[Zeta]bar) (2 \[Zeta]+\[Zeta]bar) \[CapitalDelta]))),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) 2 Sqrt[2] (-L2S \[Zeta] \[Zeta]bar (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-2 \[Zeta]bar) \[CapitalDelta])-I a S (\[Zeta]-\[Zeta]bar) sin\[Theta] (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(2 \[Zeta]-\[Zeta]bar) \[CapitalDelta])),
      -(1/(\[Zeta]bar^3 \[CapitalDelta]^3))2 I S (-I \[Zeta] K^2 (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-\[Zeta]bar) \[CapitalDelta])+\[Zeta]^2 \[CapitalDelta] (dK d\[CapitalDelta] \[Zeta]bar+(dK-\[Zeta]bar \[Omega]) \[CapitalDelta])+K (-d\[CapitalDelta]^2 \[Zeta]^2 \[Zeta]bar+\[Zeta]^2 (-d\[CapitalDelta]+\[Zeta]bar+I dK \[Zeta]bar) \[CapitalDelta]+(-2 \[Zeta]+\[Zeta]bar) \[CapitalDelta]^2)),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) 2 S (I \[Zeta] K (d\[CapitalDelta] \[Zeta] \[Zeta]bar+2 (\[Zeta]-\[Zeta]bar) \[CapitalDelta])+\[CapitalDelta] (-I dK \[Zeta]^2 \[Zeta]bar+(2 \[Zeta]-\[Zeta]bar) \[CapitalDelta])),
      (2 S \[Zeta] (\[Zeta]-\[Zeta]bar))/\[Zeta]bar^3};

dfabid\[Theta][\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,dSd\[Theta]_,dL2Sd\[Theta]_,dL1L2Sd\[Theta]_,a_]:={
      (-4 a^2 L2S sin\[Theta]^2+4 I a (dL2Sd\[Theta]-L1L2S) sin\[Theta] \[Zeta]+2 \[Zeta] (-((dL1L2Sd\[Theta]+L2S) \[Zeta])+L2S \[Zeta]bar))/\[CapitalDelta]^2,
      Sqrt[2] (-2 I a^2 K S sin\[Theta]^2 (2 \[Zeta]-\[Zeta]bar) (\[Zeta]+\[Zeta]bar)+\[Zeta] \[Zeta]bar (-I K S (\[Zeta]-\[Zeta]bar)^2+2 I dL2Sd\[Theta] K \[Zeta] \[Zeta]bar+2 dL2Sd\[Theta] \[CapitalDelta] (\[Zeta]+\[Zeta]bar))+2 a sin\[Theta] (dSd\[Theta] K \[Zeta] \[Zeta]bar (-\[Zeta]+\[Zeta]bar)+I L2S \[CapitalDelta] (\[Zeta]+\[Zeta]bar) (2 \[Zeta]+\[Zeta]bar)-K L2S \[Zeta] \[Zeta]bar (\[Zeta]+2 \[Zeta]bar)))/(\[CapitalDelta]^2 \[Zeta]bar^3),
      Sqrt[2] (-2 a^2 S sin\[Theta]^2 (2 \[Zeta]-\[Zeta]bar) (\[Zeta]+\[Zeta]bar)-\[Zeta] \[Zeta]bar (S (\[Zeta]-\[Zeta]bar)^2-2 dL2Sd\[Theta] \[Zeta] \[Zeta]bar)+2 I a sin\[Theta] \[Zeta] \[Zeta]bar (dSd\[Theta] (\[Zeta]-\[Zeta]bar)+L2S (\[Zeta]+2 \[Zeta]bar)))/(\[CapitalDelta] \[Zeta]bar^3),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) (\[Zeta] (dSd\[Theta] \[Zeta] \[Zeta]bar+2 I a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) K^2+dK \[Zeta] (I dSd\[Theta] \[Zeta] \[Zeta]bar-2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) \[CapitalDelta]+K (d\[CapitalDelta] \[Zeta] (-I dSd\[Theta] \[Zeta] \[Zeta]bar+2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar))+2 (-I dSd\[Theta] \[Zeta] \[Zeta]bar+a S sin\[Theta] (2 \[Zeta]+\[Zeta]bar)) \[CapitalDelta])),
      -(1/(\[Zeta]bar^3 \[CapitalDelta]))2 (\[Zeta] (I dSd\[Theta] \[Zeta] \[Zeta]bar-2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) K+(dSd\[Theta] \[Zeta] \[Zeta]bar+I a S sin\[Theta] (2 \[Zeta]+\[Zeta]bar)) \[CapitalDelta]),
      -((\[Zeta] (dSd\[Theta] \[Zeta] \[Zeta]bar+2 I a S sin\[Theta] (\[Zeta]+\[Zeta]bar)))/\[Zeta]bar^3)};



(* ::Subsection::Closed:: *)
(*Generic - fixed turning points*)


TeukolskySpinModeCorrectionDH[l_,m_,n_,k_,orbitCorrection_]:=Module[{a,p,e,x,h1,h2,h3,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,rg,zg,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,
    SWSH,SWSHS,R,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    lists,ir,i\[Theta],wrsub,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,
    RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,\[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1,\[CapitalDelta]tspar,\[Delta]rpar,\[Delta]zpar,\[Delta]Vrpar,\[CapitalDelta]\[Phi]spar,\[Delta]Vzpar,wrlistshort,wzlistshort,A2g,A1g,A0g,A2par,A1par,A0par,B3par,B2par,B1par},
  
  {a,p,e,x}=orbitCorrection["OrbitalElements"];
  
  (*These functions are needed for the corrections to the coordinate time and azimuthal velocities. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)    
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_]:=1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  h3[r_,z_]:=-((2 a r z)/(Sqrt[Kc] (r^2+a^2 z^2)^2));
  
  En = orbitCorrection["Eg"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzg"];
  Kc = orbitCorrection["Kg"];
  
  En1 = orbitCorrection["Es"]; (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
  Lz1 = orbitCorrection["Jzs"];
  
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Geodesic BL frequencies *)
  {\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the frequencies *)
 
  rg[wr_] := rgfun[wr,a,p,e,x,{En,Lz,Kc}]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  zg[wz_] := zgfun[wz,a,p,e,x,{En,Lz,Kc}];
  
  \[CapitalGamma]  = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega]  = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + n*\[CapitalOmega]r1 + k*\[CapitalOmega]\[Theta]1;
 
  
 {\[Lambda],\[Lambda]1,SWSH,SWSHS}=angpar[-2,l,m,a,SetPrecision[\[Omega],48],SetPrecision[\[Omega]1,48]]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
 R = RCorrection[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,e,p];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  sumPlus0 = sumPlus1 = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[Pi  ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[Pi/4]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[0   ]+k)]],32];
                
  wrlistshort=Table[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]],{ir,1,stepsr/2}];
  wzlistshort=Table[N[(iz-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]],{iz,1,steps\[Theta]/2}];
                   
  (*Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in wz"];*)
  \[Theta]list = {};(* List for functions of \[Theta] *)
  (*Print["Sampling the spin corrections to the EoM"];*)
 If[stepsr<=steps\[Theta],
   \[Delta]rpar=Table[orbitCorrection["\[Delta]rpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[Delta]zpar=Table[orbitCorrection["\[Delta]zpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[CapitalDelta]tspar=Table[orbitCorrection["\[CapitalDelta]\[Delta]tpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[CapitalDelta]\[Phi]spar=Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]par"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[Delta]Vrpar=Table[orbitCorrection["\[Delta]vrpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[Delta]Vzpar=Table[orbitCorrection["\[Delta]vzpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   ,
   \[Delta]rpar=Transpose[Table[orbitCorrection["\[Delta]rpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[Delta]zpar=Transpose[Table[orbitCorrection["\[Delta]zpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[CapitalDelta]tspar=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]tpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[CapitalDelta]\[Phi]spar=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]par"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[Delta]Vrpar=Transpose[Table[orbitCorrection["\[Delta]vrpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[Delta]Vzpar=Transpose[Table[orbitCorrection["\[Delta]vzpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   ];
   (*Print["Computing the amplitudes and fluxes"];*)  
  For[ir = 1, ir <= stepsr/2, ir++,(* integration over w_r *)
    wrsub = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = rg[wrsub];
    Ur = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];(* Geodesic radial velocity *)
    expr = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]trg"][wrsub])-m*(orbitCorrection["\[CapitalDelta]\[Phi]rg"][wrsub])+2Pi*n*(ir-1/2)/stepsr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    K1  = (rp^2+a^2)*\[Omega]1;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2*rp*\[Omega];
    dK1 = 2*rp*\[Omega]1;
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    V1  = -(2*K*K1 + 4I*(rp-1)*K1)/\[CapitalDelta] + 8*I*\[Omega]1*rp + \[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RInrp    = R[rp]["R0"]["In"]["R"]; (* radial function *)
    dRInrp   = R[rp]["R0"]["In"]["dR"];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    RInrp1   = R[rp]["R1"]["In"]["R"]; (* Linear part of the radial function *)
    dRInrp1  = R[rp]["R1"]["In"]["dR"];
    ddRInrp1 = (V*RInrp1+V1*RInrp+d\[CapitalDelta]*dRInrp1)/\[CapitalDelta];
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R[rp]["R0"]["Up"]["R"];
    dRUprp   = R[rp]["R0"]["Up"]["dR"];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    RUprp1   = R[rp]["R1"]["Up"]["R"];
    dRUprp1  = R[rp]["R1"]["Up"]["dR"];
    ddRUprp1 = (V*RUprp1+V1*RUprp+d\[CapitalDelta]*dRUprp1)/\[CapitalDelta];
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
       
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++, (* integration over w_\[Theta] *)
      w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta], Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        zp = zg[w\[Theta]];
        Uz = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]])-m*(orbitCorrection["\[CapitalDelta]\[Phi]zg"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2]; (*\[Theta] is  only defined between [0,\[Pi]), thus Sin(\[Theta]) is semi-positive.*)
        S = SWSH[ArcCos[zp]][[1]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
        S1 = SWSHS[ArcCos[zp]][[1]];  (*  Linear part of S(\[Theta](z))  *)
        dSd\[Theta] = SWSH[ArcCos[zp]][[2]]; (* First derivative of S wrt \[Theta] *)
        dS1d\[Theta] = SWSHS[ArcCos[zp]][[2]]; (* Linear part of first derivative of S wrt \[Theta] *)
        d2Sd\[Theta]2 = SWSH[ArcCos[zp]][[3]]; (* second derivative of S from Teukolsky equation *)
        d2S1d\[Theta]2 = SWSHS[ArcCos[zp]][[3]];(* Linear part of second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        L2S1 = dS1d\[Theta]-(S1 (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2) \[Omega]1))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2) \[Omega]1))/sin\[Theta]p + 
                 S1 (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) \[Omega]1 - a^2 (-1+zp^2) 2*\[Omega]*\[Omega]1);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}];,
        {zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,K1,dK1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
      vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
      Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
      Snmb = Conjugate[Snm];
      Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
      Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      
      
      rp1 =  {\[Delta]rpar[[ir,i\[Theta]]], \[Delta]rpar[[ir,-i\[Theta]]], \[Delta]rpar[[ir,-i\[Theta]]], \[Delta]rpar[[ir,i\[Theta]]]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
      zp1 =  {\[Delta]zpar[[ir,i\[Theta]]],-\[Delta]zpar[[ir,-i\[Theta]]],-\[Delta]zpar[[ir,-i\[Theta]]], \[Delta]zpar[[ir,i\[Theta]]]};
      Urp1 = {\[Delta]Vrpar[[ir,i\[Theta]]],-\[Delta]Vrpar[[ir,-i\[Theta]]], \[Delta]Vrpar[[ir,-i\[Theta]]],-\[Delta]Vrpar[[ir,i\[Theta]]]};
      Uzp1 = {\[Delta]Vzpar[[ir,i\[Theta]]], \[Delta]Vzpar[[ir,-i\[Theta]]],-\[Delta]Vzpar[[ir,-i\[Theta]]],-\[Delta]Vzpar[[ir,i\[Theta]]]};
      \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
      
      exp1 = I*{(\[Omega]*\[CapitalDelta]tspar[[ir, i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir, i\[Theta]]])+\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]+\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]],
               -(\[Omega]*\[CapitalDelta]tspar[[ir,-i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir,-i\[Theta]]])-\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]+\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]],
                (\[Omega]*\[CapitalDelta]tspar[[ir,-i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir,-i\[Theta]]])+\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]-\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]],
               -(\[Omega]*\[CapitalDelta]tspar[[ir, i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir, i\[Theta]]])-\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]-\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]]};(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      
      vn1  = -( 2*rp*rp1*En - ((rp^2+a^2)*En - a*Lz + Ur)/\[CapitalSigma]*\[CapitalSigma]1 + 
               (rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz) + Urp1)/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = ( I*2*a*zp*zp1*En - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
               -I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz)) + Uzp1)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
                          
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      
      A2g = Cmmbmb*fmbmb2;
      A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
      A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
      
      (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
      B3par = -AmbmbrS*fmbmb2;
      B2par = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
      B1par = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS fnn0;
      
      A2par = AmbmbrS*dfmbmb2dr+Ambmb\[Theta]S*dfmbmb2d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2 
              +Cmmbmb*fmbmb21;
      A1par = AmbmbrS*dfmbmb1dr+Ambmb\[Theta]S*dfmbmb1d\[Theta]+AnmbrS*dfnmb1dr+Anmb\[Theta]S*dfnmb1d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1
              +(Cmmbmb*fmbmb11+Cmnmb*fnmb11);
      A0par = AmbmbrS*dfmbmb0dr+Ambmb\[Theta]S*dfmbmb0d\[Theta]+AnmbrS*dfnmb0dr+Anmb\[Theta]S*dfnmb0d\[Theta]+AnnrS*dfnn0dr+Ann\[Theta]S*dfnn0d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0 
              +(Cmmbmb*fmbmb01+Cmnmb*fnmb01+Cmnn*fnn01);
      
      sumPlus0  += Total[\[CapitalSigma]*((A0g*RInrp - A1g*dRInrp + A2g*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Total of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*((A0g*RUprp - A1g*dRUprp + A2g*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})];
      
      sumPlus1  += Total[\[CapitalSigma]*(A0par*RInrp - (A1par+B1par)*dRInrp + (A2par+B2par)ddRInrp - B3par*dddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}+
                         \[CapitalSigma]*(A0g*RInrp1 - A1g*dRInrp1 + A2g*ddRInrp1)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*(A0par*RUprp - (A1par+B1par)*dRUprp + (A2par+B2par)ddRUprp - B3par*dddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}+
                         \[CapitalSigma]*(A0g*RUprp1 - A1g*dRUprp1 + A2g*ddRUprp1)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}];      
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/\[CapitalDelta]; (* Invariant Wronskian *)
  W1 = (RInrp1*dRUprp + RInrp*dRUprp1 - dRInrp1*RUprp - dRInrp*RUprp1)/\[CapitalDelta]; (* Linear part of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - W1/W*sumPlus0 )/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - W1/W*sumMinus0)/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  
  Association[
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]g"->\[Omega],
    "\[Omega]s"->\[Omega]1,
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  ]
]


(* ::Subsection::Closed:: *)
(*Generic - fixed constants of motion*)


TeukolskySpinModeCorrectionFC[l_,m_,n_,k_,orbitCorrection_]:=Module[{a,p,e,x,h1,h2,h3,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,rg,zg,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,
    SWSH,SWSHS,R,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    lists,ir,i\[Theta],wrsub,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,
    RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,\[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1,\[CapitalDelta]tspar,\[Delta]rpar,\[Delta]zpar,\[Delta]Vrpar,\[CapitalDelta]\[Phi]spar,\[Delta]Vzpar,wrlistshort,wzlistshort,A2g,A1g,A0g,A2par,A1par,A0par,B3par,B2par,B1par},
  
  {a,p,e,x}=orbitCorrection["OrbitalElements"];
  
  (*These functions are needed for the corrections to the coordinate time and azimuthal velocities. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)    
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_]:=1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  h3[r_,z_]:=-((2 a r z)/(Sqrt[Kc] (r^2+a^2 z^2)^2));
  
  En = orbitCorrection["Eg"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzg"];
  Kc = orbitCorrection["Kg"];
  
  En1 = 0; (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
  Lz1 = 0;
  
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Geodesic BL frequencies *)
  {\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the frequencies *)
 
  rg[wr_] := rgfun[wr,a,p,e,x,{En,Lz,Kc}]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  zg[wz_] := zgfun[wz,a,p,e,x,{En,Lz,Kc}];
  
  \[CapitalGamma]  = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega]  = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + n*\[CapitalOmega]r1 + k*\[CapitalOmega]\[Theta]1;
 
  
 {\[Lambda],\[Lambda]1,SWSH,SWSHS}=angpar[-2,l,m,a,SetPrecision[\[Omega],48],SetPrecision[\[Omega]1,48]]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
 R = RCorrection[-2,l,m,a,\[Omega],\[Omega]1,\[Lambda],\[Lambda]1,e,p];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  sumPlus0 = sumPlus1 = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[Pi  ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[Pi/4]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[0   ]+k)]],32];
                
  wrlistshort=Table[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]],{ir,1,stepsr/2}];
  wzlistshort=Table[N[(iz-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]],{iz,1,steps\[Theta]/2}];
                   
  (*Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in wz"];*)
  \[Theta]list = {};(* List for functions of \[Theta] *)
  (*Print["Sampling the spin corrections to the EoM"];*)
 If[stepsr<=steps\[Theta],
   \[Delta]rpar=Table[orbitCorrection["\[Delta]rpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[Delta]zpar=Table[orbitCorrection["\[Delta]zpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[CapitalDelta]tspar=Table[orbitCorrection["\[CapitalDelta]\[Delta]tpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[CapitalDelta]\[Phi]spar=Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]par"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[Delta]Vrpar=Table[orbitCorrection["\[Delta]vrpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   \[Delta]Vzpar=Table[orbitCorrection["\[Delta]vzpar"][wrlistshort[[ir]], wzlistshort],{ir,Length[wrlistshort]}];
   ,
   \[Delta]rpar=Transpose[Table[orbitCorrection["\[Delta]rpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[Delta]zpar=Transpose[Table[orbitCorrection["\[Delta]zpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[CapitalDelta]tspar=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]tpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[CapitalDelta]\[Phi]spar=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]par"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[Delta]Vrpar=Transpose[Table[orbitCorrection["\[Delta]vrpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   \[Delta]Vzpar=Transpose[Table[orbitCorrection["\[Delta]vzpar"][wrlistshort, wzlistshort[[iz]]],{iz,Length[wzlistshort]}]];
   ];
   (*Print["Computing the amplitudes and fluxes"];*)  
  For[ir = 1, ir <= stepsr/2, ir++,(* integration over w_r *)
    wrsub = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = rg[wrsub];
    Ur = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];(* Geodesic radial velocity *)
    expr = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]trg"][wrsub])-m*(orbitCorrection["\[CapitalDelta]\[Phi]rg"][wrsub])+2Pi*n*(ir-1/2)/stepsr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    K1  = (rp^2+a^2)*\[Omega]1;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2*rp*\[Omega];
    dK1 = 2*rp*\[Omega]1;
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    V1  = -(2*K*K1 + 4I*(rp-1)*K1)/\[CapitalDelta] + 8*I*\[Omega]1*rp + \[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RInrp    = R[rp]["R0"]["In"]["R"]; (* radial function *)
    dRInrp   = R[rp]["R0"]["In"]["dR"];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    RInrp1   = R[rp]["R1"]["In"]["R"]; (* Linear part of the radial function *)
    dRInrp1  = R[rp]["R1"]["In"]["dR"];
    ddRInrp1 = (V*RInrp1+V1*RInrp+d\[CapitalDelta]*dRInrp1)/\[CapitalDelta];
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R[rp]["R0"]["Up"]["R"];
    dRUprp   = R[rp]["R0"]["Up"]["dR"];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    RUprp1   = R[rp]["R1"]["Up"]["R"];
    dRUprp1  = R[rp]["R1"]["Up"]["dR"];
    ddRUprp1 = (V*RUprp1+V1*RUprp+d\[CapitalDelta]*dRUprp1)/\[CapitalDelta];
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
       
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++, (* integration over w_\[Theta] *)
      w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta], Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        zp = zg[w\[Theta]];
        Uz = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]])-m*(orbitCorrection["\[CapitalDelta]\[Phi]zg"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2]; (*\[Theta] is  only defined between [0,\[Pi]), thus Sin(\[Theta]) is semi-positive.*)
        S = SWSH[ArcCos[zp]][[1]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
        S1 = SWSHS[ArcCos[zp]][[1]];  (*  Linear part of S(\[Theta](z))  *)
        dSd\[Theta] = SWSH[ArcCos[zp]][[2]]; (* First derivative of S wrt \[Theta] *)
        dS1d\[Theta] = SWSHS[ArcCos[zp]][[2]]; (* Linear part of first derivative of S wrt \[Theta] *)
        d2Sd\[Theta]2 = SWSH[ArcCos[zp]][[3]]; (* second derivative of S from Teukolsky equation *)
        d2S1d\[Theta]2 = SWSHS[ArcCos[zp]][[3]];(* Linear part of second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        L2S1 = dS1d\[Theta]-(S1 (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2) \[Omega]1))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2) \[Omega]1))/sin\[Theta]p + 
                 S1 (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) \[Omega]1 - a^2 (-1+zp^2) 2*\[Omega]*\[Omega]1);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}];,
        {zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,K1,dK1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
      vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
      Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
      Snmb = Conjugate[Snm];
      Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
      Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      
      rp1 =  {\[Delta]rpar[[ir,i\[Theta]]], \[Delta]rpar[[ir,-i\[Theta]]], \[Delta]rpar[[ir,-i\[Theta]]], \[Delta]rpar[[ir,i\[Theta]]]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
      zp1 =  {\[Delta]zpar[[ir,i\[Theta]]],-\[Delta]zpar[[ir,-i\[Theta]]],-\[Delta]zpar[[ir,-i\[Theta]]], \[Delta]zpar[[ir,i\[Theta]]]};
      Urp1 = {\[Delta]Vrpar[[ir,i\[Theta]]],-\[Delta]Vrpar[[ir,-i\[Theta]]], \[Delta]Vrpar[[ir,-i\[Theta]]],-\[Delta]Vrpar[[ir,i\[Theta]]]};
      Uzp1 = {\[Delta]Vzpar[[ir,i\[Theta]]], \[Delta]Vzpar[[ir,-i\[Theta]]],-\[Delta]Vzpar[[ir,-i\[Theta]]],-\[Delta]Vzpar[[ir,i\[Theta]]]};
      \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
      
      exp1 = I*{(\[Omega]*\[CapitalDelta]tspar[[ir, i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir, i\[Theta]]])+\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]+\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]],
               -(\[Omega]*\[CapitalDelta]tspar[[ir,-i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir,-i\[Theta]]])-\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]+\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]],
                (\[Omega]*\[CapitalDelta]tspar[[ir,-i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir,-i\[Theta]]])+\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]-\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]],
               -(\[Omega]*\[CapitalDelta]tspar[[ir, i\[Theta]]]-m*\[CapitalDelta]\[Phi]spar[[ir, i\[Theta]]])-\[Omega]1*orbitCorrection["\[CapitalDelta]trg"][wrsub]-\[Omega]1*orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]]};(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      vn1  = -( 2*rp*rp1*En - ((rp^2+a^2)*En - a*Lz + Ur)/\[CapitalSigma]*\[CapitalSigma]1 + 
               (rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz) + Urp1)/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = ( I*2*a*zp*zp1*En - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
               -I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz)) + Uzp1)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
                          
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      
      A2g = Cmmbmb*fmbmb2;
      A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
      A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
      
      (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
      B3par = -AmbmbrS*fmbmb2;
      B2par = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
      B1par = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS fnn0;
      
      A2par = AmbmbrS*dfmbmb2dr+Ambmb\[Theta]S*dfmbmb2d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2 
              +Cmmbmb*fmbmb21;
      A1par = AmbmbrS*dfmbmb1dr+Ambmb\[Theta]S*dfmbmb1d\[Theta]+AnmbrS*dfnmb1dr+Anmb\[Theta]S*dfnmb1d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1
              +(Cmmbmb*fmbmb11+Cmnmb*fnmb11);
      A0par = AmbmbrS*dfmbmb0dr+Ambmb\[Theta]S*dfmbmb0d\[Theta]+AnmbrS*dfnmb0dr+Anmb\[Theta]S*dfnmb0d\[Theta]+AnnrS*dfnn0dr+Ann\[Theta]S*dfnn0d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0 
              +(Cmmbmb*fmbmb01+Cmnmb*fnmb01+Cmnn*fnn01);
      
      sumPlus0  += Total[\[CapitalSigma]*((A0g*RInrp - A1g*dRInrp + A2g*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Total of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*((A0g*RUprp - A1g*dRUprp + A2g*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})];
      
      sumPlus1  += Total[\[CapitalSigma]*(A0par*RInrp - (A1par+B1par)*dRInrp + (A2par+B2par)ddRInrp - B3par*dddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}+
                         \[CapitalSigma]*(A0g*RInrp1 - A1g*dRInrp1 + A2g*ddRInrp1)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*(A0par*RUprp - (A1par+B1par)*dRUprp + (A2par+B2par)ddRUprp - B3par*dddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}+
                         \[CapitalSigma]*(A0g*RUprp1 - A1g*dRUprp1 + A2g*ddRUprp1)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/\[CapitalDelta]; (* Invariant Wronskian *)
  W1 = (RInrp1*dRUprp + RInrp*dRUprp1 - dRInrp1*RUprp - dRInrp*RUprp1)/\[CapitalDelta]; (* Linear part of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - W1/W*sumPlus0 )/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - W1/W*sumMinus0)/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  
  Association[
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]g"->\[Omega],
    "\[Omega]s"->\[Omega]1,
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  ]
]


(* ::Section::Closed:: *)
(*Correction to the amplitudes - orthogonal component of the spin*)


(* ::Subsubsection::Closed:: *)
(*Plus mode (j=+1)*)


TeukolskySpinModeCorrectionOrtPlus[l_,m_,n_,k_,orbitCorrection_]:=Module[{a,p,e,x,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],rg,zg,\[CapitalGamma],\[Omega],\[CapitalOmega]p,
    SWSH,SWSHS,R,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    lists,ir,i\[Theta],wrsub,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,Utr,Utz,Ut,U\[Phi]r,U\[Phi]z,U\[Phi],expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,
    RUprp,dRUprp,ddRUprp,dddRUprp,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,\[Psi],Sln,Slmb,Slmbpar,Snm,Snmpar,Snmb,Snmbpar,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,Utp1,U\[Phi]p1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1,\[CapitalDelta]tsort,\[Delta]rort,\[Delta]zort,\[Delta]Vrort,\[CapitalDelta]\[Phi]sort,\[Delta]Vzort,\[Delta]Vtort,\[Delta]V\[Phi]ort,wrlist,wzlist,A2g,A1g,A0g,A2ort,A1ort,A0ort,B3ort,B2ort,B1ort},
    
  {a,p,e,x}=orbitCorrection["OrbitalElements"];
      
  En = orbitCorrection["Eg"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzg"];
  Kc = orbitCorrection["Kg"];
  
  En1 = 0; (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
  Lz1 = 0;
  
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Geodesic BL frequencies *)
  
  rg[wr_] := rgfun[wr,a,p,e,x,{En,Lz,Kc}]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  zg[wz_] := zgfun[wz,a,p,e,x,{En,Lz,Kc}];
  
  \[CapitalGamma]  = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalOmega]p = orbitCorrection["BLPrecessionFrequency"];
  \[Omega]  = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta] + \[CapitalOmega]p; (* Geodesic frequency and the linear correction *)
 
  {\[Lambda],SWSH} = angparLeading[-2,l,m,a,SetPrecision[\[Omega],48]]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  R = RLeading[-2,l,m,a,\[Omega],\[Lambda],e,p];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus0 = sumPlus1 = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[Pi  ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[Pi/4]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[0   ]+k)]],32];
                
  wrlist=Table[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]],{ir,1,stepsr}];
  wzlist=Table[N[(iz-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]],{iz,1,steps\[Theta]}];
                   
  (*Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in wz"];*)
  \[Theta]list = {};(* List for functions of \[Theta] *)
  (*Print["Sampling the spin corrections to the EoM"];*)
 If[stepsr<=steps\[Theta],
   \[Delta]rort=Table[orbitCorrection["\[Delta]rortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]zort=Table[orbitCorrection["\[Delta]zortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[CapitalDelta]tsort=Table[orbitCorrection["\[CapitalDelta]\[Delta]tortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[CapitalDelta]\[Phi]sort=Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]ortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]Vrort=Table[orbitCorrection["\[Delta]vrortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]Vzort=Table[orbitCorrection["\[Delta]vzortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]Vtort=Table[orbitCorrection["\[Delta]vtortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]V\[Phi]ort=Table[orbitCorrection["\[Delta]v\[Phi]ortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   ,
   \[Delta]rort=Transpose[Table[orbitCorrection["\[Delta]rortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]zort=Transpose[Table[orbitCorrection["\[Delta]zortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[CapitalDelta]tsort=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]tortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[CapitalDelta]\[Phi]sort=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]ortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]Vrort=Transpose[Table[orbitCorrection["\[Delta]vrortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]Vzort=Transpose[Table[orbitCorrection["\[Delta]vzortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]Vtort=Transpose[Table[orbitCorrection["\[Delta]vtortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]V\[Phi]ort=Transpose[Table[orbitCorrection["\[Delta]v\[Phi]ortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   ];
   (*Print["Computing the corrections to the amplitudes - orthogonal component in the spin"];*) 
For[ir = 1, ir <= stepsr, ir++,(* integration over w_r *)
    wrsub = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = rg[wrsub];
    Ur = Sign[Pi-Mod[wrsub,2Pi]]*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];(* Geodesic radial velocity *)
    Utr = En*(rp^2+a^2)^2/(rp^2-2rp+a^2)-2a*rp*Lz/(rp^2-2rp+a^2);
    U\[Phi]r = a/(rp^2-2rp+a^2)(En(rp^2+a^2)-a*Lz)-a*En;
    expr = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]trg"][wrsub])-m*(orbitCorrection["\[CapitalDelta]\[Phi]rg"][wrsub])+2Pi*n*(ir-1/2)/stepsr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2*rp*\[Omega];
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RInrp    = R[rp]["R0"]["In"]["R"]; (* radial function *)
    dRInrp   = R[rp]["R0"]["In"]["dR"];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R[rp]["R0"]["Up"]["R"];
    dRUprp   = R[rp]["R0"]["Up"]["dR"];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
       
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta], i\[Theta]++, (* integration over w_\[Theta] *)
      w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta], Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        zp = zg[w\[Theta]];
        Uz = dzgd\[Lambda]fun[w\[Theta],a,p,e,x,{En,Lz,Kc}];(* Polar geodesic velocity *)
        Utz = -En*a^2(1-zp^2);
        U\[Phi]z = Lz/(1-zp^2);
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]])-m*(orbitCorrection["\[CapitalDelta]\[Phi]zg"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2]; (*\[Theta] is  only defined between [0,\[Pi]), thus Sin(\[Theta]) is semi-positive.*)
        S = SWSH[ArcCos[zp]][[1]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
        dSd\[Theta] = SWSH[ArcCos[zp]][[2]]; (* First derivative of S wrt \[Theta] *)
        d2Sd\[Theta]2 = SWSH[ArcCos[zp]][[3]]; (* second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta],Utz,U\[Phi]z}];,
        {zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta],Utz,U\[Phi]z}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      
      \[Psi]=\[Psi]rfun[wrsub,a,p,e,x,{En,Lz,Kc}]+\[Psi]zfun[w\[Theta],a,p,e,x,{En,Lz,Kc}];
      
      Sln  = (Cos[\[Psi]]-I*Sin[\[Psi]])a*zp*Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]/(2Sqrt[Kc]*\[CapitalSigma]); (* Spin tensor in Kinnersley tetrad *) (*Cos[\[Psi]] and Sin[\[Psi]] comes from expanding Cos[\[Psi]+wp] and Sin[\[Psi]+wp] *)
      Snmpar  = (\[Zeta]/Sqrt[Kc])*vm*vn;
      Snmbpar = Conjugate[Snm];
      Snm = -((Cos[\[Psi]]-I*Sin[\[Psi]])(Kc-I a*rp*zp)+(I*Cos[\[Psi]]+Sin[\[Psi]])Sqrt[Kc]\[Zeta]bar)I*Snmpar/(2Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]);  
      Snmb = ((Cos[\[Psi]]-I*Sin[\[Psi]])(Kc+I a*rp*zp)+(I*Cos[\[Psi]]+Sin[\[Psi]])Sqrt[Kc]\[Zeta])I*Snmbpar/(2Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]);
      Slmbpar = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
      Slmb = -((Kc-I a*rp*zp)-Sqrt[Kc]\[Zeta]bar)I*Slmbpar/(2Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]);
      Smmb = (Cos[\[Psi]]-I*Sin[\[Psi]])I*rp*Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]/(2Sqrt[Kc]*\[CapitalSigma]);
      
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      
      rp1 =  \[Delta]rort[[ir,i\[Theta]]]; (* Corrections to the coordinates and four-velocity for each quadrant *)
      zp1 =  \[Delta]zort[[ir,i\[Theta]]];
      Urp1 = \[Delta]Vrort[[ir,i\[Theta]]];
      Uzp1 = \[Delta]Vzort[[ir,i\[Theta]]];
      Utp1 = \[Delta]Vtort[[ir,i\[Theta]]];
      U\[Phi]p1 = \[Delta]V\[Phi]ort[[ir,i\[Theta]]];
      
      Ut = Utr+Utz;
      U\[Phi] = U\[Phi]r+U\[Phi]z;
      
      \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
      
      exp1 = I*(\[Omega]*\[CapitalDelta]tsort[[ir,i\[Theta]]]-m*\[CapitalDelta]\[Phi]sort[[ir,i\[Theta]]]);(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      
      vn1  = -((Utp1*\[CapitalDelta])/(2\[CapitalSigma]^2))+(a*U\[Phi]p1(1-zp^2)*\[CapitalDelta])/(2\[CapitalSigma]^2)-Urp1/(2\[CapitalSigma])+rp1/(2\[CapitalDelta]*\[CapitalSigma]^3)((4rp*(Ut-a*U\[Phi](1-zp^2))\[CapitalDelta]^2+2(rp*(Ur-Ut)+Ut)\[CapitalDelta]*\[CapitalSigma]+2a(-1+rp)U\[Phi](1-zp^2)\[CapitalDelta]*\[CapitalSigma]))
               +zp1/(2\[CapitalDelta]*\[CapitalSigma]^3) ((4a^2zp(Ut-a*U\[Phi](1-zp^2))\[CapitalDelta]^2+2a*zp*\[CapitalDelta](a*Ur-U\[Phi]*\[CapitalDelta])\[CapitalSigma])); (* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = -((I(a^2+rp^2)U\[Phi]p1 (1-zp^2)^(1/2))/(Sqrt[2]\[Zeta]*\[CapitalSigma]))+(I*a*Utp1*(1-zp^2)^(1/2))/(Sqrt[2]\[Zeta]*\[CapitalSigma])-(Uzp1*\[Zeta]bar)/(Sqrt[2](1-zp^2)^(1/2)\[CapitalSigma])
               +rp1/(Sqrt[2](1-zp^2)^(3/2)\[Zeta]^2*\[CapitalSigma]^2)(2rp(-1+zp^2)\[Zeta](I(-a*Ut+(a^2+rp^2)U\[Phi])(-1+zp^2)-Uz*\[Zeta]*\[Zeta]bar)+I(-a*Ut+(a^2+rp^2)*U\[Phi]) (-1+zp^2)^2*\[CapitalSigma]-2I*rp*U\[Phi] (1-zp^2)^2\[Zeta]*\[CapitalSigma]+Uz(-1+zp^2)\[Zeta]^2*\[CapitalSigma])
               +zp1/(Sqrt[2](1-zp^2)^(3/2)\[Zeta]^2*\[CapitalSigma]^2)(2a^2*zp*(-1+zp^2)\[Zeta](I(-a*Ut+(a^2+rp^2)*U\[Phi])(-1+zp^2)-Uz*\[Zeta]*\[Zeta]bar)+a(-a*Ut+(a^2+rp^2)*U\[Phi])(1-zp^2)^2\[CapitalSigma]+I*a*Uz(-1+zp^2)\[Zeta]^2*\[CapitalSigma]-2Uz*zp*\[Zeta]^2\[Zeta]bar*\[CapitalSigma]-zp*\[Zeta]*(I(-a*Ut+(a^2+rp^2)U\[Phi])(-1+zp^2)-Uz*\[Zeta]*\[Zeta]bar)*\[CapitalSigma]);
      
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      
      A2g = Cmmbmb*fmbmb2;
      A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
      A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
      
      (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
      B3ort = -AmbmbrS*fmbmb2;
      B2ort = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
      B1ort = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS fnn0;
      
      A2ort = AmbmbrS*dfmbmb2dr+Ambmb\[Theta]S*dfmbmb2d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2;
      A1ort = AmbmbrS*dfmbmb1dr+Ambmb\[Theta]S*dfmbmb1d\[Theta]+AnmbrS*dfnmb1dr+Anmb\[Theta]S*dfnmb1d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1;
      A0ort = AmbmbrS*dfmbmb0dr+Ambmb\[Theta]S*dfmbmb0d\[Theta]+AnmbrS*dfnmb0dr+Anmb\[Theta]S*dfnmb0d\[Theta]+AnnrS*dfnn0dr+Ann\[Theta]S*dfnn0d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0;
      
      sumPlus0  += Total[\[CapitalSigma]*((A0g*RInrp - A1g*dRInrp + A2g*ddRInrp)*expr*exp\[Theta])]; (* Total of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*((A0g*RUprp - A1g*dRUprp + A2g*ddRUprp)*expr*exp\[Theta])];
      
      sumPlus1  += Total[\[CapitalSigma]*((A0ort*RInrp - (A1ort+B1ort)*dRInrp + (A2ort+B2ort)ddRInrp - B3ort*dddRInrp)*expr*exp\[Theta])]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*((A0ort*RUprp - (A1ort+B1ort)*dRUprp + (A2ort+B2ort)ddRUprp - B3ort*dddRUprp)*expr*exp\[Theta])];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/\[CapitalDelta]; (* Invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1)/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1)/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  
  Association[
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]g"->\[Omega],
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  ]
]


(* ::Subsubsection::Closed:: *)
(*Minus mode (j=-1)*)


TeukolskySpinModeCorrectionOrtMinus[l_,m_,n_,k_,orbitCorrection_]:=Module[{a,p,e,x,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],rg,zg,\[CapitalGamma],\[Omega],\[CapitalOmega]p,
    SWSH,SWSHS,R,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    lists,ir,i\[Theta],wrsub,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,Utr,Utz,Ut,U\[Phi]r,U\[Phi]z,U\[Phi],expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,
    RUprp,dRUprp,ddRUprp,dddRUprp,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,\[Psi],Sln,Slmb,Slmbpar,Snm,Snmpar,Snmb,Snmbpar,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,Utp1,U\[Phi]p1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1,\[CapitalDelta]tsort,\[Delta]rort,\[Delta]zort,\[Delta]Vrort,\[CapitalDelta]\[Phi]sort,\[Delta]Vzort,\[Delta]Vtort,\[Delta]V\[Phi]ort,wrlist,wzlist,A2g,A1g,A0g,A2ort,A1ort,A0ort,B3ort,B2ort,B1ort},
    
  {a,p,e,x}=orbitCorrection["OrbitalElements"];
       
  En = orbitCorrection["Eg"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzg"];
  Kc = orbitCorrection["Kg"];
  
  En1 = 0; (* Linear corrections to the constants of motion. These terms only appears in the projections to the Kinnersley tetrad of the velocity vn1 and vmb1*)
  Lz1 = 0;
  
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Geodesic BL frequencies *)
  
  rg[wr_] := rgfun[wr,a,p,e,x,{En,Lz,Kc}]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  zg[wz_] := zgfun[wz,a,p,e,x,{En,Lz,Kc}];
  
  \[CapitalGamma]  = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalOmega]p = orbitCorrection["BLPrecessionFrequency"];
  \[Omega]  = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta] - \[CapitalOmega]p; (* Geodesic frequency and the linear correction *)
 
  {\[Lambda],SWSH} = angparLeading[-2,l,m,a,SetPrecision[\[Omega],48]]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  R = RLeading[-2,l,m,a,\[Omega],\[Lambda],e,p];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus0 = sumPlus1 = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[Pi  ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]trg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]rg"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[Pi/4]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbitCorrection["\[CapitalDelta]tzg"]'[0   ]-m*orbitCorrection["\[CapitalDelta]\[Phi]zg"]'[0   ]+k)]],32];
                
  wrlist=Table[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]],{ir,1,stepsr}];
  wzlist=Table[N[(iz-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]],{iz,1,steps\[Theta]}];
                   
  (*Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in wz"];*)
  \[Theta]list = {};(* List for functions of \[Theta] *)
  (*Print["Sampling the spin corrections to the EoM"];*)
If[stepsr<=steps\[Theta],
   \[Delta]rort=Table[orbitCorrection["\[Delta]rortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]zort=Table[orbitCorrection["\[Delta]zortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[CapitalDelta]tsort=Table[orbitCorrection["\[CapitalDelta]\[Delta]tortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[CapitalDelta]\[Phi]sort=Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]ortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]Vrort=Table[orbitCorrection["\[Delta]vrortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]Vzort=Table[orbitCorrection["\[Delta]vzortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]Vtort=Table[orbitCorrection["\[Delta]vtortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   \[Delta]V\[Phi]ort=Table[orbitCorrection["\[Delta]v\[Phi]ortplus"][wrlist[[ir]], wzlist],{ir,Length[wrlist]}];
   ,
   \[Delta]rort=Transpose[Table[orbitCorrection["\[Delta]rortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]zort=Transpose[Table[orbitCorrection["\[Delta]zortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[CapitalDelta]tsort=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]tortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[CapitalDelta]\[Phi]sort=Transpose[Table[orbitCorrection["\[CapitalDelta]\[Delta]\[Phi]ortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]Vrort=Transpose[Table[orbitCorrection["\[Delta]vrortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]Vzort=Transpose[Table[orbitCorrection["\[Delta]vzortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]Vtort=Transpose[Table[orbitCorrection["\[Delta]vtortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   \[Delta]V\[Phi]ort=Transpose[Table[orbitCorrection["\[Delta]v\[Phi]ortplus"][wrlist, wzlist[[iz]]],{iz,Length[wzlist]}]];
   ];
   (*Print["Computing the corrections to the amplitudes - orthogonal component in the spin"];*)  
  For[ir = 1, ir <= stepsr, ir++,(* integration over w_r *)
    wrsub = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = rg[wrsub];
    Ur = Sign[Pi-Mod[wrsub,2Pi]]*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];(* Geodesic radial velocity *)
    Utr = En*(rp^2+a^2)^2/(rp^2-2rp+a^2)-2a*rp*Lz/(rp^2-2rp+a^2);
    U\[Phi]r = a/(rp^2-2rp+a^2)(En(rp^2+a^2)-a*Lz)-a*En;
    expr = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]trg"][wrsub])-m*(orbitCorrection["\[CapitalDelta]\[Phi]rg"][wrsub])+2Pi*n*(ir-1/2)/stepsr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2*rp*\[Omega];
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RInrp    = R[rp]["R0"]["In"]["R"]; (* radial function *)
    dRInrp   = R[rp]["R0"]["In"]["dR"];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R[rp]["R0"]["Up"]["R"];
    dRUprp   = R[rp]["R0"]["Up"]["dR"];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
       
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta], i\[Theta]++, (* integration over w_\[Theta] *)
      w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta], Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        zp = zg[w\[Theta]];
        Uz = dzgd\[Lambda]fun[w\[Theta],a,p,e,x,{En,Lz,Kc}];(* Polar geodesic velocity *)
        Utz = -En*a^2(1-zp^2);
        U\[Phi]z = Lz/(1-zp^2);
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["\[CapitalDelta]tzg"][w\[Theta]])-m*(orbitCorrection["\[CapitalDelta]\[Phi]zg"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2]; (*\[Theta] is  only defined between [0,\[Pi]), thus Sin(\[Theta]) is semi-positive.*)
        S = SWSH[ArcCos[zp]][[1]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
        dSd\[Theta] = SWSH[ArcCos[zp]][[2]]; (* First derivative of S wrt \[Theta] *)
        d2Sd\[Theta]2 = SWSH[ArcCos[zp]][[3]]; (* second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta],Utz,U\[Phi]z}];,
        {zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta],Utz,U\[Phi]z}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      (*{fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,K1,dK1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];*)
      vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      
      \[Psi]=\[Psi]rfun[wrsub,a,p,e,x,{En,Lz,Kc}]+\[Psi]zfun[w\[Theta],a,p,e,x,{En,Lz,Kc}];
      
      Sln  = (Cos[\[Psi]]+I*Sin[\[Psi]])a*zp*Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]/(2Sqrt[Kc]*\[CapitalSigma]); (* Spin tensor in Kinnersley tetrad *) (*Cos[\[Psi]] and Sin[\[Psi]] comes from expanding Cos[\[Psi]+wp] and Sin[\[Psi]+wp] *)
      Snmpar  = (\[Zeta]/Sqrt[Kc])*vm*vn;
      Snmbpar = Conjugate[Snm];
      Snm=-((Cos[\[Psi]]+I*Sin[\[Psi]])(Kc-I a*rp*zp)+(-I*Cos[\[Psi]]+Sin[\[Psi]])Sqrt[Kc]\[Zeta]bar)I*Snmpar/(2Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]);
      Snmb=((Cos[\[Psi]]+I*Sin[\[Psi]])(Kc+I a*rp*zp)+(-I*Cos[\[Psi]]+Sin[\[Psi]])Sqrt[Kc]\[Zeta])I*Snmbpar/(2Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]);
      Slmbpar = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
      Slmb=-((Cos[\[Psi]]+I*Sin[\[Psi]])(Kc-I a*rp*zp)-(-I*Cos[\[Psi]]+Sin[\[Psi]])Sqrt[Kc]\[Zeta]bar)I*Slmbpar/(2Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]);
      Smmb = (Cos[\[Psi]]+I*Sin[\[Psi]])I*rp*Sqrt[Kc+rp^2]Sqrt[Kc-a^2*zp^2]/(2Sqrt[Kc]*\[CapitalSigma]);
      
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
        
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      
      rp1 =  \[Delta]rort[[ir,i\[Theta]]]; (* Corrections to the coordinates and four-velocity for each quadrant *)
      zp1 =  \[Delta]zort[[ir,i\[Theta]]];
      Urp1 = \[Delta]Vrort[[ir,i\[Theta]]];
      Uzp1 = \[Delta]Vzort[[ir,i\[Theta]]];
      Utp1 = \[Delta]Vtort[[ir,i\[Theta]]];
      U\[Phi]p1 = \[Delta]V\[Phi]ort[[ir,i\[Theta]]];
      
      Ut = Utr+Utz;
      U\[Phi] = U\[Phi]r+U\[Phi]z;
      
      \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
      
      exp1 = I*(\[Omega]*\[CapitalDelta]tsort[[ir,i\[Theta]]]-m*\[CapitalDelta]\[Phi]sort[[ir,i\[Theta]]]);(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      
      vn1  = -((Utp1*\[CapitalDelta])/(2\[CapitalSigma]^2))+(a*U\[Phi]p1(1-zp^2)*\[CapitalDelta])/(2\[CapitalSigma]^2)-Urp1/(2\[CapitalSigma])+rp1/(2\[CapitalDelta]*\[CapitalSigma]^3)((4rp*(Ut-a*U\[Phi](1-zp^2))\[CapitalDelta]^2+2(rp*(Ur-Ut)+Ut)\[CapitalDelta]*\[CapitalSigma]+2a(-1+rp)U\[Phi](1-zp^2)\[CapitalDelta]*\[CapitalSigma]))
               +zp1/(2\[CapitalDelta]*\[CapitalSigma]^3) ((4a^2zp(Ut-a*U\[Phi](1-zp^2))\[CapitalDelta]^2+2a*zp*\[CapitalDelta](a*Ur-U\[Phi]*\[CapitalDelta])\[CapitalSigma])); (* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = -((I(a^2+rp^2)U\[Phi]p1 (1-zp^2)^(1/2))/(Sqrt[2]\[Zeta]*\[CapitalSigma]))+(I*a*Utp1*(1-zp^2)^(1/2))/(Sqrt[2]\[Zeta]*\[CapitalSigma])-(Uzp1*\[Zeta]bar)/(Sqrt[2](1-zp^2)^(1/2)\[CapitalSigma])
               +rp1/(Sqrt[2](1-zp^2)^(3/2)\[Zeta]^2*\[CapitalSigma]^2)(2rp(-1+zp^2)\[Zeta](I(-a*Ut+(a^2+rp^2)U\[Phi])(-1+zp^2)-Uz*\[Zeta]*\[Zeta]bar)+I(-a*Ut+(a^2+rp^2)*U\[Phi]) (-1+zp^2)^2*\[CapitalSigma]-2I*rp*U\[Phi] (1-zp^2)^2\[Zeta]*\[CapitalSigma]+Uz(-1+zp^2)\[Zeta]^2*\[CapitalSigma])
               +zp1/(Sqrt[2](1-zp^2)^(3/2)\[Zeta]^2*\[CapitalSigma]^2)(2a^2*zp*(-1+zp^2)\[Zeta](I(-a*Ut+(a^2+rp^2)*U\[Phi])(-1+zp^2)-Uz*\[Zeta]*\[Zeta]bar)+a(-a*Ut+(a^2+rp^2)*U\[Phi])(1-zp^2)^2\[CapitalSigma]+I*a*Uz(-1+zp^2)\[Zeta]^2*\[CapitalSigma]-2Uz*zp*\[Zeta]^2\[Zeta]bar*\[CapitalSigma]-zp*\[Zeta]*(I(-a*Ut+(a^2+rp^2)U\[Phi])(-1+zp^2)-Uz*\[Zeta]*\[Zeta]bar)*\[CapitalSigma]);
                          
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      
      A2g = Cmmbmb*fmbmb2;
      A1g = (Cmnmb*fnmb1+Cmmbmb*fmbmb1);
      A0g = (Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0);
      
      (*These coefficients are different respect to the ones define in the companion paper. Here I put together everthing, including the differentials \[Delta]r, \[Delta]z*)
      B3ort = -AmbmbrS*fmbmb2;
      B2ort = -AnmbrS*fnmb1-AmbmbrS*fmbmb1;
      B1ort = -AmbmbrS*fmbmb0-AnmbrS*fnmb0-AnnrS fnn0;
      
      A2ort = AmbmbrS*dfmbmb2dr+Ambmb\[Theta]S*dfmbmb2d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb2;
      A1ort = AmbmbrS*dfmbmb1dr+Ambmb\[Theta]S*dfmbmb1d\[Theta]+AnmbrS*dfnmb1dr+Anmb\[Theta]S*dfnmb1d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)fmbmb1+(Anmb0S+Anmbt\[Phi]S)fnmb1;
      A0ort = AmbmbrS*dfmbmb0dr+Ambmb\[Theta]S*dfmbmb0d\[Theta]+AnmbrS*dfnmb0dr+Anmb\[Theta]S*dfnmb0d\[Theta]+AnnrS*dfnn0dr+Ann\[Theta]S*dfnn0d\[Theta]+(Ambmb0S+Ambmbt\[Phi]S)*fmbmb0+(Anmb0S+Anmbt\[Phi]S)*fnmb0+(Ann0S+Annt\[Phi]S)*fnn0;
      
      sumPlus0  += Total[\[CapitalSigma]*((A0g*RInrp - A1g*dRInrp + A2g*ddRInrp)*expr*exp\[Theta])]; (* Total of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*((A0g*RUprp - A1g*dRUprp + A2g*ddRUprp)*expr*exp\[Theta])];
      
      sumPlus1  += Total[\[CapitalSigma]*((A0ort*RInrp - (A1ort+B1ort)*dRInrp + (A2ort+B2ort)ddRInrp - B3ort*dddRInrp)*expr*exp\[Theta])]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*((A0ort*RUprp - (A1ort+B1ort)*dRUprp + (A2ort+B2ort)ddRUprp - B3ort*dddRUprp)*expr*exp\[Theta])];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/\[CapitalDelta]; (* Invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1)/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1)/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  
  Association[
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]g"->\[Omega],
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  ]
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];

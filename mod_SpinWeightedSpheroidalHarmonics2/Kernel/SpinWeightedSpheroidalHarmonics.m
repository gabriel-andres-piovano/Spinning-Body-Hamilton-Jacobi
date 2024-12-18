(* ::Package:: *)

(* ::Title:: *)
(*SpinWeightedSpheoridalHarmonics package*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["SpinWeightedSpheroidalHarmonics`"];


(* ::Subsection::Closed:: *)
(*Usage messages*)


SpinWeightedSphericalHarmonicY::usage = "SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], \[Phi]] gives the spin-weighted spherical harmonic with of spin-weight s, degree l and order m.";
SpinWeightedSpheroidalHarmonicS::usage = "SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], \[Theta], \[Phi]] gives the spin-weighted oblate spheroidal harmonic with spheroidicity \[Gamma], spin-weight s, degree l and order m.";
SpinWeightedSpheroidalHarmonicSFunction::usage = "SpinWeightedSpheroidalHarmonicSFunction[s, l, m, \[Gamma], data, method] represents a function for computing numerical values of spin-weighted spheroidal harmonics.";
SpinWeightedSpheroidalEigenvalue::usage = "SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]] gives the spin-weighted oblate spheroidal eigenvalue with spheroidicity \[Gamma], spin-weight s, degree l and order m.";


(* ::Subsection::Closed:: *)
(*Error Messages*)


SpinWeightedSpheroidalEigenvalue::numterms = "Automatic determination of the number of terms to use in the SphericalExpansion method may be unreliable in certain cases. Currently using `1` terms. It is recommended to verify the result is unchanged with a different number of terms.";
SpinWeightedSpheroidalEigenvalue::optx = "Unknown options in `1`";
SpinWeightedSpheroidalEigenvalue::params = "Invalid parameters s=`1`, l=`2`, m=`3`";
SpinWeightedSpheroidalHarmonicS::numterms = "Automatic determination of the number of terms to use in the SphericalExpansion method may be unreliable in certain cases. Currently using `1` terms. It is recommended to verify the result is unchanged with a different number of terms.";
SpinWeightedSpheroidalHarmonicS::optx = "Unknown options in `1`";
SpinWeightedSpheroidalHarmonicS::normprec = "Normalisation cannot be determined for \"Leaver\" method so it will be set to 1. To obtain an accurate result either provide a higher precision input spheroidicity or use the \"SphericalExpansion\" method instead.";
SpinWeightedSpheroidalHarmonicS::prec = "Spin-weighted spheroidal harmonic cannot be computed using \"Leaver\" method with the given working precision. To obtain an accurate result either provide a higher precision input spheroidicity or use the \"SphericalExpansion\" method instead.";
SpinWeightedSpheroidalHarmonicS::params = "Invalid parameters s=`1`, l=`2`, m=`3`";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];



(* ::Section::Closed:: *)
(*Useful functions*)


(* ::Subsection::Closed:: *)
(*Series expansion coefficients*)


\[Alpha][s_,l_,m_]:=1/l (Sqrt[l^2-m^2] Sqrt[l^2-s^2])/(Sqrt[2l-1] Sqrt[2l+1]);
\[Alpha][s_,l_,m_]/;l<Abs[s]=0;
\[Alpha][0,l_,m_]:= Sqrt[l^2-m^2]/(Sqrt[2l-1] Sqrt[2l+1]);
\[Beta][s_,l_,m_]:=-((m s)/(l(l+1)));
\[Beta][s_,l_,m_]/;l<Abs[s]=0;
\[Beta][0,l_,m_]=0;
s\[Lambda]lm[s_,l_,m_][0]:=l(l+1)-s(s+1); 
s\[Lambda]lm[s_,l_,m_][i_Integer?Positive]:=-(1/d[s,l,m][0,0])(Sum[d[s,l,m][i-k,0] s\[Lambda]lm[s,l,m][k],{k,1,i-1}]+\[Alpha][s,-1+l,m] \[Alpha][s,l,m] d[s,l,m][-2+i,-2]+\[Alpha][s,l,m] (\[Beta][s,-1+l,m]+\[Beta][s,l,m]) d[s,l,m][-2+i,-1]+(-1+\[Alpha][s,l,m]^2+\[Alpha][s,1+l,m]^2+\[Beta][s,l,m]^2) d[s,l,m][-2+i,0]+\[Alpha][s,1+l,m] (\[Beta][s,l,m]+\[Beta][s,1+l,m]) d[s,l,m][-2+i,1]+\[Alpha][s,1+l,m] \[Alpha][s,2+l,m] d[s,l,m][-2+i,2]-2 s \[Alpha][s,l,m] d[s,l,m][-1+i,-1]+2 (m-s \[Beta][s,l,m]) d[s,l,m][-1+i,0]-2 s \[Alpha][s,1+l,m] d[s,l,m][-1+i,1](*+d[s,l,m][i,0] s\[Lambda]lm[s,l,m][0]-(l-s) (1+l+s) d[s,l,m][i,0]*))
d[s_,l_,m_][i_?Negative,j_]=0
d[s_,l_,m_][i_?NonNegative,j_]/;Abs[j]>i=0
d[s_,l_,m_][0,0]=1;
d[s_,l_,m_][1,1]:=-((s \[Alpha][s,1+l,m])/(1+l));
d[s_,l_,m_][1,0]=0;
d[s_,l_,m_][1,-1]:=(s \[Alpha][s,l,m])/l;
d[0,l_,m_][1,1]:=0;
d[0,l_,m_][1,-1]:=0;
d[s_,l_,m_][i_Integer?Positive,0]:=- 1/(2d[s,l,m][0,0]) Sum[d[s,l,m][i-k,j] d[s,l,m][k,j],{k,1,i-1},{j,-i,i}];
d[s_,l_,m_][i_Integer?Positive,j_Integer]:=1/(j (1+j+2 l)) (Sum[d[s,l,m][i-k,j] s\[Lambda]lm[s,l,m][k],{k,1,i-1}]+\[Alpha][s,-1+j+l,m] \[Alpha][s,j+l,m] d[s,l,m][-2+i,-2+j]+\[Alpha][s,j+l,m] (\[Beta][s,-1+j+l,m]+\[Beta][s,j+l,m]) d[s,l,m][-2+i,-1+j]+(-1+\[Alpha][s,j+l,m]^2+\[Alpha][s,1+j+l,m]^2+\[Beta][s,j+l,m]^2) d[s,l,m][-2+i,j]+\[Alpha][s,1+j+l,m] (\[Beta][s,j+l,m]+\[Beta][s,1+j+l,m]) d[s,l,m][-2+i,1+j]+\[Alpha][s,1+j+l,m] \[Alpha][s,2+j+l,m] d[s,l,m][-2+i,2+j]-2 s \[Alpha][s,j+l,m] d[s,l,m][-1+i,-1+j]+2 (m-s \[Beta][s,j+l,m]) d[s,l,m][-1+i,j]-2 s \[Alpha][s,1+j+l,m] d[s,l,m][-1+i,1+j]);
d[s_,l_,m_][i_Integer?Positive,j_Integer]/;l+j<Abs[s]:=0;

simplify[expr_] := Collect[expr, {HoldPattern[\[Alpha][__]], HoldPattern[\[Beta][__]]}, Simplify];


(* ::Subsection::Closed:: *)
(*Functions for spectral method for seed in numerical evaluation*)


kHat[s_, 0, 0, \[Gamma]_] := \[Gamma]^2/3;
kHat[s_, l_, m_, \[Gamma]_] := -l (1+l)+(2 m s^2 \[Gamma])/(l+l^2)+1/3 (1+(2 (l+l^2-3 m^2) (l+l^2-3 s^2))/(l (-3+l+8 l^2+4 l^3))) \[Gamma]^2;
k2[s_, l_, m_, \[Gamma]_] := (Sqrt[((1+l-m) (2+l-m) (1+l+m) (2+l+m) (1+l-s) (2+l-s) (1+l+s) (2+l+s))/((1+2 l) (5+2 l))] \[Gamma]^2)/((1+l) (2+l) (3+2 l));
kTilde2[s_, 0, 0, \[Gamma]_] := -2 s Sqrt[(1-s^2)/3] \[Gamma];
kTilde2[s_, l_, m_, \[Gamma]_] := -((2 s Sqrt[((1+2 l+l^2-m^2) (1+2 l+l^2-s^2))/(3+8 l+4 l^2)] \[Gamma] (2 l+l^2+m \[Gamma]))/(l (2+3 l+l^2)));


(* ::Subsection::Closed:: *)
(*Continued fraction*)


(* ::Text:: *)
(*This is preferred over Mathematica's ContinedFractionK function as we can get an error estimate on the result using this function.*)


CF[a_, b_, {n_, n0_}] := 
  Module[{A, B, ak, bk, res = Indeterminate, j = n0},
   A[n0 - 2] = 1;
   B[n0 - 2] = 0;
   ak[k_] := ak[k] = (a /. n -> k);
   bk[k_] := bk[k] = (b /. n -> k);
   A[n0 - 1] = 0(*bk[n0-1]*);
   B[n0 - 1] = 1;
   A[k_] := A[k] = bk[k] A[k - 1] + ak[k] A[k - 2];
   B[k_] := B[k] = bk[k] B[k - 1] + ak[k] B[k - 2];
   While[res =!= (res = A[j]/B[j]), j++];
   Clear[A, B, ak, bk];
   res
   ];




(* ::Section:: *)
(*SpinWeightedSpheroidalEigenvalue*)


(* ::Subsection::Closed:: *)
(*Spherical expansion method*)


Options[SWSHEigenvalueSpectral] = {"NumTerms" -> Automatic};


SWSHEigenvalueSpectral[s_, l_, m_, \[Gamma]_, OptionsPattern[]]:=
 Module[{nDown, nUp, Matrix, Eigens, lmin},
  (* FIXME: Improve the estimate of nmax. It should depend on the accuarcy sought. *)
  Switch[OptionValue["NumTerms"],
   Automatic,
    nUp = Ceiling[Abs[3/2\[Gamma]]]+5;
    Message[SpinWeightedSpheroidalEigenvalue::numterms, nUp];,
   _Integer,
    nUp = OptionValue["NumTerms"];,
   _,
    Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> {"SphericalExpansion", opts}];
    Return[$Failed];
  ];
  lmin=Max[Abs[s],Abs[m]];
  nDown = Min[nUp, l-lmin];

  Matrix=SparseArray[
	{{i_,i_}:>kHat[s,l-nDown-1+i,m,\[Gamma]],
	{i_,j_}/;j-i==-2:>k2[s,l-nDown-3+i,m,\[Gamma]],
	{i_,j_}/;j-i==-1:>kTilde2[s,l-nDown+i-2,m,\[Gamma]],
	{i_,j_}/;j-i==1:>kTilde2[s,l-nDown+i-1,m,\[Gamma]],
	{i_,j_}/;j-i==2:>k2[s,l-nDown+i+-1,m,\[Gamma]]}
  ,{nUp+nDown+1,nUp+nDown+1}];

  (* To choose the eigenvalue corrsponding to the desired l, we assume that the real
     part of eigenvalue is a monotonic function of l. *)
  Eigens=-Sort[Quiet[Eigenvalues[Matrix],Eigenvalues::arhm]];

  Eigens[[-(nDown+1)]]-s(s+1)
];


(* ::Subsection:: *)
(*Leaver's method*)


Options[SWSHEigenvalueLeaver] = {"InitialGuess" -> "SphericalExpansion", "NumInversions" -> Automatic};


SWSHEigenvalueLeaver[s_, l_, m_, \[Gamma]_, opts:OptionsPattern[]] :=
 Module[{Myprec, Nmax, nInv, \[Alpha], \[Beta], \[Alpha]n, \[Beta]n, \[Gamma]n, n, LHS, RHS, Eq, A, Aval, Avar, Aini},
  Aini = OptionValue["InitialGuess"];
  If[Aini === "SphericalExpansion",
    Aini = Quiet[SetPrecision[SWSHEigenvalueSpectral[s, l, m, N[\[Gamma]]], Precision[\[Gamma]]], SpinWeightedSpheroidalEigenvalue::numterms];
  ];
  Switch[OptionValue["NumInversions"],
   Automatic,
    nInv = l-Max[Abs[m],Abs[s]];,
   _Integer,
    nInv = OptionValue["NumInversions"];,
   _,
    Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> {"Leaver", opts}];
    Return[$Failed];
  ];
  Myprec = Max[Precision[\[Gamma]], 3];
  \[Alpha] = Abs[m+s];
  \[Beta] = Abs[m-s];
  \[Alpha]n[n_] := (-4\[Gamma](n+\[Alpha]+1)(n+\[Beta]+1)(n+(\[Alpha]+\[Beta])/2+1+s ))/((2n+\[Alpha]+\[Beta]+2)(2n+\[Alpha]+\[Beta]+3));
  \[Beta]n[n_, A_] := A + s(s+1)+\[Gamma]^2 -(n +(\[Alpha]+\[Beta])/2)(n +(\[Alpha]+\[Beta])/2+1) + If[s!=0,(8m s^2 \[Gamma])/((2n+\[Alpha]+\[Beta])(2n+\[Alpha]+\[Beta]+2)),0];
  \[Gamma]n[n_] := (4\[Gamma] n(n+\[Alpha]+\[Beta])(n+(\[Alpha]+\[Beta])/2-s))/((2n+\[Alpha]+\[Beta]-1)(2n+\[Alpha]+\[Beta]));
  RHS[Ax_] := -CF[-\[Alpha]n[n-1] \[Gamma]n[n], \[Beta]n[n,Ax], {n, nInv+1}];
  LHS[Ax_] := \[Beta]n[nInv, Ax] + ContinuedFractionK[-\[Alpha]n[nInv-n] \[Gamma]n[nInv-n+1], \[Beta]n[nInv-n, Ax], {n, 1, nInv}];
  Eq[A_?NumericQ] := LHS[A] - RHS[A];
  Aval = Avar /. Quiet[Check[FindRoot[Eq[Avar]==0, {Avar, Aini}, AccuracyGoal -> Myprec-3, WorkingPrecision -> Myprec, Method -> "Secant"], Avar -> $Failed, {Power::infy, FindRoot::nlnum}], {Power::infy, FindRoot::nlnum}];
  Clear[\[Alpha]n, \[Beta]n, \[Gamma]n, LHS, RHS, Eq];
  Aval
];


(* ::Subsection::Closed:: *)
(*SpinWeightedSpheroidalEigenvalue (uses a combination of the spectral and Leaver's method)*)


SyntaxInformation[SpinWeightedSpheroidalEigenvalue] =
 {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};


Options[SpinWeightedSpheroidalEigenvalue] = {Method -> Automatic};


SetAttributes[SpinWeightedSpheroidalEigenvalue, {NumericFunction, Listable}];


SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_, OptionsPattern[]] /; l < Abs[s] || Abs[m] > l := 
 (Message[SpinWeightedSpheroidalEigenvalue::params, s, l, m]; $Failed);


SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_, OptionsPattern[]] /; \[Gamma] == 0 :=
  l(l+1) - s(s+1);


SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_?InexactNumberQ, OptionsPattern[]] /; AllTrue[{2s, 2l, 2m}, IntegerQ] := 
 Module[{opts, Aini, \[Lambda]},
  Switch[OptionValue[Method],
    "SphericalExpansion",
      \[Lambda] = SWSHEigenvalueSpectral[s, l, m, \[Gamma]] - 2 m \[Gamma] + \[Gamma]^2,
    {"SphericalExpansion", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHEigenvalueSpectral]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> OptionValue[Method]];
      ];
      \[Lambda] = SWSHEigenvalueSpectral[s, l, m, \[Gamma], opts] - 2 m \[Gamma] + \[Gamma]^2,
    Automatic | "Leaver",
      \[Lambda] = SWSHEigenvalueLeaver[s, l, m, \[Gamma]] - 2 m \[Gamma] + \[Gamma]^2;,
    {"Leaver", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHEigenvalueLeaver]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> OptionValue[Method]];
      ];
      \[Lambda] = SWSHEigenvalueLeaver[s, l, m, \[Gamma], opts] - 2 m \[Gamma] + \[Gamma]^2,
      _,
      Message[SpinWeightedSpheroidalEigenvalue::optx, Method -> OptionValue[Method]];
      \[Lambda] = $Failed;
  ];
  \[Lambda]
];


SpinWeightedSpheroidalEigenvalue /: N[SpinWeightedSpheroidalEigenvalue[s_Integer, l_Integer, m_Integer, \[Gamma]_?NumericQ, opts:OptionsPattern[]], Nopts___] :=
  SpinWeightedSpheroidalEigenvalue[s, l, m, N[\[Gamma], Nopts], opts];


(* ::Subsection::Closed:: *)
(*Small-\[Gamma] expansion*)


SpinWeightedSpheroidalEigenvalue /: 
  Series[SpinWeightedSpheroidalEigenvalue[s_, l_, m_, \[Gamma]_], {\[Gamma]_, 0, order_}] := 
 Module[{i, j, coeffs}, Internal`InheritedBlock[{d, s\[Lambda]lm}, Block[{\[Alpha], \[Beta]},
  Do[
    Do[
      d[s, l, m][i, j] = simplify[d[s, l, m][i, j]], {j, -i, i}]; 
    s\[Lambda]lm[s, l, m][i] = simplify[s\[Lambda]lm[s, l, m][i]];
  , {i, 0, order}];
  coeffs = Table[s\[Lambda]lm[s, l, m][i], {i, 0, order}];
  SeriesData[\[Gamma], 0, coeffs, 0, order + 1, 1]
]]];


(* ::Subsection::Closed:: *)
(*Large-\[Gamma] expansion*)


(* ::Text:: *)
(*For details see "High-order asymptotics for the Spin-Weighted Spheroidal Equation at large real frequency", M. Casals, A. C. Ottewill, N. Warburton, arXiv:1810.00432*)


SWSHq[s_?NumberQ, l_?NumberQ, m_?NumberQ] := Module[{slm,z0},
	slm = Abs[m+Abs[s]] + Abs[s];
	z0  = If[EvenQ[l+m], 0, 1];

	If[l  >= slm, l+1-z0,  2l+1 - slm ]
]


SpinWeightedSpheroidalEigenvalue /: 
  Series[SpinWeightedSpheroidalEigenvalue[s_?NumberQ, l_?NumberQ, m_?NumberQ, -\[Gamma]_], {\[Gamma]_, \[Infinity], order_}]:= Series[SpinWeightedSpheroidalEigenvalue[s, l, -m, \[Gamma]], {\[Gamma], \[Infinity], order}]


SpinWeightedSpheroidalEigenvalue /: 
  Series[SpinWeightedSpheroidalEigenvalue[s_?NumberQ, l_?NumberQ, m_?NumberQ, \[Gamma]_], {\[Gamma]_, \[Infinity], order_}] :=
Module[{slm,z0,q,aFgen,AFgen,Asgen,\[Delta]gen,\[Nu]gen,RecRelgen,n,c,p,Serngen,Asum,sElm,s\[Lambda]lm,res},
  aFgen[0,_]=0;
  aFgen[0,0]=1;

  AFgen[0]=1;
  AFgen[r_]=AFgen[0]Sum[aFgen[r,k]c^-k,{k,Abs[r], order+1}];
  
  Asgen=Sum[\[Delta]gen[k]c^-k, {k, 1, order+1}];
  \[Nu]gen=1/2 (-1+SWSHq[s,l,m]-s-Abs[m+s]);
  RecRelgen=1/2  (1+n+\[Nu]gen+Abs[m+s]) (2(1+n+\[Nu]gen+s)-Abs[m-s]+Abs[m+s])AFgen[n+1]+(-1-m-m^2-2 (n+ \[Nu]gen)+4 c n-2 (n+\[Nu]gen)^2-s-m s-2 (n+\[Nu]gen) s-4 c Asgen-(1+2 n+2 \[Nu]gen+s) Abs[m+s])AFgen[n] + 1/2  (n+\[Nu]gen) (2 (n+\[Nu]gen+s)+Abs[m-s]+Abs[m+s])AFgen[n-1];

  Serngen[nn_]:=Serngen[nn]=Collect[CoefficientList[Normal[Series[RecRelgen/.n->nn,{c, \[Infinity], order+1}]],c^-1],aFgen[_,_]];

  Do[
    \[Delta]gen[p]=\[Delta]gen[p]/.Solve[Serngen[0][[p]]==0,\[Delta]gen[p]][[1]];
    Do[If[k!=0,aFgen[k,p]=aFgen[k,p]/.Solve[Serngen[k][[p]]==0,aFgen[k,p]][[1]]],{k, -p, p}]
  ,{p, 1, order+1}];
  
  Asum = 4\[Gamma] Sum[\[Delta]gen[k]\[Gamma]^-k,{k, 2, order+1}];

  sElm = -\[Gamma]^2+ 2 SWSHq[s,l,m] \[Gamma]-1/2 (SWSHq[s,l,m]^2-m^2-2s^2+1)+Asum;

  s\[Lambda]lm = sElm-s(s+1)+\[Gamma]^2-2m \[Gamma];
  
  res = SeriesData[\[Gamma], \[Infinity], Reverse[CoefficientList[Expand[s\[Lambda]lm \[Gamma]^order],\[Gamma]]], -1, order+1, 1];

  Remove[aFgen, AFgen, \[Delta]gen, Serngen, c, n];

  res
]


(* ::Section:: *)
(*SpinWeightedSpheroidalHarmonicS*)


(* ::Subsection::Closed:: *)
(*Spherical expansion method*)


Options[SWSHSSpectral] = {"NumTerms" -> Automatic};


SWSHSSpectral[s_Integer, l_Integer, m_Integer, \[Gamma]_, opts:OptionsPattern[]] :=
 Module[{lmin, nUp, nDown, A, esys,evec,eval,sign,pos},
  (* FIXME: Improve the estimate of nmax. It should depend on the accuarcy sought. *)
  Switch[OptionValue["NumTerms"],
   Automatic,
    nUp = Ceiling[Abs[3/2\[Gamma]]]+5;
    Message[SpinWeightedSpheroidalHarmonicS::numterms, nUp];,
   _Integer,
    nUp = OptionValue["NumTerms"];,
   _,
    Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> {"SphericalExpansion", opts}];
    Return[$Failed];
  ];
  lmin = Max[Abs[s],Abs[m]];
  nDown = Min[l-lmin,nUp];

  A = SparseArray[
        {{i_,i_} :> kHat[s, l-nDown-1+i, m, \[Gamma]],
         {i_,j_} /; j-i==-2 :> k2[s, l-nDown-3+i, m, \[Gamma]],
         {i_,j_} /; j-i==-1 :> kTilde2[s, l-nDown+i-2, m, \[Gamma]],
         {i_,j_} /; j-i==1 :> kTilde2[s, l-nDown+i-1, m, \[Gamma]],
         {i_,j_} /; j-i==2 :> k2[s, l-nDown+i+-1, m, \[Gamma]]},
        {nUp+nDown+1, nUp+nDown+1}];
  esys = Quiet[Eigensystem[A], Eigenvalues::arhm];
  eval = -Sort[esys[[1]]][[-(nDown+1)]];
  pos  = Position[esys[[1]], -eval][[1]];
  evec = First[esys[[2,pos]]];

  sign=Sign[evec[[Min[l-lmin+1,(nUp+nDown)/2+1]]]];
  SpinWeightedSpheroidalHarmonicSFunction[s, l, m, \[Gamma], {"SphericalExpansion", sign*evec, nDown, nUp}]
];


(* ::Subsection:: *)
(*Leaver's method*)


Options[SWSHSLeaver] = {"NumTerms" -> Automatic};


SWSHSLeaver[s_Integer, l_Integer, m_Integer, \[Gamma]_, opts:OptionsPattern[]] :=
 Module[{\[Lambda], k1, k2, \[Alpha]n, \[Beta]n, \[Gamma]n, an, n, sign, norm, anTab, nmin, nmax, nprec, normterm, prec = Precision[\[Gamma]]},
  nmin = 0;
  nmax = 1;
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, \[Gamma]];

  (* Leaver's recurrence formula, Eq. 20 of Leaver 1985 *)
  k1 = Abs[m-s]/2;
  k2 = Abs[m+s]/2;
  \[Alpha]n[n_] := -2(n+1)(n+2k1+1);
  \[Beta]n[n_] := n(n-1)+2n(k1+k2+1-2\[Gamma])-(2\[Gamma](2k1+s+1)-(k1+k2)(k1+k2+1))-(s(s+1)+\[Lambda]+2m \[Gamma]);
  \[Gamma]n[n_] := 2\[Gamma](n+k1+k2+s);
  an[nmin] = 1;
  an[nmin+1] = -\[Beta]n[nmin]/\[Alpha]n[nmin];
  an[n_] := an[n] = (-\[Beta]n[n-1]an[n-1]-\[Gamma]n[n-1]an[n-2])/\[Alpha]n[n-1];

  (* Compute coefficients until we reach the desired tolerance for computing the norm *)
  normterm[i_] := 2^i Pochhammer[i+2 (1+k1+k2), -2k2-1] Hypergeometric1F1[1+i+2 k1, i+2 (1+k1+k2), 4 \[Gamma]] Sum[an[j]an[i-nmin-j], {j, 0, i-nmin}];
  Switch[OptionValue["NumTerms"],
   Automatic,
    norm = normterm[nmin];
    While[norm != (norm += normterm[nmax]), nmax++];,
   _Integer,
    nmax = OptionValue["NumTerms"] - 1;
    norm = Sum[normterm[i], {i, nmin, nmax}];,
   _,
    Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> {"Leaver", opts}];
    Return[$Failed];
  ];
 (* This is "forward" algorithm to solve the recurrence relation. The first terms are computed with high precision, but the algorithm is unstable,
  and latter terms have few digits of precision.
 *)
 anTab = Table[an[n], {n, nmin, nmax}];

  (* Normalisation such that \[Integral]\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubsuperscriptBox[\(S\), \(lm\), \(*\)]\)\)(\[Theta],\[Phi];\[Gamma])\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubscriptBox[\(S\), \(l'm'\)]\)\)(\[Theta],\[Phi];\[Gamma])d\[CapitalOmega] = Subscript[\[Delta], ll']Subscript[\[Delta], mm'] *)
  norm = Sqrt[2\[Pi]] (2^(1+2 k1+2 k2) E^(-2 \[Gamma]) Gamma[1+2 k2] norm)^(1/2);
  If[Precision[norm] == 0.,
    Message[SpinWeightedSpheroidalHarmonicS::normprec];
    norm = 1;
  ];

  (* Overall sign such that the \[Gamma]\[Rule]0 limit is continuous *)
  sign = If[(OddQ[l]&&EvenQ[m+s])||(OddQ[l]&&OddQ[m+s]&&m>=s)||(EvenQ[l]&&OddQ[m-s]&&m<=s), -1, 1];
  
  (*We use here the (3.9) algorithm described in "Computational aspects of three-term recurrence relations" by W. Gautschi.
    This algorithm solved the recurrence relation backward, starting from a large n (called nmax) and then proceding to lower n.
    Terms with large n are computed with much higher precision compared to a forward algorithm.
    However, the Gautschi algorithm returns the first few terms of re with much lower precision, and it requires an estimate for nmax.
    In this section we combined the advantanges of the two algorithms: nmax is provided by the "while loop", and  the first nprec 
    terms (see later for definition) are given by the forward algorithmm directly. The remaining nmax-nprec are computed with the Gautschi algorithm.
  *)
  Block[{fs,rs},
  (*Number of high precision elements obtained by solving the three-terms recurrence relation directly. *)
   nprec = Length[Position[anTab,_?(Precision[#]>(prec-3)&)]];
   
   rs[n_]:= rs[n]= -((\[Gamma]n[n+1]/\[Alpha]n[n+1])/((\[Beta]n[n+1]/\[Alpha]n[n+1])+rs[n+1]));
   rs[nmax] = 0;
   Table[rs[j],{j,nmax,nprec,-1}]; 
   fs[n_] := fs[n] = rs[n-1]fs[n-1];
   Table[fs[i] = an[i],{i,nmin,nprec}]; (*The first terms of this recurrence relation are given by the forward algorithm, si  *)
   anTab = Table[fs[j],{j,nmin,nmax}];
  ];
  
  (* Clear Symbols with DownValues to prevent memory leaks *)
  Clear[\[Alpha]n, \[Beta]n, \[Gamma]n, an, normterm];

  (* Return a SpinWeightedSpheroidalHarmonicSFunction which can be evaluated for arbitratry \[Theta], \[Phi] *)
  SpinWeightedSpheroidalHarmonicSFunction[s, l, m, \[Gamma], {"Leaver", anTab, nmin, nmax}]
];


(* ::Subsection::Closed:: *)
(*SpinWeightedSpheroidalHarmonicS*)


SyntaxInformation[SpinWeightedSpheroidalHarmonicS] =
 {"ArgumentsPattern" -> {_, _, _, _, ___}};


Options[SpinWeightedSpheroidalHarmonicS] = {Method -> Automatic};


SetAttributes[SpinWeightedSpheroidalHarmonicS, {NumericFunction, Listable}];


SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_, \[Theta]_, \[Phi]_, opts:OptionsPattern[]] :=
  SpinWeightedSpheroidalHarmonicS[s, l, m, \[Gamma], opts][\[Theta], \[Phi]];


SpinWeightedSpheroidalHarmonicS[s_, l_, m_, \[Gamma]_, OptionsPattern[]] /; l < Abs[s] || Abs[m] > l := 
 (Message[SpinWeightedSpheroidalHarmonicS::params, s, l, m]; $Failed);


SpinWeightedSpheroidalHarmonicS[s_, l_, m_, \[Gamma]_, OptionsPattern[]] /; \[Gamma] == 0 :=
  SpinWeightedSpheroidalHarmonicSFunction[s, l, m, 0, {"SphericalExact", {}, 0, 0}] /; OptionValue[Method] == Automatic;


SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, OptionsPattern[]] :=
 Module[{opts, Slm},
  Switch[OptionValue[Method],
    "Eigenvalue"|"SphericalExpansion",
      Slm = SWSHSSpectral[s, l, m, \[Gamma]],
    {"Eigenvalue"|"SphericalExpansion", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHSSpectral]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> OptionValue[Method]];
      ];
      Slm = SWSHSSpectral[s, l, m, \[Gamma], opts],
    Automatic | "Leaver",
      Slm = SWSHSLeaver[s, l, m, \[Gamma]];,
    {"Leaver", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SWSHSLeaver]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> OptionValue[Method]];
      ];
      Slm = SWSHSLeaver[s, l, m, \[Gamma], opts],
    _,
     Message[SpinWeightedSpheroidalHarmonicS::optx, Method -> OptionValue[Method]];
     Slm = $Failed;
  ];
  Slm
];


SpinWeightedSpheroidalHarmonicS /: N[SpinWeightedSpheroidalHarmonicS[s_Integer, l_Integer, m_Integer, \[Gamma]_?NumericQ, \[Theta]_?NumericQ, \[Phi]_?NumericQ, opts:OptionsPattern[]], Nopts:OptionsPattern[N]] :=
  SpinWeightedSpheroidalHarmonicS[s, l, m, N[\[Gamma], Nopts], \[Theta], \[Phi], opts];


(* ::Subsection::Closed:: *)
(*Small-\[Gamma] expansion*)


SpinWeightedSpheroidalHarmonicS /: 
  Series[SpinWeightedSpheroidalHarmonicS[s_, l_, m_, \[Gamma]_, \[Theta]_, \[Phi]_], {\[Gamma]_, 0, order_}] :=
 Module[{i, j, coeffs}, Internal`InheritedBlock[{d, s\[Lambda]lm}, Block[{\[Alpha], \[Beta]},
  Do[
    Do[
      d[s, l, m][i, j] = simplify[d[s, l, m][i, j]], {j, -i, i}]; 
    s\[Lambda]lm[s, l, m][i] = simplify[s\[Lambda]lm[s, l, m][i]];
  , {i, 0, order}];
  coeffs = Table[Sum[d[s, l, m][i, j] SpinWeightedSphericalHarmonicY[s, l+j, m, \[Theta], \[Phi]], {j, -i, i}], {i, 0, order}];
  SeriesData[\[Gamma], 0, coeffs, 0, order + 1, 1]
]]];


(* ::Section:: *)
(*SpinWeightedSpheroidalHarmonicSFunction*)


SyntaxInformation[SpinWeightedSpheroidalHarmonicSFunction] =
 {"ArgumentsPattern" -> {_, _, _, _, {_, {___}, _, _}}};


SetAttributes[SpinWeightedSpheroidalHarmonicSFunction, {NumericFunction}];


(* ::Subsection::Closed:: *)
(*Output format*)


SpinWeightedSpheroidalHarmonicSFunction /:
 MakeBoxes[swshf:SpinWeightedSpheroidalHarmonicSFunction[s_, l_, m_, \[Gamma]_, {method_String, coeffs_List, nDown_Integer, nUp_Integer}], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", s}], "  ",
                  BoxForm`SummaryItem[{"l: ", l}], "  ",
                  BoxForm`SummaryItem[{"m: ", m}], "  ",
                  BoxForm`SummaryItem[{"\[Gamma]: ", \[Gamma]}]}],
                  BoxForm`SummaryItem[{"Method: ", method}]
             };
  extended = {BoxForm`SummaryItem[{"Coefficients: ", Column[coeffs]}],
              BoxForm`SummaryItem[{"NumDown: ", nDown}],
              BoxForm`SummaryItem[{"NumUp: ", nUp}]
              };
  BoxForm`ArrangeSummaryBox[
    SpinWeightedSpheroidalHarmonicSFunction,
    swshf,
    PolarPlot[{swshf[\[Theta],0],swshf[\[Theta],\[Pi]]}, {\[Theta],0,\[Pi]}, Axes->False, PlotRange->All, ImagePadding->All, PlotStyle->ColorData[97,1],
      ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]],
    summary,
    extended,
    form]
];


(* ::Subsection:: *)
(*Numerical evaluation*)


(* ::Subsubsection::Closed:: *)
(*SphericalExpansion method*)


SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {method_String, coeffs_List, nDown_Integer, nUp_Integer}][\[Theta]_?NumericQ, \[Phi]_?NumericQ] :=
  Sum[coeffs[[k+nDown+1]]SpinWeightedSphericalHarmonicY[s, l+k, m, \[Theta], 0], {k, -nDown, nUp}]Exp[I m \[Phi]] /; method == "SphericalExpansion";


(* ::Subsubsection:: *)
(*Leaver's method*)


SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {method_String, coeffs_List, nDown_Integer, nUp_Integer}][\[Theta]_?NumericQ, \[Phi]_?NumericQ] :=
 Module[{u = Cos[\[Theta]], k1 = Abs[m-s]/2, k2 = Abs[m+s]/2, oneplusu, oneminusu, res, normterm, norm,sign},
  oneplusu = 2 Cos[\[Theta]/2]^2;
  oneminusu = 2 Sin[\[Theta]/2]^2;
  
  (* Compute coefficients until we reach the desired tolerance for computing the norm *)
  normterm[i_] := 2^i Pochhammer[i+2 (1+k1+k2), -2k2-1] Hypergeometric1F1[1+i+2 k1, i+2 (1+k1+k2), 4 \[Gamma]] Sum[coeffs[[j]]coeffs[[i+2-nDown-j]], {j, 1, i-nDown+1}];
  norm = Total[Table[normterm[i], {i, nDown, nUp}], Method->"CompensatedSummation"]; (*Compensated summation improves precision here...*)
  
  (* Overall sign such that the \[Gamma]\[Rule]0 limit is continuous *)
  sign = If[(OddQ[l]&&EvenQ[m+s])||(OddQ[l]&&OddQ[m+s]&&m>=s)||(EvenQ[l]&&OddQ[m-s]&&m<=s), -1, 1];
  
  (* Normalisation such that \[Integral]\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubsuperscriptBox[\(S\), \(lm\), \(*\)]\)\)(\[Theta],\[Phi];\[Gamma])\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubscriptBox[\(S\), \(l'm'\)]\)\)(\[Theta],\[Phi];\[Gamma])d\[CapitalOmega] = Subscript[\[Delta], ll']Subscript[\[Delta], mm'] *)   
  norm = Sqrt[2\[Pi]] (2^(1+2 k1+2 k2) E^(-2 \[Gamma]) Gamma[1+2 k2] norm)^(1/2);
    
  (* Leaver's series solution, Eq. 18 of Leaver 1985 *) (*...and here.*)
  res = (sign/norm)E^(\[Gamma] u) If[k1==0, 1, oneplusu^k1] If[k2==0, 1, oneminusu^k2] Total[Table[coeffs[[i+1]]*oneplusu^i,{i,0,nUp-1}], Method->"CompensatedSummation"] Exp[I m \[Phi]];

  (* Print a warning message if the harmonic cannot be determined *)
  If[Precision[res] == 0. && Accuracy[res] < Accuracy[{\[Gamma], \[Theta], \[Phi]}], Message[SpinWeightedSpheroidalHarmonicS::prec]];

  res
] /; method == "Leaver";


(* ::Subsubsection::Closed:: *)
(*SphericalExact method*)


SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_, {"SphericalExact", {}, 0, 0}][\[Theta]_, \[Phi]_] /; \[Gamma] == 0 :=
  SpinWeightedSphericalHarmonicY[s,l,m,\[Theta],\[Phi]];


(* ::Subsubsection:: *)
(*Derivatives*)


Derivative[d1_,d2_][SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {method_String, coeffs_List, nDown_Integer, nUp_Integer}]][\[Theta]_?NumericQ, \[Phi]_?NumericQ] :=
 Module[{\[Theta]1,\[Phi]1, res},
  res = Sum[coeffs[[k+nDown+1]]D[SpinWeightedSphericalHarmonicY[s, l+k, m, \[Theta]1, 0],{\[Theta]1,d1}], {k, -nDown, nUp}]D[Exp[I m \[Phi]1],{\[Phi]1,d2}]/.{\[Theta]1->\[Theta],\[Phi]1->\[Phi]};
  Remove[\[Theta]1,\[Phi]1];
  res
] /; method == "SphericalExpansion";


Derivative[d1_,d2_][SpinWeightedSpheroidalHarmonicSFunction[s_Integer, l_Integer, m_Integer, \[Gamma]_?InexactNumberQ, {method_String, coeffs_List, nDown_Integer, nUp_Integer}]][\[Theta]_?NumericQ, \[Phi]_?NumericQ] :=
 Module[{\[Theta]1, \[Phi]1, u, k1 = Abs[m-s]/2, k2 = Abs[m+s]/2, oneplusu, oneminusu, res, normterm, norm, sign},
  u = Cos[\[Theta]1];
  oneplusu = 2 Cos[\[Theta]1/2]^2;
  oneminusu = 2 Sin[\[Theta]1/2]^2;
  
  (* Compute coefficients until we reach the desired tolerance for computing the norm *)
  normterm[i_] := 2^i Pochhammer[i+2 (1+k1+k2), -2k2-1] Hypergeometric1F1[1+i+2 k1, i+2 (1+k1+k2), 4 \[Gamma]] Sum[coeffs[[j]]coeffs[[i+2-nDown-j]], {j, 1, i-nDown+1}];
  norm = Total[Table[normterm[i], {i, nDown, nUp}], Method->"CompensatedSummation"]; (*Compensated summation improves precision here...*)
  
  (* Normalisation such that \[Integral]\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubsuperscriptBox[\(S\), \(lm\), \(*\)]\)\)(\[Theta],\[Phi];\[Gamma])\!\(
\(\*SubscriptBox[\(\[InvisiblePrefixScriptBase]\), \(s\)]\)
\(\*SubscriptBox[\(S\), \(l'm'\)]\)\)(\[Theta],\[Phi];\[Gamma])d\[CapitalOmega] = Subscript[\[Delta], ll']Subscript[\[Delta], mm'] *)   
  norm = Sqrt[2\[Pi]] (2^(1+2 k1+2 k2) E^(-2 \[Gamma]) Gamma[1+2 k2] norm)^(1/2); 
  
  (* Overall sign such that the \[Gamma]\[Rule]0 limit is continuous *)
  sign = If[(OddQ[l]&&EvenQ[m+s])||(OddQ[l]&&OddQ[m+s]&&m>=s)||(EvenQ[l]&&OddQ[m-s]&&m<=s), -1, 1]; (*...but not here, surprisingly.*)
  
  (* Leaver's series solution, Eq. 18 of Leaver 1985 *)
  res =(sign/norm) D[E^(\[Gamma] u) If[k1==0, 1, oneplusu^k1] If[k2==0, 1, oneminusu^k2] Total[Table[coeffs[[i+1]]*oneplusu^i,{i,0,nUp-1}]],{\[Theta]1,d1}] D[Exp[I m \[Phi]1],{\[Phi]1,d2}]/.{\[Theta]1->\[Theta],\[Phi]1->\[Phi]};
  Remove[\[Theta]1,\[Phi]1];
  res
] /; method == "Leaver";
  


(* ::Section::Closed:: *)
(*SpinWeightedSphericalHarmonicY*)


SyntaxInformation[SpinWeightedSphericalHarmonicY] =
 {"ArgumentsPattern" -> {_, _, _, _, _, ___}};


SetAttributes[SpinWeightedSphericalHarmonicY, {NumericFunction, Listable}];


SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] :=
  (-1)^m Sqrt[((l+m)!(l-m)!(2l+1))/(4\[Pi] (l+s)!(l-s)!)] Sum[Binomial[l-s,r] Binomial[l+s,r+s-m] (-1)^(l-r-s) If[(2 l - 2 r - s + m)==0, 1, Sin[\[Theta]/2]^(2 l - 2 r - s + m)] If[(2 r + s - m)==0, 1, Cos[\[Theta]/2]^(2 r + s - m)],{r,Max[m-s,0],Min[l-s,l+m]}]Exp[I m \[Phi]];


SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, \[Phi]_] /; Abs[m]>l = 0;


SpinWeightedSphericalHarmonicY[s_Integer, l_Integer, m_Integer, \[Theta]_, 0.] :=
  SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], 0];


SpinWeightedSphericalHarmonicY[0, l_, m_, \[Theta]_, \[Phi]_] :=
  SphericalHarmonicY[l, m, \[Theta], \[Phi]];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*End*)


End[]
EndPackage[];

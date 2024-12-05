(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["SpinningOrbit`"];


(* ::Text:: *)
(*If you make use of this package, please acknowledge "Piovano, Pantelidou, Mac Uilliam and Witzany" arXiv:2410.05769 (https://arxiv.org/abs/2410.05769v1 )*)


(* ::Text:: *)
(*IMPORTANT NOTE: the following package include all the functions to compute the spin-corrections to the orbits (trajectories and velocities), constants of motion and frequencies for the fixed turning points (or "FC") and the fixed constants of motion (or "DH") parametrizations. It includes the contributions to both the parallel and orthogonal components of the secondary spin, and the maps between the corrections to the orbits in the FC and DH parametrizations.*)


(* ::Text:: *)
(*PREREQUISITEs: none. The package is stand-alone.*)


KerrSpinOrbitCorrectionDHPar::usage = "KerrSpinOrbitCorrectionDHPar[a, p, e, x, nmax, kmax] calculates linear corrections to the orbits in the fixed turning points on average parametrization for the parallel component of the spin";


KerrSpinOrbitCorrectionFCPar::usage = "KerrSpinOrbitCorrectionFCPar[a, p, e, x, nmax, kmax] calculates linear corrections to the orbits in the fixed constants of motion parametrization for the parallel component of the spin";


KerrSpinOrbitCorrectionDHMapFC::usage = "KerrSpinOrbitCorrectionDHMapFC[a, p, e, x, nmax, kmax] calculates the map between the fixed turning points and fixed constants of motion parametrization ";


KerrSpinOrbitCorrectionOrt::usage = "KerrSpinOrbitCorrectionOrt[a, p, e, x, nmax, kmax] calculates linear corrections to the orbits for the orthogonal component of the spin ";


Begin["`Private`"];


(*ADD ME: add algorithms to automatic select radial and polar Fourier modes based on the sought precision*)


(*IMPORTANT: still to check for memory leakages*)


(*IMPORTANT2: there are several functions here that should be made available to use once the package is called (i.e. shifts turning points, shift frequencies, corrections trajectories and velocities.)*)


(*IMPORTANT3: probably easier to create subpackage. *)


(* ::Section::Closed:: *)
(*Geodesic quantities*)


(* ::Subsection::Closed:: *)
(*Constants of motion*)


EEgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,ddinc,ffinc,gginc,hhinc,\[Kappa]g,\[Epsilon]g,\[Rho]g,\[Eta]g,\[Sigma]g,\[CapitalDelta]r1g,\[CapitalDelta]r2g,EEgnum,EEgden,EEgdisc},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	ddinc[r_]:=(r^2-2r+a^2)(r^2+a^2*z1g^2);
	ffinc[r_]:=r^4+a^2(r+2)r+a^2 z1g^2 (r^2-2r+a^2);
	gginc[r_]:=2a*r;
	hhinc[r_]:=r(r-2)+(r^2-2r+a^2)z1g^2/(1-z1g^2);
	
	\[Kappa]g=ddinc[r2g]*hhinc[r1g]-ddinc[r1g]*hhinc[r2g];
	\[Epsilon]g=ddinc[r2g]*gginc[r1g]-ddinc[r1g]*gginc[r2g];
	\[Rho]g=ffinc[r2g]*hhinc[r1g]-ffinc[r1g]*hhinc[r2g];
	\[Eta]g=ffinc[r2g]*gginc[r1g]-ffinc[r1g]*gginc[r2g];
	\[Sigma]g=gginc[r2g]*hhinc[r1g]-gginc[r1g]*hhinc[r2g];
	
	\[CapitalDelta]r1g=(a^2-2 r1g+r1g^2);
	\[CapitalDelta]r2g=(a^2-2 r2g+r2g^2);
	
	EEgnum=1/(-1+z1g^2)^2 (8a^2(-1+z1g^2)(-r1g*r2g+a^2*z1g^2)(-r1g*r2g(r1g^2+r1g(-2+r2g)+(-2+r2g)r2g)+a^4 z1g^2-a^2r1g*r2g(1+z1g^2))+(r1g*r2g(r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)+a^4z1g^2(4+(-2+r1g+r2g)z1g^2)+a^2(-4 r1g*r2g+(r1g+r2g)(r1g^2+r2g^2)z1g^2))((-2+r1g)r1g(-2+r2g)r2g(r1g+r2g)+a^4*z1g^2(2+(-2+r1g+r2g) z1g^2)+a^2(-2 r1g*r2g+(r1g^3+r1g^2(-2+r2g)+r1g(-2+r2g)r2g+(-2+r2g)r2g^2)z1g^2)));
	
	EEgden=1/(-1+z1g^2)^2 (16a^2(-1+z1g^2)(-r1g*r2g+a^2*z1g^2)(-r1g*r2g(r1g^2+r1g r2g+r2g^2)+a^4*z1g^2-a^2*r1g*r2g (1+z1g^2))+(r1g*r2g(r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)+a^4*z1g^2(4+(-2+r1g+r2g)z1g^2)+a^2(-4 r1g*r2g+(r1g+r2g)(r1g^2+r2g^2)z1g^2))^2);
	
	EEgdisc=(8a^2*\[CapitalDelta]r1g(r1g+r2g)\[CapitalDelta]r2g(r1g^2+a^2 z1g^2)(-r1g*r2g+a^2*z1g^2)^3(r2g^2+a^2*z1g^2))/(-1+z1g^2)^3;
	
	Sqrt[(EEgnum-2RealSign[xg] Sqrt[EEgdisc])/EEgden]
]


Lzgfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzgzaux,ddinc,ffinc,gginc,hhinc},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	
	ddinc[r_]:=(r^2-2r+a^2) (r^2+a^2 z1g^2);
	ffinc[r_]:=r^4+a^2(r+2)r+a^2 z1g^2 (r^2-2r+a^2);
	gginc[r_]:=2a*r;
	hhinc[r_]:=r(r-2)+(r^2-2r+a^2)z1g^2/(1-z1g^2);
	
	Lzgzaux=Sqrt[(gginc[r2g]^2+hhinc[r2g]*ffinc[r2g])EEg^2-hhinc[r2g]*ddinc[r2g]];
	
	(-gginc[r2g] EEg+RealSign[xg] Lzgzaux)/hhinc[r2g]
]


KKgfun[a_,p_,e_,xg_]:=(Lzgfun[a,p,e,xg]-a EEgfun[a,p,e,xg])^2+(1-xg^2)(a^2(1-EEgfun[a,p,e,xg]^2)+Lzgfun[a,p,e,xg]^2/xg^2);


(* ::Subsubsection::Closed:: *)
(*Slightly modified version of the functions for the constants of motion *)


(* ::Text:: *)
(*These versions of the functions for the constants of motion is more useful for derivatives of the geodesics frequencies wrt to r1g r2g and z1g . The "sgn" in all the functions denotes the type of orbits. +1 is for prograde orbits, -1 for retrograde orbits. It is equivalent to the sign of x. *)


EEgfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{ddinc,ffinc,gginc,hhinc,\[Kappa]g,\[Epsilon]g,\[Rho]g,\[Eta]g,\[Sigma]g,\[CapitalDelta]r1g,\[CapitalDelta]r2g,EEgnum,EEgden,EEgdisc},
	ddinc[r_]:=(r^2-2r+a^2) (r^2+a^2 z1g^2);
	ffinc[r_]:=r^4+a^2(r+2)r+a^2 z1g^2 (r^2-2r+a^2);
	gginc[r_]:=2 a r;
	hhinc[r_]:=r(r-2)+(r^2-2r+a^2)z1g^2/(1-z1g^2);
	
	\[Kappa]g=ddinc[r2g]*hhinc[r1g]-ddinc[r1g]*hhinc[r2g];
	\[Epsilon]g=ddinc[r2g]*gginc[r1g]-ddinc[r1g]*gginc[r2g];
	\[Rho]g=ffinc[r2g]*hhinc[r1g]-ffinc[r1g]*hhinc[r2g];
	\[Eta]g=ffinc[r2g]*gginc[r1g]-ffinc[r1g]*gginc[r2g];
	\[Sigma]g=gginc[r2g]*hhinc[r1g]-gginc[r1g]*hhinc[r2g];
	
	\[CapitalDelta]r1g=(a^2-2 r1g+r1g^2);
	\[CapitalDelta]r2g=(a^2-2 r2g+r2g^2);
	
	EEgnum=1/(-1+z1g^2)^2 (8a^2(-1+z1g^2)(-r1g*r2g+a^2*z1g^2)(-r1g*r2g(r1g^2+r1g(-2+r2g)+(-2+r2g)r2g)+a^4 z1g^2-a^2r1g*r2g(1+z1g^2))+(r1g*r2g(r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)+a^4z1g^2(4+(-2+r1g+r2g)z1g^2)+a^2(-4 r1g*r2g+(r1g+r2g)(r1g^2+r2g^2)z1g^2))((-2+r1g)r1g(-2+r2g)r2g(r1g+r2g)+a^4*z1g^2(2+(-2+r1g+r2g) z1g^2)+a^2(-2 r1g*r2g+(r1g^3+r1g^2(-2+r2g)+r1g(-2+r2g)r2g+(-2+r2g)r2g^2)z1g^2)));
	
	EEgden=1/(-1+z1g^2)^2 (16a^2(-1+z1g^2)(-r1g*r2g+a^2*z1g^2)(-r1g*r2g(r1g^2+r1g r2g+r2g^2)+a^4*z1g^2-a^2*r1g*r2g (1+z1g^2))+(r1g*r2g(r1g^2(-2+r2g)+r1g(-2+r2g)r2g-2r2g^2)+a^4*z1g^2(4+(-2+r1g+r2g)z1g^2)+a^2(-4 r1g*r2g+(r1g+r2g)(r1g^2+r2g^2)z1g^2))^2);
	
	EEgdisc=(8a^2*\[CapitalDelta]r1g(r1g+r2g)\[CapitalDelta]r2g(r1g^2+a^2 z1g^2)(-r1g*r2g+a^2*z1g^2)^3(r2g^2+a^2*z1g^2))/(-1+z1g^2)^3;
	
	Sqrt[(EEgnum-2RealSign[sgn] Sqrt[EEgdisc])/EEgden]
]


Lzgfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{EEg,Lzgzaux,ddinc,ffinc,gginc,hhinc},
	EEg=EEgfunsgn[a,r1g,r2g,z1g,sgn];
	
	ddinc[r_]:=(r^2-2r+a^2) (r^2+a^2 z1g^2);
	ffinc[r_]:=r^4+a^2(r+2)r+a^2 z1g^2 (r^2-2r+a^2);
	gginc[r_]:=2a*r;
	hhinc[r_]:=r(r-2)+(r^2-2r+a^2)z1g^2/(1-z1g^2);
	
	Lzgzaux=Sqrt[(gginc[r2g]^2+hhinc[r2g]*ffinc[r2g])EEg^2-hhinc[r2g]*ddinc[r2g]];
	
	(-gginc[r2g] EEg+sgn Lzgzaux)/hhinc[r2g]
]


KKgfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=(Lzgfunsgn[a,r1g,r2g,z1g,sgn]-a EEgfunsgn[a,r1g,r2g,z1g,sgn])^2+z1g^2(a^2(1-EEgfunsgn[a,r1g,r2g,z1g,sgn]^2)+Lzgfunsgn[a,r1g,r2g,z1g,sgn]^2/(1-z1g^2));


(* ::Subsection::Closed:: *)
(*Geodesic frequencies*)


(* ::Subsubsection::Closed:: *)
(*Radial and polar geodesic frequencies*)


\[CapitalUpsilon]rgfun[a_,p_,e_,xg_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,krg},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	
	\[Pi]/(2EllipticK[krg]) Sqrt[(1-EEg^2)(r1g-r3g)(r2g-r4g)]
]


\[CapitalUpsilon]zgfun[a_,p_,e_,xg_]:=Module[{z1g,EEg,Lzg,z2g,kzg},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	(\[Pi] z2g)/(2EllipticK[kzg])
]


(* ::Text:: *)
(*These functions s are used for the derivatives of the geodesic frequencies wrt to r1g, r2g and z1g. The "sgn" in all the functions denotes the type of orbits. +1 is for prograde orbits, -1 for retrograde orbits. It is equivalent to the sign of x. *)


\[CapitalUpsilon]rgfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{EEg,Lzg,KKg,r3g,r4g,krg},
	EEg=EEgfunsgn[a,r1g,r2g,z1g,sgn];
	Lzg=Lzgfunsgn[a,r1g,r2g,z1g,sgn];
	KKg=KKgfunsgn[a,r1g,r2g,z1g,sgn];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	
	\[Pi]/(2EllipticK[krg]) Sqrt[(1-EEg^2)(r1g-r3g)(r2g-r4g)]
]


\[CapitalUpsilon]zgfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{EEg,Lzg,z2g,kzg},
	EEg=EEgfunsgn[a,r1g,r2g,z1g,sgn];
	Lzg=Lzgfunsgn[a,r1g,r2g,z1g,sgn];
	
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	(\[Pi] z2g)/(2EllipticK[kzg])
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time geodesic frequency*)


\[CapitalUpsilon]tgrfun[a_,p_,e_,xg_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,rp,rm,krg,ellK,hr,hp,hm},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr ( r3g-rp)/(r2g-rp);
	hm=hr ( r3g-rm)/(r2g-rm);
	
	(4+a^2)EEg+EEg(1/2 (4+r1g+r2g+r3g)r3g- (r1g r2g)/2+(r1g-r3g)(r2g-r4g) EllipticE[krg]/(2ellK)+1/2 (4+2/(1-EEg^2))(r2g-r3g) EllipticPi[hr,krg]/ellK)+(2EEg)/(rp-rm) (((4-a Lzg/EEg)rp-2a^2)/(r3g-rp) (1-(r2g-r3g)/(r2g-rp) EllipticPi[hp,krg]/ellK))-(2EEg)/(rp-rm) (((4-a Lzg/EEg)rm-2a^2)/(r3g-rm) (1-(r2g-r3g)/(r2g-rm) EllipticPi[hm,krg]/ellK))
]


\[CapitalUpsilon]tgzfun[a_,p_,e_,xg_]:=Module[{z1g,EEg,Lzg,KKg,z2g,kzg},	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	If[z1g==0,
		-a^2EEg,
		-a^2 EEg+(EEg(KKg-(Lzg-a EEg)^2))/((1-EEg^2)z1g^2) (1-EllipticE[kzg]/EllipticK[kzg])
	]
]


(* ::Text:: *)
(*These functions s are used for the derivatives of the geodesic frequencies wrt to r1g, r2g and z1g. The "sgn" in all the functions denotes the type of orbits. +1 is for prograde orbits, -1 for retrograde orbits. It is equivalent to the sign of x. *)


\[CapitalUpsilon]tgrfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{EEg,Lzg,KKg,r3g,r4g,rp,rm,krg,ellK,hr,hp,hm},
	EEg=EEgfunsgn[a,r1g,r2g,z1g,sgn];
	Lzg=Lzgfunsgn[a,r1g,r2g,z1g,sgn];
	KKg=KKgfunsgn[a,r1g,r2g,z1g,sgn];
	 
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	 
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr ( r3g-rp)/(r2g-rp);
	hm=hr ( r3g-rm)/(r2g-rm);
	
	  (4+a^2)EEg+EEg(1/2 (4+r1g+r2g+r3g)r3g- (r1g r2g)/2+(r1g-r3g)(r2g-r4g) EllipticE[krg]/(2ellK)+1/2 (4+2/(1-EEg^2))(r2g-r3g) EllipticPi[hr,krg]/ellK)+(2EEg)/(rp-rm) (((4-a Lzg/EEg)rp-2a^2)/(r3g-rp) (1-(r2g-r3g)/(r2g-rp) EllipticPi[hp,krg]/ellK))-(2EEg)/(rp-rm) (((4-a Lzg/EEg)rm-2a^2)/(r3g-rm) (1-(r2g-r3g)/(r2g-rm) EllipticPi[hm,krg]/ellK))
]


\[CapitalUpsilon]tgzfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{EEg,Lzg,KKg,z2g,kzg},
	EEg=EEgfunsgn[a,r1g,r2g,z1g,sgn];
	Lzg=Lzgfunsgn[a,r1g,r2g,z1g,sgn];
	KKg=KKgfunsgn[a,r1g,r2g,z1g,sgn];
	
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	
	-a^2 EEg+(EEg(KKg-(Lzg-a EEg)^2))/((1-EEg^2)z1g^2) (1-EllipticE[kzg]/EllipticK[kzg])
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal time geodesic frequency*)


\[CapitalUpsilon]\[Phi]grfun[a_,p_,e_,xg_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,rp,rm,krg,ellK,hr,hp,hm},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr ( r3g-rp)/(r2g-rp);
	hm=hr ( r3g-rm)/(r2g-rm);
	
	a/(rp-rm) ((2 EEg rp-a Lzg)/(r3g-rp) (1-(r2g-r3g)/(r2g-rp) EllipticPi[hp,krg]/ellK))-a/(rp-rm) ((2 EEg rm-a Lzg)/(r3g-rm) (1-(r2g-r3g)/(r2g-rm) EllipticPi[hm,krg]/ellK))
]


\[CapitalUpsilon]\[Phi]gzfun[a_,p_,e_,xg_]:=Module[{z1g,EEg,Lzg,z2g,kzg},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	Lzg/EllipticK[kzg] EllipticPi[z1g^2,kzg]
]


(* ::Text:: *)
(*These functions s are used for the derivatives of the geodesic frequencies wrt to r1g, r2g and z1g. The "sgn" in all the functions denotes the type of orbits. +1 is for prograde orbits, -1 for retrograde orbits. It is equivalent to the sign of x. *)


\[CapitalUpsilon]\[Phi]grfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{EEg,Lzg,KKg,r3g,r4g,rp,rm,krg,ellK,hr,hp,hm},
	EEg=EEgfunsgn[a,r1g,r2g,z1g,sgn];
	Lzg=Lzgfunsgn[a,r1g,r2g,z1g,sgn];
	KKg=KKgfunsgn[a,r1g,r2g,z1g,sgn];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	rp=1+Sqrt[1-a^2];
	rm=1-Sqrt[1-a^2];
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	hr=(r1g-r2g)/(r1g-r3g);
	hp=hr ( r3g-rp)/(r2g-rp);
	hm=hr ( r3g-rm)/(r2g-rm);
	
	a/(rp-rm) ((2 EEg rp-a Lzg)/(r3g-rp) (1-(r2g-r3g)/(r2g-rp) EllipticPi[hp,krg]/ellK))-a/(rp-rm) ((2 EEg rm-a Lzg)/(r3g-rm) (1-(r2g-r3g)/(r2g-rm) EllipticPi[hm,krg]/ellK))
]


\[CapitalUpsilon]\[Phi]gzfunsgn[a_,r1g_,r2g_,z1g_,sgn_]:=Module[{EEg,Lzg,KKg,z2g,kzg},
	EEg=EEgfunsgn[a,r1g,r2g,z1g,sgn];
	Lzg=Lzgfunsgn[a,r1g,r2g,z1g,sgn];
	
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	Lzg/EllipticK[kzg] EllipticPi[z1g^2,kzg]
]


(* ::Subsubsection::Closed:: *)
(*Spin precession frequency*)


\[CapitalUpsilon]pfun[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,krg,z2g,kzg,\[Psi]freq},
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	
	\[Psi]freq=(Sqrt[KKg]((r3g^2+a^2)EEg-a Lzg))/(r3g^2+KKg)+((r2g-r3g)((KKg-a^2)EEg+a Lzg))/((r2g^2+KKg)(r3g^2+KKg)) 1/(2I) ((r2g+I Sqrt[KKg])(r3g+I Sqrt[KKg])EllipticPi[(r1g-r2g)/(r1g-r3g) (r3g-I Sqrt[KKg])/(r2g-I Sqrt[KKg]),krg] 1/EllipticK[krg]-(r2g-I Sqrt[KKg])(r3g-I Sqrt[KKg])EllipticPi[(r1g-r2g)/(r1g-r3g) (-r3g-I Sqrt[KKg])/(-r2g-I Sqrt[KKg]),krg] 1/EllipticK[krg])-EEg Sqrt[KKg]+((KKg-a^2)EEg+a Lzg)/Sqrt[KKg] EllipticPi[(a^2 z1g^2)/KKg,kzg] 1/EllipticK[kzg];
	Re[\[Psi]freq]
]


(* ::Subsection::Closed:: *)
(*Homogeneous angles*)


(* ::Text:: *)
(*Inversion of the integrals in A12 and A13 in Phys. Rev. D 100, 104030*)


\[Chi]rfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,jSN,Wrg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi] ellK wr,krg];
	Wrg = jSN^2/(r1g-r3g);
	
	If[Mod[wr,2Pi]<Pi,
		ArcSin[((r1g+r2g-2r3g)Wrg-1)/(1-(r1g-r2g)Wrg)]+2\[Pi]*Floor[wr/(2\[Pi])]
		,
		\[Pi]-ArcSin[((r1g+r2g-2r3g)Wrg-1)/(1-(r1g-r2g)Wrg)]+2\[Pi]*Floor[wr/(2\[Pi])]
	 ]
]


Sin\[Chi]rfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,jSN,Wrg},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi] ellK wr,krg];
	Wrg = jSN^2/(r1g-r3g);
	
	((r1g+r2g-2r3g)Wrg-1)/(1-(r1g-r2g)Wrg)
]


(* ::Text:: *)
(*The \[Pi]/2 term is added to match the convention in KerrGeodesics package of the Toolkit*)


\[Chi]zfun[wz_,a_,p_,e_,xg_,const_]:=Module[{z1g,EEg,Lzg,z2g,kzg,ellK},
	EEg=const[[1]];
	Lzg=const[[2]];
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	JacobiAmplitude[2/\[Pi] ellK(wz+\[Pi]/2),kzg]
]


(* ::Subsubsection::Closed:: *)
(*Derivatives of \[Chi]z with respect to z1g*)


d\[Chi]zdz1gfun[wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,kzg,ellK,\[Alpha]z},
	EEg=const[[1]];
	Lzg=const[[2]];
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	kzg=a^2(1-EEg^2)z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	\[Alpha]z=(\[Pi]+2wz)*ellK/\[Pi] ;
	
	1/(2(1-kzg)kzg*\[Pi])(kzg*\[Pi]*JacobiCD[(2wz)ellK/\[Pi],kzg] JacobiCN[\[Alpha]z,kzg]+JacobiDN[\[Alpha]z,kzg]((\[Pi]+2 wz)EllipticE[kzg]-\[Pi]*JacobiEpsilon[\[Alpha]z,kzg]))2a^2*z1g(1-EEg^2)(-((Lzg ^2z1g^2)/((1-z1g^2)^2 z2g^4))+1/z2g^2)
]


(* ::Subsubsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Geodesics radial and polar trajectories and derivatives*)


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*Derivatives trajectories wrt to wr and wz*)


drgdwrfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,jSN},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	jSN=JacobiSN[1/\[Pi] ellK wr,krg];
	
	(2 ellK (r1g-r2g) (r1g-r3g) (r2g-r3g) JacobiCN[(ellK wr)/\[Pi],krg] JacobiDN[(ellK wr)/\[Pi],krg] jSN)/(\[Pi] (-r1g+r3g+(r1g-r2g) jSN^2)^2)
]


dzgdwzfun[wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,z1g,z2g,kzg,ellK},
	EEg=const[[1]];
	Lzg=const[[2]];
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	2/\[Pi] ellK z1g JacobiCN[2/\[Pi] ellK(wz+\[Pi]/2),kzg] JacobiDN[2/\[Pi] ellK(wz+\[Pi]/2),kzg]
]


(* ::Subsubsection::Closed:: *)
(*Derivatives radial trajectories wrt to constant of motions*)


drgdEgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,\[Alpha]r,jSN,drgdkrg,drgdr3g,dkrgdr3g,dkrgdr4g,dr3gdEg,dr4gdEg},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	\[Alpha]r= wr /\[Pi] ellK;
	jSN=JacobiSN[\[Alpha]r,krg];
	
	drgdkrg=(jSN (r1g-r2g) (r1g-r3g) (r2g-r3g) JacobiCN[(ellK wr)/\[Pi],krg] JacobiDN[(ellK wr)/\[Pi],krg] (wr EllipticE[krg]+jSN krg \[Pi] JacobiCD[(ellK wr)/\[Pi],krg]-\[Pi] JacobiEpsilon[(ellK wr)/\[Pi],krg]))/((1-krg) krg \[Pi] (-r1g+r3g+jSN^2 (r1g-r2g))^2);
	drgdr3g=-(((r1g-r2g)^2 jSN^2 (1-jSN^2))/(-r1g+r3g+(r1g-r2g) jSN^2)^2);
	
	dkrgdr3g =((r1g-r2g) (r1g-r4g))/((r1g-r3g)^2 (r2g-r4g));
	dkrgdr4g =-(((r1g-r2g) (r2g-r3g))/((r1g-r3g) (r2g-r4g)^2));
	dr3gdEg=(2 EEg)/(1-EEg^2)^2+1/(2 Sqrt[-((a^2 (KKg-(-a EEg+Lzg)^2))/((1-EEg^2) r1g r2g))+(-(1/(1-EEg^2))+(r1g+r2g)/2)^2]) (-((2 a^3 (-a EEg+Lzg))/((1-EEg^2) r1g r2g))-(2 a^2 EEg (KKg-(-a EEg+Lzg)^2))/((1-EEg^2)^2r1g r2g)-(4 EEg (-(1/(1-EEg^2))+(r1g+r2g)/2))/(1-EEg^2)^2);
	dr4gdEg=(a^2 (2 (-a^2 EEg+a (1+EEg^2) Lzg+EEg (KKg-Lzg^2)) r3g+(1-EEg^2) (-KKg+(-a EEg+Lzg)^2) dr3gdEg))/((1-EEg^2)^2 r1g r2g r3g^2);
	
	(drgdr3g +drgdkrg dkrgdr3g) dr3gdEg+drgdkrg dkrgdr4g dr4gdEg
]


drgdLzgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,\[Alpha]r,jSN,drgdkrg,drgdr3g,dkrgdr3g,dkrgdr4g,dr3gdLzg,dr4gdLzg},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	\[Alpha]r= wr /\[Pi] ellK;
	jSN=JacobiSN[\[Alpha]r,krg];
	
	drgdkrg=(jSN (r1g-r2g) (r1g-r3g) (r2g-r3g) JacobiCN[(ellK wr)/\[Pi],krg] JacobiDN[(ellK wr)/\[Pi],krg] (wr EllipticE[krg]+jSN krg \[Pi] JacobiCD[(ellK wr)/\[Pi],krg]-\[Pi] JacobiEpsilon[(ellK wr)/\[Pi],krg]))/((1-krg) krg \[Pi] (-r1g+r3g+jSN^2 (r1g-r2g))^2);
	drgdr3g=-(((r1g-r2g)^2 jSN^2 (1-jSN^2))/(-r1g+r3g+(r1g-r2g) jSN^2)^2);
	
	dkrgdr3g =((r1g-r2g) (r1g-r4g))/((r1g-r3g)^2 (r2g-r4g));
	dkrgdr4g =-(((r1g-r2g) (r2g-r3g))/((r1g-r3g) (r2g-r4g)^2));
	
	dr3gdLzg=(a^2 (-a EEg+Lzg))/((1-EEg^2) r1g r2g Sqrt[(a^2 (KKg-(-a EEg+Lzg)^2))/((-1+EEg^2) r1g r2g)+1/4 (2/(-1+EEg^2)+r1g+r2g)^2]);
	dr4gdLzg=(a^2 (2 (-a EEg+Lzg) r3g+(KKg-(-a EEg+Lzg)^2) dr3gdLzg))/((-1+EEg^2) r1g r2g r3g^2);
	
	(drgdr3g +drgdkrg dkrgdr3g) dr3gdLzg+drgdkrg dkrgdr4g dr4gdLzg
]


drgdKgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,\[Alpha]r,jSN,drgdkrg,drgdr3g,dkrgdr3g,dkrgdr4g,dr3gdKg,dr4gdKg},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+\[Sqrt](((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)));
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	\[Alpha]r= wr /\[Pi] ellK;
	jSN=JacobiSN[\[Alpha]r,krg];
	
	drgdkrg=(jSN (r1g-r2g) (r1g-r3g) (r2g-r3g) JacobiCN[(ellK wr)/\[Pi],krg] JacobiDN[(ellK wr)/\[Pi],krg] (wr EllipticE[krg]+jSN krg \[Pi] JacobiCD[(ellK wr)/\[Pi],krg]-\[Pi] JacobiEpsilon[(ellK wr)/\[Pi],krg]))/((1-krg) krg \[Pi] (-r1g+r3g+jSN^2 (r1g-r2g))^2);
	drgdr3g=-(((r1g-r2g)^2 jSN^2 (1-jSN^2))/(-r1g+r3g+(r1g-r2g) jSN^2)^2);
	
	dkrgdr3g =((r1g-r2g) (r1g-r4g))/((r1g-r3g)^2 (r2g-r4g));
	dkrgdr4g =-(((r1g-r2g) (r2g-r3g))/((r1g-r3g) (r2g-r4g)^2));
	
	dr3gdKg=-(a^2/(2 (1-EEg^2) r1g r2g Sqrt[(a^2 (KKg-(-a EEg+Lzg)^2))/((-1+EEg^2) r1g r2g)+1/4 (2/(-1+EEg^2)+r1g+r2g)^2]));
	dr4gdKg=(a^2 (r3g-(KKg-(-a EEg+Lzg)^2) dr3gdKg))/((1-EEg^2) r1g r2g r3g^2);
	
	(drgdr3g +drgdkrg dkrgdr3g)dr3gdKg+drgdkrg dkrgdr4g dr4gdKg
]


(* ::Subsubsection::Closed:: *)
(*Derivatives polar trajectories wrt to constant of motions*)


dzgdEgfun[wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,kzg,ellK,\[Alpha]z},
	EEg=const[[1]];
	Lzg=const[[2]];
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	\[Alpha]z=(\[Pi]+2 wz) /\[Pi] ellK;
	
	(EEg z1g (z2g^2-a^2 (1-EEg^2)) JacobiCN[\[Alpha]z,kzg] JacobiDN[\[Alpha]z,kzg] ((\[Pi]+2 wz) EllipticE[kzg]-\[Pi] JacobiEpsilon[\[Alpha]z,kzg]+kzg \[Pi] JacobiCD[\[Alpha]z,kzg] JacobiSN[\[Alpha]z,kzg]))/((1-EEg^2) (-1+kzg) \[Pi] z2g^2)
]


dzgdLzgfun[wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,kzg,ellK,\[Alpha]z},
	EEg=const[[1]];
	Lzg=const[[2]];
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	\[Alpha]z=(\[Pi]+2 wz) /\[Pi] ellK;
	
	(Lzg z1g JacobiCN[\[Alpha]z,kzg] JacobiDN[\[Alpha]z,kzg] ((\[Pi]+2 wz) EllipticE[kzg]-\[Pi] JacobiEpsilon[\[Alpha]z,kzg]+kzg \[Pi] JacobiCD[\[Alpha]z,kzg] JacobiSN[\[Alpha]z,kzg]))/((-1+kzg) \[Pi] (1-z1g^2) z2g^2)
]


(* ::Subsubsection::Closed:: *)
(*Derivatives radial trajectories wrt to r1g and r2g*)


drgdr1gfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,\[Alpha]r,jSN,drgdr1gexplicit,drgdkrg,drgdr3g,dkrgdr1g,dkrgdr3g,dkrgdr4g,dr3gdr1g,dr4gdr1g},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	\[Alpha]r= wr /\[Pi] ellK;
	jSN=JacobiSN[\[Alpha]r,krg];
	
	drgdr1gexplicit=((r2g-r3g)^2 jSN^2)/(-r1g+r3g+(r1g-r2g)jSN^2)^2;
	
	drgdkrg=(jSN (r1g-r2g) (r1g-r3g) (r2g-r3g) JacobiCN[(ellK wr)/\[Pi],krg] JacobiDN[(ellK wr)/\[Pi],krg] (wr EllipticE[krg]+jSN krg \[Pi] JacobiCD[(ellK wr)/\[Pi],krg]-\[Pi] JacobiEpsilon[(ellK wr)/\[Pi],krg]))/((1-krg) krg \[Pi] (-r1g+r3g+jSN^2 (r1g-r2g))^2);
	drgdr3g=-(((r1g-r2g)^2 jSN^2 (1-jSN^2))/(-r1g+r3g+(r1g-r2g) jSN^2)^2);
	
	dkrgdr1g =((r2g-r3g) (r3g-r4g))/((r1g-r3g)^2 (r2g-r4g));
	dkrgdr3g =((r1g-r2g) (r1g-r4g))/((r1g-r3g)^2 (r2g-r4g));
	dkrgdr4g =-(((r1g-r2g) (r2g-r3g))/((r1g-r3g) (r2g-r4g)^2));
	
	dr3gdr1g=-(1/2)+1/(2 Sqrt[-((a^2 (-(-a EEg+Lzg)^2+KKg))/((1-EEg^2) r1g r2g))+1/4 (2/(-1+EEg^2)+r1g+r2g)^2]) (1/(-1+EEg^2)+(r1g+r2g)/2+(a^2 (-(-a EEg+Lzg)^2+KKg))/(r1g^2r2g (1-EEg^2)));
	
	dr4gdr1g=(-(dr3gdr1g/r3g^2) (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g(1-EEg^2))-1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g^2 r2g(1-EEg^2)));
	
	drgdr1gexplicit+drgdr3g dr3gdr1g+drgdkrg(dkrgdr1g +dkrgdr3g dr3gdr1g+dkrgdr4g dr4gdr1g)
]


drgdr2gfun[wr_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,krg,ellK,\[Alpha]r,jSN,drgdr2gexplicit,drgdkrg,drgdr3g,dkrgdr2g,dkrgdr3g,dkrgdr4g,dr3gdr2g,dr4gdr2g},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	krg=((r1g-r2g)(r3g-r4g))/((r1g-r3g)(r2g-r4g));
	ellK=EllipticK[krg];
	\[Alpha]r= wr /\[Pi] ellK;
	jSN=JacobiSN[\[Alpha]r,krg];
	
	drgdr2gexplicit=((r1g-r3g)^2 (1-jSN^2))/(-r1g+r3g+(r1g-r2g)jSN^2)^2;
	
	drgdkrg=(jSN (r1g-r2g) (r1g-r3g) (r2g-r3g) JacobiCN[(ellK wr)/\[Pi],krg] JacobiDN[(ellK wr)/\[Pi],krg] (wr EllipticE[krg]+jSN krg \[Pi] JacobiCD[(ellK wr)/\[Pi],krg]-\[Pi] JacobiEpsilon[(ellK wr)/\[Pi],krg]))/((1-krg) krg \[Pi] (-r1g+r3g+jSN^2 (r1g-r2g))^2);
	drgdr3g=-(((r1g-r2g)^2 jSN^2 (1-jSN^2))/(-r1g+r3g+(r1g-r2g) jSN^2)^2);
	
	dkrgdr2g =-(((r1g-r4g) (r3g-r4g))/((r1g-r3g) (r2g-r4g)^2));
	dkrgdr3g =((r1g-r2g) (r1g-r4g))/((r1g-r3g)^2 (r2g-r4g));
	dkrgdr4g =-(((r1g-r2g) (r2g-r3g))/((r1g-r3g) (r2g-r4g)^2));
	
	dr3gdr2g=-(1/2)+1/(2 Sqrt[-((a^2 (-(-a EEg+Lzg)^2+KKg))/((1-EEg^2) r1g r2g))+1/4 (2/(-1+EEg^2)+r1g+r2g)^2]) (1/(-1+EEg^2)+(r1g+r2g)/2+(a^2 (-(-a EEg+Lzg)^2+KKg))/(r1g r2g^2 (1-EEg^2)));
	dr4gdr2g=(-(dr3gdr2g/r3g^2) (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g(1-EEg^2))-1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g^2(1-EEg^2)));
	
	drgdr2gexplicit+drgdr3g dr3gdr2g+drgdkrg(dkrgdr2g +dkrgdr3g dr3gdr2g+dkrgdr4g dr4gdr2g)
]


(* ::Subsection::Closed:: *)
(*Expansion coordinate time and azimuthal trajectories as Fourier series, purely oscillating part*)


\[CapitalDelta]trgfun[wr_,coeff_,\[CapitalUpsilon]rg_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2Sin[n wr]coeff[[n+dim+1]]/(n \[CapitalUpsilon]rg),{n,1,dim}]
];


\[CapitalDelta]tzgfun[wz_,coeff_,\[CapitalUpsilon]zg_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2Sin[k wz]coeff[[k+dim+1]]/(k \[CapitalUpsilon]zg),{k,1,dim}]
];


\[CapitalDelta]\[Phi]rgfun[wr_,coeff_,\[CapitalUpsilon]rg_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2Sin[n wr]coeff[[n+dim+1]]/(n \[CapitalUpsilon]rg),{n,1,dim}]
];


\[CapitalDelta]\[Phi]zgfun[wz_,coeff_,\[CapitalUpsilon]zg_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2Sin[k wz]coeff[[k+dim+1]]/(k \[CapitalUpsilon]zg),{k,1,dim}]
];


(* ::Subsection::Closed:: *)
(*Geodesics velocities*)


(* ::Subsubsection::Closed:: *)
(*Coordinate-time geodesic velocity*)


Vtrgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	EEg((rg^2+a^2)^2/(rg^2-2rg+a^2))-(2a*rg)/(rg^2-2rg+a^2)Lzg
]


Vtzgfun[wz_,a_,p_,e_,xg_,const_]:=Module[{z1g,EEg,Lzg,KKg,zg},
	{EEg,Lzg,KKg}=const;
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	-EEg*a^2(1-zg^2)
]


dVtrgdrgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg,\[CapitalDelta]rg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	(2a*Lzg(rg^2-a^2))/\[CapitalDelta]rg^2+(4 EEg*rg(a^2+rg^2))/\[CapitalDelta]rg-(2EEg(rg-1)(a^2+rg^2)^2)/\[CapitalDelta]rg^2
]


dVtzgdzgfun[wz_,a_,p_,e_,xg_,const_]:=Module[{z1g,EEg,Lzg,KKg,zg},
	{EEg,Lzg,KKg}=const;
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	2a^2*EEg*zg
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal geodesic velocity*)


V\[Phi]rgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	a/(rg^2-2rg+a^2)(EEg(rg^2+a^2)-a*Lzg)-a*EEg
]


V\[Phi]zgfun[wz_,a_,p_,e_,xg_,const_]:=Module[{z1g,EEg,Lzg,KKg,\[Chi]z,zg},
	{EEg,Lzg,KKg}=const;
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	Lzg/(1-zg^2)
]


dV\[Phi]rgdrgfun[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg,\[CapitalDelta]rg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	
	(2a*EEg*rg)/\[CapitalDelta]rg-(2a*(rg-1)(-a*Lzg+EEg(a^2+rg^2)))/\[CapitalDelta]rg^2
]


dV\[Phi]zgdzgfun[wz_,a_,p_,e_,xg_,const_]:=Module[{z1g,EEg,Lzg,KKg,zg},
	{EEg,Lzg,KKg}=const;
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	(2Lzg*zg)/(1-zg^2)^2
]


(* ::Subsubsection::Closed:: *)
(*Geodesic coordinate time and azimuthal Fourier coefficients of the velocities*)


Vtrgcoeff[nmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,wrlist,Vtrglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	Vtrglist=Vtrgfun[wrlist,a,p,e,xg,const];
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . Vtrglist)/stepsr,10^(-16)]
]


Vtzgcoeff[kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsz,wzlist,Vtzglist,ExpjkTable},
	stepsz=4*kmax;
	wzlist=Table[i,{i,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	Vtzglist=Vtzgfun[wzlist,a,p,e,xg,const];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	Chop[(Vtzglist . ExpjkTable)/stepsz,10^(-16)]
]


V\[Phi]rgcoeff[nmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,wrlist,V\[Phi]rglist,ExpniTable},
	stepsr=4*nmax;
	wrlist=Table[i,{i,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	V\[Phi]rglist=V\[Phi]rgfun[wrlist,a,p,e,xg,const];
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	Chop[(ExpniTable . V\[Phi]rglist)/stepsr,10^(-16)]
]


V\[Phi]zgcoeff[kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsz,wzlist,V\[Phi]zglist,ExpjkTable},
	stepsz=4*kmax;
	wzlist=Table[i,{i,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	V\[Phi]zglist=V\[Phi]zgfun[wzlist,a,p,e,xg,const];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	Chop[(V\[Phi]zglist . ExpjkTable)/stepsz,10^(-16)]	
]


(* ::Subsubsection::Closed:: *)
(*Expansion coordinate time and azimuthal velocities as Fourier series*)


dtrgd\[Lambda][wr_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	coeff[[dim+1]]+Sum[2Cos[n wr]coeff[[n+dim+1]],{n,1,dim}]
];


dtzgd\[Lambda][wz_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	coeff[[dim+1]]+Sum[2Cos[k wz]coeff[[k+dim+1]],{k,1,dim}]
];


d\[Phi]rgd\[Lambda][wr_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	coeff[[dim+1]]+Sum[2Cos[n wr]coeff[[n+dim+1]],{n,1,dim}]
];


d\[Phi]zgd\[Lambda][wz_,coeff_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	coeff[[dim+1]]+Sum[2Cos[k wz]coeff[[k+dim+1]],{k,1,dim}]
];


(* ::Section::Closed:: *)
(*Spin coefficients and projections Marck tetrad over Christoffel symbols*)


(* ::Subsection::Closed:: *)
(*Spin coefficients - radial part*)


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalS]^r_{12}*)


\[ScriptCapitalS]12rfun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrgdvrgdrg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrgdvrgdrg=(2rg EEg Rrg-rg \[CapitalDelta]rg -(rg-1)(KKg+rg^2));
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	vrg vzg (2 a^2rg zg Rrg)/(Sqrt[KKg]\[CapitalDelta]rg^2 \[CapitalSigma]^2)+vrgdvrgdrg (rg Rrg(KKg-a^2zg^2))/(Sqrt[KKg](KKg+rg^2)\[CapitalDelta]rg^2 \[CapitalSigma])+vrg^2/(Sqrt[KKg]\[CapitalDelta]rg^2) (Rrg/\[CapitalSigma]-(2rg^2 Rrg)/\[CapitalSigma]^2-EEg (KKg-rg^2)/(KKg+rg^2))
]


\[ScriptCapitalS]12rfuneq[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,rg,\[CapitalDelta]rg,Rrg,Yrg,vrgdvrgdrg,vrg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	vrgdvrgdrg=(2rg EEg Rrg-rg \[CapitalDelta]rg -(rg-1)(KKg+rg^2));
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	
	vrgdvrgdrg (Rrg Sqrt[KKg])/((KKg+rg^2)\[CapitalDelta]rg^2 rg)+vrg^2/(Sqrt[KKg]\[CapitalDelta]rg^2) (Rrg/rg^2-(2 Rrg)/rg^2-EEg (KKg-rg^2)/(KKg+rg^2))
]


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalS]^r_{23}*)


\[ScriptCapitalS]23rfun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrgdvrgdrg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrgdvrgdrg=(2rg EEg Rrg-rg \[CapitalDelta]rg -(rg-1)(KKg+rg^2));
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	vrg vzg (a Rrg(rg^2(KKg-a^2zg^2)-a^2zg^2(KKg+rg^2)))/(Sqrt[KKg]\[CapitalDelta]rg^2 \[CapitalSigma]^2 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])-vrgdvrgdrg (a zg Rrg(KKg-a^2zg^2))/(Sqrt[KKg]\[CapitalDelta]rg^2 \[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vrg^2 (a rg zg(2 KKg Rrg-a rg^2(Lzg-a EEg)-a^4 zg^4 EEg-(3a^2 rg^2 EEg-a^3 (Lzg-a EEg))zg^2))/(Sqrt[KKg]\[CapitalDelta]rg^2 \[CapitalSigma]^2 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])
]


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalS]^r_{31}*)


\[ScriptCapitalS]31rfun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrgdvrgdrg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrgdvrgdrg=(2rg EEg Rrg-rg \[CapitalDelta]rg -(rg-1)(KKg+rg^2));
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	vrg (a zg rg Rrg(rg(a Lzg-rg EEg)-a(Lzg-a EEg))Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta]rg^3 \[CapitalSigma] Sqrt[KKg+rg^2])-vrg vzg^2 (a^3zg^3 Sqrt[KKg+rg^2])/(KKg(1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg-a^2zg^2])-vrg (a^3zg^3Zzg(Lzg \[CapitalSigma] -2rg Zzg)Sqrt[KKg+rg^2])/(KKg \[CapitalDelta]rg^2\[CapitalSigma]^2(1-zg^2)Sqrt[KKg-a^2zg^2])-vrg vrgdvrgdrg (a rg zg (KKg-a^2zg^2))/(KKg \[CapitalDelta]rg^2\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vrg^2vzg (KKg a rg)/(KKg \[CapitalDelta]rg^2\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vrg^3 (a zg rg(\[CapitalSigma](rg-1)+rg \[CapitalDelta]rg)Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta]rg^3\[CapitalSigma]^2 Sqrt[KKg+rg^2])
]


(* ::Text:: *)
(*The spin coefficient is rescaled by the radial velocity vrg. This is convenient for computation of the radial velocities*)


\[ScriptCapitalS]31rfunres[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrgdvrgdrg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrgdvrgdrg=(2rg EEg Rrg-rg \[CapitalDelta]rg -(rg-1)(KKg+rg^2));
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	(a zg rg Rrg(rg(a Lzg-rg EEg)-a(Lzg-a EEg))Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta]rg^3 \[CapitalSigma] Sqrt[KKg+rg^2])-vzg^2 (a^3zg^3 Sqrt[KKg+rg^2])/(KKg(1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg-a^2zg^2])-(a^3zg^3Zzg(Lzg \[CapitalSigma] -2rg Zzg)Sqrt[KKg+rg^2])/(KKg \[CapitalDelta]rg^2\[CapitalSigma]^2(1-zg^2)Sqrt[KKg-a^2zg^2])-vrgdvrgdrg (a rg zg (KKg-a^2zg^2))/(KKg \[CapitalDelta]rg^2\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vrg vzg (KKg a rg)/(KKg \[CapitalDelta]rg^2\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vrg^2 (a zg rg(\[CapitalSigma](rg-1)+rg \[CapitalDelta]rg)Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta]rg^3\[CapitalSigma]^2 Sqrt[KKg+rg^2])
]


(* ::Subsection::Closed:: *)
(*Spin coefficients - polar part*)


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalS]^z_{12}*)


\[ScriptCapitalS]12zfun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vzgdvzgdzg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vzgdvzgdzg=(1-EEg^2)a^2zg(2zg^2-z1g^2)-z2g^2zg;
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	vzg vrg (2 a rg zg Zzg)/(Sqrt[KKg](1-zg^2) \[CapitalSigma]^2)+vzgdvzgdzg (a zg(KKg+rg^2) Zzg)/(Sqrt[KKg](1-zg^2)(KKg-a^2zg^2) \[CapitalSigma])-vzg^2/(Sqrt[KKg](1-zg^2)) ((Rrg(rg^2-a^2 zg^2))/\[CapitalSigma]^2+(2a^2EEg zg^2(KKg+rg^2))/(\[CapitalSigma] (KKg-a^2zg^2)))
]


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalS]^z_{23}*)


\[ScriptCapitalS]23zfun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vzgdvzgdzg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vzgdvzgdzg=(1-EEg^2)a^2zg(2zg^2-z1g^2)-z2g^2zg;
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	vzg vrg ((rg^2(KKg-a^2 zg^2)-a^2zg^2(KKg+rg^2))Zzg)/(Sqrt[KKg](1-zg^2) \[CapitalSigma]^2 Sqrt[(KKg+rg^2)] Sqrt[(KKg-a^2zg^2)])+vzgdvzgdzg (rg Sqrt[(KKg+rg^2)] Zzg)/(Sqrt[KKg](1-zg^2)\[CapitalSigma] Sqrt[(KKg-a^2zg^2)] )+vzg^2((a rg zg(a^2 EEg(1-zg^2)(rg^2-a^2zg^2)-2a KKg Zzg-EEg \[CapitalSigma]^2-a Lzg(rg^2-a^2zg^2)))/(Sqrt[KKg](1-zg^2) \[CapitalSigma]^2 Sqrt[(KKg-a^2zg^2)] Sqrt[(KKg+rg^2)]))
]


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalS]^z_{31}*)


\[ScriptCapitalS]31zfun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vzgdvzgdzg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vzgdvzgdzg=(1-EEg^2)a^2zg(2zg^2-z1g^2)-z2g^2zg;
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	-vzg((a Lzg rg zg^2 Zzg Sqrt[KKg+rg^2])/(KKg (1-zg^2)^2\[CapitalSigma] Sqrt[KKg-a^2zg^2])+vrg^2  (a rg^3 Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2]))-vzg vzgdvzgdzg (a zg rg(KKg+rg^2))/(KKg(1-zg^2)\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vzg^2vrg (a zg KKg)/(KKg(1-zg^2)\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])-vzg((rg^3Rrg(2rg Zzg-Lzg \[CapitalSigma])Sqrt[KKg-a^2zg^2])/(KKg(1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2]))-vzg^3 (a rg zg^2(\[CapitalSigma]-a^2(1-zg^2))Sqrt[KKg+rg^2])/(KKg (1-zg^2)^2\[CapitalSigma]^2 Sqrt[KKg-a^2zg^2])
]


(* ::Text:: *)
(*The spin coefficient is rescaled by the polar velocity v2zg. This is convenient for computation of the polar velocities*)


\[ScriptCapitalS]31zfunres[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vzgdvzgdzg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vzgdvzgdzg=(1-EEg^2)a^2zg(2zg^2-z1g^2)-z2g^2zg;
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	-((a Lzg rg zg^2 Zzg Sqrt[KKg+rg^2])/(KKg (1-zg^2)^2\[CapitalSigma] Sqrt[KKg-a^2zg^2])+vrg^2  (a rg^3 Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2]))-vzgdvzgdzg (a zg rg(KKg+rg^2))/(KKg(1-zg^2)\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vzg vrg (a zg KKg)/(KKg(1-zg^2)\[CapitalSigma] Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])-((rg^3Rrg(2rg Zzg-Lzg \[CapitalSigma])Sqrt[KKg-a^2zg^2])/(KKg(1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2]))-vzg^2 (a rg zg^2(\[CapitalSigma]-a^2(1-zg^2))Sqrt[KKg+rg^2])/(KKg (1-zg^2)^2\[CapitalSigma]^2 Sqrt[KKg-a^2zg^2])
]


(* ::Subsection::Closed:: *)
(*Projection Christoffel symbols over Marck tetrad - coordinate time part*)


(* ::Subsubsection::Closed:: *)
(*\[CapitalGamma]^t_{12}*)


\[CapitalGamma]t12fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrg,vzg,Lzred,\[ScriptCapitalQ]1,\[ScriptCapitalQ]2},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg-a EEg(1-zg^2));
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	Lzred=Lzg-a EEg;
	\[ScriptCapitalQ]1=a^2Lzred \[CapitalSigma]+a^3 zg^4 Rrg+rg^3(rg-2)Zzg+KKg(Lzg \[CapitalSigma]-2rg Zzg);
	\[ScriptCapitalQ]2=2 a rg(KKg+rg^2)Zzg-(rg^2+a^2)\[CapitalSigma](EEg KKg+a Lzred);
	
	-2\[ScriptCapitalQ]1((a rg (Rrg(KKg-a^2zg^2)(a^2(rg^2-a^2zg^2)+2rg^4+rg^2\[CapitalSigma])+2a^3zg^2(KKg+rg^2)\[CapitalDelta]rg Zzg))/(Sqrt[KKg]\[CapitalDelta]rg^2\[CapitalSigma]^4(KKg+rg^2)(KKg-a^2zg^2)))-2\[ScriptCapitalQ]2 (rg(Rrg(rg^2+a^2)(rg^2-a^2zg^2)(KKg-a^2zg^2)+2a^3zg^2(KKg+rg^2)\[CapitalDelta]rg Zzg))/(Sqrt[KKg]\[CapitalDelta]rg^2\[CapitalSigma]^4(KKg+rg^2)(KKg-a^2zg^2))-vrg vzg ( 8a^2rg^2 zg)/(Sqrt[KKg]\[CapitalDelta]rg \[CapitalSigma]^3)-vrg^2 (2rg(rg^2-a^2)(KKg-a^2zg^2))/(Sqrt[KKg]\[CapitalDelta]rg^2\[CapitalSigma]^2(KKg+rg^2))
]


\[CapitalGamma]t12funeq[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,rg,\[CapitalDelta]rg,Rrg,Yrg,vrg,Lzred,\[ScriptCapitalQ]1,\[ScriptCapitalQ]2},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	rg=rgfun[wr,a,p,e,xg,const];
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	Lzred=Lzg-a EEg;
	\[ScriptCapitalQ]1=a^2Lzred rg^2+rg^3(rg-2)Lzred+KKg(Lzg rg^2-2rg Lzred);
	\[ScriptCapitalQ]2=2 a rg(KKg+rg^2)Lzred-(rg^2+a^2)rg^2(EEg KKg+a Lzred);
	
	-2\[ScriptCapitalQ]1 (a (Rrg(a^2+2rg^2+rg^2)))/(Sqrt[KKg]\[CapitalDelta]rg^2rg^5(KKg+rg^2))-2\[ScriptCapitalQ]2 (Rrg(rg^2+a^2))/(Sqrt[KKg]\[CapitalDelta]rg^2rg^5(KKg+rg^2))-vrg^2 (2(rg^2-a^2)KKg)/(Sqrt[KKg]\[CapitalDelta]rg^2rg^3(KKg+rg^2))
]


(* ::Subsubsection::Closed:: *)
(*\[CapitalGamma]^t_{23}*)


\[CapitalGamma]t23fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrg,vzg,Lzred,\[ScriptCapitalQ]1,\[ScriptCapitalQ]2},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg-a EEg(1-zg^2));
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	Lzred=Lzg-a EEg;
	\[ScriptCapitalQ]1=a^2Lzred \[CapitalSigma]+a^3 zg^4 Rrg+rg^3(rg-2)Zzg+KKg(Lzg \[CapitalSigma]-2rg Zzg);
	\[ScriptCapitalQ]2=2 a rg(KKg+rg^2)Zzg-(rg^2+a^2)\[CapitalSigma](EEg KKg+a Lzred);
	
	(2a^2zg \[ScriptCapitalQ]1)/(Sqrt[KKg]\[CapitalDelta]rg \[CapitalSigma]^4 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2]) (Rrg/\[CapitalDelta]rg (3rg^4-a^4zg^2+a^2rg^2(1+zg^2))-2a rg^2Zzg)+(2a zg \[ScriptCapitalQ]2)/(Sqrt[KKg]\[CapitalDelta]rg \[CapitalSigma]^4 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2]) (Rrg/\[CapitalDelta]rg (rg^2+a^2)(rg^2-a^2zg^2)-2a rg^2Zzg)+vrg vzg (4rg(2a^3 rg^2zg^2-a KKg(rg^2-a^2zg^2)))/(Sqrt[KKg]\[CapitalDelta]rg \[CapitalSigma]^3 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vrg^2 (2a zg (rg^2-a^2)Sqrt[KKg-a^2zg^2])/(Sqrt[KKg]\[CapitalDelta]rg^2\[CapitalSigma]^2 Sqrt[KKg+rg^2])
]


(* ::Subsubsection::Closed:: *)
(*\[CapitalGamma]^t_{31}*)


\[CapitalGamma]t31fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	vrg (4a^2 rg zg Zzg)/(\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])-vzg (4a rg^2 Rrg)/(\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])
]


(* ::Subsection::Closed:: *)
(*Projection Christoffel symbols over Marck tetrad - azimuthal part*)


(* ::Subsubsection::Closed:: *)
(*\[CapitalGamma]^\[Phi]_{12}*)


\[CapitalGamma]\[Phi]12fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrg,vzg,Lzred,\[ScriptCapitalQ]1,\[ScriptCapitalQ]2},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg-a EEg(1-zg^2));
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	Lzred=Lzg-a EEg;
	\[ScriptCapitalQ]1=a^2Lzred \[CapitalSigma]+a^3 zg^4 Rrg+rg^3(rg-2)Zzg+KKg(Lzg \[CapitalSigma]-2rg Zzg);
	\[ScriptCapitalQ]2=2 a rg(KKg+rg^2)Zzg-(rg^2+a^2)\[CapitalSigma](EEg KKg+a Lzred);
	
	-2\[ScriptCapitalQ]1((a zg^2 Zzg (2a^2rg(1-zg^2)+\[CapitalSigma]^2))/(Sqrt[KKg](1-zg^2)^2\[CapitalDelta]rg \[CapitalSigma]^4(KKg-a^2zg^2)))
	-2\[ScriptCapitalQ]1((rg Rrg(rg^2(rg^2+a^2)+rg \[CapitalSigma](rg-\[CapitalSigma])-a^4zg^2(1-zg^2)))/(Sqrt[KKg](1-zg^2)\[CapitalDelta]rg^2\[CapitalSigma]^4(KKg+rg^2)))-2\[ScriptCapitalQ]2 (a rg(Rrg(rg^2-a^2zg^2)(1-zg^2)(KKg-a^2zg^2)+2a zg^2(KKg+rg^2)\[CapitalDelta]rg Zzg))/(Sqrt[KKg](1-zg^2)\[CapitalDelta]rg^2\[CapitalSigma]^4(KKg+rg^2)(KKg-a^2zg^2))+vrg vzg ( 4a rg zg(\[CapitalSigma]-2rg))/(Sqrt[KKg](1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^3)-vrg^2 (2a rg(rg-1)(KKg-a^2zg^2))/(Sqrt[KKg]\[CapitalDelta]rg^2\[CapitalSigma]^2(KKg+rg^2))-vzg^2 (2 a zg^2(KKg+rg^2))/(Sqrt[KKg](1-zg^2)^2\[CapitalSigma]^2(KKg-a^2zg^2) )
]


\[CapitalGamma]\[Phi]12funeq[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,rg,\[CapitalDelta]rg,Rrg,Yrg,vrg,Lzred,\[ScriptCapitalQ]1,\[ScriptCapitalQ]2},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	Lzred=Lzg-a EEg;
	\[ScriptCapitalQ]1=a^2Lzred rg^2+rg^3(rg-2)Lzred+KKg(Lzg rg^2-2rg Lzred);
	\[ScriptCapitalQ]2=2 a rg(KKg+rg^2)Lzred-(rg^2+a^2)rg^2(EEg KKg+a Lzred);
	
	-2\[ScriptCapitalQ]1 ( Rrg((rg^2+a^2)+rg (rg-rg^2)))/(Sqrt[KKg]\[CapitalDelta]rg^2rg^5(KKg+rg^2))-2\[ScriptCapitalQ]2 (a Rrg)/(Sqrt[KKg]\[CapitalDelta]rg^2rg^5(KKg+rg^2))-vrg^2 (2a (rg-1)KKg)/(Sqrt[KKg]\[CapitalDelta]rg^2rg^3(KKg+rg^2))
]


(* ::Subsubsection::Closed:: *)
(*\[CapitalGamma]^\[Phi]_{23}*)


\[CapitalGamma]\[Phi]23fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrg,vzg,Lzred,\[ScriptCapitalQ]1,\[ScriptCapitalQ]2},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg+a EEg(-1+zg^2));
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	Lzred=Lzg-a EEg;
	\[ScriptCapitalQ]1=a^2Lzred \[CapitalSigma]+a^3 zg^4 Rrg+rg^3(rg-2)Zzg+KKg(Lzg \[CapitalSigma]-2rg Zzg);
	\[ScriptCapitalQ]2=2 a rg(KKg+rg^2)Zzg-(rg^2+a^2)\[CapitalSigma](EEg KKg+a Lzred);
	
	(2zg \[ScriptCapitalQ]1)/(Sqrt[KKg](1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^4 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2]) (a Rrg/\[CapitalDelta]rg (a^2(rg^2-a^2zg^2)+2rg^4+\[CapitalSigma](a^2zg^2-rg \[CapitalSigma]))-2a^2 rg^2Zzg-(rg Zzg \[CapitalSigma]^2)/(1-zg^2))
	-(2a zg \[ScriptCapitalQ]2)/(Sqrt[KKg](1-zg^2)\[CapitalDelta]rg^2 \[CapitalSigma]^4 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2]) (2 Lzg rg^4+(rg^2+a^2)(3EEg rg^2-a Lzred)a zg^2-a Rrg(3rg^2+a^2zg^4)-4rg^3Zzg)+vrg vzg (2(\[CapitalSigma]-2rg)(KKg(rg^2-a^2zg^2)-2a^2zg^2rg^2))/(Sqrt[KKg](1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^3 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])+vrg^2 (2a^2 zg (rg-1)Sqrt[KKg-a^2zg^2])/(Sqrt[KKg]\[CapitalDelta]rg^2\[CapitalSigma]^2 Sqrt[KKg+rg^2])-vzg^2 (2rg zg Sqrt[KKg+rg^2])/(Sqrt[KKg](1-zg^2)^2\[CapitalSigma]^2 Sqrt[KKg-a^2zg^2])
]


(* ::Subsubsection::Closed:: *)
(*\[CapitalGamma]^\[Phi]_{31}*)


\[CapitalGamma]\[Phi]31fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,Rrg,\[CapitalSigma],Zzg,Yrg,Yzg,vrg,vzg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	\[CapitalSigma]=rg^2+a^2zg^2;
	Zzg=(Lzg-a EEg(1-zg^2));
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	vrg=Sign[Pi-Mod[wr,2Pi]]Sqrt[(r1g-rg)(rg-r2g)]Yrg;
	vzg=z1g Cos[\[Chi]z]Yzg;
	
	vzg (2rg Rrg(\[CapitalSigma]-2rg))/((1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])-vrg (2a zg Zzg(\[CapitalSigma]-2rg))/((1-zg^2)\[CapitalDelta]rg \[CapitalSigma]^2 Sqrt[KKg+rg^2] Sqrt[KKg-a^2zg^2])
]


(* ::Section::Closed:: *)
(*Shifts to the constants of motion, turning points and their averages*)


(* ::Subsection::Closed:: *)
(*Radial and polar corrections to the turning points*)


(* ::Subsubsection::Closed:: *)
(*Parallel component of the spin - fixed constants of motion*)


r1sparFCfun[wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[CapitalDelta]r1g,Rr1g,Yr1g,zg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[CapitalDelta]r1g=(a^2-2 r1g+r1g^2);
	Rr1g=(-a Lzg+EEg (a^2+r1g^2));
	Yr1g=Sqrt[(1-EEg^2) (-r1g+r3g) (-r1g+r4g)];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Rr1g/(KKg+r1g^2) ((2 \[CapitalDelta]r1g Sqrt[KKg])/( (r1g-r2g)Yr1g^2 )-(r1g(KKg-a^2 zg^2))/(Sqrt[KKg](r1g^2+a^2zg^2)))+(2 a \[CapitalDelta]r1g)/((r1g-r2g)Yr1g^2 ) RealSign[-a EEg+Lzg]
]


r2sparFCfun[wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[CapitalDelta]r2g,Rr2g,Yr2g,zg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[CapitalDelta]r2g=(a^2-2 r2g+r2g^2);
	Rr2g=(-a Lzg+EEg (a^2+r2g^2));
	Yr2g=Sqrt[(1-EEg^2) (-r2g+r3g) (-r2g+r4g)];
	zg=zgfun[wz,a,p,e,xg,const];
	
	-(Rr2g/(KKg+r2g^2))((2 \[CapitalDelta]r2g Sqrt[KKg])/( (r1g-r2g)Yr2g^2 )+(r2g(KKg-a^2 zg^2))/(Sqrt[KKg](r2g^2+a^2zg^2)))-(2 a \[CapitalDelta]r2g)/((r1g-r2g)Yr2g^2 ) RealSign[-a EEg+Lzg]
]


z1sparFCfun[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,Zz1g,Yz1g,rg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	Zz1g=(Lzg+a EEg(-1+z1g^2));
	Yz1g=Sqrt[-a^2 (1-EEg^2) z1g^2+z2g^2];
	rg=rgfun[wr,a,p,e,xg,const];
	
	If[z1g==0,
		0
		,
		(a Zz1g)/(KKg-a^2z1g^2) ((Sqrt[KKg](1-z1g^2))/( z1g Yz1g^2 )-(z1g(KKg+rg^2))/(Sqrt[KKg](rg^2+a^2z1g^2)))-(a (1-z1g^2))/(z1g Yz1g^2 ) RealSign[-a EEg+Lzg]
	]
]


(* ::Subsubsection::Closed:: *)
(*Parallel component of the spin - fixed turning points (on average)*)


r1sparDHfun[wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[CapitalDelta]r1g,Rr1g,Yr1g,zg,dVrgdEg,dVrgdLzg,dVrgdKg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[CapitalDelta]r1g=(a^2-2 r1g+r1g^2);
	Rr1g=(-a Lzg+EEg (a^2+r1g^2));
	Yr1g=Sqrt[(1-EEg^2) (-r1g+r3g) (-r1g+r4g)];
	zg=zgfun[wz,a,p,e,xg,const];
	
	dVrgdEg=2(r1g^2+a^2)Rr1g;
	dVrgdLzg=-2a Rr1g;
	dVrgdKg=-\[CapitalDelta]r1g;
	
	Rr1g/(KKg+r1g^2) ((2 \[CapitalDelta]r1g Sqrt[KKg])/( (r1g-r2g)Yr1g^2 )-(r1g(KKg-a^2 zg^2))/(Sqrt[KKg](r1g^2+a^2zg^2)))+(2 a \[CapitalDelta]r1g)/((r1g-r2g)Yr1g^2 ) RealSign[-a EEg+Lzg]+1/((r1g-r2g)Yr1g^2) (dVrgdEg \[Delta]const[[1]]+dVrgdLzg \[Delta]const[[2]]+dVrgdKg \[Delta]const[[3]])
]


r2sparDHfun[wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[CapitalDelta]r2g,Rr2g,Yr2g,zg,dVrgdEg,dVrgdLzg,dVrgdKg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[CapitalDelta]r2g=(a^2-2 r2g+r2g^2);
	Rr2g=(-a Lzg+EEg (a^2+r2g^2));
	Yr2g=Sqrt[(1-EEg^2) (-r2g+r3g) (-r2g+r4g)];
	zg=zgfun[wz,a,p,e,xg,const];
	
	dVrgdEg=2(r2g^2+a^2)Rr2g;
	dVrgdLzg=-2a Rr2g;
	dVrgdKg=-\[CapitalDelta]r2g;
	
	-(Rr2g/(KKg+r2g^2))((2 \[CapitalDelta]r2g Sqrt[KKg])/( (r1g-r2g)Yr2g^2 )+(r2g(KKg-a^2 zg^2))/(Sqrt[KKg](r2g^2+a^2zg^2)))-(2 a \[CapitalDelta]r2g)/((r1g-r2g)Yr2g^2 ) RealSign[-a EEg+Lzg]-1/((r1g-r2g)Yr2g^2) (dVrgdEg \[Delta]const[[1]]+dVrgdLzg \[Delta]const[[2]]+dVrgdKg \[Delta]const[[3]])
]


z1sparDHfun[wr_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,Zz1g,Yz1g,rg,dVzgdEg,dVzgdLzg,dVzgdKg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	Zz1g=(Lzg-a EEg (1-z1g^2));
	Yz1g=Sqrt[-a^2 (1-EEg^2) z1g^2+z2g^2];
	rg=rgfun[wr,a,p,e,xg,const];
	
	dVzgdEg=2 a (1-z1g^2) Zz1g;
	dVzgdLzg=-2 Zz1g;
	dVzgdKg=1-z1g^2;
	
	If[z1g==0,
		0
		,
		(a Zz1g)/(KKg-a^2z1g^2) ((Sqrt[KKg](1-z1g^2))/( z1g Yz1g^2 )-(z1g(KKg+rg^2))/(Sqrt[KKg](rg^2+a^2z1g^2)))-(a (1-z1g^2))/(z1g Yz1g^2 ) RealSign[-a EEg+Lzg]+1/(2z1g Yz1g^2) (dVzgdEg \[Delta]const[[1]]+dVzgdLzg \[Delta]const[[2]]+dVzgdKg \[Delta]const[[3]])
	]
]


(* ::Subsubsection::Closed:: *)
(*Orthogonal component of the spin*)


r1sortCos\[Psi]fun[wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,zg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	(a (EEg(r1g^2+a^2)-a Lzg) zg Sqrt[KKg-a^2 zg^2])/(Sqrt[KKg (KKg+r1g^2)] (r1g^2+a^2 zg^2))
]


r2sortCos\[Psi]fun[wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,zg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	(a (EEg(r2g^2+a^2)-a Lzg)zg Sqrt[KKg-a^2zg^2] )/(Sqrt[KKg (KKg+r2g^2)] (r2g^2+a^2 zg^2))
]


z1sortCos\[Psi]fun[wr_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,z1g,rg},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	-((rg Sqrt[KKg+rg^2] (Lzg-a EEg (1-z1g^2)))/(Sqrt[KKg] Sqrt[KKg-a^2 z1g^2] (rg^2+a^2 z1g^2)))
]


(* ::Subsection::Closed:: *)
(*Averages shifts turning points - fixed constants of motion*)


r1savgFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]r1g,\[CapitalDelta]r1g,Rr1g,Yr1g,ellK},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]r1g=-a^2z1g^2/r1g^2;
	
	\[CapitalDelta]r1g=(a^2-2 r1g+r1g^2);
	Rr1g=(-a Lzg+EEg (a^2+r1g^2));
	Yr1g=Sqrt[(1-EEg^2) (-r1g+r3g) (-r1g+r4g)];
	
	ellK=EllipticK[kzg];
	
	Rr1g (r1g^2ellK -(KKg+r1g^2)EllipticPi[\[Gamma]r1g,kzg] )/(Sqrt[KKg]r1g(KKg+r1g^2)ellK)+(2 \[CapitalDelta]r1g (Sqrt[KKg]Rr1g+a(KKg+r1g^2)RealSign[Lzg-a EEg]))/((KKg+r1g^2)(r1g-r2g)Yr1g^2)
]


r2savgFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]r2g,\[CapitalDelta]r2g,Rr2g,Yr2g,ellK},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]r2g=-a^2z1g^2/r2g^2;
	
	\[CapitalDelta]r2g=(a^2-2 r2g+r2g^2);
	Rr2g=(-a Lzg+EEg (a^2+r2g^2));
	Yr2g=Sqrt[(1-EEg^2) (-r2g+r3g) (-r2g+r4g)];
	
	ellK=EllipticK[kzg];
	
	Rr2g (r2g^2ellK -(KKg+r2g^2)EllipticPi[\[Gamma]r2g,kzg] )/(Sqrt[KKg]r2g(KKg+r2g^2)ellK)-(2 \[CapitalDelta]r2g (Sqrt[KKg]Rr2g+a(KKg+r2g^2)RealSign[Lzg-a EEg]))/((KKg+r2g^2)(r1g-r2g)Yr2g^2)
]


z1savgFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,\[CapitalDelta]r2g,Zz1g,Yz1g,ellK,arg,argcc},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
	
	Zz1g=(Lzg+a EEg(-1+z1g^2));
	Yz1g=Sqrt[-a^2 (1-EEg^2) z1g^2+z2g^2];
	
	ellK=EllipticK[krg];
	
	arg=(r3g-I a z1g)/(r2g-I a z1g) (r1g-r2g)/(r1g-r3g);
	argcc=(r3g+I a z1g)/(r2g+I a z1g) (r1g-r2g)/(r1g-r3g);
	
	If[(z1g==0||a==0),
		0,
		(-a z1g Zz1g (r3g^2+KKg))/(Sqrt[KKg](r3g^2+a^2z1g^2)(KKg-a^2z1g^2))+(a (1-z1g^2)(Sqrt[KKg]Zz1g-(KKg-a^2z1g^2)RealSign[Lzg-a EEg]))/(z1g(KKg-a^2z1g^2)Yz1g^2)+(Zz1g(r2g-r3g))/(2Sqrt[KKg]ellK(r2g^2+a^2z1g^2)(r3g^2+a^2z1g^2)) Re[(I a^2z1g^2+a z1g (r2g+r3g)-I r2g r3g )EllipticPi[arg,krg]+(-I a^2z1g^2+a z1g (r2g+r3g)+I r2g r3g )EllipticPi[argcc,krg]]
	]
]


(* ::Subsection::Closed:: *)
(*Shifts constants of motion in the fixed turning points (on average)*)


\[Delta]constmotionfun[a_,p_,e_,xg_]:=Module[{r1g,r1d,r2g,r2d,z1g,z1d,sgn,dEgdr1g,dEgdr2g,dEgdz1g,dLzgdr1g,dLzgdr2g,dLzgdz1g,dKgdr1g,dKgdr2g,dKgdz1g,r1savgsub,r2savgsub,z1savgsub},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	sgn=RealSign[xg];
	
	dEgdr1g=D[EEgfunsgn[a,r1d,r2g,z1g,sgn],r1d]/.r1d->r1g;
	dEgdr2g=D[EEgfunsgn[a,r1g,r2d,z1g,sgn],r2d]/.r2d->r2g;
	dEgdz1g=D[EEgfunsgn[a,r1g,r2g,z1d,sgn],z1d]/.z1d->z1g;
	
	dLzgdr1g=D[Lzgfunsgn[a,r1d,r2g,z1g,sgn],r1d]/.r1d->r1g;
	dLzgdr2g=D[Lzgfunsgn[a,r1g,r2d,z1g,sgn],r2d]/.r2d->r2g;
	dLzgdz1g=D[Lzgfunsgn[a,r1g,r2g,z1d,sgn],z1d]/.z1d->z1g;
	
	dKgdr1g=D[KKgfunsgn[a,r1d,r2g,z1g,sgn],r1d]/.r1d->r1g;
	dKgdr2g=D[KKgfunsgn[a,r1g,r2d,z1g,sgn],r2d]/.r2d->r2g;
	dKgdz1g=D[KKgfunsgn[a,r1g,r2g,z1d,sgn],z1d]/.z1d->z1g;
	
	r1savgsub=r1savgFC[a,p,e,xg];
	r2savgsub=r2savgFC[a,p,e,xg];
	z1savgsub=z1savgFC[a,p,e,xg];
	
	-{dEgdr1g r1savgsub+dEgdr2g r2savgsub+dEgdz1g z1savgsub,dLzgdr1g r1savgsub+dLzgdr2g r2savgsub+dLzgdz1g z1savgsub,dKgdr1g r1savgsub+dKgdr2g r2savgsub+dKgdz1g z1savgsub}
]


(* ::Section::Closed:: *)
(*Shift frequencies - fixed constants of motion*)


(* ::Subsection::Closed:: *)
(*Shift radial frequency*)


(* ::Subsubsection::Closed:: *)
(*Part 1*)


\[ScriptCapitalY]1rFC[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,\[Gamma]r1g,\[Gamma]r2g,ellK,ellPi,ellPir1g,ellPir2g,\[CapitalSigma]z1g,Rrg,Rr1g,Rr2g,Yrg,Yr1g,Yr2g,\[CapitalDelta]rg,\[CapitalDelta]r1g,\[CapitalDelta]r2g,rg,r1savgres,r2savgres},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	\[Gamma]r1g=-a^2z1g^2/r1g^2;
	\[Gamma]r2g=-a^2z1g^2/r2g^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	ellPir1g=EllipticPi[\[Gamma]r1g,kzg];
	ellPir2g=EllipticPi[\[Gamma]r2g,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Rr1g=(-a Lzg+EEg (a^2+r1g^2));
	Rr2g=(-a Lzg+EEg (a^2+r2g^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yr1g=Sqrt[(1-EEg^2) (-r1g+r3g) (-r1g+r4g)];
	Yr2g=Sqrt[(1-EEg^2) (-r2g+r3g) (-r2g+r4g)];
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	\[CapitalDelta]r1g=(a^2-2 r1g+r1g^2);
	\[CapitalDelta]r2g=(a^2-2 r2g+r2g^2);
	4(Rrg(ellK rg^2-(KKg+rg^2)ellPi)(KKg(rg-1)+2a EEg (Lzg-a EEg)rg +a^2 rg+rg^2(2(1-EEg^2)rg-3)))/(Sqrt[KKg]z2g rg(KKg+rg^2)(r1g-rg) (rg-r2g) Yrg^3)+4 ellK/z2g (\[CapitalDelta]rg(Sqrt[KKg]Rrg+a(KKg+rg^2)RealSign[-a EEg+Lzg]))/((r1g-rg) (rg-r2g)(KKg+rg^2)Yrg^3)-4 ellK/z2g r1savgFC[a,p,e,xg]/(2(r1g-rg)Yrg)+4 ellK/z2g r2savgFC[a,p,e,xg]/(2(rg-r2g)Yrg)
]


(* ::Subsubsection::Closed:: *)
(*Part 2*)


\[ScriptCapitalY]2r[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPi,\[CapitalSigma]z1g,Rrg,Yrg,rg,r1savgres,r2savgres},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	4/(Sqrt[KKg]Yrg z2g)(-EEg (KKg-rg^2)/(KKg+rg^2) ellK)+(4Rrg )/(Sqrt[KKg]Yrg z2g) (ellPi/rg^2 -(z2g^2 EllipticE[kzg])/(\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2))+ellK/\[CapitalSigma]z1g-((1-EEg^2)(3rg^2+2a^2z1g^2))/(\[CapitalSigma]z1g ((1-EEg^2)rg^2+z2g^2)) ellPi-((a^2 z1g^2 +2rg^2)z2g^2)/(rg^2\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2))ellPi)
]


(* ::Subsubsection::Closed:: *)
(*Part 3*)


\[ScriptCapitalY]3rFC[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,Yrg,rg,r1savgres,r2savgres},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	r1savgres=r1savgFC[a,p,e,xg];
	r2savgres=r2savgFC[a,p,e,xg];
	
	1/Yrg^2 ((1-EEg^2) (2 rg-r3g-r4g))/(2Yrg) ((r1savgres+r2savgres)/2+((r1savgres-r2savgres)/2)Sin[\[Chi]r])
]


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]rsfunFC[a_,p_,e_,xg_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];

	((\[CapitalUpsilon]rg^2\[CapitalUpsilon]zg)/(2\[Pi])^2 NIntegrate[(\[ScriptCapitalY]1rFC[a,p,e,xg,\[Chi]r]+\[ScriptCapitalY]2r[a,p,e,xg,\[Chi]r]),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]rg^2/(2\[Pi]) NIntegrate[\[ScriptCapitalY]3rFC[a,p,e,xg,\[Chi]r],{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec])
]


(* ::Subsection::Closed:: *)
(*Shift polar frequency*)


(* ::Subsubsection::Closed:: *)
(*Part 1*)


\[ScriptCapitalY]1zFC[a_,p_,e_,xg_,\[Chi]z_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,ellK,Zzg,zg,Yzg,argzg,argzgcc,Lzred,aux},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
	
	ellK=EllipticK[krg];
	
	zg=z1g Sin[\[Chi]z];
	
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	argzg=(r3g-I a zg)/(r2g-I a zg) (r1g-r2g)/(r1g-r3g);
	argzgcc=(r3g+I a zg)/(r2g+I a zg) (r1g-r2g)/(r1g-r3g);
	Lzred=Lzg-a EEg;
	aux=a^2(1-EEg^2)+Lzg^2;
	
	If[z1g==0,
		(a RealSign[Lzg-a EEg])/aux^(1/2) (Lzg /(aux Lzred)+1/(3 r1g^2 r2g  r3g^2) (-2(r1g-r3g) (r2g r3g+r1g (r2g+r3g)) EllipticE[krg]/EllipticK[krg]+r1g ((r2g-r3g) r3g+r1g (2 r2g+r3g))))
		,
		Sec[\[Chi]z]^2/Yzg^3 ((a(1-zg^2))/z1g^2 ((Sqrt[KKg] Zzg)/(KKg-a^2zg^2)-RealSign[-a EEg+Lzg])+(Sin[\[Chi]z]Zzg(a^2(1-EEg^2)(2 zg^2-z1g^2)-z2g^2))/(Sqrt[KKg](r2g^2+a^2 zg^2)(r3g^2+a^2 zg^2)) ((r2g-r3g) /(2ellK z1g) Re[EllipticPi[argzgcc,krg](-I r2g r3g-a(r2g+r3g)zg+I a^2 zg^2)+EllipticPi[argzg,krg](I r2g r3g-a(r2g+r3g)zg-I a^2 zg^2)]+ a (KKg+r3g^2)/(KKg-a^2 zg^2)(r2g^2+a^2 zg^2)Sin[\[Chi]z]))-z1savgFC[a,p,e,xg]/(z1g Yzg Cos[\[Chi]z]^2)
	]
]


(* ::Subsubsection::Closed:: *)
(*Part 2*)


\[ScriptCapitalY]2z[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPi,\[CapitalSigma]z1g,Rrg,Yrg,rg,r1savgres,r2savgres},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	4/(Sqrt[KKg]Yrg z2g) (-2EEg EllipticPi[(a^2z1g^2)/KKg,kzg] +1/rg^2 (Rrg+2rg^2 EEg)ellPi-( z2g^2Rrg EllipticE[kzg])/(\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2))+(Rrg ellK)/\[CapitalSigma]z1g-((1-EEg^2)Rrg ellPi)/(\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2)) (3rg^2+2a^2z1g^2)-(z2g^2Rrg (2rg^2+a^2 z1g^2) ellPi)/(\[CapitalSigma]z1g rg^2((1-EEg^2)rg^2+z2g^2)))
]


(* ::Subsubsection::Closed:: *)
(*Part 3*)


\[ScriptCapitalY]3zFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,Yz1g,\[Chi]z,kzg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	
	z2g=Sqrt[a^2 (1-EEg^2)+Lzg^2/(1-z1g^2)];
	Yz1g=Sqrt[-a^2 (1-EEg^2) z1g^2+z2g^2];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	
	If[z1g==0,
		0,
		- 4z1savgFC[a,p,e,xg](z2g^2 EllipticE[kzg]-Yz1g^2EllipticK[kzg])/(z1g z2g Yz1g^2)
	]
]


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]zsfunFC[a_,p_,e_,xg_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];

	\[CapitalUpsilon]zg^2/(2\[Pi])(If[a==0,0,NIntegrate[\[ScriptCapitalY]1zFC[a,p,e,xg,\[Chi]z],{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]]+\[CapitalUpsilon]rg/(2\[Pi]) NIntegrate[\[ScriptCapitalY]2z[a,p,e,xg,\[Chi]r],{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]+\[ScriptCapitalY]3zFC[a,p,e,xg])]


(* ::Subsection::Closed:: *)
(*Shift "coordinate time frequency"*)


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]rs/\[CapitalUpsilon]rg \[CapitalUpsilon]^t_{rg}\[RightAngleBracket] *)


\[CapitalUpsilon]tsr0funFC[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]rsfunFC[a,p,e,xg,prec]/\[CapitalUpsilon]rgfun[a,p,e,xg] \[CapitalUpsilon]tgrfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]zs/\[CapitalUpsilon]zg \[CapitalUpsilon]^t_{zg}\[RightAngleBracket] *)


\[CapitalUpsilon]tsz0funFC[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]zsfunFC[a,p,e,xg,prec]/\[CapitalUpsilon]zgfun[a,p,e,xg] \[CapitalUpsilon]tgzfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*Radial part  \[LeftAngleBracket]\[PartialD]V^t_{rg}/\[PartialD]rg rs\[RightAngleBracket]*)


\[CapitalUpsilon]tsr2funFC[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,Yrg,rg,r1savgres,r2savgres,\[CapitalDelta]rg,Rrg,dVtrgdrg,\[Chi]r},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	r1savgres=r1savgFC[a,p,e,xg];
	r2savgres=r2savgFC[a,p,e,xg];
	\[CapitalDelta]rg=rg^2-2rg+a^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	dVtrgdrg=2/\[CapitalDelta]rg (EEg rg(a^2+rg^2)-((rg-1)(a^2+rg^2)Rrg)/\[CapitalDelta]rg+ rg Rrg);
	
	\[CapitalUpsilon]rgfun[a,p,e,xg]/(2\[Pi]) NIntegrate[1/Yrg dVtrgdrg((r1savgres+r2savgres)/2+((r1savgres-r2savgres)/2)Sin[\[Chi]r]),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


(* ::Subsubsection::Closed:: *)
(*Polar part  \[LeftAngleBracket]\[PartialD]V^t_{zg}/\[PartialD]zg zs\[RightAngleBracket]*)


\[CapitalUpsilon]tsz2funFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,kzg,ellK},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	
	z2g=Sqrt[a^2 (1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	If[(z1g==0||a==0),
		0,
		-z1savgFC[a,p,e,xg] z2g/ellK 2 (EEg z2g (EllipticE[kzg]-ellK))/((1-EEg^2) z1g)
	]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^t_{yg}\[RightAngleBracket] - radial part*)


\[CapitalUpsilon]tsr1funFC[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,\[CapitalDelta]rg,Rrg,rg,Vtrg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	\[CapitalDelta]rg=rg^2-2rg+a^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	Vtrg=a Lzg+(Rrg(a^2+rg^2) )/\[CapitalDelta]rg;
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	(\[CapitalUpsilon]rg \[CapitalUpsilon]zg)/(2\[Pi])^2 NIntegrate[(\[ScriptCapitalY]1r[a,p,e,xg,\[Chi]r]+\[ScriptCapitalY]2r[a,p,e,xg,\[Chi]r])Vtrg,{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]rg/(2\[Pi]) NIntegrate[\[ScriptCapitalY]3r[a,p,e,xg,\[Chi]r]Vtrg,{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^t_{yg}\[RightAngleBracket]  - polar part*)


\[CapitalUpsilon]tsz11funFC[a_,p_,e_,xg_,\[Chi]z_]:=Module[{z1g,zg,Vtzg},
	z1g=Sqrt[1-xg^2];
	zg=z1g Sin[\[Chi]z];
	Vtzg=-a^2EEgfun[a,p,e,xg](1-zg^2);
	
	\[ScriptCapitalY]1z[a,p,e,xg,\[Chi]z]Vtzg
]


\[CapitalUpsilon]tsz12funFC[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPi,\[CapitalSigma]z1g,Rrg,Yrg,rg,Lzred},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Lzred=Lzg-a EEg;
	
	(2 EEg)/(Sqrt[KKg] Yrg z2g) (-2 (a Lzred-3 EEg rg^2) ellK+(2 (a^2+rg^2) Rrg (z2g^2 (EllipticE[kzg]-ellK)-(1-EEg^2) rg^2 ellK))/(\[CapitalSigma]z1g((1-EEg^2) rg^2+z2g^2))+4 EEg (KKg ellK+(a^2-KKg) EllipticPi[(a^2 z1g^2)/KKg,kzg]))+1/(Sqrt[KKg] Yrg z2g) (4 EEg )/rg^2 ((a^2+rg^2) (a Lzred-3 EEg  rg^2)-(Rrg ((-1+EEg^2) rg^2 (3 a^2  rg^2+rg^4+2 a^4 z1g^2)+a^2 (-a^2 z1g^2+rg^2 (-2+z1g^2))z2g^2))/(\[CapitalSigma]z1g  ((1-EEg^2) rg^2+z2g^2))) EllipticPi[\[Gamma]rg,kzg]
]


\[CapitalUpsilon]tsz13funFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,Yz1g,\[Chi]z,kzg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	
	z2g=Sqrt[a^2 (1-EEg^2)+Lzg^2/(1-z1g^2)];
	Yz1g=Sqrt[-a^2 (1-EEg^2) z1g^2+z2g^2];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	
	If[(z1g==0||a==0),
		0,
		- 4z1savgFC[a,p,e,xg](-EEg)(z2g^2(a^2 (1-EEg^2) (1+z1g^2) -2 z2g^2) EllipticE[kzg]-Yz1g^2 (a^2 (1-EEg^2)-2 z2g^2) EllipticK[kzg])/((1-EEg^2)z1g z2g Yz1g^2)
	]
]


\[CapitalUpsilon]tsz1funFC[a_,p_,e_,xg_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	If[a==0,
		0,
		\[CapitalUpsilon]zg /(2\[Pi]) (NIntegrate[\[CapitalUpsilon]tsz11fun[a,p,e,xg,\[Chi]z],{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]rg/(2\[Pi]) NIntegrate[\[CapitalUpsilon]tsz12fun[a,p,e,xg,\[Chi]r],{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]tsz13fun[a,p,e,xg])
	]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[CapitalSigma]/2 \[CapitalGamma]^t_{12}\[RightAngleBracket]*)


\[CapitalGamma]t1avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,ellK,ellPi,Lzred,Rrg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r,int2,int1,int0},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	ellK= EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	Lzred=Lzg-a EEg;
	
	int2=((1-EEg^2) rg^4ellK-(1-EEg^2) rg^2 (3 rg^2 +2 a^2 z1g^2) ellPi-z2g^2 (rg^2 EllipticE[kzg]-rg^2 ellK+(2 rg^2 +a^2 z1g^2)ellPi)) (8Rrg^2)/(Sqrt[KKg] rg z2g ((1-EEg^2) rg^2+z2g^2) \[CapitalDelta]rg (rg^2+a^2 z1g^2));
	int1=(4Rrg ellPi)/(Sqrt[KKg] rg z2g \[CapitalDelta]rg^2) (-3a^3 Lzred+4a Lzred rg+a(8 a EEg-Lzg) rg^2+EEg rg^3(-12 +5 rg));
	int0=-((4rg Rrg ellK)/(Sqrt[KKg] (KKg+rg^2) z2g \[CapitalDelta]rg^2))(-a^3Lzred+2EEg KKg(a^2 -2 rg)+(2 a^2 EEg+2 EEg KKg+a Lzg) rg^2 +rg^3EEg(rg-4));
	
	NIntegrate[1/Yrg (int2+int1+int0),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


\[CapitalGamma]t2avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	NIntegrate[((a^2-rg^2) (r1g-r2g)^2 Yrg Cos[\[Chi]r]^2 (rg^2 EllipticK[kzg]-(KKg+rg^2) EllipticPi[\[Gamma]rg,kzg]))/(Sqrt[KKg] rg \[CapitalDelta]rg^2(KKg+rg^2) z2g),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


\[CapitalGamma]tavg[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]rgfun[a,p,e,xg]/(2\[Pi]) \[CapitalUpsilon]zgfun[a,p,e,xg]/(2\[Pi]) (\[CapitalGamma]t1avg[a,p,e,xg,prec]+\[CapitalGamma]t2avg[a,p,e,xg,prec])


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]tsfunFC[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]tsr0funFC[a,p,e,xg,prec]-\[CapitalUpsilon]tsr1funFC[a,p,e,xg,prec]+\[CapitalUpsilon]tsr2funFC[a,p,e,xg,prec]+\[CapitalUpsilon]tsz0funFC[a,p,e,xg,prec]-\[CapitalUpsilon]tsz1funFC[a,p,e,xg,prec]+\[CapitalUpsilon]tsz2funFC[a,p,e,xg]+\[CapitalGamma]tavg[a,p,e,xg,prec]


(* ::Subsection::Closed:: *)
(*Shift azimuthal frequency*)


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]rs/\[CapitalUpsilon]rg \[CapitalUpsilon]^\[Phi]_{rg}\[RightAngleBracket] *)


\[CapitalUpsilon]\[Phi]sr0funFC[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]rsfunFC[a,p,e,xg,prec]/\[CapitalUpsilon]rgfun[a,p,e,xg] \[CapitalUpsilon]\[Phi]grfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]zs/\[CapitalUpsilon]zg \[CapitalUpsilon]^\[Phi]_{zg}\[RightAngleBracket] *)


\[CapitalUpsilon]\[Phi]sz0funFC[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]zsfunFC[a,p,e,xg,prec]/\[CapitalUpsilon]zgfun[a,p,e,xg] \[CapitalUpsilon]\[Phi]gzfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*Radial part   \[LeftAngleBracket]\[PartialD]V^\[Phi]_{rg}/\[PartialD]rg rs\[RightAngleBracket]*)


\[CapitalUpsilon]\[Phi]sr2funFC[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,Yrg,rg,r1savgres,r2savgres,\[CapitalDelta]rg,Rrg,dV\[Phi]rgdrg,\[Chi]r},
	If[a==0,
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		z1g=Sqrt[1-xg^2];
	
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
	
		r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
		r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
		rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
		
		Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
		
		r1savgres=r1savgFC[a,p,e,xg];
		r2savgres=r2savgFC[a,p,e,xg];
		\[CapitalDelta]rg=rg^2-2rg+a^2;
		Rrg=(-a Lzg+EEg (a^2+rg^2));
		dV\[Phi]rgdrg=(2 a EEg rg)/\[CapitalDelta]rg-(a (-2+2 rg) Rrg)/\[CapitalDelta]rg^2;
		
		\[CapitalUpsilon]rgfun[a,p,e,xg]/(2\[Pi]) NIntegrate[1/Yrg dV\[Phi]rgdrg((r1savgres+r2savgres)/2+((r1savgres-r2savgres)/2)Sin[\[Chi]r]),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
	]
]


(* ::Subsubsection::Closed:: *)
(*Polar part  \[LeftAngleBracket]\[PartialD]V^\[Phi]_{zg}/\[PartialD]zg /\[PartialD]zg zs\[RightAngleBracket]*)


\[CapitalUpsilon]\[Phi]sz2funFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,kzg,ellK},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	
	z2g=Sqrt[a^2 (1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	ellK=EllipticK[kzg];
	
	If[z1g==0,
		0,
		1/z1g 1/ellK (z1savgFC[a,p,e,xg]Lzg)/((1-z1g^2)(z2g^2-a^2(1-EEg^2))) ((z2g^2z1g^2-a^2(1-EEg^2))EllipticPi[z1g^2,kzg]+(a^2(1-EEg^2)-z2g^2)ellK+z2g^2EllipticE[kzg])
	]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^\[Phi]_{yg}\[RightAngleBracket] - radial part*)


\[CapitalUpsilon]\[Phi]sr1funFC[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,\[CapitalDelta]rg,Rrg,rg,V\[Phi]rg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r},
	If[a==0,
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		z1g=Sqrt[1-xg^2];
	
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
	
		\[CapitalDelta]rg=rg^2-2rg+a^2;
		Rrg=(-a Lzg+EEg (a^2+rg^2));
	
		rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
		V\[Phi]rg=(a Rrg)/\[CapitalDelta]rg-a EEg;
		\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
		\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
		(\[CapitalUpsilon]rg \[CapitalUpsilon]zg)/(2\[Pi])^2 NIntegrate[(\[ScriptCapitalY]1r[a,p,e,xg,\[Chi]r]+\[ScriptCapitalY]2r[a,p,e,xg,\[Chi]r])V\[Phi]rg,{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]rg/(2\[Pi]) NIntegrate[\[ScriptCapitalY]3r[a,p,e,xg,\[Chi]r]V\[Phi]rg,{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
	]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^\[Phi]_{yg}\[RightAngleBracket] - polar part*)


\[CapitalUpsilon]\[Phi]sz11funFC[a_,p_,e_,xg_,\[Chi]z_]:=Module[{z1g,zg,V\[Phi]zg},
	z1g=Sqrt[1-xg^2];
	zg=z1g Sin[\[Chi]z];
	V\[Phi]zg=Lzgfun[a,p,e,xg] 1/(1-zg^2);
	
	\[ScriptCapitalY]1z[a,p,e,xg,\[Chi]z]V\[Phi]zg
]


\[CapitalUpsilon]\[Phi]sz12funFC[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPiz1g,\[CapitalSigma]z1g,Rrg,Yrg,rg,Lzred},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPiz1g=EllipticPi[z1g^2,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Lzred=Lzg-a EEg;
	
	(2 Lzg)/(Sqrt[KKg] Yrg z2g) 1/(a^2+rg^2) (-(2 a^2 Rrg (z2g^2 EllipticE[kzg]-((1-EEg^2) rg^2+z2g^2) ellK))/(\[CapitalSigma]z1g((1-EEg^2) rg^2+z2g^2))-(4 rg^2 Rrg ellPiz1g)/(a^2+rg^2)+2 (Rrg+2 EEg rg^2) ellPiz1g+(4 EEg(a^2+rg^2) (KKg ellPiz1g-a^2EllipticPi[(a^2 z1g^2)/KKg,kzg]))/(a^2-KKg))+(2 Lzg )/(Sqrt[KKg] Yrg z2g) (2 a^2)/(rg^2 (a^2+rg^2)) ((Rrg+2 EEg rg^2)-(Rrg (5 (1-EEg^2) rg^6+4 rg^4 z2g^2+a^4 z1g^2 (2 (1-EEg^2) rg^2+z2g^2)+a^2 ((1-EEg^2) rg^4 (3+4 z1g^2)+rg^2 (2+3 z1g^2) z2g^2)))/((a^2+rg^2)\[CapitalSigma]z1g((1-EEg^2) rg^2+z2g^2))) EllipticPi[\[Gamma]rg,kzg]
]


\[CapitalUpsilon]\[Phi]sz13funFC[a_,p_,e_,xg_]:=Module[{r1g,r2g,z1g,EEg,Lzg,z2g,Yz1g,\[Chi]z,kzg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	
	z2g=Sqrt[a^2 (1-EEg^2)+Lzg^2/(1-z1g^2)];
	Yz1g=Sqrt[-a^2 (1-EEg^2) z1g^2+z2g^2];
	
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	
	If[(z1g==0||a==0),
		0,
		- 4z1savgFC[a,p,e,xg]Lzg a^2 (1-EEg^2) (z2g^2 EllipticE[kzg]-Yz1g^2 EllipticPi[z1g^2,kzg])/(z1g z2g (a^2 (1-EEg^2)-z2g^2 ) Yz1g^2)
	]
]


\[CapitalUpsilon]\[Phi]sz1funFC[a_,p_,e_,xg_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	\[CapitalUpsilon]zg/(2\[Pi]) (If[a==0,0,NIntegrate[\[CapitalUpsilon]\[Phi]sz11funFC[a,p,e,xg,\[Chi]z],{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]]+\[CapitalUpsilon]rg/(2\[Pi]) NIntegrate[\[CapitalUpsilon]\[Phi]sz12funFC[a,p,e,xg,\[Chi]r],{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]\[Phi]sz13funFC[a,p,e,xg])
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[CapitalSigma]/2 \[CapitalGamma]^\[Phi]_{12}\[RightAngleBracket]*)


\[CapitalGamma]\[Phi]1avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,ellK,ellPiz1g,ellPirg,ellE,Lzred,Rrg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r,int5,int4,int3,int2,int1,int0},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	ellK= EllipticK[kzg];
	ellPiz1g=EllipticPi[z1g^2,kzg];
	ellPirg=EllipticPi[\[Gamma]rg,kzg];
	ellE=EllipticE[kzg];
	Lzred=Lzg-a EEg;
	
	int5=(8a Rrg^2 )/(Sqrt[KKg] rg (a^2+rg^2) z2g ((1-EEg^2)rg^2+z2g^2) \[CapitalDelta]rg (rg^2+a^2z1g^2)) ((1-EEg^2) rg^4 ellK-(1-EEg^2) rg^2 (3 rg^2+2 a^2 z1g^2)ellPirg-z2g^2 (rg^2 (ellE-ellK)+(2 rg^2+a^2 z1g^2) ellPirg));
	
	int4=(4a Rrg ellPirg)/(Sqrt[KKg] rg (a^2+rg^2)^2 z2g \[CapitalDelta]rg^2) (-3a^5Lzred+4a^3Lzred rg+a^3(9 a EEg-2  Lzg) rg^2-12 a^2 EEg rg^3+a(9 a EEg+Lzg) rg^4+EEg rg^5(-8+3rg));
	
	int3=1/(Sqrt[KKg] (-a^2+KKg) (a^2+rg^2) (1-z1g^2) z2g (a^2 (1-EEg^2)-z2g^2)) (2 a Lzg^2 (KKg+rg^2) (a^2 (1-EEg^2) (2ellPiz1g-ellK)-z2g^2 (ellE-ellK+(1+kzg)ellPiz1g)));
	
	int2=(2 Lzg )/(Sqrt[KKg] z2g) (((a^2+KKg) (EEg( a^2 -2 KKg)-a Lzred))/(KKg-a^2)^2-(2 a^3 Lzg)/(rg^2+a^2)^2+(a Lzg)/(rg^2+a^2)) ellPiz1g;
	
	int1=(4 a (EEg KKg+a Lzred)^2 EllipticPi[(a^2 z1g^2)/KKg,kzg])/((KKg-a^2)^2 Sqrt[KKg] z2g);
	
	int0=-((4a Rrg ellK)/(Sqrt[KKg] (KKg+rg^2) z2g \[CapitalDelta]rg^2))(EEg KKg(rg^2+a^2)-(2 EEg KKg+a Lzred) rg+rg^2(Lzg a-EEg rg));
	
	NIntegrate[1/Yrg (int5+int4+int3+int2+int1+int0),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


\[CapitalGamma]\[Phi]2avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r},
	If[a==0,
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		z1g=Sqrt[1-xg^2];
		
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
		
		r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
		r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
		z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
		kzg=a^2(1-EEg^2) z1g^2/z2g^2;
		\[Gamma]rg=-a^2z1g^2/rg^2;
		
		rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
		
		\[CapitalDelta]rg=(a^2-2 rg+rg^2);
		Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
		
		NIntegrate[( (-1+rg) rg a Yrg(r1g-rg)(rg-r2g))/(Sqrt[KKg] \[CapitalDelta]rg^2 (KKg+rg^2)) 4/(rg^2 z2g) ((rg^2+KKg)EllipticPi[\[Gamma]rg,kzg]-rg^2EllipticK[kzg]),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
	]
]


\[CapitalGamma]\[Phi]3avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,\[Gamma]rg,zg,Yzg,\[Chi]z},
	z1g=Sqrt[1-xg^2];
	If[(a==0||z1g==0),
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
		
		r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
		r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
		z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
		krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
		zg=z1g Sin[\[Chi]z];
		Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
		
		NIntegrate[(2 zg z1g^2 Cos[\[Chi]z]^2 Yzg )/(Sqrt[1-EEg^2] Sqrt[KKg] Sqrt[r1g-r3g] Sqrt[r2g-r4g](1-zg^2)^2 (r3g^2+a^2 zg^2)) ((2 a (KKg+r3g^2) zg EllipticK[krg])/  (KKg-a^2 zg^2)-(r2g-r3g) / (r2g^2+a^2 zg^2)  Re[(-I r2g r3g+a (r2g+r3g) zg+I a^2 zg^2) EllipticPi[((r1g-r2g) (r3g-I a zg))/((r1g-r3g) (r2g-I a zg)),krg]-(-I r2g r3g-a (r2g+r3g) zg+I  a^2 zg^2) EllipticPi[((r1g-r2g) (r3g+I a zg))/((r1g-r3g) (r2g+I a zg)),krg]]),{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]
	]
]


\[CapitalGamma]\[Phi]avg[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]rgfun[a,p,e,xg]/(2\[Pi]) \[CapitalUpsilon]zgfun[a,p,e,xg]/(2\[Pi]) (\[CapitalGamma]\[Phi]1avg[a,p,e,xg,prec]+\[CapitalGamma]\[Phi]2avg[a,p,e,xg,prec]+\[CapitalGamma]\[Phi]3avg[a,p,e,xg,prec])


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]\[Phi]sfunFC[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]\[Phi]sr0funFC[a,p,e,xg,prec]-\[CapitalUpsilon]\[Phi]sr1funFC[a,p,e,xg,prec]+\[CapitalUpsilon]\[Phi]sr2funFC[a,p,e,xg,prec]+\[CapitalUpsilon]\[Phi]sz0funFC[a,p,e,xg,prec]-\[CapitalUpsilon]\[Phi]sz1funFC[a,p,e,xg,prec]+\[CapitalUpsilon]\[Phi]sz2funFC[a,p,e,xg]+\[CapitalGamma]\[Phi]avg[a,p,e,xg,prec]


(* ::Subsection::Closed:: *)
(*Function for gauge shift between DH and FC frequencies*)


freqdevrad[a_,p_,e_,xg_]:=Module[{sgn,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]z,\[Chi]r,r1g,r2g,z1g,d\[CapitalUpsilon]rdr1g,d\[CapitalUpsilon]rdr2g,d\[CapitalUpsilon]rdz1g,r1p,r2p,z1p,const},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	sgn=RealSign[xg];
	
	const={EEgfun[a,p,e,xg],Lzgfun[a,p,e,xg],KKgfun[a,p,e,xg]};
	
	{d\[CapitalUpsilon]rdr1g,d\[CapitalUpsilon]rdr2g,d\[CapitalUpsilon]rdz1g}={D[\[CapitalUpsilon]rgfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]rgfunsgn[a,r1g,r2p,z1g,sgn],r2p],D[\[CapitalUpsilon]rgfunsgn[a,r1g,r2g,z1p,sgn],z1p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
	
	If[(a==0)||(z1g==0),
		d\[CapitalUpsilon]rdr1g r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]rdr2g r2savgFC[a,p,e,xg]
		,
		d\[CapitalUpsilon]rdr1g*r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]rdr2g*r2savgFC[a,p,e,xg]+d\[CapitalUpsilon]rdz1g*z1savgFC[a,p,e,xg]
	]
]


freqdevpol[a_,p_,e_,xg_]:=Module[{sgn,\[Chi]z,\[Chi]r,r1g,r2g,z1g,d\[CapitalUpsilon]zdr1g,d\[CapitalUpsilon]zdr2g,d\[CapitalUpsilon]zdz1g,r1p,r2p,z1p,const},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	sgn=RealSign[xg];
	
	const={EEgfun[a,p,e,xg],Lzgfun[a,p,e,xg],KKgfun[a,p,e,xg]};
	
	If[(a==0)||(z1g==0),
		{d\[CapitalUpsilon]zdr1g,d\[CapitalUpsilon]zdr2g}={D[\[CapitalUpsilon]zgfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]zgfunsgn[a,r1g,r2p,z1g,sgn],r2p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
		,
		{d\[CapitalUpsilon]zdr1g,d\[CapitalUpsilon]zdr2g,d\[CapitalUpsilon]zdz1g}={D[\[CapitalUpsilon]zgfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]zgfunsgn[a,r1g,r2p,z1g,sgn],r2p],D[\[CapitalUpsilon]zgfunsgn[a,r1g,r2g,z1p,sgn],z1p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
	];
	
	If[(a==0)||(z1g==0)
		,
		d\[CapitalUpsilon]zdr1g r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]zdr2g r2savgFC[a,p,e,xg]
		,
		d\[CapitalUpsilon]zdr1g*r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]zdr2g*r2savgFC[a,p,e,xg]+d\[CapitalUpsilon]zdz1g*z1savgFC[a,p,e,xg]
	]
]


freqdevtime[a_,p_,e_,xg_]:=Module[{sgn,\[Chi]z,\[Chi]r,r1g,r2g,z1g,d\[CapitalUpsilon]trdr1g,d\[CapitalUpsilon]trdr2g,d\[CapitalUpsilon]trdz1g,d\[CapitalUpsilon]tzdr1g,d\[CapitalUpsilon]tzdr2g,d\[CapitalUpsilon]tzdz1g,r1p,r2p,z1p,const},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	sgn=RealSign[xg];
	
	const={EEgfun[a,p,e,xg],Lzgfun[a,p,e,xg],KKgfun[a,p,e,xg]};
	
	If[(a==0)||(z1g==0),
		{d\[CapitalUpsilon]trdr1g,d\[CapitalUpsilon]trdr2g}={D[\[CapitalUpsilon]tgrfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]tgrfunsgn[a,r1g,r2p,z1g,sgn],r2p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
		{d\[CapitalUpsilon]tzdr1g,d\[CapitalUpsilon]tzdr2g}={D[-a^2 EEgfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[-a^2 EEgfunsgn[a,r1g,r2p,z1g,sgn],r2p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
		,
		{d\[CapitalUpsilon]trdr1g,d\[CapitalUpsilon]trdr2g,d\[CapitalUpsilon]trdz1g}={D[\[CapitalUpsilon]tgrfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]tgrfunsgn[a,r1g,r2p,z1g,sgn],r2p],D[\[CapitalUpsilon]tgrfunsgn[a,r1g,r2g,z1p,sgn],z1p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
		{d\[CapitalUpsilon]tzdr1g,d\[CapitalUpsilon]tzdr2g,d\[CapitalUpsilon]tzdz1g}={D[\[CapitalUpsilon]tgzfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]tgzfunsgn[a,r1g,r2p,z1g,sgn],r2p],D[\[CapitalUpsilon]tgzfunsgn[a,r1g,r2g,z1p,sgn],z1p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
	];
	
	If[(a==0)||(z1g==0)
		,
		(d\[CapitalUpsilon]trdr1g r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]trdr2g r2savgFC[a,p,e,xg]+d\[CapitalUpsilon]tzdr1g r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]tzdr2g r2savgFC[a,p,e,xg])
		,
		(d\[CapitalUpsilon]trdr1g+d\[CapitalUpsilon]tzdr1g)*r1savgFC[a,p,e,xg] +(d\[CapitalUpsilon]trdr2g+d\[CapitalUpsilon]tzdr2g)*r2savgFC[a,p,e,xg]+(d\[CapitalUpsilon]trdz1g+d\[CapitalUpsilon]tzdz1g)*z1savgFC[a,p,e,xg]
	]
]


freqdev\[Phi][a_,p_,e_,xg_]:=Module[{sgn,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]z,\[Chi]r,r1g,r2g,z1g,d\[CapitalUpsilon]\[Phi]rdr1g,d\[CapitalUpsilon]\[Phi]rdr2g,d\[CapitalUpsilon]\[Phi]rdz1g,d\[CapitalUpsilon]\[Phi]zdr1g,d\[CapitalUpsilon]\[Phi]zdr2g,d\[CapitalUpsilon]\[Phi]zdz1g,r1p,r2p,z1p,const},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	sgn=RealSign[xg];
	
	const={EEgfun[a,p,e,xg],Lzgfun[a,p,e,xg],KKgfun[a,p,e,xg]};
	
	If[(a==0)||(z1g==0),
		{d\[CapitalUpsilon]\[Phi]rdr1g,d\[CapitalUpsilon]\[Phi]rdr2g}={D[\[CapitalUpsilon]\[Phi]grfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]\[Phi]grfunsgn[a,r1g,r2p,z1g,sgn],r2p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
		{d\[CapitalUpsilon]\[Phi]zdr1g,d\[CapitalUpsilon]\[Phi]zdr2g}={D[\[CapitalUpsilon]\[Phi]gzfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]\[Phi]gzfunsgn[a,r1g,r2p,z1g,sgn],r2p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
		,
		{d\[CapitalUpsilon]\[Phi]rdr1g,d\[CapitalUpsilon]\[Phi]rdr2g,d\[CapitalUpsilon]\[Phi]rdz1g}={D[\[CapitalUpsilon]\[Phi]grfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]\[Phi]grfunsgn[a,r1g,r2p,z1g,sgn],r2p],D[\[CapitalUpsilon]\[Phi]grfunsgn[a,r1g,r2g,z1p,sgn],z1p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
		{d\[CapitalUpsilon]\[Phi]zdr1g,d\[CapitalUpsilon]\[Phi]zdr2g,d\[CapitalUpsilon]\[Phi]zdz1g}={D[\[CapitalUpsilon]\[Phi]gzfunsgn[a,r1p,r2g,z1g,sgn],r1p],D[\[CapitalUpsilon]\[Phi]gzfunsgn[a,r1g,r2p,z1g,sgn],r2p],D[\[CapitalUpsilon]\[Phi]gzfunsgn[a,r1g,r2g,z1p,sgn],z1p]}/.{r1p->r1g,r2p->r2g,z1p->z1g};
	];
	
	If[(a==0)||(z1g==0),
		d\[CapitalUpsilon]\[Phi]rdr1g r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]\[Phi]rdr2g r2savgFC[a,p,e,xg]+d\[CapitalUpsilon]\[Phi]zdr1g r1savgFC[a,p,e,xg]+d\[CapitalUpsilon]\[Phi]zdr2g r2savgFC[a,p,e,xg]
		,
		(d\[CapitalUpsilon]\[Phi]rdr1g+d\[CapitalUpsilon]\[Phi]zdr1g)r1savgFC[a,p,e,xg]+(d\[CapitalUpsilon]\[Phi]rdr2g+d\[CapitalUpsilon]\[Phi]zdr2g)*r2savgFC[a,p,e,xg]+(d\[CapitalUpsilon]\[Phi]rdz1g+d\[CapitalUpsilon]\[Phi]zdz1g)*z1savgFC[a,p,e,xg]
	]
]


(* ::Section::Closed:: *)
(*Shift frequencies - fixed turning points (on average)*)


(* ::Subsection::Closed:: *)
(*Shift radial frequency*)


(* ::Subsubsection::Closed:: *)
(*Part 1*)


\[ScriptCapitalY]1rDH[a_,p_,e_,xg_,\[Delta]const_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,\[Gamma]r1g,\[Gamma]r2g,ellK,ellPi,ellPir1g,ellPir2g,\[CapitalSigma]z1g,Rrg,Rr1g,Rr2g,Yrg,Yr1g,Yr2g,\[CapitalDelta]rg,\[CapitalDelta]r1g,\[CapitalDelta]r2g,rg,r1savgres,r2savgres,dVrgdEg,dVrgdLzg,dVrgdKg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	\[Gamma]r1g=-a^2z1g^2/r1g^2;
	\[Gamma]r2g=-a^2z1g^2/r2g^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	ellPir1g=EllipticPi[\[Gamma]r1g,kzg];
	ellPir2g=EllipticPi[\[Gamma]r2g,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Rr1g=(-a Lzg+EEg (a^2+r1g^2));
	Rr2g=(-a Lzg+EEg (a^2+r2g^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yr1g=Sqrt[(1-EEg^2) (-r1g+r3g) (-r1g+r4g)];
	Yr2g=Sqrt[(1-EEg^2) (-r2g+r3g) (-r2g+r4g)];
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	\[CapitalDelta]r1g=(a^2-2 r1g+r1g^2);
	\[CapitalDelta]r2g=(a^2-2 r2g+r2g^2);
	
	dVrgdEg=2(rg^2+a^2)Rrg;
	dVrgdLzg=-2a Rrg;
	dVrgdKg=-\[CapitalDelta]rg;
	4(Rrg(ellK rg^2-(KKg+rg^2)ellPi)(KKg(rg-1)+2a EEg (Lzg-a EEg)rg +a^2 rg+rg^2(2(1-EEg^2)rg-3)))/(Sqrt[KKg]z2g rg(KKg+rg^2)(r1g-rg)(rg-r2g)Yrg^3)+4(ellK/z2g)\[CapitalDelta]rg (Sqrt[KKg]Rrg+a(KKg+rg^2)RealSign[-a EEg+Lzg])/((r1g-rg)(rg-r2g)(KKg+rg^2)Yrg^3)+2(ellK/z2g)(\[Delta]const[[1]]dVrgdEg+\[Delta]const[[2]]dVrgdLzg+\[Delta]const[[3]]dVrgdKg)/((-rg+r1g)(rg-r2g) Yrg^3)
]


(* ::Subsubsection::Closed:: *)
(*Part 2*)


\[ScriptCapitalY]2r[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPi,\[CapitalSigma]z1g,Rrg,Yrg,rg,r1savgres,r2savgres},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	4/(Sqrt[KKg]Yrg z2g) (-EEg (KKg-rg^2)/(KKg+rg^2) ellK)+(4Rrg )/(Sqrt[KKg]Yrg z2g) (ellPi /rg^2 -(z2g^2 EllipticE[kzg])/(\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2))+ellK/\[CapitalSigma]z1g-((1-EEg^2)(3rg^2+2a^2z1g^2))/(\[CapitalSigma]z1g ((1-EEg^2)rg^2+z2g^2)) ellPi-((a^2 z1g^2 +2rg^2)z2g^2)/(rg^2\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2)) ellPi)
]


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]rsfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];

	((\[CapitalUpsilon]rg^2\[CapitalUpsilon]zg)/(2\[Pi])^2 NIntegrate[(\[ScriptCapitalY]1rDH[a,p,e,xg,\[Delta]const,\[Chi]r]+\[ScriptCapitalY]2r[a,p,e,xg,\[Chi]r]),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec])
]


(* ::Subsection::Closed:: *)
(*Shift polar frequency*)


(* ::Subsubsection::Closed:: *)
(*Part 1*)


\[ScriptCapitalY]1zDH[a_,p_,e_,xg_,\[Delta]const_,\[Chi]z_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,ellK,Zzg,zg,Yzg,argzg,argzgcc,Lzred,aux,dVzgdEg,dVzgdLzg,dVzgdKg},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
	
	ellK=EllipticK[krg];
	
	zg=z1g Sin[\[Chi]z];
	
	Zzg=(Lzg+a EEg(-1+zg^2));
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	
	argzg=(r3g-I a zg)/(r2g-I a zg) (r1g-r2g)/(r1g-r3g);
	argzgcc=(r3g+I a zg)/(r2g+I a zg) (r1g-r2g)/(r1g-r3g);
	Lzred=Lzg-a EEg;
	aux=a^2(1-EEg^2)+Lzg^2;
	
	If[z1g==0,
		(a RealSign[Lzg-a EEg])/aux^(1/2) (Lzg /(aux Lzred)+1/(3 r1g^2 r2g  r3g^2) (-2(r1g-r3g) (r2g r3g+r1g (r2g+r3g)) EllipticE[krg]/EllipticK[krg]+r1g ((r2g-r3g) r3g+r1g (2 r2g+r3g))))+RealSign[Lzg]/(2*z2g*z2g^2)(2a(-2 a EEg+Lzg)\[Delta]const[[1]]+2a*EEg*\[Delta]const[[2]]+\[Delta]const[[3]])
		,
		Sec[\[Chi]z]^2/Yzg^3((a(1-zg^2))/z1g^2 ((Sqrt[KKg] Zzg)/(KKg-a^2zg^2)-RealSign[-a EEg+Lzg])+(Sin[\[Chi]z]Zzg(a^2(1-EEg^2)(2 zg^2-z1g^2)-z2g^2))/(Sqrt[KKg](r2g^2+a^2 zg^2)(r3g^2+a^2 zg^2)) ((r2g-r3g) /(2ellK z1g) Re[EllipticPi[argzgcc,krg](-I r2g r3g-a(r2g+r3g)zg+I a^2 zg^2)+EllipticPi[argzg,krg](I r2g r3g-a(r2g+r3g)zg-I a^2 zg^2)]+ a (KKg+r3g^2)/(KKg-a^2 zg^2)(r2g^2+a^2 zg^2)Sin[\[Chi]z]))+(\[Delta]const[[1]]dVzgdEg+\[Delta]const[[2]]dVzgdLzg+\[Delta]const[[3]]dVzgdKg)/(2z1g^2*Cos[\[Chi]z]^2*Yzg^3)
	]
]


(* ::Subsubsection::Closed:: *)
(*Part 2*)


\[ScriptCapitalY]2z[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPi,\[CapitalSigma]z1g,Rrg,Yrg,rg,r1savgres,r2savgres},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	4/(Sqrt[KKg]Yrg z2g) (-2EEg EllipticPi[(a^2z1g^2)/KKg,kzg] +1/rg^2 (Rrg+2rg^2 EEg)ellPi-( z2g^2Rrg EllipticE[kzg])/(\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2))+(Rrg ellK)/\[CapitalSigma]z1g-((1-EEg^2)Rrg ellPi)/(\[CapitalSigma]z1g((1-EEg^2)rg^2+z2g^2)) (3rg^2+2a^2z1g^2)-(z2g^2Rrg (2rg^2+a^2 z1g^2) ellPi)/(\[CapitalSigma]z1g rg^2((1-EEg^2)rg^2+z2g^2)))
]


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]rsfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];

	\[CapitalUpsilon]zg^2/(2\[Pi])(NIntegrate[\[ScriptCapitalY]1zDH[a,p,e,xg,\[Delta]const,\[Chi]z],{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]rg/(2\[Pi])NIntegrate[\[ScriptCapitalY]2z[a,p,e,xg,\[Chi]r],{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec])
]		


(* ::Subsection::Closed:: *)
(*Shift "coordinate time frequency"*)


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]rs/\[CapitalUpsilon]rg \[CapitalUpsilon]^t_{rg}\[RightAngleBracket] *)


\[CapitalUpsilon]tsr0funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=\[CapitalUpsilon]rsfunDH[a,p,e,xg,\[Delta]const,prec]/\[CapitalUpsilon]rgfun[a,p,e,xg] \[CapitalUpsilon]tgrfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]zs/\[CapitalUpsilon]zg \[CapitalUpsilon]^t_{zg}\[RightAngleBracket] *)


\[CapitalUpsilon]tsz0funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=\[CapitalUpsilon]zsfunDH[a,p,e,xg,\[Delta]const,prec]/\[CapitalUpsilon]zgfun[a,p,e,xg] \[CapitalUpsilon]tgzfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*Radial part  \[LeftAngleBracket]\[PartialD]V^t_{rg}/\[PartialD]Cig Cis\[RightAngleBracket]*)


\[CapitalUpsilon]tsrCisfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,\[CapitalDelta]rg,rg,Yrg,dVtrgdEg,dVtrgdLzg,\[CapitalUpsilon]rg,\[Chi]r},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];

	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];

	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	\[CapitalDelta]rg=rg^2-2rg+a^2;

	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];

	dVtrgdEg=(a^2+rg^2)^2/\[CapitalDelta]rg;
	dVtrgdLzg=-2a*rg/\[CapitalDelta]rg;

	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	
	\[CapitalUpsilon]rg/(2\[Pi])NIntegrate[1/Yrg(\[Delta]const[[1]]dVtrgdEg+\[Delta]const[[2]]dVtrgdLzg),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


(* ::Subsubsection::Closed:: *)
(*Polar part  \[LeftAngleBracket]\[PartialD]V^t_{zg}/\[PartialD]Cig Cis\[RightAngleBracket]*)


\[CapitalUpsilon]tszCisfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{z1g,zg,EEg,Lzg,z2g,dVtzgdEg,Yzg,\[CapitalUpsilon]zg,\[Chi]z},
	If[a==0,
		0,
		z1g=Sqrt[1-xg^2];
		zg=z1g Sin[\[Chi]z];

		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];

		z2g=Sqrt[a^2 (1-EEg^2)+Lzg^2/(1-z1g^2)];

		dVtzgdEg=-a^2 (1-zg^2);
		Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];

		\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];

		\[CapitalUpsilon]zg/(2\[Pi])NIntegrate[(1/Yzg)\[Delta]const[[1]]dVtzgdEg,{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]
		]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^t_{yg}\[RightAngleBracket] - radial part*)


\[CapitalUpsilon]tsr1funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,\[CapitalDelta]rg,Rrg,rg,Vtrg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	\[CapitalDelta]rg=rg^2-2rg+a^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	Vtrg=a Lzg+(Rrg(a^2+rg^2) )/\[CapitalDelta]rg;
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	(\[CapitalUpsilon]rg*\[CapitalUpsilon]zg)/(2\[Pi])^2 NIntegrate[(\[ScriptCapitalY]1rDH[a,p,e,xg,\[Delta]const,\[Chi]r]+\[ScriptCapitalY]2r[a,p,e,xg,\[Chi]r])Vtrg,{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^t_{yg}\[RightAngleBracket]  - polar part*)


\[CapitalUpsilon]tsz11funDH[a_,p_,e_,xg_,\[Delta]const_,\[Chi]z_]:=Module[{z1g,zg,Vtzg},
	z1g=Sqrt[1-xg^2];
	zg=z1g Sin[\[Chi]z];
	Vtzg=-a^2EEgfun[a,p,e,xg](1-zg^2);
	
	\[ScriptCapitalY]1zDH[a,p,e,xg,\[Delta]const,\[Chi]z]Vtzg
]


\[CapitalUpsilon]tsz12funDH[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPi,\[CapitalSigma]z1g,Rrg,Yrg,rg,Lzred},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Lzred=Lzg-a EEg;
	
	(2 EEg)/(Sqrt[KKg] Yrg z2g) (-2 (a Lzred-3 EEg rg^2) ellK+(2 (a^2+rg^2) Rrg (z2g^2 (EllipticE[kzg]-ellK)-(1-EEg^2) rg^2 ellK))/(\[CapitalSigma]z1g((1-EEg^2) rg^2+z2g^2))+4 EEg (KKg ellK+(a^2-KKg) EllipticPi[(a^2 z1g^2)/KKg,kzg]))+1/(Sqrt[KKg] Yrg z2g) (4 EEg )/rg^2 ((a^2+rg^2) (a Lzred-3 EEg  rg^2)-(Rrg ((-1+EEg^2) rg^2 (3 a^2  rg^2+rg^4+2 a^4 z1g^2)+a^2 (-a^2 z1g^2+rg^2 (-2+z1g^2))z2g^2))/(\[CapitalSigma]z1g  ((1-EEg^2) rg^2+z2g^2))) EllipticPi[\[Gamma]rg,kzg]
]


\[CapitalUpsilon]tsz1funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	If[a==0,
		0,
		\[CapitalUpsilon]zg/(2\[Pi])(NIntegrate[\[CapitalUpsilon]tsz11funDH[a,p,e,xg,\[Delta]const,\[Chi]z],{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]rg/(2\[Pi]) NIntegrate[\[CapitalUpsilon]tsz12funDH[a,p,e,xg,\[Chi]r],{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec])
	]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[CapitalSigma]/2 \[CapitalGamma]^t_{12}\[RightAngleBracket]*)


\[CapitalGamma]t1avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,ellK,ellPi,Lzred,Rrg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r,int2,int1,int0},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	ellK= EllipticK[kzg];
	ellPi=EllipticPi[\[Gamma]rg,kzg];
	Lzred=Lzg-a EEg;
	
	int2=((1-EEg^2) rg^4ellK-(1-EEg^2) rg^2 (3 rg^2 +2 a^2 z1g^2) ellPi-z2g^2 (rg^2 EllipticE[kzg]-rg^2 ellK+(2 rg^2 +a^2 z1g^2)ellPi)) (8Rrg^2)/(Sqrt[KKg] rg z2g ((1-EEg^2) rg^2+z2g^2) \[CapitalDelta]rg (rg^2+a^2 z1g^2));
	int1=(4Rrg ellPi)/(Sqrt[KKg] rg z2g \[CapitalDelta]rg^2) (-3a^3 Lzred+4a Lzred rg+a(8 a EEg-Lzg) rg^2+EEg rg^3(-12 +5 rg));
	int0=-((4rg Rrg ellK)/(Sqrt[KKg] (KKg+rg^2) z2g \[CapitalDelta]rg^2))(-a^3Lzred+2EEg KKg(a^2 -2 rg)+(2 a^2 EEg+2 EEg KKg+a Lzg) rg^2 +rg^3EEg(rg-4));
	
	NIntegrate[1/Yrg (int2+int1+int0),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


\[CapitalGamma]t2avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	NIntegrate[((a^2-rg^2) (r1g-r2g)^2 Yrg Cos[\[Chi]r]^2 (rg^2 EllipticK[kzg]-(KKg+rg^2) EllipticPi[\[Gamma]rg,kzg]))/(Sqrt[KKg] rg \[CapitalDelta]rg^2(KKg+rg^2) z2g),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


\[CapitalGamma]tavg[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]rgfun[a,p,e,xg]/(2\[Pi]) \[CapitalUpsilon]zgfun[a,p,e,xg]/(2\[Pi]) (\[CapitalGamma]t1avg[a,p,e,xg,prec]+\[CapitalGamma]t2avg[a,p,e,xg,prec])


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]tsfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=\[CapitalUpsilon]tsr0funDH[a,p,e,xg,\[Delta]const,prec]-\[CapitalUpsilon]tsr1funDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalUpsilon]tsz0funDH[a,p,e,xg,\[Delta]const,prec]-\[CapitalUpsilon]tsz1funDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalUpsilon]tsrCisfunDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalUpsilon]tszCisfunDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalGamma]tavg[a,p,e,xg,prec]


(* ::Subsection::Closed:: *)
(*Shift azimuthal frequency*)


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]rs/\[CapitalUpsilon]rg \[CapitalUpsilon]^\[Phi]_{rg}\[RightAngleBracket] *)


\[CapitalUpsilon]\[Phi]sr0funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=\[CapitalUpsilon]rsfunDH[a,p,e,xg,\[Delta]const,prec]/\[CapitalUpsilon]rgfun[a,p,e,xg] \[CapitalUpsilon]\[Phi]grfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*\[LeftAngleBracket]\[CapitalUpsilon]zs/\[CapitalUpsilon]zg \[CapitalUpsilon]^\[Phi]_{zg}\[RightAngleBracket] *)


\[CapitalUpsilon]\[Phi]sz0funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=\[CapitalUpsilon]zsfunDH[a,p,e,xg,\[Delta]const,prec]/\[CapitalUpsilon]zgfun[a,p,e,xg] \[CapitalUpsilon]\[Phi]gzfun[a,p,e,xg]


(* ::Subsubsection::Closed:: *)
(*Radial part  \[LeftAngleBracket]\[PartialD]V^\[Phi]_{rg}/\[PartialD]Cig Cis\[RightAngleBracket]*)


\[CapitalUpsilon]\[Phi]srCisfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,\[CapitalDelta]rg,rg,Yrg,dV\[Phi]rgdEg,dV\[Phi]rgdLzg,\[CapitalUpsilon]rg,\[Chi]r},
	If[a==0,
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		z1g=Sqrt[1-xg^2];
	
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
	
		r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
		r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
		rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
		
		Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
		\[CapitalDelta]rg=rg^2-2rg+a^2;

		Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];

		dV\[Phi]rgdEg=2a*rg/\[CapitalDelta]rg;
		dV\[Phi]rgdLzg=-a^2/\[CapitalDelta]rg;

		\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];

		\[CapitalUpsilon]rg/(2\[Pi])NIntegrate[1/Yrg(\[Delta]const[[1]]dV\[Phi]rgdEg+\[Delta]const[[2]]dV\[Phi]rgdLzg),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
	]
]


(* ::Subsubsection::Closed:: *)
(*Polar part   \[LeftAngleBracket]\[PartialD]V^\[Phi]_{zg}/\[PartialD]Cig Cis\[RightAngleBracket]*)


\[CapitalUpsilon]\[Phi]szCisfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{z1g,zg,EEg,Lzg,z2g,dV\[Phi]zgdLzg,Yzg,\[CapitalUpsilon]zg,\[Chi]z},
	z1g=Sqrt[1-xg^2];
	zg=z1g Sin[\[Chi]z];

	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];

	z2g=Sqrt[a^2 (1-EEg^2)+Lzg^2/(1-z1g^2)];

	dV\[Phi]zgdLzg=1/(1-zg^2);
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];

	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];

	\[CapitalUpsilon]zg/(2\[Pi])NIntegrate[1/Yzg*\[Delta]const[[2]]dV\[Phi]zgdLzg,{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^\[Phi]_{yg}\[RightAngleBracket] - radial part*)


\[CapitalUpsilon]\[Phi]sr1funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,\[CapitalDelta]rg,Rrg,rg,V\[Phi]rg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r},
	If[a==0,
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		z1g=Sqrt[1-xg^2];
	
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
	
		\[CapitalDelta]rg=rg^2-2rg+a^2;
		Rrg=(-a Lzg+EEg (a^2+rg^2));
	
		rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
		V\[Phi]rg=(a Rrg)/\[CapitalDelta]rg-a EEg;
		\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
		\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
		(\[CapitalUpsilon]rg \[CapitalUpsilon]zg)/(2\[Pi])^2 NIntegrate[(\[ScriptCapitalY]1rDH[a,p,e,xg,\[Delta]const,\[Chi]r]+\[ScriptCapitalY]2r[a,p,e,xg,\[Chi]r])V\[Phi]rg,{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
	]
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[Delta]Yys/Yyg  V^\[Phi]_{yg}\[RightAngleBracket] - polar part*)


\[CapitalUpsilon]\[Phi]sz11funDH[a_,p_,e_,xg_,\[Delta]const_,\[Chi]z_]:=Module[{z1g,zg,V\[Phi]zg},
	z1g=Sqrt[1-xg^2];
	zg=z1g Sin[\[Chi]z];
	V\[Phi]zg=Lzgfun[a,p,e,xg] 1/(1-zg^2);
	
	\[ScriptCapitalY]1zDH[a,p,e,xg,\[Delta]const,\[Chi]z]V\[Phi]zg
]


\[CapitalUpsilon]\[Phi]sz12funDH[a_,p_,e_,xg_,\[Chi]r_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,kzg,\[Gamma]rg,ellK,ellPiz1g,\[CapitalSigma]z1g,Rrg,Yrg,rg,Lzred},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	ellK=EllipticK[kzg];
	ellPiz1g=EllipticPi[z1g^2,kzg];
	
	\[CapitalSigma]z1g=rg^2+a^2z1g^2;
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Lzred=Lzg-a EEg;
	
	(2 Lzg)/(Sqrt[KKg] Yrg z2g) 1/(a^2+rg^2) (-(2 a^2 Rrg (z2g^2 EllipticE[kzg]-((1-EEg^2) rg^2+z2g^2) ellK))/(\[CapitalSigma]z1g((1-EEg^2) rg^2+z2g^2))-(4 rg^2 Rrg ellPiz1g)/(a^2+rg^2)+2 (Rrg+2 EEg rg^2) ellPiz1g+(4 EEg(a^2+rg^2) (KKg ellPiz1g-a^2EllipticPi[(a^2 z1g^2)/KKg,kzg]))/(a^2-KKg))+(2 Lzg )/(Sqrt[KKg] Yrg z2g) (2 a^2)/(rg^2 (a^2+rg^2)) ((Rrg+2 EEg rg^2)-(Rrg (5 (1-EEg^2) rg^6+4 rg^4 z2g^2+a^4 z1g^2 (2 (1-EEg^2) rg^2+z2g^2)+a^2 ((1-EEg^2) rg^4 (3+4 z1g^2)+rg^2 (2+3 z1g^2) z2g^2)))/((a^2+rg^2)\[CapitalSigma]z1g((1-EEg^2) rg^2+z2g^2))) EllipticPi[\[Gamma]rg,kzg]
]


\[CapitalUpsilon]\[Phi]sz1funDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=Module[{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	\[CapitalUpsilon]zg/(2\[Pi])(NIntegrate[\[CapitalUpsilon]\[Phi]sz11funDH[a,p,e,xg,\[Delta]const,\[Chi]z],{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]+\[CapitalUpsilon]rg/(2\[Pi])NIntegrate[\[CapitalUpsilon]\[Phi]sz12funDH[a,p,e,xg,\[Chi]r],{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec])
]


(* ::Subsubsection::Closed:: *)
(*Semi-analytic expressions  \[LeftAngleBracket]\[CapitalSigma]/2 \[CapitalGamma]^\[Phi]_{12}\[RightAngleBracket]*)


\[CapitalGamma]\[Phi]1avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,ellK,ellPiz1g,ellPirg,ellE,Lzred,Rrg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r,int5,int4,int3,int2,int1,int0},
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	
	EEg=EEgfun[a,p,e,xg];
	Lzg=Lzgfun[a,p,e,xg];
	KKg=KKgfun[a,p,e,xg];
	
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
	kzg=a^2(1-EEg^2) z1g^2/z2g^2;
	\[Gamma]rg=-a^2z1g^2/rg^2;
	
	rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	ellK= EllipticK[kzg];
	ellPiz1g=EllipticPi[z1g^2,kzg];
	ellPirg=EllipticPi[\[Gamma]rg,kzg];
	ellE=EllipticE[kzg];
	Lzred=Lzg-a EEg;
	
	int5=(8a Rrg^2 )/(Sqrt[KKg] rg (a^2+rg^2) z2g ((1-EEg^2)rg^2+z2g^2) \[CapitalDelta]rg (rg^2+a^2z1g^2)) ((1-EEg^2) rg^4 ellK-(1-EEg^2) rg^2 (3 rg^2+2 a^2 z1g^2)ellPirg-z2g^2 (rg^2 (ellE-ellK)+(2 rg^2+a^2 z1g^2) ellPirg));
	
	int4=(4a Rrg ellPirg)/(Sqrt[KKg] rg (a^2+rg^2)^2 z2g \[CapitalDelta]rg^2) (-3a^5Lzred+4a^3Lzred rg+a^3(9 a EEg-2  Lzg) rg^2-12 a^2 EEg rg^3+a(9 a EEg+Lzg) rg^4+EEg rg^5(-8+3rg));
	
	int3=1/(Sqrt[KKg] (-a^2+KKg) (a^2+rg^2) (1-z1g^2) z2g (a^2 (1-EEg^2)-z2g^2)) (2 a Lzg^2 (KKg+rg^2) (a^2 (1-EEg^2) (2ellPiz1g-ellK)-z2g^2 (ellE-ellK+(1+kzg)ellPiz1g)));
	
	int2=(2 Lzg )/(Sqrt[KKg] z2g) (((a^2+KKg) (EEg( a^2 -2 KKg)-a Lzred))/(KKg-a^2)^2-(2 a^3 Lzg)/(rg^2+a^2)^2+(a Lzg)/(rg^2+a^2)) ellPiz1g;
	
	int1=(4 a (EEg KKg+a Lzred)^2 EllipticPi[(a^2 z1g^2)/KKg,kzg])/((KKg-a^2)^2 Sqrt[KKg] z2g);
	
	int0=-((4a Rrg ellK)/(Sqrt[KKg] (KKg+rg^2) z2g \[CapitalDelta]rg^2))(EEg KKg(rg^2+a^2)-(2 EEg KKg+a Lzred) rg+rg^2(Lzg a-EEg rg));
	
	NIntegrate[1/Yrg (int5+int4+int3+int2+int1+int0),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
]


\[CapitalGamma]\[Phi]2avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,kzg,\[Gamma]rg,zg,rg,\[CapitalDelta]rg,Yrg,\[Chi]r},
	If[a==0,
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		z1g=Sqrt[1-xg^2];
		
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
		
		r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
		r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
		z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
		kzg=a^2(1-EEg^2) z1g^2/z2g^2;
		\[Gamma]rg=-a^2z1g^2/rg^2;
		
		rg=(r1g+r2g)/2+1/2 (r1g-r2g) Sin[\[Chi]r];
		
		\[CapitalDelta]rg=(a^2-2 rg+rg^2);
		Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
		
		NIntegrate[( (-1+rg) rg a Yrg(r1g-rg)(rg-r2g))/(Sqrt[KKg] \[CapitalDelta]rg^2 (KKg+rg^2)) 4/(rg^2 z2g) ((rg^2+KKg)EllipticPi[\[Gamma]rg,kzg]-rg^2EllipticK[kzg]),{\[Chi]r,0,2\[Pi]},WorkingPrecision->prec]
	]
]


\[CapitalGamma]\[Phi]3avg[a_,p_,e_,xg_,prec_]:=Module[{r1g,r2g,z1g,EEg,Lzg,KKg,r3g,r4g,z2g,krg,\[Gamma]rg,zg,Yzg,\[Chi]z},
	z1g=Sqrt[1-xg^2];
	If[(a==0||z1g==0),
		0,
		r1g=p/(1-e);
		r2g=p/(1+e);
		
		EEg=EEgfun[a,p,e,xg];
		Lzg=Lzgfun[a,p,e,xg];
		KKg=KKgfun[a,p,e,xg];
		
		r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
		r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
		z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
		krg=((r1g-r2g) (r3g-r4g))/((r1g-r3g) (r2g-r4g));
		zg=z1g Sin[\[Chi]z];
		Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
		
		NIntegrate[(2 zg z1g^2 Cos[\[Chi]z]^2 Yzg )/(Sqrt[1-EEg^2] Sqrt[KKg] Sqrt[r1g-r3g] Sqrt[r2g-r4g](1-zg^2)^2 (r3g^2+a^2 zg^2)) ((2 a (KKg+r3g^2) zg EllipticK[krg])/  (KKg-a^2 zg^2)-(r2g-r3g) / (r2g^2+a^2 zg^2)  Re[(-I r2g r3g+a (r2g+r3g) zg+I a^2 zg^2) EllipticPi[((r1g-r2g) (r3g-I a zg))/((r1g-r3g) (r2g-I a zg)),krg]-(-I r2g r3g-a (r2g+r3g) zg+I  a^2 zg^2) EllipticPi[((r1g-r2g) (r3g+I a zg))/((r1g-r3g) (r2g+I a zg)),krg]]),{\[Chi]z,0,2\[Pi]},WorkingPrecision->prec]
	]
]


\[CapitalGamma]\[Phi]avg[a_,p_,e_,xg_,prec_]:=\[CapitalUpsilon]rgfun[a,p,e,xg]/(2\[Pi]) \[CapitalUpsilon]zgfun[a,p,e,xg]/(2\[Pi]) (\[CapitalGamma]\[Phi]1avg[a,p,e,xg,prec]+\[CapitalGamma]\[Phi]2avg[a,p,e,xg,prec]+\[CapitalGamma]\[Phi]3avg[a,p,e,xg,prec])


(* ::Subsubsection::Closed:: *)
(*Full expression*)


\[CapitalUpsilon]\[Phi]sfunDH[a_,p_,e_,xg_,\[Delta]const_,prec_]:=\[CapitalUpsilon]\[Phi]sr0funDH[a,p,e,xg,\[Delta]const,prec]-\[CapitalUpsilon]\[Phi]sr1funDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalUpsilon]\[Phi]sz0funDH[a,p,e,xg,\[Delta]const,prec]-\[CapitalUpsilon]\[Phi]sz1funDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalUpsilon]\[Phi]srCisfunDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalUpsilon]\[Phi]szCisfunDH[a,p,e,xg,\[Delta]const,prec]+\[CapitalGamma]\[Phi]avg[a,p,e,xg,prec]


(* ::Section::Closed:: *)
(*\[Xi]r and \[Xi]z shifts for near-identity transformation*)


(* ::Subsection::Closed:: *)
(*Functions for \[Xi]r - parallel component of the spin *)


\[Xi]rpar[wr_,wz_,\[CapitalUpsilon]rg_,\[CapitalUpsilon]zg_,coeff_]:=Module[{dimr,dimz},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	
	-2\[CapitalUpsilon]rg Sum[Sin[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]]/(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg),{n,1,dimr},{k,1,dimz}]-2\[CapitalUpsilon]rg Sum[Sin[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]]/(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg),{n,1,dimr},{k,-dimz,-1}]-2Sum[Sin[n wr]coeff[[n+dimr+1,dimz+1]]/n,{n,1,dimr}]-2(\[CapitalUpsilon]rg/\[CapitalUpsilon]zg) Sum[Sin[k wz]coeff[[dimr+1,k+dimz+1]]/k,{k,1,dimz}]
];


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


\[Delta]YroverYrgsparFC[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],Rrg,Rr1g,Rr2g,r1spar,r2spar,fraux2,fraux1,fraux0,fraux,hr},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg rg^2);
	Rr1g=(a^2 EEg-a Lzg+EEg r1g^2);
	Rr2g=(a^2 EEg-a Lzg+EEg r2g^2);
	
	r1spar=r1sparFCfun[wz,a,p,e,xg,const];
	r2spar=r2sparFCfun[wz,a,p,e,xg,const];
	
	fraux2=Yrg/Sqrt[KKg] (EEg-(2 EEg KKg)/(KKg+rg^2)-(2 rg^2Rrg)/\[CapitalSigma]^2+Rrg/\[CapitalSigma]);
	fraux1=(2 a^2 rg Rrg zg z1g Yzg Cos[\[Chi]z])/(Sqrt[KKg] \[CapitalSigma]^2Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(-rg+r1g) (rg-r2g)]  );
	fraux0=(Sqrt[KKg] \[CapitalDelta] Rrg)/((KKg+rg^2) (-rg+r1g) (rg-r2g) Yrg)+(a \[CapitalDelta]  RealSign[Lzg-a EEg])/((-rg+r1g) (rg-r2g)Yrg)+(rg Rrg (KKg (-1+rg)+a (a-2 a EEg^2+2 EEg Lzg) rg+rg^2 (-3+2 (1-EEg^2) rg)) (KKg-a^2 zg^2))/(Sqrt[KKg]\[CapitalSigma] (KKg+rg^2) (rg-r1g) (rg-r2g) Yrg )-(r1spar-r2spar +(r1spar+r2spar)Sin\[Chi]rfun[wr,a,p,e,xg,const]) ((r1g-r2g)Yrg)/(4(r1g-rg)(rg-r2g));
	fraux=((r1spar+r2spar) /2+(r1spar-r2spar )/2 Sin\[Chi]rfun[wr,a,p,e,xg,const]) ((1-EEg^2)(2 rg-r3g-r4g))/(2Yrg);
	hr=z1g Cos[\[Chi]z]Yzg (a^2zg)/Sqrt[KKg] ((r1g Rr1g)/ (r1g^2+a^2 zg^2)^2+(r2g Rr2g)/(r2g^2+a^2 zg^2)^2+((r1g Rr1g)/(r1g^2+a^2 zg^2)^2-(r2g Rr2g )/ (r2g^2+a^2 zg^2)^2)Sin\[Chi]rfun[wr,a,p,e,xg,const]);
	
	1/Yrg (fraux0+fraux1+fraux2+fraux-hr/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(-rg+r1g) (rg-r2g)]))
]


\[Delta]YroverYrgsparFCcoeff[nmax_,kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[\[Delta]YroverYrgsparFC[wrlist[[i]],wzlist,a,p,e,xg,const],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fixed turning points (on average)*)


\[Delta]YroverYrgsparDH[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]z,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],Rrg,Rr1g,Rr2g,r1spar,r2spar,dVrgdEg,dVrgdLzg,dVrgdKg,fraux2,fraux1,fraux0,fraux,fauxdC,hr},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg rg^2);
	Rr1g=(a^2 EEg-a Lzg+EEg r1g^2);
	Rr2g=(a^2 EEg-a Lzg+EEg r2g^2);
	
	r1spar=r1sparDHfun[wz,a,p,e,xg,const,\[Delta]const];
	r2spar=r2sparDHfun[wz,a,p,e,xg,const,\[Delta]const];
	
	dVrgdEg=2(rg^2+a^2)Rrg;
	dVrgdLzg=-2a Rrg;
	dVrgdKg=-\[CapitalDelta];
	
	fraux2=Yrg/Sqrt[KKg] (EEg-(2 EEg KKg)/(KKg+rg^2)-(2 rg^2Rrg)/\[CapitalSigma]^2+Rrg/\[CapitalSigma]);
	fraux1=(2 a^2 rg Rrg zg z1g Yzg Cos[\[Chi]z])/(Sqrt[KKg] \[CapitalSigma]^2Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(-rg+r1g) (rg-r2g)]  );
	fraux0=(Sqrt[KKg] \[CapitalDelta] Rrg)/((KKg+rg^2) (-rg+r1g) (rg-r2g) Yrg)+(a \[CapitalDelta]  RealSign[Lzg-a EEg])/((-rg+r1g) (rg-r2g)Yrg)+(rg Rrg (KKg (-1+rg)+a (a-2 a EEg^2+2 EEg Lzg) rg+rg^2 (-3+2 (1-EEg^2) rg)) (KKg-a^2 zg^2))/(Sqrt[KKg]\[CapitalSigma] (KKg+rg^2) (rg-r1g) (rg-r2g) Yrg )-(r1spar-r2spar +(r1spar+r2spar)Sin\[Chi]rfun[wr,a,p,e,xg,const]) ((r1g-r2g)Yrg)/(4(r1g-rg)(rg-r2g));
	fraux=((r1spar+r2spar) /2+(r1spar-r2spar )/2 Sin\[Chi]rfun[wr,a,p,e,xg,const]) ((1-EEg^2)(2 rg-r3g-r4g))/(2Yrg);
	fauxdC=(\[Delta]const[[1]]dVrgdEg+\[Delta]const[[2]]dVrgdLzg+\[Delta]const[[3]]dVrgdKg)/(2(-rg+r1g) (rg-r2g) Yrg);
	hr=z1g Cos[\[Chi]z]Yzg (a^2zg)/Sqrt[KKg] ((r1g Rr1g)/ (r1g^2+a^2 zg^2)^2+(r2g Rr2g)/(r2g^2+a^2 zg^2)^2+((r1g Rr1g)/(r1g^2+a^2 zg^2)^2-(r2g Rr2g )/ (r2g^2+a^2 zg^2)^2)Sin\[Chi]rfun[wr,a,p,e,xg,const]);
	
	1/Yrg (fraux0+fraux1+fraux2+fraux+fauxdC-hr/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(-rg+r1g) (rg-r2g)]))
]


\[Delta]YroverYrgsparDHcoeff[nmax_,kmax_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[\[Delta]YroverYrgsparDH[wrlist[[i]],wzlist,a,p,e,xg,const,\[Delta]const],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsection::Closed:: *)
(*Functions for \[Xi]r - parallel component of the spin*)


\[Xi]zpar[wr_,wz_,\[CapitalUpsilon]rg_,\[CapitalUpsilon]zg_,coeff_]:=Module[{dimr,dimz},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	
	-2\[CapitalUpsilon]zg Sum[Sin[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]]/(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg),{n,1,dimr},{k,1,dimz}]-2\[CapitalUpsilon]zg Sum[Sin[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]]/(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg),{n,1,dimr},{k,-dimz,-1}]-2(\[CapitalUpsilon]zg/\[CapitalUpsilon]rg)Sum[Sin[n wr]coeff[[n+dimr+1,dimz+1]]/n,{n,1,dimr}]-2 Sum[Sin[k wz]coeff[[dimr+1,k+dimz+1]]/k,{k,1,dimz}]
];


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


\[Delta]YzoverYzgsparFC[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],Rrg,Rr1g,Rr2g,Zzg,z1s,fzaux2,fzaux1,fzaux0,fzaux,hz},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg rg^2);
	Rr1g=(a^2 EEg-a Lzg+EEg r1g^2);
	Rr2g=(a^2 EEg-a Lzg+EEg r2g^2);
	Zzg=(Lzg-a EEg (1-zg^2));
	
	z1s=z1sparFCfun[wr,a,p,e,xg,const];
	
	fzaux2=(Yzg (-KKg rg^2 Rrg+a^2 (KKg+rg^2)(Rrg-2EEg rg^2) zg^2- (a^2 EEg-a Lzg+2 EEg KKg+3 EEg rg^2) a^4zg^4))/(Sqrt[KKg] (KKg-a^2 zg^2) \[CapitalSigma]^2);
	fzaux1=(2 a rg)/(Sqrt[KKg] \[CapitalSigma]^2 ) Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)] Yrg Zzg Tan[\[Chi]z];
	
	(*fzaux0 is the only term that requires special care for z1g->0*)
	fzaux0=If[z1g==0,
	RealSign[Lzg]((a^2(1-EEg^2)+a EEg Lzg)/(-a EEg+Lzg)^2-((a^2 (1-EEg^2)+Lzg^2) ((-a EEg+Lzg)^2+rg^2 ) )/((-a EEg+Lzg)^2 rg^2 )) (a Tan[\[Chi]z]^2)/z2g-RealSign[Lzg]a (Lzg/((a EEg-Lzg) (a^2 (1-EEg^2)+Lzg^2))-1/rg^2) z2g/Cos[\[Chi]z]^2
	,
	(a Sqrt[KKg](1-zg^2)Zzg )/((KKg-a^2 zg^2) z1g^2Cos[\[Chi]z]^2Yzg )+(a (KKg+rg^2)Zzg (a^2 (1-EEg^2) (2 zg^2-z1g^2)-z2g^2) )/(Sqrt[KKg] (KKg-a^2zg^2) \[CapitalSigma] Yzg) Tan[\[Chi]z]^2-(a (1-zg^2) RealSign[-a EEg+Lzg] )/(z1g^2  Cos[\[Chi]z]^2Yzg)-z1s Yzg/(z1g Cos[\[Chi]z]^2)];
	fzaux=-z1s Sin[\[Chi]z](a^2 (1-EEg^2) zg)/Yzg;
	hz=Sign[\[Pi]-Mod[wr,2\[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg (2 a rg (Lzg-a EEg (1-z1g^2)))/(Sqrt[KKg] (rg^2+a^2 z1g^2)^2) Tan[\[Chi]z];
	
	1/Yzg (fzaux0+fzaux1+fzaux2+fzaux-hz)
]


\[Delta]YzoverYzgsparFCcoeff[nmax_,kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[-2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[\[Delta]YzoverYzgsparFC[wrlist[[i]],wzlist,a,p,e,xg,const],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fixed turning points (on average)*)


\[Delta]YzoverYzgsparDH[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]r,\[Chi]z,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],Rrg,Rr1g,Rr2g,Zzg,Lzred,z1s,dVzgdEg,dVzgdLzg,dVzgdKg,fzaux2,fzaux1,fzaux0,fzaux,fauxdC,hz},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg rg^2);
	Rr1g=(a^2 EEg-a Lzg+EEg r1g^2);
	Rr2g=(a^2 EEg-a Lzg+EEg r2g^2);
	Zzg=(Lzg-a EEg (1-zg^2));
	Lzred=Lzg-a EEg;
	
	z1s=z1sparDHfun[wr,a,p,e,xg,const,\[Delta]const];
	
	dVzgdEg=2 a (1-zg^2) Zzg;
	dVzgdLzg=-2 Zzg;
	dVzgdKg=1-zg^2;
	
	If[z1g==0,
		fzaux2=-z2g Rrg/(Lzred rg^2);
		fzaux0=RealSign[Lzg]((a^2(1-EEg^2)+a EEg Lzg)/Lzred^2-((a^2 (1-EEg^2)+Lzg^2) (Lzred^2+rg^2 ) )/(Lzred^2 rg^2 )) (a Tan[\[Chi]z]^2)/z2g+a RealSign[Lzg](( (a^2 (-1+EEg^2)-a EEg Lzg) rg^2+(Lzred^2+rg^2) z2g^2 )/(Lzred^2rg^2z2g Cos[\[Chi]z]^2))+RealSign[Lzg] (2 a (-2 a EEg+Lzg) \[Delta]const[[1]]+2 a EEg \[Delta]const[[2]]+\[Delta]const[[3]])/(2 z2g);
	
		1/Yzg (fzaux0+fzaux2)
		,
	
		fzaux2=(Yzg (-KKg rg^2 Rrg+a^2 (KKg+rg^2)(Rrg-2EEg rg^2) zg^2- (a^2 EEg-a Lzg+2 EEg KKg+3 EEg rg^2) a^4zg^4))/(Sqrt[KKg] (KKg-a^2 zg^2) \[CapitalSigma]^2);
		fzaux1=(2 a rg)/(Sqrt[KKg] \[CapitalSigma]^2 ) Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)] Yrg Zzg Tan[\[Chi]z];
		fzaux0=(a Sqrt[KKg]  (1-zg^2)Zzg )/((KKg-a^2 zg^2) z1g^2Cos[\[Chi]z]^2Yzg )+(a (KKg+rg^2)Zzg (a^2 (1-EEg^2) (2 zg^2-z1g^2)-z2g^2) )/(Sqrt[KKg] (KKg-a^2zg^2) \[CapitalSigma] Yzg) Tan[\[Chi]z]^2-(a (1-zg^2) RealSign[-a EEg+Lzg] )/(z1g^2  Cos[\[Chi]z]^2Yzg)-z1s Yzg/(z1g Cos[\[Chi]z]^2);
	
		fzaux=-z1s Sin[\[Chi]z] (a^2 (1-EEg^2) zg)/Yzg;
		fauxdC=(\[Delta]const[[1]]dVzgdEg+\[Delta]const[[2]]dVzgdLzg+\[Delta]const[[3]]dVzgdKg)/(2z1g^2Cos[\[Chi]z]^2Yzg);
	
		hz=Sign[\[Pi]-Mod[wr,2\[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg (2 a rg (Lzg-a EEg (1-z1g^2)))/(Sqrt[KKg] (rg^2+a^2 z1g^2)^2) Tan[\[Chi]z];
	
		1/Yzg (fzaux0+fzaux1+fzaux2+fzaux+fauxdC-hz)
	]
]


\[Delta]YzoverYzgsparDHcoeff[nmax_,kmax_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[\[Delta]YzoverYzgsparDH[wrlist[[i]],wzlist,a,p,e,xg,const,\[Delta]const],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsection::Closed:: *)
(*Functions for \[Xi]r - orthogonal component of the spin*)


(* ::Subsubsection::Closed:: *)
(*Function for \[Xi]r - Sin\[Psi]p*)


\[Delta]YroverYrgsortSin\[Psi][wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]z,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],Rrg,Zzg,r1sortCos\[Psi],r2sortCos\[Psi],fraux3,fraux2,fraux1,fraux0,fraux,d\[Psi]d\[Lambda] ,hrSin\[Psi]},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg rg^2);
	Zzg=(Lzg-a EEg (1-zg^2));
	
	r1sortCos\[Psi]=r1sortCos\[Psi]fun[wz,a,p,e,xg,const];
	r2sortCos\[Psi]=r2sortCos\[Psi]fun[wz,a,p,e,xg,const];
	
	fraux3=(a rg zg Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg) (rg-r2g)] Yrg^2  Sqrt[(KKg-a^2 zg^2)])/(KKg \[CapitalDelta] Sqrt[(KKg+rg^2)] \[CapitalSigma]^2) (rg^2 (-3+2 rg)+a^2  (rg-zg^2+rg zg^2));
	fraux2=(a rg Yrg  z1g Cos[\[Chi]z]Yzg)/(Sqrt[(KKg+rg^2) (KKg-a^2 zg^2)]\[CapitalSigma]);
	fraux1=1/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg) (rg-r2g)]) ((a zg rg Rrg (rg(a Lzg - rg EEg)-a (Lzg - a EEg))Sqrt[KKg-a^2 zg^2])/(KKg \[CapitalDelta] \[CapitalSigma] Sqrt[KKg+rg^2])-(a^3zg^3z1g^2Cos[\[Chi]z]^2Yzg^2\[CapitalDelta] Sqrt[KKg+rg^2] )/(KKg(1-zg^2) \[CapitalSigma]^2 Sqrt[KKg-a^2 zg^2])-(a^3zg^3Zzg(\[CapitalSigma] Lzg-2rg Zzg)Sqrt[KKg+rg^2])/(KKg(1-zg^2) \[CapitalSigma]^2 Sqrt[KKg-a^2 zg^2])-(2 rg EEg Rrg - rg \[CapitalDelta]-(rg-1)(KKg+rg^2)) (a rg zg Sqrt[KKg-a^2 zg^2])/(KKg \[CapitalSigma] Sqrt[KKg+rg^2]));
	
	d\[Psi]d\[Lambda]=Sqrt[KKg] (Rrg/(KKg+rg^2)+(a Zzg)/(KKg-a^2 zg^2));
	hrSin\[Psi]=-d\[Psi]d\[Lambda]((r1sortCos\[Psi]+r2sortCos\[Psi])/2+(r1sortCos\[Psi]-r2sortCos\[Psi])/2 Sin\[Chi]rfun[wr,a,p,e,xg,const]);
	
	1/Yrg (fraux1+fraux2+fraux3-hrSin\[Psi]/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(-rg+r1g) (rg-r2g)]))
]


(* ::Subsubsection::Closed:: *)
(*Function for \[Xi]r - Cos\[Psi]p*)


\[Delta]YroverYrgsortCos\[Psi][wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,z1g,z2g,\[Chi]z,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],\[CapitalSigma]r1g,\[CapitalSigma]r2g,Rrg,Zzg,Rr1g,Rr2g,r1sortCos\[Psi],r2sortCos\[Psi],fraux2,fraux1,fraux0,fraux,hrCos\[Psi]},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	\[CapitalSigma]r1g=r1g^2+a^2 zg^2;
	\[CapitalSigma]r2g=r2g^2+a^2 zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg rg^2);
	Rr1g=(a^2 EEg-a Lzg+EEg r1g^2);
	Rr2g=(a^2 EEg-a Lzg+EEg r2g^2);
	
	r1sortCos\[Psi]=r1sortCos\[Psi]fun[wz,a,p,e,xg,const];
	r2sortCos\[Psi]=r2sortCos\[Psi]fun[wz,a,p,e,xg,const];
	
	fraux2=-(1/(Sqrt[KKg]\[CapitalSigma]^2)) (a rg zg  Yrg)/(Sqrt[(KKg+rg^2) ] Sqrt[(KKg-a^2 zg^2)]) (a rg^2(Lzg-a EEg) -2KKg Rrg+a^2 zg^2Rrg+EEg a^2zg^2(a^2zg^2+2rg^2));
	fraux1=(a Rrg (KKg rg^2-a^2 (KKg+2 rg^2) zg^2) z1g Yzg Cos[\[Chi]z])/(Sqrt[KKg]\[CapitalSigma]^2 Sqrt[(KKg-a^2 zg^2)] Sqrt[(KKg+rg^2)] Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg) (rg-r2g)]);
	fraux0=-(r1sortCos\[Psi]-r2sortCos\[Psi]+(r1sortCos\[Psi]+r2sortCos\[Psi])Sin\[Chi]rfun[wr,a,p,e,xg,const]) ((r1g-r2g)Yrg)/(4(r1g-rg)(rg-r2g))+Sqrt[KKg-a^2 zg^2]/(Sqrt[KKg (KKg+rg^2)](rg-r1g) (rg-r2g)Yrg  \[CapitalSigma]) (a zg Rrg (KKg (1-rg)+rg (-a^2 (1-2 EEg^2)-2 a EEg Lzg+rg (3-2rg (1-EEg^2)))));
	fraux=((r1sortCos\[Psi]+r2sortCos\[Psi])/2+(r1sortCos\[Psi]-r2sortCos\[Psi])/2 Sin\[Chi]rfun[wr,a,p,e,xg,const]) ((1-EEg^2)(2 rg-r3g-r4g))/(2Yrg);
	
	hrCos\[Psi]=(a z1g Cos[\[Chi]z]Yzg)/(2Sqrt[KKg] Sqrt[(KKg-a^2 zg^2)]) ((Rr1g (KKg r1g^2-a^2 (KKg+2 r1g^2) zg^2))/(Sqrt[(KKg+r1g^2) ] \[CapitalSigma]r1g^2)+(Rr2g(KKg r2g^2-a^2 (KKg+2 r2g^2) zg^2))/(Sqrt[(KKg+r2g^2)] \[CapitalSigma]r2g^2)+((Rr1g (KKg r1g^2-a^2 (KKg+2 r1g^2) zg^2))/(Sqrt[(KKg+r1g^2) ] \[CapitalSigma]r1g^2)-(Rr2g(KKg r2g^2-a^2 (KKg+2 r2g^2) zg^2))/(Sqrt[(KKg+r2g^2)] \[CapitalSigma]r2g^2))Sin\[Chi]rfun[wr,a,p,e,xg,const]);
	
	1/Yrg (fraux0+fraux1+fraux2+fraux-hrCos\[Psi]/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(-rg+r1g) (rg-r2g)]))
]


(* ::Subsubsection::Closed:: *)
(*Functions for sampling*)


\[Delta]YroverYrgsortpluscoeff[nmax_,kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable,\[Psi]pbar},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	\[Psi]pbar[wr_,wz_]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	
	Chop[(ExpniTable . (Table[1/2 (\[Delta]YroverYrgsortCos\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const]+I*\[Delta]YroverYrgsortSin\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const])(Cos[\[Psi]pbar[wrlist[[i]],wzlist]]-I*Sin[\[Psi]pbar[wrlist[[i]],wzlist]]),{i,stepsr}]) . ExpjkTable) 1/(stepsr*stepsz),10^(-16)]
]


\[Delta]YroverYrgsortminuscoeff[nmax_,kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable,\[Psi]pbar},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	\[Psi]pbar[wr_,wz_]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	
	Chop[(ExpniTable . (Table[1/2 (\[Delta]YroverYrgsortCos\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const]-I*\[Delta]YroverYrgsortSin\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const])(Cos[\[Psi]pbar[wrlist[[i]],wzlist]]+I*Sin[\[Psi]pbar[wrlist[[i]],wzlist]]),{i,stepsr}]) . ExpjkTable) 1/(stepsr*stepsz),10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Final result*)


(* ::Text:: *)
(*Term proportional to e^(-i wp)   (j=+1)*)


\[Xi]rortplus[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	 \[CapitalUpsilon]rg Sum[Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg+\[CapitalUpsilon]p)),{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Text:: *)
(*Term proportional to e^(i wp)   (j=-1)*)


\[Xi]rortminus[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	 \[CapitalUpsilon]rg Sum[Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg-\[CapitalUpsilon]p)),{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Subsection::Closed:: *)
(*Functions for \[Xi]z - orthogonal component of the spin*)


(* ::Subsubsection::Closed:: *)
(*Function for \[Xi]z - Sin\[Psi]p*)


\[Delta]YzoverYzgsortSin\[Psi][wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,r3g,r4g,z2g,\[Chi]z,EEg,Lzg,KKg,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],Rrg,Zzg,Zz1g,z1sortCos\[Psi],fzaux3,fzaux2,fzaux1,hzSin\[Psi]},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2)(-rg+r3g)(-rg+r4g)];
	Yzg=Sqrt[-a^2(1-EEg^2)zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg*rg^2);
	Zzg=(Lzg-a*EEg(1-zg^2));
	Zz1g=(Lzg-a*EEg(1-z1g^2));
	
	z1sortCos\[Psi]=z1sortCos\[Psi]fun[wr,a,p,e,xg,const];
	
	If[z1g==0,
		0
		,
		fzaux3=-((a rg Sqrt[KKg+rg^2] zg^2 (rg^2-a^2 (1-2zg^2)) z1g Cos[\[Chi]z]Yzg^2)/(KKg (1-zg^2) Sqrt[KKg-a^2 zg^2]\[CapitalSigma]^2));
	
		fzaux2=(a Yrg zg  Yzg Sign[\[Pi]-Mod[wr,2 \[Pi]]] Sqrt[ (r1g-rg) (rg-r2g)])/(Sqrt[(KKg+rg^2)] Sqrt[KKg -a^2 zg^2] \[CapitalSigma]);
		fzaux1=-(1/(z1g*Cos[\[Chi]z]))((a Lzg rg zg^2 Zzg Sqrt[KKg+rg^2])/(KKg(1-zg^2)\[CapitalSigma] Sqrt[KKg-a^2zg^2])+(r1g-rg) (rg-r2g)(1-zg^2) (a rg^3Yrg^2 Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta] \[CapitalSigma]^2 Sqrt[KKg+rg^2])+(rg^3Rrg(2rg Zzg-Lzg \[CapitalSigma])Sqrt[KKg-a^2zg^2])/(KKg \[CapitalDelta] \[CapitalSigma]^2 Sqrt[KKg+rg^2])+((1-EEg^2)a^2zg(2zg^2-z1g^2)-zg z2g^2) (a zg rg Sqrt[KKg+rg^2])/(KKg \[CapitalSigma] Sqrt[KKg-a^2zg^2]));
	
		hzSin\[Psi]=((a(Lzg-a EEg)+EEg KKg) rg \[CapitalSigma] Zz1g)/((KKg-a^2 zg^2) Sqrt[(KKg+rg^2) (KKg-a^2 z1g^2)] (rg^2+a^2 z1g^2));
	
		1/Yzg(fzaux1+fzaux2+fzaux3-hzSin\[Psi]/(z1g*Cos[\[Chi]z]))
	]
]


(* ::Subsubsection::Closed:: *)
(*Function for \[Xi]z - Cos\[Psi]p*)


\[Delta]YzoverYzgsortCos\[Psi][wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{r1g,r2g,z1g,r3g,r4g,z2g,\[Chi]z,EEg,Lzg,KKg,rg,zg,Yrg,Yzg,\[CapitalDelta],\[CapitalSigma],Rrg,Zzg,Zz1g,z1sortCos\[Psi],fzaux2,fzaux1,fzaux0,fzaux,hzCos\[Psi]},
	{EEg,Lzg,KKg}=const;
	r1g=p/(1-e);
	r2g=p/(1+e);
	z1g=Sqrt[1-xg^2];
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Yrg=Sqrt[(1-EEg^2)(-rg+r3g)(-rg+r4g)];
	Yzg=Sqrt[-a^2(1-EEg^2)zg^2+z2g^2];
	\[CapitalDelta]=(a^2-2rg+rg^2);
	\[CapitalSigma]=rg^2+a^2 zg^2;
	Rrg=(a^2*EEg-a*Lzg+EEg*rg^2);
	Zzg=(Lzg-a*EEg(1-zg^2));
	Zz1g=(Lzg-a*EEg(1-z1g^2));
	z1sortCos\[Psi]=z1sortCos\[Psi]fun[wr,a,p,e,xg,const];
	
	If[z1g==0,
		0
		,
		fzaux2=-(a rg zg  (EEg rg^2 \[CapitalSigma]+a rg^2 Zzg+Rrg a^2zg^2+2a KKg Zzg) Yzg )/(Sqrt[KKg (KKg+rg^2)] Sqrt[KKg-a^2zg^2] \[CapitalSigma]^2);
	
		fzaux1=( Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg) (rg-r2g)] Yrg (KKg rg^2-a^2zg^2 (KKg+2rg^2)) Zzg)/(Sqrt[KKg (KKg+rg^2)] Sqrt[KKg-a^2zg^2] \[CapitalSigma]^2 z1g Cos[\[Chi]z]);
		fzaux0=-((rg zg Zzg Sqrt[KKg+rg^2] (-a^2 (1-EEg^2) (2 zg^2-z1g^2)+z2g^2))/(Sqrt[KKg] Sqrt[KKg-a^2 zg^2] \[CapitalSigma] z1g^2 Cos[\[Chi]z]^2Yzg))-z1sortCos\[Psi] Sin[\[Chi]z] Yzg/(z1g Cos[\[Chi]z]^2);
		fzaux=-z1sortCos\[Psi] (a^2 (1-EEg^2) zg)/Yzg;
	
		hzCos\[Psi]=(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg) (rg-r2g)] Yrg Zz1g)(-2 rg^2a^2 z1g^2+KKg (rg^2-a^2 z1g^2))/(Sqrt[KKg (KKg+rg^2) (KKg-a^2 z1g^2)](rg^2+a^2z1g^2)^2);
	
		1/Yzg(fzaux0+fzaux1+fzaux2+fzaux-hzCos\[Psi]/(z1g*Cos[\[Chi]z]))
	]
]


(* ::Subsubsection::Closed:: *)
(*Functions for sampling*)


\[Delta]YzoverYzgsortpluscoeff[nmax_,kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable,\[Psi]pbar},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	\[Psi]pbar[wr_,wz_]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	
	Chop[(ExpniTable . (Table[1/2 (\[Delta]YzoverYzgsortCos\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const]+I*\[Delta]YzoverYzgsortSin\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const])(Cos[\[Psi]pbar[wrlist[[i]],wzlist]]-I*Sin[\[Psi]pbar[wrlist[[i]],wzlist]]),{i,stepsr}]) . ExpjkTable) 1/(stepsr*stepsz),10^(-16)]
]


\[Delta]YzoverYzgsortminuscoeff[nmax_,kmax_,a_,p_,e_,xg_,const_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable,\[Psi]pbar},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	\[Psi]pbar[wr_,wz_]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	
	Chop[(ExpniTable . (Table[1/2 (\[Delta]YzoverYzgsortCos\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const]-I*\[Delta]YzoverYzgsortSin\[Psi][wrlist[[i]],wzlist,a,p,e,xg,const])(Cos[\[Psi]pbar[wrlist[[i]],wzlist]]+I*Sin[\[Psi]pbar[wrlist[[i]],wzlist]]),{i,stepsr}]) . ExpjkTable) 1/(stepsr*stepsz),10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Final result*)


(* ::Text:: *)
(*Term proportional to e^(-i wp)   (j=+1)*)


\[Xi]zortplus[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	\[CapitalUpsilon]zg Sum[Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg+\[CapitalUpsilon]p)),{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Text:: *)
(*Term proportional to e^(i wp)   (j=-1)*)


\[Xi]zortminus[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	\[CapitalUpsilon]zg Sum[Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg-\[CapitalUpsilon]p)),{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Section::Closed:: *)
(*Spin corrections to the orbit  - parallel component of the spin*)


(* ::Text:: *)
(*General function for Fourier series expansion of the spin corrections to the velocities - parallel component in the spin *)


dvspard\[Lambda]Fourier[wr_,wz_,coeff_]:=Module[{dimr,dimz},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	
	coeff[[dimr+1,dimz+1]]+2 Sum[Cos[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]],{n,1,dimr},{k,1,dimz}]+2 Sum[Cos[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]],{n,1,dimr},{k,-dimz,-1}]+2Sum[Cos[n wr]coeff[[n+dimr+1,dimz+1]],{n,1,dimr}]+2Sum[Cos[k wz]coeff[[dimr+1,k+dimz+1]],{k,1,dimz}]
];


(* ::Text:: *)
(*General function for Fourier series expansion of the spin corrections to the trajectories - parallel component in the spin (identical functional structure as above, but different name for convenience)*)


\[Delta]ysparFourier[wr_,wz_,coeff_]:=Module[{dimr,dimz},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	
	coeff[[dimr+1,dimz+1]]+2 Sum[Cos[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]],{n,1,dimr},{k,1,dimz}]+2 Sum[Cos[n wr+k wz]coeff[[n+dimr+1,k+dimz+1]],{n,1,dimr},{k,-dimz,-1}]+2Sum[Cos[n wr]coeff[[n+dimr+1,dimz+1]],{n,1,dimr}]+2Sum[Cos[k wz]coeff[[dimr+1,k+dimz+1]],{k,1,dimz}]
];


(* ::Subsection::Closed:: *)
(*Radial and polar trajectories*)


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


\[Delta]rparFCfun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r1spar,r2spar,rspar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r,\[Chi]z},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	r1spar=r1sparFCfun[wz,a,p,e,xg,const];
	r2spar=r2sparFCfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	rspar=1/2 (r1spar+r2spar)+1/2 (r1spar-r2spar)Sin\[Chi]rfun[wr,a,p,e,xg,const];
	rspar-drgdwrfun[wr,a,p,e,xg,const]Re[\[Xi]rpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]
]


\[Delta]zparDHfun[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,z1spar,zspar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r,\[Chi]z},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	z1spar=z1sparDHfun[wr,a,p,e,xg,const,\[Delta]const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	zspar=z1spar Sin[\[Chi]z];
	
	If[z1g==0,
		0,
		zspar- dzgdwzfun[wz,a,p,e,xg,const]Re[\[Xi]zpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]
	]
]


(* ::Subsubsection::Closed:: *)
(*Fixed turning points (on average)*)


\[Delta]rparDHfun[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_,coeff_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r1spar,r2spar,rspar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r,\[Chi]z},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	
	r1spar=r1sparDHfun[wz,a,p,e,xg,const,\[Delta]const];
	r2spar=r2sparDHfun[wz,a,p,e,xg,const,\[Delta]const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	rspar=1/2 (r1spar+r2spar)+1/2 (r1spar-r2spar)Sin\[Chi]rfun[wr,a,p,e,xg,const];
	rspar-drgdwrfun[wr,a,p,e,xg,const]Re[\[Xi]rpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]
]


\[Delta]zparFCfun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,z1spar,zspar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[Chi]r,\[Chi]z},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	z1spar=z1sparFCfun[wr,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	zspar=z1spar Sin[\[Chi]z];
	
	If[z1g==0,
		0,
		zspar-dzgdwzfun[wz,a,p,e,xg,const]Re[\[Xi]zpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]
	]
]


(* ::Subsection::Closed:: *)
(*Radial and polar velocities *)


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


drsparFCd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,\[Chi]r,\[Chi]z,r1spar,r2spar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,rg,rspar,\[CapitalDelta]rg,Rrg,Yrg,dVrgdrg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	r1spar=r1sparFCfun[wz,a,p,e,xg,const];
	r2spar=r2sparFCfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	rg=rgfun[wr,a,p,e,xg,const];
	rspar=1/2 (r1spar+r2spar)+1/2 (r1spar-r2spar)Sin\[Chi]rfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	dVrgdrg=2(1-rg)(KKg+rg^2)-2rg*\[CapitalDelta]rg+4EEg*rg*Rrg;
	
	1/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg) (\[CapitalDelta]rg^2 \[ScriptCapitalS]12rfun[wr,wz,a,p,e,xg,const]+a*\[CapitalDelta]rg*RealSign[Lzg-a EEg]+\[CapitalDelta]rg Sqrt[KKg] Rrg/(KKg+rg^2)+1/2 rspar dVrgdrg)-(drgdwrfun[wr,a,p,e,xg,const]dVrgdrg)/(2Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg) Re[\[Xi]rpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]
]


dzsparFCd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,z2g,\[Chi]r,\[Chi]z,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,zg,z1spar,zspar,Zzg,Yzg,dVzgdzg},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	zg=zgfun[wz,a,p,e,xg,const];
	z1spar=z1sparFCfun[wr,a,p,e,xg,const];
	zspar=z1spar Sin[\[Chi]z];
	
	Zzg=(Lzg-a EEg (1-zg^2));
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	dVzgdzg=-2 zg (KKg+a^2 (1-2 zg^2))-2 Zzg 2 a EEg zg;
	
	If[z1g==0,
		0
		,
		1/(z1g Cos[\[Chi]z]Yzg) ((1-zg^2)\[ScriptCapitalS]12zfun[wr,wz,a,p,e,xg,const]-a(1-zg^2)RealSign[Lzg-a EEg]+(1-zg^2)Sqrt[KKg]a Zzg/(KKg-a^2zg^2)+1/2 dVzgdzg zspar)-dzgdwzfun[wz,a,p,e,xg,const]/(2z1g Cos[\[Chi]z]Yzg) Re[\[Xi]zpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]dVzgdzg
	]
]


(* ::Subsubsection::Closed:: *)
(*Fixed turning points (on average)*)


drsparDHd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_,coeff_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,\[Chi]r,\[Chi]z,r1spar,r2spar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,rg,rspar,\[CapitalDelta]rg,Rrg,Yrg,dVrgdrg,dVrgdEg,dVrgdLzg,dVrgdKg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	r1spar=r1sparDHfun[wz,a,p,e,xg,const,\[Delta]const];
	r2spar=r2sparDHfun[wz,a,p,e,xg,const,\[Delta]const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	rg=rgfun[wr,a,p,e,xg,const];
	rspar=1/2 (r1spar+r2spar)+1/2 (r1spar-r2spar)Sin\[Chi]rfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	dVrgdrg=2(1-rg) (KKg+rg^2)-2 rg \[CapitalDelta]rg+4 EEg rg Rrg;
	dVrgdEg=2(rg^2+a^2)Rrg;
	dVrgdLzg=-2a Rrg;
	dVrgdKg=-\[CapitalDelta]rg;
	
	1/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg) (\[CapitalDelta]rg^2 \[ScriptCapitalS]12rfun[wr,wz,a,p,e,xg,const]+a \[CapitalDelta]rg RealSign[Lzg-a EEg]+\[CapitalDelta]rg Sqrt[KKg] Rrg/(KKg+rg^2)+1/2 rspar dVrgdrg+1/2 (\[Delta]const[[1]]dVrgdEg +\[Delta]const[[2]]dVrgdLzg+\[Delta]const[[3]]dVrgdKg))-(drgdwrfun[wr,a,p,e,xg,const]dVrgdrg)/(2Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg) Re[\[Xi]rpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]
]


dzsparDHd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,z2g,\[Chi]r,\[Chi]z,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,zg,z1spar,zspar,Zzg,Yzg,dVzgdzg,dVzgdEg,dVzgdLzg,dVzgdKg},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	
	zg=zgfun[wz,a,p,e,xg,const];
	z1spar=z1sparDHfun[wr,a,p,e,xg,const,\[Delta]const];
	zspar=z1spar Sin[\[Chi]z];
	
	Zzg=(Lzg-a EEg (1-zg^2));
	Yzg=Sqrt[-a^2 (1-EEg^2) zg^2+z2g^2];
	dVzgdzg=-2 zg (KKg+a^2 (1-2 zg^2))-2 Zzg 2 a EEg zg;
	dVzgdEg=2 a (1-zg^2) Zzg;
	dVzgdLzg=-2 Zzg;
	dVzgdKg=1-zg^2;
	
	If[z1g==0,
		0
		,
		1/(z1g Cos[\[Chi]z]Yzg) ((1-zg^2) \[ScriptCapitalS]12zfun[wr,wz,a,p,e,xg,const]-a (1-zg^2) RealSign[Lzg-a EEg]+(1-zg^2)Sqrt[KKg]a Zzg/(KKg-a^2zg^2)+1/2 dVzgdzg zspar+1/2 (\[Delta]const[[1]]dVzgdEg +\[Delta]const[[2]]dVzgdLzg+\[Delta]const[[3]]dVzgdKg))- dzgdwzfun[wz,a,p,e,xg,const]/(2z1g Cos[\[Chi]z]Yzg) Re[\[Xi]zpar[wr,wz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,coeff]]dVzgdzg
	]
]


(* ::Subsection::Closed:: *)
(*Coordinate time velocity*)


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


dtsparFCd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,\[Chi]r,\[Chi],rg,zg,\[CapitalDelta]rg,\[CapitalSigma],dVtrgdrg,dVtzgdzg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2zg^2;
	
	dVtrgdrg=(2 a Lzg (rg^2-a^2))/\[CapitalDelta]rg^2+(4 EEg rg (a^2+rg^2))/\[CapitalDelta]rg-(2EEg (rg-1)  (a^2+rg^2)^2)/\[CapitalDelta]rg^2;
	dVtzgdzg=2 a^2 EEg zg;
	
	dVtrgdrg \[Delta]rparFCfun[wr,wz,a,p,e,xg,const,coeffrad]+dVtzgdzg \[Delta]zparFCfun[wr,wz,a,p,e,xg,const,coeffpol] -\[CapitalSigma]/2 \[CapitalGamma]t12fun[wr,wz,a,p,e,xg,const]
]


dtsparFCd\[Lambda]coefffun[nmax_,kmax_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[dtsparFCd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,xg,const,coeffrad,coeffpol],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fixed turning points (on average)*)


dtsparDHd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,\[Chi]r,\[Chi],rg,zg,\[CapitalDelta]rg,\[CapitalSigma],Rrg,Zzg,dVtrgdrg,dVtrgdEg,dVtrgdLzg,dVtzgdzg,dVtzgdEg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2zg^2;
	Rrg=(a^2 EEg-a Lzg+EEg rg^2);
	Zzg=(Lzg-a EEg (1-zg^2));
	
	dVtrgdrg=(2 a Lzg (rg^2-a^2))/\[CapitalDelta]rg^2+(4 EEg rg (a^2+rg^2))/\[CapitalDelta]rg-(2EEg (rg-1)  (a^2+rg^2)^2)/\[CapitalDelta]rg^2;
	dVtrgdEg=(a^2+rg^2)^2/\[CapitalDelta]rg;
	dVtrgdLzg=-((2 a rg)/\[CapitalDelta]rg);
	dVtzgdzg=2 a^2 EEg zg;
	dVtzgdEg=-a^2 (1-zg^2);
	
	dVtrgdrg \[Delta]rparDHfun[wr,wz,a,p,e,xg,const,\[Delta]const,coeffrad]+dVtzgdzg \[Delta]zparDHfun[wr,wz,a,p,e,xg,const,\[Delta]const,coeffpol]+(\[Delta]const[[1]]dVtrgdEg+\[Delta]const[[2]]dVtrgdLzg)+(\[Delta]const[[1]]dVtzgdEg)-\[CapitalSigma]/2 \[CapitalGamma]t12fun[wr,wz,a,p,e,xg,const]
]


dtsparDHd\[Lambda]coefffun[nmax_,kmax_,a_,p_,e_,xg_,const_,\[Delta]const_,coeffrad_,coeffpol_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[dtsparDHd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,xg,const,\[Delta]const,coeffrad,coeffpol],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
	]


(* ::Subsection::Closed:: *)
(*Azimuthal velocity*)


(* ::Subsubsection::Closed:: *)
(*Fixed constants of motion*)


d\[Phi]sparFCd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,\[Chi]r,\[Chi]z,rg,zg,\[CapitalDelta]rg,\[CapitalSigma],dV\[Phi]rgdrg,dV\[Phi]zgdzg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2zg^2;
	
	dV\[Phi]rgdrg=(2 a EEg rg)/\[CapitalDelta]rg-(2a (rg-1) (-a Lzg+EEg (a^2+rg^2)))/\[CapitalDelta]rg^2;
	dV\[Phi]zgdzg=(2 Lzg zg)/(1-zg^2)^2;
	
	dV\[Phi]rgdrg \[Delta]rparFCfun[wr,wz,a,p,e,xg,const,coeffrad]+dV\[Phi]zgdzg \[Delta]zparFCfun[wr,wz,a,p,e,xg,const,coeffpol]-\[CapitalSigma]/2 \[CapitalGamma]\[Phi]12fun[wr,wz,a,p,e,xg,const]
]


d\[Phi]sparFCd\[Lambda]coefffun[nmax_,kmax_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[d\[Phi]sparFCd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,xg,const,coeffrad,coeffpol],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Fixed turning points (on average)*)


d\[Phi]sparDHd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,\[Delta]const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,rg,zg,\[CapitalDelta]rg,\[CapitalSigma],dV\[Phi]rgdrg,dV\[Phi]rgdEg,dV\[Phi]rgdLzg,dV\[Phi]zgdzg,dV\[Phi]zgdLzg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2zg^2;
	
	dV\[Phi]rgdrg=(2 a EEg rg)/\[CapitalDelta]rg-(2a (rg-1) (-a Lzg+EEg (a^2+rg^2)))/\[CapitalDelta]rg^2;
	dV\[Phi]rgdEg=(2 a rg)/\[CapitalDelta]rg;
	dV\[Phi]rgdLzg=-(a^2/\[CapitalDelta]rg);
	dV\[Phi]zgdzg=(2 Lzg zg)/(1-zg^2)^2;
	dV\[Phi]zgdLzg=1/(1-zg^2);
	
	dV\[Phi]rgdrg \[Delta]rparDHfun[wr,wz,a,p,e,xg,const,\[Delta]const,coeffrad]+dV\[Phi]zgdzg \[Delta]zparDHfun[wr,wz,a,p,e,xg,const,\[Delta]const,coeffpol]+(\[Delta]const[[1]]dV\[Phi]rgdEg+\[Delta]const[[2]]dV\[Phi]rgdLzg)+(\[Delta]const[[2]]dV\[Phi]zgdLzg)-\[CapitalSigma]/2 \[CapitalGamma]\[Phi]12fun[wr,wz,a,p,e,xg,const]
]


d\[Phi]sparDHd\[Lambda]coefffun[nmax_,kmax_,a_,p_,e_,xg_,const_,\[Delta]const_,coeffrad_,coeffpol_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[d\[Phi]sparDHd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,xg,const,\[Delta]const,coeffrad,coeffpol],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsection::Closed:: *)
(*Coordinate time and azimuthal trajectories*)


(* ::Text:: *)
(*The functional structure for Subscript[\[CapitalDelta]t, s\[DoubleVerticalBar]] and Subscript[\[CapitalDelta]\[Phi], s\[DoubleVerticalBar]] is the same for both parametrizations. So I create just one Mathematica function that I will specialize later on.*)


\[CapitalDelta]spar[wr_,wz_,{\[CapitalUpsilon]rg_,\[CapitalUpsilon]zg_},{\[CapitalUpsilon]rs_,\[CapitalUpsilon]zs_},\[Delta]coeff_,coeffVtr_,coeffVtz_]:=Module[{dimr,dimz},
	dimr=(Dimensions[\[Delta]coeff][[1]]-1)/2;
	dimz=(Dimensions[\[Delta]coeff][[2]]-1)/2;
	
	2Sum[Sin[n wr+k wz]\[Delta]coeff[[n+dimr+1,k+dimz+1]]/(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg),{n,1,dimr},{k,1,dimz}]+2Sum[Sin[n wr+k wz]\[Delta]coeff[[n+dimr+1,k+dimz+1]]/(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg),{n,1,dimr},{k,-dimz,-1}]+2Sum[Sin[n wr]\[Delta]coeff[[n+dimr+1,dimz+1]]/(n*\[CapitalUpsilon]rg),{n,1,dimr}]+2Sum[Sin[k wz]\[Delta]coeff[[dimr+1,k+dimz+1]]/(k*\[CapitalUpsilon]zg),{k,1,dimz}]-2(\[CapitalUpsilon]rs/\[CapitalUpsilon]rg^2)Sum[Sin[n wr]coeffVtr[[n+dimr+1]]/n,{n,1,dimr}]-2(\[CapitalUpsilon]zs/\[CapitalUpsilon]zg^2)Sum[Sin[k wz]coeffVtz[[k+dimz+1]]/k,{k,1,dimz}]
]


(* ::Section::Closed:: *)
(*Spin corrections to the orbit  - orthogonal component of the spin*)


(* ::Subsection::Closed:: *)
(*Radial velocity*)


(* ::Subsubsection::Closed:: *)
(*Radial velocity orthogonal component of the spin - Sin\[Psi]p*)


drsortSin\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg,\[CapitalDelta]rg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	
	\[CapitalDelta]rg^2\[ScriptCapitalS]31rfunres[wr,wz,a,p,e,xg,const]
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity orthogonal component of the spin - Cos\[Psi]p*)


drsortCos\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,r1sort,r2sort,rg,rsort,\[CapitalDelta]rg,Rrg,Yrg,dVrgdrg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	r1sort=r1sortCos\[Psi]fun[wz,a,p,e,xg,const];
	r2sort=r2sortCos\[Psi]fun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	rsort=1/2 (r1sort+r2sort)+1/2 (r1sort-r2sort)Sin\[Chi]rfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	dVrgdrg=2(1-rg) (KKg+rg^2)-2 rg \[CapitalDelta]rg+4 EEg rg Rrg;
	
	1/(Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg) (\[CapitalDelta]rg^2 \[ScriptCapitalS]23rfun[wr,wz,a,p,e,xg,const]+1/2 rsort dVrgdrg )
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity orthogonal component of the spin - term proportional to  e^(-i wp)   (j=+1)*)


drsortplusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,rg,\[CapitalDelta]rg,Rrg,dVrgdrg,Yrg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	dVrgdrg=2(1-rg) (KKg+rg^2)-2 rg \[CapitalDelta]rg+4 EEg rg Rrg;
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[Abs[xg]==1,
		0
		,
		1/2(drsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]+I*drsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])-(dVrgdrg*drgdwrfun[wr,a,p,e,xg,const])/(2Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg)\[Xi]rortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]
	]
]


(* ::Subsubsection::Closed:: *)
(*Radial velocity orthogonal component of the spin - term proportional to e^(i wp)   (j=-1)*)


drsortminusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,rg,\[CapitalDelta]rg,Rrg,dVrgdrg,Yrg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g (a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2));
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	Rrg=(-a Lzg+EEg (a^2+rg^2));
	dVrgdrg=2(1-rg) (KKg+rg^2)-2 rg \[CapitalDelta]rg+4 EEg rg Rrg;
	Yrg=Sqrt[(1-EEg^2) (-rg+r3g) (-rg+r4g)];
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[Abs[xg]==1,
		0,
		1/2(drsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]-I*drsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])-(dVrgdrg*drgdwrfun[wr,a,p,e,xg,const])/(2Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg)\[Xi]rortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]
	]
]


(* ::Subsection::Closed:: *)
(*Polar velocity*)


(* ::Subsubsection::Closed:: *)
(*Polar velocity orthogonal component of the spin -  Sin\[Psi]p*)


dzsortSin\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,zg},
	{EEg,Lzg,KKg}=const;
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	(1-zg^2)\[ScriptCapitalS]31zfunres[wr,wz,a,p,e,xg,const]
]


(* ::Subsubsection::Closed:: *)
(*Polar velocity orthogonal component of the spin - Cos\[Psi]p*)


dzsortCos\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,z1g,z2g,\[Chi]z,rg,zg,zsort,\[CapitalDelta],\[CapitalSigma],Zzg,Rrg,Yzg,dVzgdzg},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	zsort=z1sortCos\[Psi]fun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]=(a^2-2rg+rg^2);
	\[CapitalSigma]=rg^2+a^2*zg^2;
	Rrg=(a^2*EEg-a*Lzg+EEg*rg^2);
	Zzg=(Lzg-a*EEg(1-zg^2));
	Yzg=Sqrt[-a^2(1-EEg^2)zg^2+z2g^2];
	dVzgdzg=-2zg(KKg+a^2(1-2 zg^2))-2Zzg*2a*EEg*zg;
	
	If[z1g==0,
		((-a EEg+Lzg)Sign[Pi-Mod[wr,2Pi]]Sqrt[Rrg^2-\[CapitalDelta](KKg+rg^2)])/(rg^2*Sqrt[(-a EEg+Lzg)^2+rg^2])
		,
		1/(z1g Cos[\[Chi]z]Yzg)((1-zg^2)\[ScriptCapitalS]23zfun[wr,wz,a,p,e,xg,const]+1/2dVzgdzg*zsort)
	]
]


(* ::Subsubsection::Closed:: *)
(*Polar velocity orthogonal component of the spin - term proportional to e^(-i wp)   (j=+1)*)


dzsortplusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,z2g,\[Chi]z,zg,Zzg,dVzgdzg,Yzg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	\[Chi]z=zgfun[wz,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Zzg=(Lzg-a*EEg(1-zg^2));
	dVzgdzg=-2zg(KKg+a^2(1-2 zg^2))-2Zzg*2a*EEg*zg;
	Yzg=Sqrt[-a^2(1-EEg^2)zg^2+z2g^2];
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		1/2*(dzsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]+I*dzsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])
		,
		1/2*(dzsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]+I*dzsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])-dzgdwzfun[wz,a,p,e,xg,const]/(2z1g*Cos[\[Chi]z]Yzg) \[Xi]zortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]*dVzgdzg
	]
]


(* ::Subsubsection::Closed:: *)
(*Polar velocity orthogonal component of the spin - term proportional to e^(i wp)   (j=-1)*)


dzsortminusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,z2g,\[Chi]z,zg,Zzg,dVzgdzg,Yzg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=zgfun[wz,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	Zzg=(Lzg-a*EEg(1-zg^2));
	dVzgdzg=-2zg(KKg+a^2(1-2 zg^2))-2Zzg*2a*EEg*zg;
	Yzg=Sqrt[-a^2(1-EEg^2)zg^2+z2g^2];
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		1/2 (dzsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]-I*dzsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])
		,
		1/2 (dzsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]-I*dzsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])-dzgdwzfun[wz,a,p,e,xg,const]/(2z1g*Cos[\[Chi]z]Yzg) \[Xi]zortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]*dVzgdzg
	]
]


(* ::Subsection::Closed:: *)
(*Coordinate time velocity *)


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity orthogonal component of the spin -  Sin\[Psi]p*)


dtsortSin\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg,zg,\[CapitalSigma]},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalSigma]=rg^2+a^2zg^2;
	
	 -\[CapitalSigma]/2 \[CapitalGamma]t31fun[wr,wz,a,p,e,xg,const]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity parallel component of the spin - Cos\[Psi]p*)


dtsortCos\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg,zg,\[CapitalDelta]rg,\[CapitalSigma],r1sortCos\[Psi],r2sortCos\[Psi],rsort,zsort,dVtrgdrg,dVtzgdzg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	\[CapitalSigma]=rg^2+a^2zg^2;
	
	r1sortCos\[Psi]=r1sortCos\[Psi]fun[wz,a,p,e,xg,const];
	r2sortCos\[Psi]=r2sortCos\[Psi]fun[wz,a,p,e,xg,const];
	rsort=1/2 (r1sortCos\[Psi]+r2sortCos\[Psi])+1/2 (r1sortCos\[Psi]-r2sortCos\[Psi])Sin\[Chi]rfun[wr,a,p,e,xg,const];
	zsort=z1sortCos\[Psi]fun[wr,a,p,e,xg,const];
	
	dVtrgdrg=(2 a Lzg (rg^2-a^2))/\[CapitalDelta]rg^2+(4 EEg rg (a^2+rg^2))/\[CapitalDelta]rg-(2EEg(rg-1)(a^2+rg^2)^2)/\[CapitalDelta]rg^2;
	dVtzgdzg=2a^2*EEg*zg;
	
	dVtrgdrg*rsort+dVtzgdzg*zsort -\[CapitalSigma]/2*\[CapitalGamma]t23fun[wr,wz,a,p,e,xg,const]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity orthogonal component of the spin - term proportional to e^(-i wp)   (j=+1)*)


dtsortplusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,z1g,\[Chi]z,rg,zg,\[CapitalDelta]rg,dVtrgdrg,dVtzgdzg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	dVtrgdrg=(2a*Lzg(rg^2-a^2))/\[CapitalDelta]rg^2+(4 EEg*rg(a^2+rg^2))/\[CapitalDelta]rg-(2EEg(rg-1)(a^2+rg^2)^2)/\[CapitalDelta]rg^2;
	dVtzgdzg=2a^2*EEg*zg;
	
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		0
		,
		1/2(dtsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]+I*dtsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])-drgdwrfun[wr,a,p,e,xg,const]*\[Xi]rortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffrad]*dVtrgdrg-dzgdwzfun[wz,a,p,e,xg,const]*\[Xi]zortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffpol]*dVtzgdzg
	]
]


dtsortplusd\[Lambda]coefffun[nmax_,kmax_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[dtsortplusd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,xg,const,coeffrad,coeffpol],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)]
]


(* ::Subsubsection::Closed:: *)
(*Coordinate time velocity orthogonal component of the spin - term proportional to e^(i wp)   (j=-1)*)


dtsortminusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,z1g,\[Chi]z,rg,zg,\[CapitalDelta]rg,dVtrgdrg,dVtzgdzg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	dVtrgdrg=(2a*Lzg*(rg^2-a^2))/\[CapitalDelta]rg^2+(4EEg*rg*(a^2+rg^2))/\[CapitalDelta]rg-(2EEg(rg-1)(a^2+rg^2)^2)/\[CapitalDelta]rg^2;
	dVtzgdzg=2a^2*EEg*zg;
	
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		0
		,
		1/2(dtsortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]-I*dtsortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])-drgdwrfun[wr,a,p,e,xg,const]*\[Xi]rortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffrad]dVtrgdrg-dzgdwzfun[wz,a,p,e,xg,const]*\[Xi]zortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffpol]*dVtzgdzg
	]
]


dtsortminusd\[Lambda]coefffun[nmax_,kmax_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{stepsr,stepsz,wrlist,wzlist,ExpniTable,ExpjkTable},
	stepsr=4*nmax;
	stepsz=4*kmax;
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,xg}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,xg}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	Chop[(ExpniTable . (Table[dtsortminusd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,xg,const,coeffrad,coeffpol],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^-(16)]
]


(* ::Subsection::Closed:: *)
(*Azimuthal velocity *)


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity orthogonal component of the spin -  Sin\[Psi]p*)


d\[Phi]sortSin\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg,zg,\[CapitalSigma]},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalSigma]=rg^2+a^2zg^2;
	
	-(\[CapitalSigma]/2)\[CapitalGamma]\[Phi]31fun[wr,wz,a,p,e,xg,const]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity parallel component of the spin - Cos\[Psi]p*)


d\[Phi]sortCos\[Psi]d\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_]:=Module[{EEg,Lzg,KKg,rg,zg,\[CapitalDelta]rg,\[CapitalSigma],r1sortCos\[Psi],r2sortCos\[Psi],rsort,zsort,dV\[Phi]rgdrg,dV\[Phi]zgdzg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	\[CapitalSigma]=rg^2+a^2zg^2;
	
	r1sortCos\[Psi]=r1sortCos\[Psi]fun[wz,a,p,e,xg,const];
	r2sortCos\[Psi]=r2sortCos\[Psi]fun[wz,a,p,e,xg,const];
	rsort=1/2(r1sortCos\[Psi]+r2sortCos\[Psi])+1/2(r1sortCos\[Psi]-r2sortCos\[Psi])Sin\[Chi]rfun[wr,a,p,e,xg,const];
	zsort=z1sortCos\[Psi]fun[wr,a,p,e,xg,const];
	
	dV\[Phi]rgdrg=(2a*EEg*rg)/\[CapitalDelta]rg-(2a(rg-1)(-a Lzg+EEg(a^2+rg^2)))/\[CapitalDelta]rg^2;
	dV\[Phi]zgdzg=(2Lzg*zg)/(1-zg^2)^2;
	
	dV\[Phi]rgdrg*rsort+dV\[Phi]zgdzg*zsort-\[CapitalSigma]/2*\[CapitalGamma]\[Phi]23fun[wr,wz,a,p,e,xg,const]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity orthogonal component of the spin - term proportional to e^(-i wp)   (j=+1)*)


d\[Phi]sortplusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,z1g,\[Chi]z,rg,zg,\[CapitalDelta]rg,dV\[Phi]rgdrg,dV\[Phi]zgdzg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	dV\[Phi]rgdrg=(2a*EEg*rg)/\[CapitalDelta]rg-(2a(rg-1)(-a*Lzg+EEg(a^2+rg^2)))/\[CapitalDelta]rg^2;
	dV\[Phi]zgdzg=(2Lzg*zg)/(1-zg^2)^2;
	
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		0
		,
		1/2(d\[Phi]sortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]+I*d\[Phi]sortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])-drgdwrfun[wr,a,p,e,xg,const]*\[Xi]rortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffrad]*dV\[Phi]rgdrg-dzgdwzfun[wz,a,p,e,xg,const]*\[Xi]zortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffpol]*dV\[Phi]zgdzg
	]
]


(* ::Subsubsection::Closed:: *)
(*Azimuthal velocity orthogonal component of the spin - term proportional to e^(i wp)   (j=-1)*)


d\[Phi]sortminusd\[Lambda]fun[wr_,wz_,a_,p_,e_,xg_,const_,coeffrad_,coeffpol_]:=Module[{EEg,Lzg,KKg,z1g,\[Chi]z,rg,zg,\[CapitalDelta]rg,dV\[Phi]rgdrg,dV\[Phi]zgdzg,\[Psi]pbar,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	rg=rgfun[wr,a,p,e,xg,const];
	zg=zgfun[wz,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	dV\[Phi]rgdrg=(2a*EEg*rg)/\[CapitalDelta]rg-(2a(rg-1)(-a*Lzg+EEg(a^2+rg^2)))/\[CapitalDelta]rg^2;
	dV\[Phi]zgdzg=(2Lzg*zg)/(1-zg^2)^2;
	
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		0
		,
		1/2(d\[Phi]sortCos\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const]-I*d\[Phi]sortSin\[Psi]d\[Lambda]fun[wr,wz,a,p,e,xg,const])(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])-drgdwrfun[wr,a,p,e,xg,const]*\[Xi]rortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffrad]*dV\[Phi]rgdrg-dzgdwzfun[wz,a,p,e,xg,const]*\[Xi]zortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeffpol]*dV\[Phi]zgdzg
	]
]


(* ::Subsection::Closed:: *)
(*Radial trajectory*)


(* ::Subsubsection::Closed:: *)
(* Term proportional to e^(-i wp)   (j=+1)*)


\[Delta]rortplusfun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,r1sortCos\[Psi],r2sortCos\[Psi],\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p,rsort,\[Psi]pbar},
	{EEg,Lzg,KKg}=const;
	
	r1sortCos\[Psi]=r1sortCos\[Psi]fun[wz,a,p,e,xg,const];
	r2sortCos\[Psi]=r2sortCos\[Psi]fun[wz,a,p,e,xg,const];
	rsort=1/2 (r1sortCos\[Psi]+r2sortCos\[Psi])+1/2 (r1sortCos\[Psi]-r2sortCos\[Psi])Sin\[Chi]rfun[wr,a,p,e,xg,const];
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[Abs[xg]==1,
		0,
		1/2 rsort(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])-drgdwrfun[wr,a,p,e,xg,const]\[Xi]rortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]
	]
]


(* ::Subsubsection::Closed:: *)
(*Term proportional to e^(i wp)   (j=-1)*)


\[Delta]rortminusfun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,r1sortCos\[Psi],r2sortCos\[Psi],\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p,rsort,\[Psi]pbar},
	{EEg,Lzg,KKg}=const;
	
	r1sortCos\[Psi]=r1sortCos\[Psi]fun[wz,a,p,e,xg,const];
	r2sortCos\[Psi]=r2sortCos\[Psi]fun[wz,a,p,e,xg,const];
	rsort=1/2 (r1sortCos\[Psi]+r2sortCos\[Psi])+1/2 (r1sortCos\[Psi]-r2sortCos\[Psi])Sin\[Chi]rfun[wr,a,p,e,xg,const];
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[Abs[xg]==1,
		0,
		1/2 rsort(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])-drgdwrfun[wr,a,p,e,xg,const]\[Xi]rortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]
	]
]


(* ::Subsection::Closed:: *)
(*Polar trajectory*)


(* ::Subsubsection::Closed:: *)
(* Term proportional to e^(-i wp)   (j=+1)*)


\[Delta]zortplusfun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,\[Chi]z,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p,zsort,\[Psi]pbar},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	zsort=z1sortCos\[Psi]fun[wr,a,p,e,xg,const];
	
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		zsort(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])
		,
		zsort(Cos[\[Psi]pbar[wr,wz]]-I*Sin[\[Psi]pbar[wr,wz]])-dzgdwzfun[wz,a,p,e,xg,const]\[Xi]zortplus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]
	]
]


(* ::Subsubsection::Closed:: *)
(* Term proportional to e^(i wp)   (j=-1)*)


\[Delta]zortminusfun[wr_,wz_,a_,p_,e_,xg_,const_,coeff_]:=Module[{EEg,Lzg,KKg,z1g,\[Chi]z,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p,zsort,\[Psi]pbar},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	zsort=z1sortCos\[Psi]fun[wr,a,p,e,xg,const];
	
	\[Psi]pbar[wr,wz]:=\[Psi]rfun[wr,a,p,e,xg,const]+\[Psi]zfun[wz,a,p,e,xg,const];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,xg];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,xg];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,xg];
	
	If[z1g==0,
		zsort(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])
		,
		zsort(Cos[\[Psi]pbar[wr,wz]]+I*Sin[\[Psi]pbar[wr,wz]])-dzgdwzfun[wz,a,p,e,xg,const]\[Xi]zortminus[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},coeff]
	]
]


(* ::Subsection::Closed:: *)
(*Coordinate trajectory*)


(* ::Text:: *)
(*Term proportional to e^(-i wp)   (j=+1)*)


\[CapitalDelta]tsortplusfun[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	 -Sum[If[n==0&&k==0,0,Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg+\[CapitalUpsilon]p))],{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Text:: *)
(*Term proportional to e^(i wp)   (j=-1)*)


\[CapitalDelta]tsortminusfun[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	 -Sum[If[n==0&&k==0,0,Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg-\[CapitalUpsilon]p))],{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Subsection::Closed:: *)
(*Azimuthal trajectory*)


(* ::Text:: *)
(*Term proportional to e^(-i wp)   (j=+1)*)


\[CapitalDelta]\[Phi]sortplusfun[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	 -Sum[Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg+\[CapitalUpsilon]p)),{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Text:: *)
(*Term proportional to e^(i wp)   (j=-1)*)


\[CapitalDelta]\[Phi]sortminusfun[wr_,wz_,\[CapitalUpsilon]g_,coeff_]:=Module[{dimr,dimz,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},
	dimr=(Dimensions[coeff][[1]]-1)/2;
	dimz=(Dimensions[coeff][[2]]-1)/2;
	{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p}=\[CapitalUpsilon]g;
	
	-Sum[Exp[-I(n*wr+k*wz)]coeff[[n+dimr+1,k+dimz+1]]/(I(n*\[CapitalUpsilon]rg+k*\[CapitalUpsilon]zg-\[CapitalUpsilon]p)),{n,-dimr,dimr},{k,-dimz,dimz}]
];


(* ::Section::Closed:: *)
(*Generic orbits - parallel component of the spin*)


(* ::Subsection::Closed:: *)
(*Generic - fixed turning points on average*)


KerrSpinOrbitCorrectionDHPar[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&), kmax_?EvenQ]:=Module[{EEg,Lzg,KKg,constmot,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg,\[Delta]constmot,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,stepsr,stepsz,\[CapitalDelta]tspar,\[CapitalDelta]\[Phi]spar,\[Delta]rpar,\[Delta]zpar,
	ExpniTable,ExpjkTable,\[Xi]rparcoeff,\[Xi]zparcoeff,dtrgd\[Lambda]coeff,dtzgd\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,d\[Phi]zgd\[Lambda]coeff,wrlist,wzlist,dtspard\[Lambda]coeff,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]trg,\[CapitalDelta]tzg,\[CapitalDelta]\[Phi]rg,\[CapitalDelta]\[Phi]zg,dtspard\[Lambda],drspard\[Lambda],dzspard\[Lambda],d\[Phi]spard\[Lambda]},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	KKg=KKgfun[a,p,e,x];
	constmot={EEg,Lzg,KKg};
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	
	(* spin correction constants of motion *)
	{\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg}=\[Delta]constmotionfun[a,p,e,x];
	\[Delta]constmot={\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg};
	  
	(* steps for numerical integration *)
	stepsr=4*nmax;
	stepsz=4*kmax;
	(* matrices of discrete Fourier transform *)
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,x}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
  
	(* Fourier coeffients geodesic functions *)
	dtrgd\[Lambda]coeff=Chop[(ExpniTable . Vtrgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	dtzgd\[Lambda]coeff=Chop[(Vtzgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];	
	d\[Phi]rgd\[Lambda]coeff=Chop[(ExpniTable . V\[Phi]rgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	d\[Phi]zgd\[Lambda]coeff=Chop[(V\[Phi]zgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];
	
	Print["Calculating \!\(\*SubscriptBox[\(\[Xi]\), \(r\)]\), \!\(\*SubscriptBox[\(\[Xi]\), \(z\)]\) Fourier coefficients"];
	\[Xi]rparcoeff=Chop[(ExpniTable . (Table[\[Delta]YroverYrgsparDH[wrlist[[i]],wzlist,a,p,e,x,constmot,\[Delta]constmot],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	\[Xi]zparcoeff=Chop[(ExpniTable . (Table[\[Delta]YzoverYzgsparDH[wrlist[[i]],wzlist,a,p,e,x,constmot,\[Delta]constmot],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[CapitalUpsilon]rg*\[Xi]rparcoeff[[nmax+1,kmax+1]];
	\[CapitalUpsilon]zs=\[CapitalUpsilon]zg*\[Xi]zparcoeff[[nmax+1,kmax+1]];
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(dt\), \(s\)]\)/d\[Lambda]"];
	dtspard\[Lambda]coeff=Chop[(ExpniTable . (Table[dtsparDHd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot,\[Delta]constmot,\[Xi]rparcoeff,\[Xi]zparcoeff],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	\[CapitalUpsilon]ts=dtspard\[Lambda]coeff[[nmax+1,kmax+1]];
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(d\[Phi]\), \(s\)]\)/d\[Lambda]"];
	d\[Phi]spard\[Lambda]coeff=Chop[(ExpniTable . (Table[d\[Phi]sparDHd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot,\[Delta]constmot,\[Xi]rparcoeff,\[Xi]zparcoeff],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	\[CapitalUpsilon]\[Phi]s=d\[Phi]spard\[Lambda]coeff[[nmax+1,kmax+1]];
	
	\[CapitalDelta]trg[wr_]:=\[CapitalDelta]trgfun[wr,dtrgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]tzg[wz_]:=\[CapitalDelta]tzgfun[wz,dtzgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[CapitalDelta]\[Phi]rgfun[wr,d\[Phi]rgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]\[Phi]zg[wz_]:=\[CapitalDelta]\[Phi]zgfun[wz,d\[Phi]zgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	
	\[CapitalDelta]tspar[wr_,wz_]:=\[CapitalDelta]spar[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},{\[CapitalUpsilon]rs,\[CapitalUpsilon]zs},dtspard\[Lambda]coeff,dtrgd\[Lambda]coeff,dtzgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]spar[wr_,wz_]:=\[CapitalDelta]spar[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},{\[CapitalUpsilon]rs,\[CapitalUpsilon]zs},d\[Phi]spard\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,d\[Phi]zgd\[Lambda]coeff];
	\[Delta]rpar[wr_,wz_]:=\[Delta]rparDHfun[wr,wz,a,p,e,x,constmot,\[Delta]constmot,\[Xi]rparcoeff];
	\[Delta]zpar[wr_,wz_]:=\[Delta]zparDHfun[wr,wz,a,p,e,x,constmot,\[Delta]constmot,\[Xi]zparcoeff];
	
	dtspard\[Lambda][wr_,wz_]:=dtsparDHd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Delta]constmot,\[Xi]rparcoeff,\[Xi]zparcoeff];
	d\[Phi]spard\[Lambda][wr_,wz_]:=d\[Phi]sparDHd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Delta]constmot,\[Xi]rparcoeff,\[Xi]zparcoeff];
	drspard\[Lambda][wr_,wz_]:=drsparDHd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Delta]constmot,\[Xi]rparcoeff];
	dzspard\[Lambda][wr_,wz_]:=dzsparDHd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Delta]constmot,\[Xi]zparcoeff];
			
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Kg"->KKg,
	 "Es"->\[Delta]EEg,
	 "Jzs"->\[Delta]Lzg,
	 "Ks"->\[Delta]KKg,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]zs/\[CapitalUpsilon]tg-\[CapitalUpsilon]zg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectories*)
	 "\[CapitalDelta]trg"->Function[{wr},\[CapitalDelta]trg[wr]],
	 "\[CapitalDelta]tzg"->Function[{wz},\[CapitalDelta]tzg[wz]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[CapitalDelta]\[Phi]rg[wr]],
	 "\[CapitalDelta]\[Phi]zg"->Function[{wz},\[CapitalDelta]\[Phi]zg[wz]],
	 (*Keys corrections trajectories*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr,wz},\[CapitalDelta]tspar[wr,wz]],
	 "\[Delta]rpar"->Function[{wr,wz},\[Delta]rpar[wr,wz]],
	 "\[Delta]zpar"->Function[{wr,wz},\[Delta]zpar[wr,wz]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr,wz},\[CapitalDelta]\[Phi]spar[wr,wz]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr,wz},dtspard\[Lambda][wr,wz]],
	 "\[Delta]vrpar"->Function[{wr,wz},drspard\[Lambda][wr,wz]],
	 "\[Delta]vzpar"->Function[{wr,wz},dzspard\[Lambda][wr,wz]],
	 "\[Delta]v\[Phi]par"->Function[{wr,wz},d\[Phi]spard\[Lambda][wr,wz]]
	|>
]


(* ::Subsection::Closed:: *)
(*Generic - fixed constants of motion*)


KerrSpinOrbitCorrectionFCPar[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&), kmax_?EvenQ]:=Module[{EEg,Lzg,KKg,constmot,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg,\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s,stepsr,stepsz,\[CapitalDelta]tspar,\[CapitalDelta]\[Phi]spar,\[Delta]rpar,\[Delta]zpar,
	ExpniTable,ExpjkTable,\[Xi]rparcoeff,\[Xi]zparcoeff,dtrgd\[Lambda]coeff,dtzgd\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,d\[Phi]zgd\[Lambda]coeff,wrlist,wzlist,dtspard\[Lambda]coeff,d\[Phi]spard\[Lambda]coeff,\[CapitalDelta]trg,\[CapitalDelta]tzg,\[CapitalDelta]\[Phi]rg,\[CapitalDelta]\[Phi]zg,dtspard\[Lambda],drspard\[Lambda],dzspard\[Lambda],d\[Phi]spard\[Lambda]},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	KKg=KKgfun[a,p,e,x];
	constmot={EEg,Lzg,KKg};
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	
	(* spin correction constants of motion *)
	{\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg}=\[Delta]constmotionfun[a,p,e,x];
	\[Delta]constmot={\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg};
	  
	(* steps for numerical integration *)
	stepsr=4*nmax;
	stepsz=4*kmax;
	(* matrices of discrete Fourier transform *)
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,x}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
  
	(* Fourier coeffients geodesic functions *)
	dtrgd\[Lambda]coeff=Chop[(ExpniTable . Vtrgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	dtzgd\[Lambda]coeff=Chop[(Vtzgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];	
	d\[Phi]rgd\[Lambda]coeff=Chop[(ExpniTable . V\[Phi]rgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	d\[Phi]zgd\[Lambda]coeff=Chop[(V\[Phi]zgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];
	
	Print["Calculating \!\(\*SubscriptBox[\(\[Xi]\), \(r\)]\), \!\(\*SubscriptBox[\(\[Xi]\), \(z\)]\) Fourier coefficients"];
	\[Xi]rparcoeff=Chop[(ExpniTable . (Table[\[Delta]YroverYrgsparFC[wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	\[Xi]zparcoeff=Chop[(ExpniTable . (Table[\[Delta]YzoverYzgsparFC[wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	
	(* spin correction radial and polar frequencies *)
	\[CapitalUpsilon]rs=\[CapitalUpsilon]rg*\[Xi]rparcoeff[[nmax+1,kmax+1]];
	\[CapitalUpsilon]zs=\[CapitalUpsilon]zg*\[Xi]zparcoeff[[nmax+1,kmax+1]];
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(dt\), \(s\)]\)/d\[Lambda]"];
	dtspard\[Lambda]coeff=Chop[(ExpniTable . (Table[dtsparFCd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot,\[Xi]rparcoeff,\[Xi]zparcoeff],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	\[CapitalUpsilon]ts=dtspard\[Lambda]coeff[[nmax+1,kmax+1]];
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(d\[Phi]\), \(s\)]\)/d\[Lambda]"];
	d\[Phi]spard\[Lambda]coeff=Chop[(ExpniTable . (Table[d\[Phi]sparFCd\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot,\[Xi]rparcoeff,\[Xi]zparcoeff],{i,stepsr}]) . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	\[CapitalUpsilon]\[Phi]s=d\[Phi]spard\[Lambda]coeff[[nmax+1,kmax+1]];
	
	\[CapitalDelta]trg[wr_]:=\[CapitalDelta]trgfun[wr,dtrgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]tzg[wz_]:=\[CapitalDelta]tzgfun[wz,dtzgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[CapitalDelta]\[Phi]rgfun[wr,d\[Phi]rgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]\[Phi]zg[wz_]:=\[CapitalDelta]\[Phi]zgfun[wz,d\[Phi]zgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	
	\[CapitalDelta]tspar[wr_,wz_]:=\[CapitalDelta]spar[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},{\[CapitalUpsilon]rs,\[CapitalUpsilon]zs},dtspard\[Lambda]coeff,dtrgd\[Lambda]coeff,dtzgd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]spar[wr_,wz_]:=\[CapitalDelta]spar[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg},{\[CapitalUpsilon]rs,\[CapitalUpsilon]zs},d\[Phi]spard\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,d\[Phi]zgd\[Lambda]coeff];
	\[Delta]rpar[wr_,wz_]:=\[Delta]rparFCfun[wr,wz,a,p,e,x,constmot,\[Xi]rparcoeff];
	\[Delta]zpar[wr_,wz_]:=\[Delta]zparFCfun[wr,wz,a,p,e,x,constmot,\[Xi]zparcoeff];
	
	dtspard\[Lambda][wr_,wz_]:=dtsparFCd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rparcoeff,\[Xi]zparcoeff];
	d\[Phi]spard\[Lambda][wr_,wz_]:=d\[Phi]sparFCd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rparcoeff,\[Xi]zparcoeff];
	drspard\[Lambda][wr_,wz_]:=drsparFCd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rparcoeff];
	dzspard\[Lambda][wr_,wz_]:=dzsparFCd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]zparcoeff];
			
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Kg"->KKg,
	 "Es"->\[Delta]EEg,
	 "Jzs"->\[Delta]Lzg,
	 "Ks"->\[Delta]KKg,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoFrequenciesCorrection"->{\[CapitalUpsilon]ts,\[CapitalUpsilon]rs,\[CapitalUpsilon]zs,\[CapitalUpsilon]\[Phi]s},
	 "BLFrequenciesCorrection"->{\[CapitalUpsilon]rs/\[CapitalUpsilon]tg-\[CapitalUpsilon]rg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]zs/\[CapitalUpsilon]tg-\[CapitalUpsilon]zg/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts,\[CapitalUpsilon]\[Phi]s/\[CapitalUpsilon]tg-\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg^2*\[CapitalUpsilon]ts},
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectories*)
	 "\[CapitalDelta]trg"->Function[{wr},\[CapitalDelta]trg[wr]],
	 "\[CapitalDelta]tzg"->Function[{wz},\[CapitalDelta]tzg[wz]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[CapitalDelta]\[Phi]rg[wr]],
	 "\[CapitalDelta]\[Phi]zg"->Function[{wz},\[CapitalDelta]\[Phi]zg[wz]],
	 (*Keys corrections trajectories*)
	 "\[CapitalDelta]\[Delta]tpar"->Function[{wr,wz},\[CapitalDelta]tspar[wr,wz]],
	 "\[Delta]rpar"->Function[{wr,wz},\[Delta]rpar[wr,wz]],
	 "\[Delta]zpar"->Function[{wr,wz},\[Delta]zpar[wr,wz]],
	 "\[CapitalDelta]\[Delta]\[Phi]par"->Function[{wr,wz},\[CapitalDelta]\[Phi]spar[wr,wz]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtpar"->Function[{wr,wz},dtspard\[Lambda][wr,wz]],
	 "\[Delta]vrpar"->Function[{wr,wz},drspard\[Lambda][wr,wz]],
	 "\[Delta]vzpar"->Function[{wr,wz},dzspard\[Lambda][wr,wz]],
	 "\[Delta]v\[Phi]par"->Function[{wr,wz},d\[Phi]spard\[Lambda][wr,wz]]
	|>
]


(* ::Subsection::Closed:: *)
(*Generic - map between parametrizations*)


(* ::Subsubsection::Closed:: *)
(*Functions for mapping radial and polar trajectories in the two parametrizations*)


\[Eta]rfun[wr_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,rspar,\[Chi]r,r1savgsub,r2savgsub},
	{EEg,Lzg,KKg}=const;
	
	r1savgsub=r1savgFC[a,p,e,xg];
	r2savgsub=r2savgFC[a,p,e,xg];
	
	-(drgdr1gfun[wr,a,p,e,xg,const]r1savgsub+drgdr2gfun[wr,a,p,e,xg,const]r2savgsub)+(\[Delta]const[[1]]drgdEgfun[wr,a,p,e,xg,const]+\[Delta]const[[2]]drgdLzgfun[wr,a,p,e,xg,const]+\[Delta]const[[3]]drgdKgfun[wr,a,p,e,xg,const])
]


\[Eta]zfun[wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,z1g,z1savgsub,\[Chi]z},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	z1savgsub=z1savgFC[a,p,e,xg];
	
	- z1g*Cos[\[Chi]z]z1savgsub d\[Chi]zdz1gfun[wz,a,p,e,xg,const]+(\[Delta]const[[1]]dzgdEgfun[wz,a,p,e,xg,const]+\[Delta]const[[2]]dzgdLzgfun[wz,a,p,e,xg,const])-Sin[\[Chi]z]z1savgsub
]


(* ::Subsubsection::Closed:: *)
(*Functions for mapping radial and polar velocities in the two parametrizations*)


drsd\[Lambda]FCmapDHfun[wr_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,r1g,r2g,r3g,r4g,rg,\[CapitalDelta]rg,Rrg,dVrgdrg,dVrgdEg,dVrgdLzg,dVrgdKg,Yrg},
	{EEg,Lzg,KKg}=const;
	
	r1g=p/(1-e);
	r2g=p/(1+e);
	r3g=1/(1-EEg^2)-(r1g+r2g)/2+Sqrt[(((r1g+r2g)/2-1/(1-EEg^2))^2-(a^2(KKg-(Lzg-a EEg)^2))/(r1g r2g (1-EEg^2)))];
	r4g=1/r3g(a^2(KKg-(Lzg-a EEg)^2))/(r1g*r2g*(1-EEg^2));
	
	rg=rgfun[wr,a,p,e,xg,const];
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	Rrg=(-a*Lzg+EEg(a^2+rg^2));
	dVrgdrg=-(-2+2rg)(KKg+rg^2)-2rg*\[CapitalDelta]rg+4EEg*rg*Rrg;
	dVrgdEg=2(rg^2+a^2)Rrg;
	dVrgdLzg=-2a*Rrg;
	dVrgdKg=-\[CapitalDelta]rg;
	Yrg=Sqrt[(1-EEg^2)(-rg+r3g)(-rg+r4g)];
	
	1/(2Sign[\[Pi]-Mod[wr,2 \[Pi]]]Sqrt[(r1g-rg)(rg-r2g)]Yrg)(\[Eta]rfun[wr,a,p,e,xg,const,\[Delta]const]dVrgdrg+\[Delta]const[[1]]dVrgdEg+\[Delta]const[[2]]dVrgdLzg+\[Delta]const[[3]]dVrgdKg)
]


dzsd\[Lambda]FCmapDHfun[wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,z1g,z2g,rspar,\[Chi]z,zg,Zzg,Yzg,dVzgdzg,dVzgdEg,dVzgdLzg,dVzgdKg},
	{EEg,Lzg,KKg}=const;
	
	z1g=Sqrt[1-xg^2];
	z2g=Sqrt[a^2(1-EEg^2)+Lzg^2/(1-z1g^2)];
	
	\[Chi]z=\[Chi]zfun[wz,a,p,e,xg,const];
	
	zg=zgfun[wz,a,p,e,xg,const];
	Zzg=(Lzg-a*EEg(1-zg^2));
	Yzg=Sqrt[-a^2(1-EEg^2)zg^2+z2g^2];
	dVzgdzg=-2 zg (KKg+a^2 (1-2 zg^2))-2 Zzg 2 a EEg zg;
	dVzgdEg=2a(1-zg^2)Zzg;
	dVzgdLzg=-2Zzg;
	dVzgdKg=1-zg^2;
	
	1/(2z1g*Cos[\[Chi]z]Yzg)(\[Eta]zfun[wz,a,p,e,xg,const,\[Delta]const]dVzgdzg+\[Delta]const[[1]]dVzgdEg+\[Delta]const[[2]]dVzgdLzg+\[Delta]const[[3]]dVzgdKg)
]


(* ::Subsubsection::Closed:: *)
(*Functions for mapping coordinate time and azimuthal velocities in the two parametrizations*)


\[Eta]trfun[wr_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,rg,\[CapitalDelta]rg,dVtrgdrg,dVtrgdEg,dVtrgdLzg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2rg+rg^2);
	dVtrgdrg=(2a*Lzg(rg^2-a^2))/\[CapitalDelta]rg^2+(4EEg*rg(a^2+rg^2))/\[CapitalDelta]rg-(2EEg(rg-1)(a^2+rg^2)^2)/\[CapitalDelta]rg^2;
	dVtrgdEg=(a^2+rg^2)^2/\[CapitalDelta]rg;
	dVtrgdLzg=-((2a*rg)/\[CapitalDelta]rg);
	
	\[Eta]rfun[wr,a,p,e,xg,const,\[Delta]const]dVtrgdrg+\[Delta]const[[1]]dVtrgdEg+\[Delta]const[[2]]dVtrgdLzg
]


\[Eta]tzfun[wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,zg,dVtzgdzg,dVtzgdEg},
	{EEg,Lzg,KKg}=const;
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	dVtzgdzg=2 a^2 EEg zg;
	dVtzgdEg=-a^2 (1-zg^2);
	
	\[Eta]zfun[wz,a,p,e,xg,const,\[Delta]const]dVtzgdzg+\[Delta]const[[1]]dVtzgdEg
]


\[Eta]\[Phi]rfun[wr_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,rg,\[CapitalDelta]rg,dV\[Phi]rgdrg,dV\[Phi]rgdEg,dV\[Phi]rgdLzg},
	{EEg,Lzg,KKg}=const;
	
	rg=rgfun[wr,a,p,e,xg,const];
	
	\[CapitalDelta]rg=(a^2-2 rg+rg^2);
	dV\[Phi]rgdrg=(2 a EEg rg)/\[CapitalDelta]rg-(2a (rg-1) (-a Lzg+EEg (a^2+rg^2)))/\[CapitalDelta]rg^2;
	dV\[Phi]rgdEg=(2 a rg)/\[CapitalDelta]rg;
	dV\[Phi]rgdLzg=-(a^2/\[CapitalDelta]rg);
	
	\[Eta]rfun[wr,a,p,e,xg,const,\[Delta]const]dV\[Phi]rgdrg+\[Delta]const[[1]]dV\[Phi]rgdEg+\[Delta]const[[2]]dV\[Phi]rgdLzg
]


\[Eta]\[Phi]zfun[wz_,a_,p_,e_,xg_,const_,\[Delta]const_]:=Module[{EEg,Lzg,KKg,zg,dV\[Phi]zgdzg,dV\[Phi]zgdLzg},
	{EEg,Lzg,KKg}=const;
	
	zg=zgfun[wz,a,p,e,xg,const];
	
	dV\[Phi]zgdzg=(2 Lzg zg)/(1-zg^2)^2;
	dV\[Phi]zgdLzg=1/(1-zg^2);
	
	\[Eta]zfun[wz,a,p,e,xg,const,\[Delta]const]dV\[Phi]zgdzg+\[Delta]const[[2]]dV\[Phi]zgdLzg
]


(* ::Subsubsection::Closed:: *)
(*Functions for mapping coordinate time and azimuthal trajectories in the two parametrizations*)


\[Eta]tr\[CapitalDelta]tsfun[wr_,\[CapitalUpsilon]rg_,dev_,coeff_,coeffVtr_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2/(n*\[CapitalUpsilon]rg)Sin[n wr]coeff[[n+dim+1]],{n,1,dim}]+dev*Sum[2/(n*\[CapitalUpsilon]rg^2)Sin[n wr]coeffVtr[[n+dim+1]],{n,1,dim}]
];


\[Eta]tz\[CapitalDelta]tsfun[wz_,\[CapitalUpsilon]zg_,dev_,coeff_,coeffVtz_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2/(k*\[CapitalUpsilon]zg)Sin[k wz]coeff[[k+dim+1]],{k,1,dim}]+dev*Sum[2/(k*\[CapitalUpsilon]zg^2)Sin[k wz]coeffVtz[[k+dim+1]],{k,1,dim}]
];


\[Eta]\[Phi]r\[CapitalDelta]\[Phi]sfun[wr_,\[CapitalUpsilon]rg_,dev_,coeff_,coeffV\[Phi]r_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2/(n*\[CapitalUpsilon]rg)Sin[n wr]coeff[[n+dim+1]],{n,1,dim}]+dev*Sum[2/(n*\[CapitalUpsilon]rg^2)Sin[n wr]coeffV\[Phi]r[[n+dim+1]],{n,1,dim}]
];


\[Eta]\[Phi]z\[CapitalDelta]\[Phi]sfun[wz_,\[CapitalUpsilon]zg_,dev_,coeff_,coeffV\[Phi]z_]:=Module[{dim},
	dim=(Length[coeff]-1)/2;
	Sum[2/(k*\[CapitalUpsilon]zg)Sin[k wz]coeff[[k+dim+1]],{k,1,dim}]+dev*Sum[2/(k*\[CapitalUpsilon]zg^2)Sin[k wz]coeffV\[Phi]z[[k+dim+1]],{k,1,dim}]
];


(* ::Subsubsection::Closed:: *)
(*Wrapper*)


KerrSpinOrbitCorrectionDHMapFC[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&), kmax_?EvenQ]:=Module[{EEg,Lzg,KKg,constmot,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg,\[Delta]constmot,devtime,devrad,devpol,devphi,stepsr,stepsz,
	ExpniTable,ExpjkTable,dtrgd\[Lambda]coeff,dtzgd\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,d\[Phi]zgd\[Lambda]coeff,wrlist,wzlist,\[Eta]trcoeff,\[Eta]tzcoeff,\[Eta]\[Phi]rcoeff,\[Eta]\[Phi]zcoeff,\[Eta]r,\[Eta]z,\[Eta]tr\[CapitalDelta]ts,\[Eta]tz\[CapitalDelta]ts,\[Eta]\[Phi]r\[CapitalDelta]\[Phi]s,\[Eta]\[Phi]z\[CapitalDelta]\[Phi]s,drsd\[Lambda]FCmapDH,dzsd\[Lambda]FCmapDH,
	\[Eta]tr,\[Eta]tz,\[Eta]\[Phi]r,\[Eta]\[Phi]z,\[CapitalDelta]trg,\[CapitalDelta]tzg,\[CapitalDelta]\[Phi]rg,\[CapitalDelta]\[Phi]zg},
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	KKg=KKgfun[a,p,e,x];
	constmot={EEg,Lzg,KKg};
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	
	(* spin correction constants of motion *)
	{\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg}=\[Delta]constmotionfun[a,p,e,x];
	\[Delta]constmot={\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg};
	
	devtime=freqdevtime[a,p,e,x];
	devrad=freqdevrad[a,p,e,x];
	devpol=freqdevpol[a,p,e,x];
	devphi=freqdev\[Phi][a,p,e,x];
	  
	(* steps for numerical integration *)
	stepsr=4*nmax;
	stepsz=4*kmax;
	(* matrices of discrete Fourier transform *)
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,x}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
  
	(* Fourier coeffients geodesic functions *)
	dtrgd\[Lambda]coeff=Chop[(ExpniTable . Vtrgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	dtzgd\[Lambda]coeff=Chop[(Vtzgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];	
	d\[Phi]rgd\[Lambda]coeff=Chop[(ExpniTable . V\[Phi]rgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	d\[Phi]zgd\[Lambda]coeff=Chop[(V\[Phi]zgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];
	
	\[CapitalDelta]trg[wr_]:=\[CapitalDelta]trgfun[wr,dtrgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]tzg[wz_]:=\[CapitalDelta]tzgfun[wz,dtzgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[CapitalDelta]\[Phi]rgfun[wr,d\[Phi]rgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]\[Phi]zg[wz_]:=\[CapitalDelta]\[Phi]zgfun[wz,d\[Phi]zgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	
	(* Fourier coeffients map functions *)
	\[Eta]trcoeff=Chop[(ExpniTable . \[Eta]trfun[wrlist,a,p,e,x,constmot,\[Delta]constmot])/stepsr,10^(-16)];
	\[Eta]tzcoeff=Chop[(\[Eta]tzfun[wzlist,a,p,e,x,constmot,\[Delta]constmot] . ExpjkTable)/stepsz,10^(-16)];	
	\[Eta]\[Phi]rcoeff=Chop[(ExpniTable . \[Eta]\[Phi]rfun[wrlist,a,p,e,x,constmot,\[Delta]constmot])/stepsr,10^(-16)];
	\[Eta]\[Phi]zcoeff=Chop[(\[Eta]\[Phi]zfun[wzlist,a,p,e,x,constmot,\[Delta]constmot] . ExpjkTable)/stepsz,10^(-16)];
		
	\[Eta]r[wr_]:=\[Eta]rfun[wr,a,p,e,x,constmot,\[Delta]constmot];
	\[Eta]z[wz_]:=\[Eta]zfun[wz,a,p,e,x,constmot,\[Delta]constmot];
	\[Eta]tr\[CapitalDelta]ts[wr_]:=\[Eta]tr\[CapitalDelta]tsfun[wr,\[CapitalUpsilon]rg,devrad,\[Eta]trcoeff,dtrgd\[Lambda]coeff];
	\[Eta]tz\[CapitalDelta]ts[wz_]:=\[Eta]tz\[CapitalDelta]tsfun[wz,\[CapitalUpsilon]zg,devpol,\[Eta]tzcoeff,dtzgd\[Lambda]coeff];
	\[Eta]\[Phi]r\[CapitalDelta]\[Phi]s[wr_]:=\[Eta]\[Phi]r\[CapitalDelta]\[Phi]sfun[wr,\[CapitalUpsilon]rg,devrad,\[Eta]\[Phi]rcoeff,d\[Phi]rgd\[Lambda]coeff];
	\[Eta]\[Phi]z\[CapitalDelta]\[Phi]s[wz_]:=\[Eta]\[Phi]z\[CapitalDelta]\[Phi]sfun[wz,\[CapitalUpsilon]zg,devpol,\[Eta]\[Phi]zcoeff,d\[Phi]zgd\[Lambda]coeff];
			
	drsd\[Lambda]FCmapDH[wr_]:=drsd\[Lambda]FCmapDHfun[wr,a,p,e,x,constmot,\[Delta]constmot];
	dzsd\[Lambda]FCmapDH[wz_]:=dzsd\[Lambda]FCmapDHfun[wz,a,p,e,x,constmot,\[Delta]constmot];
	\[Eta]tr[wr_]:=\[Eta]trfun[wr,a,p,e,x,constmot,\[Delta]constmot];
	\[Eta]tz[wz_]:=\[Eta]tzfun[wz,a,p,e,x,constmot,\[Delta]constmot];
	\[Eta]\[Phi]r[wr_]:=\[Eta]\[Phi]rfun[wr,a,p,e,x,constmot,\[Delta]constmot];
	\[Eta]\[Phi]z[wz_]:=\[Eta]\[Phi]zfun[wz,a,p,e,x,constmot,\[Delta]constmot];
			
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Kg"->KKg,
	 "Es"->\[Delta]EEg,
	 "Jzs"->\[Delta]Lzg,
	 "Ks"->\[Delta]KKg,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectories*)
	 "\[CapitalDelta]trg"->Function[{wr},\[CapitalDelta]trg[wr]],
	 "\[CapitalDelta]tzg"->Function[{wz},\[CapitalDelta]tzg[wz]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[CapitalDelta]\[Phi]rg[wr]],
	 "\[CapitalDelta]\[Phi]zg"->Function[{wz},\[CapitalDelta]\[Phi]zg[wz]],
	 (*Keys maps trajectories*)
	 "\[Eta]r"->Function[{wr},\[Eta]r[wr]],
	 "\[Eta]z"->Function[{wz},\[Eta]z[wz]],
	 "\[Eta]tr\[CapitalDelta]\[Delta]t"->Function[{wr},\[Eta]tr\[CapitalDelta]ts[wr]],
	 "\[Eta]tz\[CapitalDelta]\[Delta]t"->Function[{wz},\[Eta]tz\[CapitalDelta]ts[wz]],
	 "\[Eta]\[Phi]r\[CapitalDelta]\[Delta]\[Phi]"->Function[{wr},\[Eta]\[Phi]r\[CapitalDelta]\[Phi]s[wr]],
	 "\[Eta]\[Phi]z\[CapitalDelta]\[Delta]\[Phi]"->Function[{wz},\[Eta]\[Phi]z\[CapitalDelta]\[Phi]s[wz]],
	 (*Keys maps velocities*)
	 "d\[Delta]rd\[Lambda]FCmapDH"->Function[{wr},drsd\[Lambda]FCmapDH[wr]],
	 "d\[Delta]zd\[Lambda]FCmapDH"->Function[{wz},dzsd\[Lambda]FCmapDH[wz]],
	 "\[Eta]tr"->Function[{wr},\[Eta]tr[wr]],
	 "\[Eta]tz"->Function[{wz},\[Eta]tz[wz]],
	 "\[Eta]\[Phi]r"->Function[{wr},\[Eta]\[Phi]r[wr]],
	 "\[Eta]\[Phi]z"->Function[{wz},\[Eta]\[Phi]z[wz]],
	 (*Keys maps frequencies*)
	 "freqdevtime"->devtime,
	 "freqdevpol"->devpol,
	 "freqdevrad"->devrad,
	 "freqdevphi"->devphi
	|>
]


(* ::Section::Closed:: *)
(*Generic orbits - orthogonal component of the spin*)


(* ::Subsection::Closed:: *)
(*Generic*)


KerrSpinOrbitCorrectionOrt[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&), kmax_?EvenQ]:=Module[{EEg,Lzg,KKg,constmot,\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g,\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg,\[Delta]constmot,\[CapitalUpsilon]p,stepsr,stepsz,
	ExpniTable,ExpjkTable,wrlist,wzlist,dtrgd\[Lambda]coeff,dtzgd\[Lambda]coeff,d\[Phi]rgd\[Lambda]coeff,d\[Phi]zgd\[Lambda]coeff,\[Psi]pbar,\[Psi]pbarlist,Aplus,Aminus,\[Delta]YroverYrgsortSin\[Psi]list,\[Delta]YroverYrgsortCos\[Psi]list,\[Delta]YzoverYzgsortSin\[Psi]list,\[Delta]YzoverYzgsortCos\[Psi]list,
	\[Xi]rortpluscoeff,\[Xi]rortminuscoeff,\[Xi]zortpluscoeff,\[Xi]zortminuscoeff,dtsortSin\[Psi]d\[Lambda]list,dtsortCos\[Psi]d\[Lambda]list,d\[Phi]sortSin\[Psi]d\[Lambda]list,d\[Phi]sortCos\[Psi]d\[Lambda]list,\[CapitalDelta]trg,\[CapitalDelta]tzg,dVtrgdwrlist,dVtzgdwzlist,dV\[Phi]rgdwrlist,dV\[Phi]zgdwzlist,\[CapitalDelta]\[Phi]rg,\[CapitalDelta]\[Phi]zg,
	dtsortplusd\[Lambda]list,dtsortminusd\[Lambda]list,d\[Phi]sortplusd\[Lambda]list,d\[Phi]sortminusd\[Lambda]list,\[Xi]rortpluslist,\[Xi]rortminuslist,\[Xi]zortpluslist,\[Xi]zortminuslist,dtsortplusd\[Lambda]coeff,dtsortminusd\[Lambda]coeff,d\[Phi]sortplusd\[Lambda]coeff,d\[Phi]sortminusd\[Lambda]coeff,
	\[Delta]rortplus,\[Delta]rortminus,\[Delta]rort,\[Delta]zortplus,\[Delta]zortminus,\[Delta]zort,drsortplusd\[Lambda],drsortminusd\[Lambda],drsortd\[Lambda],dzsortplusd\[Lambda],dzsortminusd\[Lambda],dzsortd\[Lambda],dtsortplusd\[Lambda],dtsortminusd\[Lambda],dtsortd\[Lambda],d\[Phi]sortplusd\[Lambda],d\[Phi]sortminusd\[Lambda],d\[Phi]sortd\[Lambda],
	\[CapitalDelta]tsortplus,\[CapitalDelta]tsortminus,\[CapitalDelta]tsort,\[CapitalDelta]\[Phi]sortplus,\[CapitalDelta]\[Phi]sortminus,\[CapitalDelta]\[Phi]sort
   },
	(* geodesic constants of motion *)
	EEg=EEgfun[a,p,e,x];
	Lzg=Lzgfun[a,p,e,x];
	KKg=KKgfun[a,p,e,x];
	constmot={EEg,Lzg,KKg};
	
	(* geodesic frequencies *)
	\[CapitalUpsilon]tg=\[CapitalUpsilon]tgrfun[a,p,e,x]+\[CapitalUpsilon]tgzfun[a,p,e,x];
	\[CapitalUpsilon]rg=\[CapitalUpsilon]rgfun[a,p,e,x];
	\[CapitalUpsilon]zg=\[CapitalUpsilon]zgfun[a,p,e,x];
	\[CapitalUpsilon]\[Phi]g=\[CapitalUpsilon]\[Phi]grfun[a,p,e,x]+\[CapitalUpsilon]\[Phi]gzfun[a,p,e,x];
	\[CapitalUpsilon]p=\[CapitalUpsilon]pfun[a,p,e,x];
	
	(* spin correction constants of motion *)
	{\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg}=\[Delta]constmotionfun[a,p,e,x];
	\[Delta]constmot={\[Delta]EEg,\[Delta]Lzg,\[Delta]KKg};
	  
	(* steps for numerical integration *)
	stepsr=4*nmax;
	stepsz=4*kmax;
	(* matrices of discrete Fourier transform *)
	ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
	ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/stepsz],Precision[{a,p,e,x}]],{j,1,stepsz},{k,-kmax,kmax}];
	
	wrlist=Table[wr,{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)}];
	wzlist=Table[wz,{wz,2Pi/(stepsz)/2,2Pi,2Pi/(stepsz)}];
	
	(* Fourier coeffients geodesic functions *)
	dtrgd\[Lambda]coeff=Chop[(ExpniTable . Vtrgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	dtzgd\[Lambda]coeff=Chop[(Vtzgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];	
	d\[Phi]rgd\[Lambda]coeff=Chop[(ExpniTable . V\[Phi]rgfun[wrlist,a,p,e,x,constmot])/stepsr,10^(-16)];
	d\[Phi]zgd\[Lambda]coeff=Chop[(V\[Phi]zgfun[wzlist,a,p,e,x,constmot] . ExpjkTable)/stepsz,10^(-16)];
	
	\[Psi]pbar[wr_,wz_]:=\[Psi]rfun[wr,a,p,e,x,constmot]+\[Psi]zfun[wz,a,p,e,x,constmot];
	
	\[Psi]pbarlist=Table[\[Psi]pbar[wrlist[[i]],wzlist],{i,stepsr}];
	Aplus=Cos[\[Psi]pbarlist]+I*Sin[\[Psi]pbarlist];
	Aminus=Cos[\[Psi]pbarlist]-I*Sin[\[Psi]pbarlist];
	
	(* Fourier coeffients function for near-identiy transformation *)
	Print["Calculating \!\(\*SubscriptBox[\(\[Xi]\), \(r\)]\), \!\(\*SubscriptBox[\(\[Xi]\), \(z\)]\) Fourier coefficients"];
	\[Delta]YroverYrgsortSin\[Psi]list=Table[\[Delta]YroverYrgsortSin\[Psi][wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	\[Delta]YroverYrgsortCos\[Psi]list=Table[\[Delta]YroverYrgsortCos\[Psi][wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	\[Delta]YzoverYzgsortSin\[Psi]list=Table[\[Delta]YzoverYzgsortSin\[Psi][wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	\[Delta]YzoverYzgsortCos\[Psi]list=Table[\[Delta]YzoverYzgsortCos\[Psi][wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	
	\[Xi]rortpluscoeff=Chop[(ExpniTable . (1/2(\[Delta]YroverYrgsortCos\[Psi]list+I*\[Delta]YroverYrgsortSin\[Psi]list)Aminus) . ExpjkTable)1/(stepsr*stepsz),10^(-16)];
	\[Xi]rortminuscoeff=Chop[(ExpniTable . (1/2(\[Delta]YroverYrgsortCos\[Psi]list-I*\[Delta]YroverYrgsortSin\[Psi]list)Aplus) . ExpjkTable)1/(stepsr*stepsz),10^(-16)];	
	\[Xi]zortpluscoeff=Chop[(ExpniTable . (1/2(\[Delta]YzoverYzgsortCos\[Psi]list+I*\[Delta]YzoverYzgsortSin\[Psi]list)Aminus) . ExpjkTable)1/(stepsr*stepsz),10^(-16)];  
	\[Xi]zortminuscoeff=Chop[(ExpniTable . (1/2(\[Delta]YzoverYzgsortCos\[Psi]list-I*\[Delta]YzoverYzgsortSin\[Psi]list)Aplus) . ExpjkTable)1/(stepsr*stepsz),10^(-16)];
	
	(* Here I am sampling the near-identiy transformation functions*)
	Print["Sampling \!\(\*SubscriptBox[\(\[Xi]\), \(r\)]\) function"];
	\[Xi]rortpluslist=Table[\[Xi]rortplus[wrlist[[i]],wzlist,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},\[Xi]rortpluscoeff],{i,stepsr}];
	\[Xi]rortminuslist=Table[\[Xi]rortminus[wrlist[[i]],wzlist,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},\[Xi]rortminuscoeff],{i,stepsr}];
	Print["Sampling \!\(\*SubscriptBox[\(\[Xi]\), \(z\)]\) function"];
	\[Xi]zortpluslist=Table[\[Xi]zortplus[wrlist[[i]],wzlist,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},\[Xi]zortpluscoeff],{i,stepsr}];
	\[Xi]zortminuslist=Table[\[Xi]zortminus[wrlist[[i]],wzlist,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},\[Xi]zortminuscoeff],{i,stepsr}];
		
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(dt\), \(s\)]\)/d\[Lambda]"];
	dtsortSin\[Psi]d\[Lambda]list=Table[dtsortSin\[Psi]d\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	dtsortCos\[Psi]d\[Lambda]list=Table[dtsortCos\[Psi]d\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	
	dVtrgdwrlist=dVtrgdrgfun[wrlist,a,p,e,x,constmot]*drgdwrfun[wrlist,a,p,e,x,constmot];
	dVtzgdwzlist=dVtzgdzgfun[wzlist,a,p,e,x,constmot]*dzgdwzfun[wzlist,a,p,e,x,constmot];
	dtsortplusd\[Lambda]list=1/2(dtsortCos\[Psi]d\[Lambda]list+I*dtsortSin\[Psi]d\[Lambda]list)Aminus-Table[dVtrgdwrlist[[i]]*\[Xi]rortpluslist[[i,;;]]+dVtzgdwzlist*\[Xi]zortpluslist[[i,;;]],{i,stepsr}];
	dtsortminusd\[Lambda]list=1/2(dtsortCos\[Psi]d\[Lambda]list-I*dtsortSin\[Psi]d\[Lambda]list)Aplus-Table[dVtrgdwrlist[[i]]*\[Xi]rortminuslist[[i,;;]]+dVtzgdwzlist*\[Xi]zortminuslist[[i,;;]],{i,stepsr}];
	
	dtsortplusd\[Lambda]coeff=Chop[(ExpniTable . dtsortplusd\[Lambda]list . ExpjkTable)/(stepsr*stepsz),10^(-16)];
	dtsortminusd\[Lambda]coeff=Chop[(ExpniTable . dtsortminusd\[Lambda]list . ExpjkTable)/(stepsr*stepsz),10^-(16)];
	
	Print["Calculating Fourier coefficients of \!\(\*SubscriptBox[\(d\[Phi]\), \(s\)]\)/d\[Lambda]"];
	d\[Phi]sortSin\[Psi]d\[Lambda]list=Table[d\[Phi]sortSin\[Psi]d\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	d\[Phi]sortCos\[Psi]d\[Lambda]list=Table[d\[Phi]sortCos\[Psi]d\[Lambda]fun[wrlist[[i]],wzlist,a,p,e,x,constmot],{i,stepsr}];
	
	dV\[Phi]rgdwrlist=dV\[Phi]rgdrgfun[wrlist,a,p,e,x,constmot]*drgdwrfun[wrlist,a,p,e,x,constmot];
	dV\[Phi]zgdwzlist=dV\[Phi]zgdzgfun[wzlist,a,p,e,x,constmot]*dzgdwzfun[wzlist,a,p,e,x,constmot];
	d\[Phi]sortplusd\[Lambda]list=1/2(d\[Phi]sortCos\[Psi]d\[Lambda]list+I*d\[Phi]sortSin\[Psi]d\[Lambda]list)Aminus-Table[dV\[Phi]rgdwrlist[[i]]*\[Xi]rortpluslist[[i,;;]]+dV\[Phi]zgdwzlist*\[Xi]zortpluslist[[i,;;]],{i,stepsr}];
	d\[Phi]sortminusd\[Lambda]list=1/2(d\[Phi]sortCos\[Psi]d\[Lambda]list-I*d\[Phi]sortSin\[Psi]d\[Lambda]list)Aplus-Table[dV\[Phi]rgdwrlist[[i]]*\[Xi]rortminuslist[[i,;;]]+dV\[Phi]zgdwzlist*\[Xi]zortminuslist[[i,;;]],{i,stepsr}];
	  
	d\[Phi]sortplusd\[Lambda]coeff=Chop[(ExpniTable . d\[Phi]sortplusd\[Lambda]list . ExpjkTable)/(stepsr*stepsz),10^-(16)];
	d\[Phi]sortminusd\[Lambda]coeff=Chop[(ExpniTable . d\[Phi]sortminusd\[Lambda]list . ExpjkTable)/(stepsr*stepsz),10^-(16)];
	
	\[CapitalDelta]trg[wr_]:=\[CapitalDelta]trgfun[wr,dtrgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]tzg[wz_]:=\[CapitalDelta]tzgfun[wz,dtzgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	\[CapitalDelta]\[Phi]rg[wr_]:=\[CapitalDelta]\[Phi]rgfun[wr,d\[Phi]rgd\[Lambda]coeff,\[CapitalUpsilon]rg];
	\[CapitalDelta]\[Phi]zg[wz_]:=\[CapitalDelta]\[Phi]zgfun[wz,d\[Phi]zgd\[Lambda]coeff,\[CapitalUpsilon]zg];
	
	\[Delta]rortplus[wr_,wz_]:=\[Delta]rortplusfun[wr,wz,a,p,e,x,constmot,\[Xi]rortpluscoeff];
	\[Delta]rortminus[wr_,wz_]:=\[Delta]rortminusfun[wr,wz,a,p,e,x,constmot,\[Xi]rortminuscoeff];
	\[Delta]rort[wr_,wz_,wp_]:=Re[\[Delta]rortplus[wr,wz]Exp[-I*wp]+\[Delta]rortminus[wr,wz]Exp[I*wp]];
	
	\[Delta]zortplus[wr_,wz_]:=\[Delta]zortplusfun[wr,wz,a,p,e,x,constmot,\[Xi]zortpluscoeff];
	\[Delta]zortminus[wr_,wz_]:=\[Delta]zortminusfun[wr,wz,a,p,e,x,constmot,\[Xi]zortminuscoeff];
	\[Delta]zort[wr_,wz_,wp_]:=Re[\[Delta]zortplus[wr,wz]Exp[-I*wp]+\[Delta]zortminus[wr,wz]Exp[I*wp]];
	
	drsortplusd\[Lambda][wr_,wz_]:=drsortplusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rortpluscoeff];
	drsortminusd\[Lambda][wr_,wz_]:=drsortminusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rortminuscoeff];
	drsortd\[Lambda][wr_,wz_,wp_]:=Re[drsortplusd\[Lambda][wr,wz]Exp[-I*wp]+drsortminusd\[Lambda][wr,wz]Exp[I*wp]];
	
	dzsortplusd\[Lambda][wr_,wz_]:=dzsortplusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]zortpluscoeff];
	dzsortminusd\[Lambda][wr_,wz_]:=dzsortminusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]zortminuscoeff];
	dzsortd\[Lambda][wr_,wz_,wp_]:=Re[dzsortplusd\[Lambda][wr,wz]Exp[-I*wp]+dzsortminusd\[Lambda][wr,wz]Exp[I*wp]];
	
	dtsortplusd\[Lambda][wr_,wz_]:=dtsortplusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rortpluscoeff,\[Xi]zortpluscoeff];
	dtsortminusd\[Lambda][wr_,wz_]:=dtsortminusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rortminuscoeff,\[Xi]zortminuscoeff];
	dtsortd\[Lambda][wr_,wz_,wp_]:=Re[dtsortplusd\[Lambda][wr,wz]Exp[-I*wp]+dtsortminusd\[Lambda][wr,wz]Exp[I*wp]];
	
	d\[Phi]sortplusd\[Lambda][wr_,wz_]:=d\[Phi]sortplusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rortpluscoeff,\[Xi]zortpluscoeff];
	d\[Phi]sortminusd\[Lambda][wr_,wz_]:=d\[Phi]sortminusd\[Lambda]fun[wr,wz,a,p,e,x,constmot,\[Xi]rortminuscoeff,\[Xi]zortminuscoeff];
	d\[Phi]sortd\[Lambda][wr_,wz_,wp_]:=Re[d\[Phi]sortplusd\[Lambda][wr,wz]Exp[-I*wp]+d\[Phi]sortminusd\[Lambda][wr,wz]Exp[I*wp]];
	
	\[CapitalDelta]tsortplus[wr_,wz_]:=\[CapitalDelta]tsortplusfun[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},dtsortplusd\[Lambda]coeff];
	\[CapitalDelta]tsortminus[wr_,wz_]:=\[CapitalDelta]tsortminusfun[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},dtsortminusd\[Lambda]coeff];
	\[CapitalDelta]tsort[wr_,wz_,wp_]:=Re[Exp[-I*wp]\[CapitalDelta]tsortplus[wr,wz]+Exp[I*wp]\[CapitalDelta]tsortminus[wr,wz]];
	
	\[CapitalDelta]\[Phi]sortplus[wr_,wz_]:=\[CapitalDelta]\[Phi]sortplusfun[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},d\[Phi]sortplusd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]sortminus[wr_,wz_]:=\[CapitalDelta]\[Phi]sortminusfun[wr,wz,{\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]p},d\[Phi]sortminusd\[Lambda]coeff];
	\[CapitalDelta]\[Phi]sort[wr_,wz_,wp_]:=Re[Exp[-I*wp]\[CapitalDelta]\[Phi]sortplus[wr,wz]+Exp[I*wp]\[CapitalDelta]\[Phi]sortminus[wr,wz]];
				
	<|
	 "OrbitalElements"->{a,p,e,x},
	 "Eg"->EEg,
	 "Lzg"->Lzg,
	 "Kg"->KKg,
	 "Es"->\[Delta]EEg,
	 "Jzs"->\[Delta]Lzg,
	 "Ks"->\[Delta]KKg,
	 "MinoFrequenciesGeo"->{\[CapitalUpsilon]tg,\[CapitalUpsilon]rg,\[CapitalUpsilon]zg,\[CapitalUpsilon]\[Phi]g},
	 "BLFrequenciesGeo"->{\[CapitalUpsilon]rg/\[CapitalUpsilon]tg,\[CapitalUpsilon]zg/\[CapitalUpsilon]tg,\[CapitalUpsilon]\[Phi]g/\[CapitalUpsilon]tg},
	 "MinoPrecessionFrequency"->\[CapitalUpsilon]p,
	 "BLPrecessionFrequency"->\[CapitalUpsilon]p/\[CapitalUpsilon]tg,
	 (*Keys purely oscillatory part geodesic coordinate time and azimuthal trajectories*)
	 "\[CapitalDelta]trg"->Function[{wr},\[CapitalDelta]trg[wr]],
	 "\[CapitalDelta]tzg"->Function[{wz},\[CapitalDelta]tzg[wz]],
	 "\[CapitalDelta]\[Phi]rg"->Function[{wr},\[CapitalDelta]\[Phi]rg[wr]],
	 "\[CapitalDelta]\[Phi]zg"->Function[{wz},\[CapitalDelta]\[Phi]zg[wz]],
	 (*Keys corrections trajectories*)
	 "\[CapitalDelta]\[Delta]tort"->Function[{wr,wz,wp},\[CapitalDelta]tsort[wr,wz,wp]],
	 "\[CapitalDelta]\[Delta]tortplus"->Function[{wr,wz},\[CapitalDelta]tsortplus[wr,wz]],
	 "\[CapitalDelta]\[Delta]tortminus"->Function[{wr,wz},\[CapitalDelta]tsortminus[wr,wz]],
	 "\[Delta]rort"->Function[{wr,wz,wp},\[Delta]rort[wr,wz,wp]],
	 "\[Delta]rortplus"->Function[{wr,wz},\[Delta]rortplus[wr,wz]],
	 "\[Delta]rortminus"->Function[{wr,wz},\[Delta]rortminus[wr,wz]],
	 "\[Delta]zort"->Function[{wr,wz,wp},\[Delta]zort[wr,wz,wp]],
	 "\[Delta]zortplus"->Function[{wr,wz},\[Delta]zortplus[wr,wz]],
	 "\[Delta]zortminus"->Function[{wr,wz},\[Delta]zortminus[wr,wz]],
	 "\[CapitalDelta]\[Delta]\[Phi]ort"->Function[{wr,wz,wp},\[CapitalDelta]\[Phi]sort[wr,wz,wp]],
	 "\[CapitalDelta]\[Delta]\[Phi]ortplus"->Function[{wr,wz},\[CapitalDelta]\[Phi]sortplus[wr,wz]],
	 "\[CapitalDelta]\[Delta]\[Phi]ortminus"->Function[{wr,wz},\[CapitalDelta]\[Phi]sortminus[wr,wz]],
	 (*Keys corrections velocities*)
	 "\[Delta]vtort"->Function[{wr,wz,wp},dtsortd\[Lambda][wr,wz,wp]],
	 "\[Delta]vtortplus"->Function[{wr,wz},dtsortplusd\[Lambda][wr,wz]],
	 "\[Delta]vtortminus"->Function[{wr,wz},dtsortminusd\[Lambda][wr,wz]],
	 "\[Delta]vrort"->Function[{wr,wz,wp},drsortd\[Lambda][wr,wz,wp]],
	 "\[Delta]vrortplus"->Function[{wr,wz},drsortplusd\[Lambda][wr,wz]],
	 "\[Delta]vrortminus"->Function[{wr,wz},drsortminusd\[Lambda][wr,wz]],
	 "\[Delta]vzort"->Function[{wr,wz,wp},dzsortd\[Lambda][wr,wz,wp]],
	 "\[Delta]vzortplus"->Function[{wr,wz},dzsortplusd\[Lambda][wr,wz]],
	 "\[Delta]vzortminus"->Function[{wr,wz},dzsortminusd\[Lambda][wr,wz]],
	 "\[Delta]v\[Phi]ort"->Function[{wr,wz,wp},d\[Phi]sortd\[Lambda][wr,wz,wp]],
	 "\[Delta]v\[Phi]ortplus"->Function[{wr,wz},d\[Phi]sortplusd\[Lambda][wr,wz]],
	 "\[Delta]v\[Phi]ortminus"->Function[{wr,wz},d\[Phi]sortminusd\[Lambda][wr,wz]]
	|>
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];

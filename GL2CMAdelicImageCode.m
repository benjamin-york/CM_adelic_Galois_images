///////////////////////////////////

//README: This code assumes the user already has access to the ell-adic-galois-images code from Andrew Sutherland for CM elliptic curves.
//This code is an implementation of Algorithm 8.4 of our paper which takes as input a CM elliptic curve E/Q with j(E) not 0 or 1728
//and returns the Galois image of E at an adelic level of definition

//The function GL2CMAdelicImage() returns the image of E at a level of definition defined as in Algorithm 8.4.
//The optional paratemeter MinimalImageLevel is set to false by defult. If it is set to true, 
//then GL2CMAdelicImage() instead returns the image of E at the minimal adelic level of defintion

//The function GL2CMImageModN() takes as input E and integer N > 1 and returns the image of E modulo N

///////////////////////////////////

//A list of all Simplest CM curves defined over Q
//implimented as a function for easier access
SIMPLESTCMCURVESOVERQ := function(dfvec)
	if [] eq dfvec then
		print("This is not a CM elliptic curve");
    elif [-3, 1] eq dfvec then
        SC := [[0, -432], [0, 16], [0, -48], [0, 1296], [0, -3888], [0, 144]];
    elif [-3, 2] eq dfvec then 
        SC := [[-135, -594], [-15, 22]];
    elif [-3, 3] eq dfvec then
        SC := [[-4320, -109296], [-480, 4048]];
    elif [-4, 1] eq dfvec then
        SC := [[-1, 0], [4, 0], [-4, 0], [1, 0], [-2, 0], [8, 0], [-8, 0], [-2, 0]];
    elif [-4, 2] eq dfvec then
        SC := [[-11, -14], [-11, 14], [-44, -112], [-44, 112]];
    elif [-7, 1] eq dfvec then 
        SC := [[-1715, 33614], [-35, -98]];
    elif [-7, 2] eq dfvec then
        SC := [[-29155, 1915998], [-595,-5586]];
    elif [-8, 1] eq dfvec then 
        SC := [[-1080, -12096], [-270, 1512], [-1080, 12096], [-270, -1512]];
    elif [-11, 1] eq dfvec then 
        SC := [[-1149984, -487018224], [-9504, 365904]];
    elif [-19, 1] eq dfvec then 
        SC := [[-219488, -39617584], [-608, 5776]];
    elif [-43, 1] eq dfvec then 
        SC := [[-25442240, -49394836848], [-13760, 621264]];
    elif [-67, 1] eq dfvec then 
        SC := [[-529342880, -4687634371504], [-117920, 15585808]];
    elif [-163, 1] eq dfvec then 
        SC := [[-924354639680, -342062961763303088], [-34790720, 78984748304]];
    else
        print("This is not the discriminant of a CM curve over Q");
    end if;
    return SC;
end function;



//Given a CM elliptic curve E/Q with CM by O in K,
//Returns the discriminant of K and the conductor of O
CMdf:=function(E)
	yesno,disc:=HasComplexMultiplication(E);
	if yesno eq false then
	  	print "The elliptic curve has no CM";
		return yesno,[];
	else 
		x,y:=SquarefreeFactorization(disc);
		if x eq -1 then
			d:=-4;
			f:=y div 2;
		elif x eq -2 then
			d:=-8;
			f:=1;
		else
			d:=x;
			f:=y;
		end if;
		return yesno,[d,f];
	end if;
end function;

//Given dfvec = [d,f],
//Returns values delta and phi with which matrix groups are defined
deltaphi:=function(dfvec)
	d:=dfvec[1];
	f:=dfvec[2];
	df2:=Integers()!d*f^2;
	if (df2 mod 4) eq 0 then
		delt:=Integers()!(d*f^2/4);
		phi:=0;
	else
		delt:=Integers()!((d-1)*f^2/4);
		phi:=f;
	end if;
	return [delt,phi];
end function;


//Given an elliptic curve E over Q
//Returns a short Weierstrass model y^2 = x^3 + Ax + B for E
//such that if d^4 | A and d^6 | B, then |d| = 1
ShortWeierstrassForm := function(E)
	E := WeierstrassModel(E);
	Coeffs := Coefficients(E);
    A := Integers()!Coeffs[4];
	B := Integers()!Coeffs[5];
    if A eq 0 then
        PrimesDividing := PrimeDivisors(B);
        for p in PrimesDividing do
            InfiniteLoopPreventer := 1;
            while (B/p^6) in Integers() and InfiniteLoopPreventer le 10^3 do
				B := Integers()!(B/p^6);
                InfiniteLoopPreventer +:= 1;
            end while;
        end for;
    elif B eq 0 then
        PrimesDividing := PrimeDivisors(A);
        for p in PrimesDividing do
            InfiniteLoopPreventer := 1;
            while (A/p^4) in Integers() and InfiniteLoopPreventer le 10^3 do
				A := Integers()!(A/p^4);
                InfiniteLoopPreventer +:= 1;
            end while;
        end for;
    else
        PrimesDividing := PrimeDivisors(A*B);
        for p in PrimesDividing do
            InfiniteLoopPreventer := 1;
            while (A/p^4) in Integers() and (B/p^6) in Integers() and InfiniteLoopPreventer le 10^3 do
				A := Integers()!(A/p^4);
                B := Integers()!(B/p^6);
                InfiniteLoopPreventer +:= 1;
            end while;
        end for;
	end if;
	return EllipticCurve([A,B]);
end function;

//Given dfvec = [d,f] and matrix mat,
//Returns whether mat is a Cartan matrix
IsCartanMatrix:=function(dfvec,mat)
	dphi:=deltaphi(dfvec);
	delta:=dphi[1];
	phi:=dphi[2];
	A:=mat[1][1];
	B:=mat[1][2];
	C:=mat[2][1];
	D:=mat[2][2];
	yesno:=false;
	if (C eq B*delta) and (A-D eq B*phi) then
		yesno:=true;
	end if;
	return yesno;
end function;

//Given dfvec = [d,f] and list of matrices matvec with coprime characteristics
//Returns matrix lift of the matrices in matvac using CRT
CRTLift:=function(dfvec,matvec)
	dphi:=deltaphi(dfvec);
	delta:=dphi[1];
	phi:=dphi[2];
	N:=1;
	plist:=[];
	coe11:=[];
	coe12:=[];
	coe21:=[];
	coe22:=[];
	yesno:=1;
	for i:=1 to #matvec do
		if IsCartanMatrix(dfvec,matvec[i]) then
			yesno:=yesno*1;
		else
			yesno:=yesno*0;
		end if;
	end for;
	for i:=1 to #matvec do
		m:=matvec[i];
		char:=Characteristic(CoefficientRing(m));
		plist:=Append(plist,char);
		N:=N*char;
		coe11:=Append(coe11,Integers()!(m[1][1]));
		coe12:=Append(coe12,Integers()!(m[1][2]));
		coe21:=Append(coe21,Integers()!(m[2][1]));
		coe22:=Append(coe22,Integers()!(m[2][2]));
	end for;
	GL:=GL(2,Integers(N));
	if yesno eq 1 then
		a:=CRT(coe22,plist);
		b:=CRT(coe12,plist);
 		M:=GL![a+b*phi,b,delta*b,a];
	else
		a:=CRT(coe22,plist);
		b:=CRT(coe12,plist);
 		M:=GL![-a,b,-delta*b+a*phi,a];
	end if;
	return M;
end function;

//Given dfvec = [d,f] and a matrix group G, a subgroup of the normalizer of Cartan
//Return G so that the first generator is contained in the normalizer, if such exists
//and all other generators are contained in the Cartan
function CleanCartanSubgroup(dfvec,G)
	//G:=CleanGroup(G);
	gens:=[g: g in Generators(G)];
	newgens:=[];
	done:=0;
	tau:=G![1,0,0,1];
	for i:=1 to #gens do
		yesnoi:=IsCartanMatrix(dfvec,gens[i]);
		if yesnoi eq false then
			newgens:=Append(newgens,gens[i]*tau);
		else
			newgens:=Append(newgens,gens[i]);
		end if;
		if (done eq 0) and (yesnoi eq false) then
			tau:=gens[i];
			done:=1;
		end if;
	end for;
	newgens2:=[tau];
	for g in newgens do
		if g ne tau then
			newgens2:=Append(newgens2,g);
		end if;
	end for;
	H:=sub<G|newgens2>;
	return H;
end function;


//Given dfvec = [d,f], and matrix groups H and K,
//Returns a lift of HK mod Char(H)*Char(K)
//Char(H) and Char(K) are assumed to be relatively prime.
CRTSubgroupLift := function(dfvec, H, K)
	H := CleanCartanSubgroup(dfvec,H);
	K := CleanCartanSubgroup(dfvec,K);
	RH:=CoefficientRing(H);
	RK:=CoefficientRing(K);
	charH:=Characteristic(RH);
	charK:=Characteristic(RK);
	GL2:=GL(2,Integers(charH*charK));
	gensH:=[g: g in Generators(H)];
	gensK:=[g: g in Generators(K)];
	IdH:=H![1,0,0,1];
	IdK:=K![1,0,0,1];
	if #gensH eq 0 then
		gensH := [IdH];
	end if;
	if #gensK eq 0 then
		gensK := [IdK];
	end if;
	gensHK:=[];
	if IsCartanMatrix(dfvec,gensH[1]) or IsCartanMatrix(dfvec,gensK[1]) then
		gensHK:=Append(gensHK,CRTLift(dfvec,<gensH[1],IdK>));
		gensHK:=Append(gensHK,CRTLift(dfvec,<IdH,gensK[1]>));
	else
		gensHK:=Append(gensHK,CRTLift(dfvec,<gensH[1],gensK[1]>));	
	end if;
	for i in [2..#gensH] do
		gensHK:=Append(gensHK,CRTLift(dfvec,<gensH[i],IdK>));
	end for;
	for j in [2..#gensK] do
		gensHK:=Append(gensHK,CRTLift(dfvec,<IdH,gensK[j]>));
	end for;
	HK:=sub<GL2|gensHK>;
	return HK; 
end function;


//Takes as input a discriminant dfvec and Normalizer subgroup N
//returns Cartan subgroup of N
GL2CartanFromNormalizer := function(dfvec,N)
	gens := [g : g in Generators(N)];
	cartangens:=[];
	for i:=1 to #gens do
		if IsCartanMatrix(dfvec,gens[i]) then
			cartangens:=Append(cartangens,gens[i]);
		end if;
	end for;
	C:=sub<N|cartangens>;
	return C;
end function;

//Takes as input discriminant dfvec and positive integer n
//returns the Cartan subgroup and Normalizer of Cartan at level n
//GL2CartanNormalizer is a function created by Drew Sutherland
CartanSubgroup := function(dfvec,n)
	Nn:=GL2CartanNormalizer(dfvec[1]*dfvec[2]^2,n);
	Cn:=GL2CartanFromNormalizer(dfvec,Nn);
	return Cn, Nn;
end function;

//////////////////////////////////////////////////////////////////////////////
//Explicit representations of matrices from classification of ell-adic images (see Section 2 of paper)
C1 := function(dfvec)
    dphi := deltaphi(dfvec);
    delta := dphi[1];
    phi := dphi[2];
    return [1, 0, -phi, -1];
end function;

Cm1 := function(dfvec)
    dphi := deltaphi(dfvec);
    delta := dphi[1];
    phi := dphi[2];
    return [-1, 0, phi, 1];
end function;

G21D0List := function(dfvec)
    dphi := deltaphi(dfvec);
    delta := dphi[1];
    phi := dphi[2];
    return [[5,0,0,5], [1,1,delta,1]];
end function;

G22D0List := function(dfvec)
    dphi := deltaphi(dfvec);
    delta := dphi[1];
    phi := dphi[2];
    return [[5,0,0,5], [-1,-1,-1*delta,-1]];
end function;

G23D0List := function(dfvec)
    dphi := deltaphi(dfvec);
    delta := dphi[1];
    phi := dphi[2];
    return [ [3,0,0,3], [1,1,delta,1] ];
end function;

G24D0List := function(dfvec)
    dphi := deltaphi(dfvec);
    delta := dphi[1];
    phi := dphi[2];
    return [ [3,0,0,3], [-1,-1,-1*delta,-1] ];
end function;

G21DPList := function(dfvec, ell)
    dphi := deltaphi(dfvec);
    delta := dphi[1];
    phi := dphi[2];
    return [[Integers(ell)!(a^2 + b*phi/2),Integers(ell)!b,Integers(ell)!(delta*b),Integers(ell)!(a^2 - b*phi/2)] : a, b in [0..(ell-1)] | (a mod ell) ne 0];
end function;
//////////////////////////////////////////////////////////////////////////////


//Takes as input a Cartan subgroup Cell and a subgroup CG of Cell
//and returns the image of CG with a reduced number of generators
//without intervention, CG may sometimes have thousands of generators
CGR := function(Cell, CG)
	A_ell, mappytoA_ell := AbelianGroup(Cell);
    AG := mappytoA_ell(CG);
    return AG @@ mappytoA_ell; //this simplifies generators for CG
end function;

//Given a CM elliptic curve E/Q with j not 0, 1728 and prime ell
//Returns the ell-adic image of E at level ell if ell ne 2 or 16 if ell eq 2
//Also returns intersection of ell-adic image with Cartan, and image of complex conjugation
//This utilizes the classification of models realizing all possible ell-adic images 
//for curves over Q previously proved by the authors and Gonzalez-Jimenez
//Significant time savings achieved for discriminants -43, -67, and -163
GL2CMEllAdicImageFast := function(E, ell)
    yesno, dfvec := CMdf(E);
    Coeffs := SIMPLESTCMCURVESOVERQ(dfvec);
    if ell eq 2 then
        Cell, Nell := CartanSubgroup(dfvec, ell^4);
        if dfvec eq [-4, 2] then
            if IsIsomorphic(E,EllipticCurve(Coeffs[1])) then
                CG := sub<Cell | G22D0List(dfvec)>;
                ccmat := Cm1(dfvec);
                G := sub<Nell | CG, ccmat>;
            elif IsIsomorphic(E,EllipticCurve(Coeffs[2])) then
                CG := sub<Cell | G22D0List(dfvec)>;
                ccmat := C1(dfvec);
                G := sub<Nell | CG, ccmat>;
            elif IsIsomorphic(E,EllipticCurve(Coeffs[3])) then
                CG := sub<Cell | G21D0List(dfvec)>;
                ccmat := Cm1(dfvec);
                G := sub<Nell | CG, ccmat>;
            elif IsIsomorphic(E,EllipticCurve(Coeffs[4])) then
                CG := sub<Cell | G21D0List(dfvec)>;
                ccmat := C1(dfvec);
                G := sub<Nell | CG, ccmat>;
            else
                G := Nell;
                CG := Cell;
                ccmat := C1(dfvec);
            end if;
        elif dfvec eq [-8, 1] then
            if IsIsomorphic(E,EllipticCurve(Coeffs[1])) then
                CG := sub<Cell | G24D0List(dfvec)>;
                ccmat := Cm1(dfvec);
                G := sub<Nell | CG, ccmat>;
            elif IsIsomorphic(E,EllipticCurve(Coeffs[2])) then
                CG := sub<Cell | G24D0List(dfvec)>;
                ccmat := C1(dfvec);
                G := sub<Nell | CG, ccmat>;
            elif IsIsomorphic(E,EllipticCurve(Coeffs[3])) then
                CG := sub<Cell | G23D0List(dfvec)>;
                ccmat := C1(dfvec);
                G := sub<Nell | CG, ccmat>;
            elif IsIsomorphic(E,EllipticCurve(Coeffs[4])) then
                CG := sub<Cell | G23D0List(dfvec)>;
                ccmat := Cm1(dfvec);
                G := sub<Nell | CG, ccmat>;
            else
                G := Nell;
                CG := Cell;
                ccmat := C1(dfvec);
            end if;
        elif dfvec eq [-4, 1] then
            print("The case of j = 1728 is not yet supported");
        else
            G := Nell;
            CG := Cell;
            ccmat := C1(dfvec);
        end if;
    else
        Cell, Nell := CartanSubgroup(dfvec, ell);
        if ell eq 3 then
            if dfvec eq [-3, 2] or dfvec eq [-3, 3] then
                if IsIsomorphic(E,EllipticCurve(Coeffs[1])) then
                    CG := sub<Cell | G21DPList(dfvec, ell)>;
					CG := CGR(Cell, CG);
                    ccmat := Cm1(dfvec);
                    G := sub<Nell | CG, ccmat>;
                elif IsIsomorphic(E,EllipticCurve(Coeffs[2])) then
                    CG := sub<Cell | G21DPList(dfvec, ell)>;
					CG := CGR(Cell, CG);
                    ccmat := C1(dfvec);
                    G := sub<Nell | CG, ccmat>;
                else
                    G := Nell;
                    CG := Cell;
                    ccmat := C1(dfvec);
                end if;
            elif dfvec eq [-3,1] then
                print("The case of j = 0 is not yet supported");
            else
                G := Nell;
                CG := Cell;
                ccmat := C1(dfvec);
            end if;

        elif ell eq 7 then
            if dfvec eq [-7, 1] or dfvec eq [-7, 2] then
                if IsIsomorphic(E,EllipticCurve(Coeffs[1])) then
                    CG := sub<Cell | G21DPList(dfvec, ell)>;
					CG := CGR(Cell, CG);
                    ccmat := C1(dfvec);
                    G := sub<Nell | CG, ccmat>;
                elif IsIsomorphic(E,EllipticCurve(Coeffs[2])) then
                    CG := sub<Cell | G21DPList(dfvec, ell)>;
					CG := CGR(Cell, CG);
                    ccmat := Cm1(dfvec);
                    G := sub<Nell | CG, ccmat>;
                else
                    G := Nell;
                    CG := Cell;
                    ccmat := C1(dfvec);
                end if;
            else
                G := Nell; //ell adic image always maximal all other cases
                CG := Cell;
                ccmat := C1(dfvec);
            end if;

        elif ell in [11, 19, 43, 67, 163] then
            if dfvec in [[-11, 1], [-19, 1], [-43, 1], [-67, 1], [-163, 1]] then
                if IsIsomorphic(E,EllipticCurve(Coeffs[1])) then
                    CG := sub<Cell | G21DPList(dfvec, ell)>;
					CG := CGR(Cell, CG);
                    ccmat := Cm1(dfvec);
                    G := sub<Nell | CG, ccmat>;
                elif IsIsomorphic(E,EllipticCurve(Coeffs[2])) then
                    CG := sub<Cell | G21DPList(dfvec, ell)>;
					CG := CGR(Cell, CG);
                    ccmat := C1(dfvec);
                    G := sub<Nell | CG, ccmat>;
                else
                    G := Nell;
                    CG := Cell;
                    ccmat := C1(dfvec);
                end if;
            else
                G := Nell; //ell adic image always maximal in this case
                CG := Cell;
                ccmat := C1(dfvec);
            end if;
        else
            G := Nell;
            CG := Cell;
            ccmat := C1(dfvec);
        end if;
    end if;
    return G, CG, ccmat;
end function;



//Given a square-free integer N,
//Return N^dagger, the absolute value of the discriminant of the order of Q(sqrt(N))
Ndagger := function(N)
	if N mod 4 eq 1 then
		Nd := AbsoluteValue(N);
	else
		Nd := 4*AbsoluteValue(N);
	end if;
	return Nd;
end function;

//Given integer d and positive integer N,
//Returns powers of zeta_N for which sqrt(d) is fixed
ExponentsFixingSqrtd:=function(d,N)
	C<z>:=CyclotomicField(N);
	P<x>:=PolynomialRing(C);
	list:=[];
	roo:=Roots(x^2-d);
	if #roo eq 0 then
		print "Sqrt(",d,") is not present in Q(zeta_",N,")";
	else
		roo1:=roo[1][1];
		A:=Automorphisms(C);
		exps:=[j : j in [1..N] | GCD(j,N) eq 1];
		for phi in A do
			expophi:=1;
			for j in exps do
				if phi(z) eq z^j then
					expophi:=j;
				end if;
			end for;
			if phi(roo1) eq roo1 then
				list:=Append(list,expophi);
			end if;
		end for;
	end if;
	return Sort(list);
end function;

//Given a discriminant dfvec, integer d, and positive integer N,
//Returns Cartan subgroup mod N that fixes sqrt(d)
CartanSubgroupModNThatFixes := function(dfvec,d,N)
	C:=CartanSubgroup(dfvec,N);
	dets := ExponentsFixingSqrtd(d,N);
	UnitGrp, mappy := UnitGroup(Integers(N));
	phi := hom<C -> UnitGrp | [<g, Inverse(mappy)(Determinant(g))> : g in GeneratorsSequence(C)]>;
	UnitSub := sub<UnitGrp | [Inverse(mappy)(m) : m in dets]>;
	return UnitSub @@ phi;
end function;


//Given a CM elliptic curve E, 
//Returns the prime ell for which E has maximal ell-adic index
//This function MaximalPrimeEll is specific to curves defined over Q!!
MaximalPrimeEll := function(E)
	bool, dfvec := CMdf(E);
	primes:=PrimeDivisors(dfvec[1]);
	ell:=primes[1];
	return ell;	
end function;

//Given CM elliptic curve E with j not 0,1728
//Returns a simplest CM curve E', and the twist sending E to E'
SimplestCMCurveV2 := function(E)
	yesno, dfvec := CMdf(E);
	SimplestCurves := SIMPLESTCMCURVESOVERQ(dfvec);
	ell := MaximalPrimeEll(E);
	E := ShortWeierstrassForm(E);
	Coeffs := Coefficients(E);
	A := Integers()!Coeffs[4];
	B := Integers()!Coeffs[5];

	if A*B eq 0 then
		print("We do not currently support curves with j-invariant 0 or 1728");
	else
		TwistsToTest := #SimplestCurves;
		i := 1;
		twistfinding := true;

		while i le TwistsToTest and twistfinding do
			Ai := Integers()!SimplestCurves[i][1];
			Bi := Integers()!SimplestCurves[i][2];
			TwistFactorInitial := (B/Bi)/(A/Ai);
			if TwistFactorInitial in Integers() then
				TwistFactor := Integers()!TwistFactorInitial;
			else
				TwistFactor := Integers()!(TwistFactorInitial * Denominator(TwistFactorInitial)^2);
			end if;

			if GCD(ell, Ndagger(TwistFactor)) eq 1 then
				Enew := EllipticCurve([Ai, Bi]);
				twistfinding := false;
			else
				i +:= 1;
			end if;
		end while;
		testbool := IsIsomorphic(E, QuadraticTwist(Enew, TwistFactor));
		if not testbool then
			print("Something has gone wrong with SimplestCMCurveV2");
		end if;
	end if;
	return MinimalModel(Enew), TwistFactor;
end function;



//Given CM elliptic curve E, discriminant dfvec, prime ell, and square-free integer N,
//Returns the Cartan image of E at a level for which the adelic image is defined
//Assumes GCD(ell,N) == 1
QuadraticTwistCartanImage := function(E,dfvec,ell,N)
	M := Ndagger(N);
	if ell eq 2 then
		C_ellM := CartanSubgroup(dfvec,16*M);
		C_ell := CartanSubgroup(dfvec,16);
		level:=16*M;
	else
		C_ellM := CartanSubgroup(dfvec,ell*M);
		C_ell := CartanSubgroup(dfvec,ell);
		level:=ell*M;
	end if;
	C_M := CartanSubgroup(dfvec,M);

	Idell := sub<C_ell|C_ell![1,0,0,1]>;
	IdM := sub<C_M|C_M![1,0,0,1]>;
	G_ell := CRTSubgroupLift(dfvec,C_ell,IdM);
	G_M := CRTSubgroupLift(dfvec,Idell,C_M);
	if ell eq 2 then 
		N_ell1, C_ell1 := GL2CMEllAdicImageFast(E,ell); //returns image at level 16
	else
		N_ell1, C_ell1 := GL2CMEllAdicImageFast(E,ell);
        C_ell1 := ChangeRing(C_ell1, Integers(ell)); //returns image at level ell
	end if;
	
	C_M1 := CartanSubgroupModNThatFixes(dfvec,N,M);
	G_ell1xG_M := CRTSubgroupLift(dfvec,C_ell1,C_M);
	G_ellxG_M1 := CRTSubgroupLift(dfvec,C_ell,C_M1);
	G_ell1xG_M1 := CRTSubgroupLift(dfvec,C_ell1,C_M1);

	G:=C_ellM;
	AG,mappytoAG:=AbelianGroup(G);

    H2 := G_ell1xG_M1;
	AH2 := mappytoAG(H2);

    G_ell1xG_M := mappytoAG(G_ell1xG_M) @@ mappytoAG; //cleans up generators
    G_ellxG_M1 := mappytoAG(G_ellxG_M1) @@ mappytoAG; //cleans up generators

	Q,mapptoQ:= AG / AH2;
	gens := [g : g in Q];

	CandList := [];

    i := 1;
    GensToTest := #gens;
    loopbool := true;

	while i le GensToTest and loopbool do
        h := gens[i];
        sub_AG := sub<AG | AH2, Inverse(mapptoQ)(h)>;
        sub := sub_AG @@ mappytoAG;
		if (sub meet G_ell1xG_M) eq (sub meet G_ellxG_M1) then
			if #sub eq (#G div 2) then
				subell:=ChangeRing(sub,CoefficientRing(C_ell));
				subM:=ChangeRing(sub,CoefficientRing(C_M));
				if (subell eq C_ell) and (subM eq C_M) then
					CandList := Append(CandList,sub);
                    loopbool := false;
				end if;
			end if;
		end if;
        i +:= 1;
	end while;
	if #CandList ne 1 then
		print("There is not one and only one candidate subgroup");
	end if;
	HG := CandList[1];
	return HG, level;
end function;


//Given curve E with Galois image G at level N, an even integer
//Return Galois image of E at level of pullback, 
//i.e. the smallest level for which the image is index 2 in the normalizer
//This is only necessary when the level is even for j not 0
MinimalAdelicImage := function(E, G)
	bool, dfvec := CMdf(E);
	testbool := true;
    N := Characteristic(CoefficientRing(G));
	while testbool do
		N := N/2;
		if N notin Integers() then
			testbool := false;
			N := Integers()!(N*2);	
		else 
			N := Integers()!N;
			Cc, Nc := CartanSubgroup(dfvec,N);
			Gtemp := ChangeRing(G,Integers(N));
			if #Gtemp eq #Nc then
				testbool := false;
				N := Integers()!(N*2);
			else
				G := Gtemp;
			end if;
		end if;
	end while;
	return G;
end function;


//Given a CM elliptic curve E, j not 0, 1728, 
//Returns Galois image of G at level M, 
//where M is an adelic level of definition
GL2CMAdelicImage := function(E : MinimalImageLevel := false)
	E_simple, N := SimplestCMCurveV2(E); 
	//The Simplest CM curve for E, and twist factor N
	ell := MaximalPrimeEll(E_simple);
	//Prime where E_simple has maximal index
	NG, CG, CpxConj := GL2CMEllAdicImageFast(E_simple,ell); 
	if N eq 1 then
		G := NG;
	else
        TwistSign := Sign(N);
		for i in [1..4] do
			CpxConj[i] := TwistSign * CpxConj[i];
		end for;
		bool, dfvec := CMdf(E);
		Cartan, Level := QuadraticTwistCartanImage(E_simple,dfvec,ell,N);
		GL2 := GL(2,Integers(Level));
		G := sub<GL2 | Cartan, CpxConj>;
	end if;
    if MinimalImageLevel then
        G := MinimalAdelicImage(E,G);
    end if;
	return G;
end function;


//Given normalizer of cartan H and Galois image K,
//Return K lifted to be a subgroup of H
NormalizerProjectionLift := function(H,K)
    R := BaseRing(K);
    parentK := GL(2,R);
    phi := hom<H -> parentK | [ChangeRing(g,R) : g in GeneratorsSequence(H)]>;
    return K @@ phi;
end function;



//Given CM elliptic curve E, j not 0, 1728, and level N,
//Return the Galois image of E mod N
GL2CMImageModN := function(E,N)
	ell := MaximalPrimeEll(E);
	Es, twistfact := SimplestCMCurveV2(E);

	test_level := ell*Ndagger(twistfact);
	bool, dfvec := CMdf(E);

	//Easy test to see if image is full normalizer
	if GCD(N,test_level) ne 1 then

		Gtemp := GL2CMAdelicImage(E : MinimalImageLevel := true);

		char := Characteristic(CoefficientRing(Gtemp));
		
		//if N divides char, project G to level N
		if (char mod N) eq 0 then
			G := ChangeRing(Gtemp,Integers(N));

		//if char divides N, lift G to level N
		elif (N mod char) eq 0 and char ne N then
			Cn, Nn := CartanSubgroup(dfvec, N);
			G := NormalizerProjectionLift(Nn, Gtemp);

		//else, image is full normalizer
		else
			C, G := CartanSubgroup(dfvec,N);

		end if;
	else
		C,G := CartanSubgroup(dfvec,N);
	end if;
	return G;
end function;

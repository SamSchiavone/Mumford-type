load "Shimura/Shimura_models.m"; 

S<t, Y, T> := PolynomialRing(QQ, 3);

Q, E := GenericShimuraCurve239();
f := Evaluate(E, [t, Y^2, Y, 1, T]);


// fiber at 0
f_0 := Evaluate(f, [t^2, t*Y, t^3*(-T-1/3*Y)-t^2]) div t^8;
Factorization(Evaluate(f1, [0, Y, T]));
// we get (Y^3-T)^2

f_0 := -1/3*Evaluate(S!(1/t^2*Evaluate(f_0, [t^2, Y, Y^3])), [0, Y, 0]);
// y^2 = x*(x^8 - 4*x^6 + 6*x^4 - 44/27*x^2 + 1)
//Does it have CM by Q(i)+Q(zeta_9+zeta_9^(-1), i) ?


// fiber at 1
f_1 := Evaluate(f, [t^3+1, Y-1, t*T]) div t^3;
ell := Evaluate(f_1, [0,Y,T]);
Factorization((ell/3-T^3) div Y^3);
//Get an elliptic curve isomorphic to y^2=x^3-1 as a component
f_1 := Evaluate(f_1, [t^3, t^3*Y, t^4*T]) div t^12;
pic := Evaluate(f_1, [0, Y, T]);
pic;
//Get the Picard curve y^3=x^4+x as another component. It has CM by Z[zeta_9]


// fiber at infinity
f1 := Evaluate(f, [t^3, t*Y, t^5*T]);
Terms(f1, t)[16] div t^15;
// -2Y^6 + 15Y^4T - 14Y^3 + 6YT + 3T^3 - 3

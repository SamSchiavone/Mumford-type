load "Shimura/Shimura_models.m"; 

f := GenericShimuraCurve237();
S<t, x> := Parent(f);
S0<x0> := PolynomialRing(QQ);

// fiber at 0
f_0 := -64/567*Evaluate(Evaluate(f, [t^2, t*x]) div t^11, [0, x0]);
// y^2 = x0^9 + 16/3*x0^7 + 32/3*x0^5 - 256/21*x0^3 + 256/81*x0


// fiber at 1
ell := &*[l[1] : l in Factorization(Evaluate(f, [1, x]))];
ell := 11/7*Evaluate(S0!(x0^4*Evaluate(ell, [1, 1/x0])), x0+1);
// y^2 = x^3 - 45/28*x + 27/28

f_1 := Evaluate(f, [t+1, x]);
f_1 := S!((x-1)^10*Evaluate(f_1, [t, 2/(x-1)]));
f_1 := -1/1152*Evaluate(Evaluate(f_1, [t^7, t^2*x]) div t^14, [0, x0]);
// y^2 = x^7 - 3


// fiber at infinity
f_inf := Evaluate(f, [t^3, t*x]);
Coeffs := Terms(f_inf, t);
f_inf := Evaluate(Coeffs[#Coeffs] div t^16, [0, x0]);
// y^2 = x^10 - 84*x^7 + 84*x^4 - 28*x

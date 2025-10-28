///////////////////////////////////////////////////
///////////// Discriminant of C_{7,t} /////////////
///////////////////////////////////////////////////

R<t> := PolynomialRing(Rationals());
S<x> := PolynomialRing(FieldOfFractions(R), 2);
f := t*((t - 27/16)*x^10 - 567/64*x^9 - 189/4*t*x^8 + (-84*t^2 - 189/4*t)*x^7 - 189*t^2*x^6 - 189/2*t^2*x^5 + 84*t^3*x^4 + 108*t^3*x^3 - 28*t^4*x);

disc := R!Discriminant(f, 1);
factors := Factorization(disc);
factors;
Factorization(Integers()!(disc/&*[fact[1]^fact[2] : fact in factors]));

///////////////////////////////////////////////////
///////////// Discriminant of C_{9,t} /////////////
///////////////////////////////////////////////////

R<t> := PolynomialRing(Rationals());
S<X, Y, Z, T> := PolynomialRing(FieldOfFractions(R), 4);
Q := X*Z - Y^2;
E := 3*T^3 + t*(t - 1)*(3*T*(5*X^2 + 6*X*Y + 2*t*Y*Z + 3*t*Z^2) + (-2*t + 9)*X^3 + 22*t*X^2*Y + 21*t*X^2*Z + (-14*t^2 + 18*t)*X*Y*Z + t^2*X*Z^2 + 6*t^2*Y*Z^2 + (-3*t^3 + 6*t^2)*Z^3);
  
invs, wgt, norm := InvariantsGenus4Curves(Q,E); // norm is the normalization factor due to the determinants of the transformations
disc := R!DiscriminantFromInvariantsGenus4(invs : norm := norm); // this is the exact discriminant, 
factors := Factorization(disc);
factors;
Factorization(Integers()!(disc/&*[fact[1]^fact[2] : fact in factors]));

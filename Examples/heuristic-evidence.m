// Heuristic evidence for fibers being of Mumford-type

AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/Genus-4/magma/spec");
load "~/github/Mumford-type/Shimura/Shimura_models.m";
load "~/github/monodromy/monodromy.m";

// hyperelliptic 2,3,7 family
f := GetShimuraCurve237(QQ!2);
C := HyperellipticCurve(f);
polys := [];
  for p in PrimesInInterval(11,100) do
  Append(~polys, <p, LPolynomial(Curve(Reduction(C,p)))>);
end for;

// Noot criterion: Abelian varieties with ell-adic...
// Proposition 3.6
slopes := [];
for pair in polys do
  p, poly := Explode(pair);
  Append(~slopes, <p^2 mod 7, LowerSlopes(NewtonPolygon(poly,p))>);
end for;

// computing Heuristic monodromy group using Zywina's package
// https://github.com/davidzywina/monodromy
S := [el[2] : el in polys];
S[1..5];
flag,T:=ComputeRootDatumOfMT(S);
flag;
T;

// non-hyperelliptic 2,3,9 family
Q, E := GetShimuraCurve239(QQ!2);
C := Curve(Proj(Parent(Q)), [Q,E]);

polys := [];
  for p in PrimesInInterval(11,100) do
  Append(~polys, <p, LPolynomial(Curve(Reduction(C,p)))>);
end for;

// Noot criterion: Abelian varieties with ell-adic...
// Proposition 3.6
slopes := [];
for pair in polys do
  p, poly := Explode(pair);
  Append(~slopes, <p^2 mod 9, LowerSlopes(NewtonPolygon(poly,p))>);
end for;

// computing Heuristic monodromy group using Zywina's package
// https://github.com/davidzywina/monodromy
S := [el[2] : el in polys];
flag,T:=ComputeRootDatumOfMT(S);
flag;
T;

// at singular fibers
// 2,3,7 family

R<x> := PolynomialRing(QQ);
// fiber above 0 (see Semistable/237.m)
f := x^9 + 16/3*x^7 + 32/3*x^5 - 256/21*x^3 + 256/81*x;
C := HyperellipticCurve(f);

polys := [];
  for p in PrimesInInterval(11,100) do
  Append(~polys, <p, LPolynomial(Curve(Reduction(C,p)))>);
end for;

slopes := [];
for pair in polys do
  p, poly := Explode(pair);
  Append(~slopes, <p^2 mod 7, LowerSlopes(NewtonPolygon(poly,p))>);
end for;
slopes;

S := [el[2] : el in polys];
flag,T:=ComputeRootDatumOfMT(S);
flag;
T;


// 2,3,9 family
// fiber above 0 (see Semistable/239.m)
f := x*(x^8 - 4*x^6 + 6*x^4 - 44/27*x^2 + 1);
C := HyperellipticCurve(f);

polys := [];
  for p in PrimesInInterval(11,100) do
  Append(~polys, <p, LPolynomial(Curve(Reduction(C,p)))>);
end for;

slopes := [];
for pair in polys do
  p, poly := Explode(pair);
  Append(~slopes, <p^2 mod 7, LowerSlopes(NewtonPolygon(poly,p))>);
end for;
slopes;

S := [el[2] : el in polys];
flag,T:=ComputeRootDatumOfMT(S);
flag;
T;

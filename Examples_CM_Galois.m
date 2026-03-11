load "GL2CMAdelicImageCode.m";

//Computational Examples

E := EllipticCurve([1, -1, 0, -107, 552]); //7-simplest CM curve with disc -7
time GL2CMAdelicImage(E); //the adelic image of E
/*
MatrixGroup(2, IntegerRing(7))
Generators:
    [4 4]
    [6 0]

    [1 0]
    [6 6]
Time: 0.000
*/


E := EllipticCurve([0, 0, 1, -57772164980, -5344733777551611]); //163-simplest CM curve with disc -163
time GL2CMAdelicImage(E);
/*
MatrixGroup(2, IntegerRing(163))
Generators:
    [ 20  83]
    [ 20 100]

    [162   0]
    [  1   1]
Time: 2.172
*/


E := EllipticCurve([1, -1, 1, -965, -13940]); //CM curve with disc -7
time GL2CMAdelicImage(E);
/*
MatrixGroup(2, IntegerRing(21))
Generators:
    [ 4 19]
    [ 4  6]

    [20  0]
    [ 1  1]
Time: 0.031
*/


E := EllipticCurve([0, 0, 0, -1100, -14000]); //CM curve with disc -16
time GL2CMAdelicImage(E);
/*
MatrixGroup(2, IntegerRing(80))
Generators:
    [27 16]
    [16 27]

    [57 25]
    [60 57]

    [ 5 56]
    [16  5]

    [31 48]
    [48 31]

    [79  0]
    [ 0  1]
Time: 0.172
*/
time GL2CMAdelicImage(E : MinimalImageLevel := true); //returns adelic image at minimal level of definition
/*
MatrixGroup(2, IntegerRing(20)) of order 2^7
Generators:
    [7 4]
    [4 7]

    [ 3 16]
    [16  3]

    [17  0]
    [ 0 17]

    [ 5  7]
    [12  5]

    [19  0]
    [ 0  1]
Time: 0.188
*/

E := EllipticCurve([0, 0, 0, -29155, -1915998]); //CM curve with disc -28 and quadratic entanglement at level 14
time GL2CMAdelicImage(E);
/*
MatrixGroup(2, IntegerRing(28))
Generators:
    [20 21]
    [21 20]

    [23  4]
    [ 0 23]

    [15 14]
    [14 15]

    [27  0]
    [ 0  1]
Time: 0.016
*/

time GL2CMAdelicImage(E : MinimalImageLevel := true);
/*
MatrixGroup(2, IntegerRing(14)) of order 2^2 * 3 * 7
Generators:
    [12  3]
    [ 7 12]

    [13  0]
    [ 0  1]
Time: 0.047
*/


//In testing, we found that the higher the conductor of the curve, the larger the adelic level of definition, 
// and so the longer the computation took.
E := EllipticCurve([0, 0, 0, -485100, 129654000]); //CM curve with disc -16 and maximal conductor among curves of its CM class on the LMFDB
//time GL2CMAdelicImage(E);
/*
MatrixGroup(2, IntegerRing(1680))
Generators:
    [  43 1344]
    [1344   43]

    [ 645 1232]
    [ 112  645]

    [1133   85]
    [1340 1133]

    [ 527 1384]
    [1184  527]

    [ 691  890]
    [1480  691]

    [1393 1320]
    [1440 1393]

    [   1    0]
    [   0 1679]
Time: 649.875
*/

E := EllipticCurve([1, -1, 1, -965, -13940]); //CM curve with disc -7 and adelic level of def. 21
//If 21 divides N, the image returned is index 2 inside the normalizer
//If 21 does not divide N, the image is the full normalizer of Cartan
for N in [2..100] do
    GL2CMImageModN(E, N);
end for;
//Time: 5.375


//Error Examples

//Our code does not currently support curves with j-invariants 0 or 1728
E := EllipticCurve([1, 0]); //CM curve with j(E) = 1728
G := GL2CMAdelicImage(E);


E := EllipticCurve([0, 1]); //CM curve with j(E) = 0
G := GL2CMAdelicImage(E);

//Our code also doesn't support the images of non-CM curves (for obvious reasons)
E := EllipticCurve([1, 1]); //non-CM curve with j-invariant 6912/31

G := GL2CMAdelicImage(E);

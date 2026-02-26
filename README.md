# CM_adelic_Galois_images

This repository contains `Magma` code implementing an algorithm to compute the adelic Galois image of an elliptic curve E/Q with CM and j(E) not 0, 1728. This algorithm is the implementation of Algorithm 8.4 in the paper [*The Image of the Adelic Galois Representation of an Elliptic Curve with Complex Multiplication*] by Àlvaro Lozano-Robledo and Benjamin York.

The main code of interest is the function `GL2CMAdelicImage()` which takes as input an elliptic curve E/Q with CM and j(E) not 0, 1728, and outputs a Galois image modulo M, where M is an adelic level of definition for E.

We also implement the function `GL2CMImageModN()` which takes as input an elliptic curve E/Q with CM and j(E) not 0, 1728 and a positive integer N > 1, and outputs the Galois image modulo N of E.


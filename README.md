# CM_adelic_Galois_images

This repository contains `Magma` code implementing an algorithm to compute the adelic Galois image of an elliptic curve $E/\mathbb{Q}$ with CM and $j(E) \neq 0, 1728$. This algorithm is an implementation of Algorithm 8.4 in the paper [*The Image of the Adelic Galois Representation of an Elliptic Curve with Complex Multiplication*] by Àlvaro Lozano-Robledo and Benjamin York.

The main function is `GL2CMAdelicImage()` which takes as input an elliptic curve $E/\mathbb{Q}$ with CM and $j(E) \neq 0, 1728$, and outputs a Galois image of $E$ modulo $M$, where $M$ is an adelic level of definition for $E$. This function also has the optional parameter `MinimalImageLevel` which is set to `false` by default. If set to `true` then `GL2CMAdelicImage()` returns the Galois image of $E$ modulo $M_E$, where $M_E$ is the minimal adelic level of definition for $E$.

We also implement the function `GL2CMImageModN()` which takes as input an elliptic curve $E/\mathbb{Q}$ with CM and $j(E) \neq 0, 1728$ and a positive integer $N > 1$, and outputs the Galois image of $E$ modulo $N$.

To test the speed our code, we ran the function `GL2CMAdelicImage()` on all CM elliptic curves defined over $\mathbb{Q}$ with conductors $\leq 100000$ (this data set was taken from the [LMFDB](https://www.lmfdb.org/EllipticCurve/Q/) and consisted of 420 curves). The computation was performed on an HP EliteBook with an IntelCore i7-6600U processor and 16GB of RAM, which computed the images in about 30 minutes.


# Factorizer
 
This project is still a work in progress, but I implemented some driver programs to help debug what I have done.
 
`testeuclid` finds the formal derivatives of polynomials, the quotient and remainder after dividing the first polynomial inputted by the second, and the gcd of these polynomials over finite fields.
 
`testberlekamp` takes a polynomial as input and finds its Berlekamp matrix, Berlekamp subalgebra and its factors. It does this in two different ways. The first one sometimes finds trivial and reducible factors, and the second one trys to find all the irreducible factors by applying the first method recursively to reducible factors. This is still pretty buggy though (eg /test/berlekamp/poly6.txt).

`testlift` lifts roots of polynomials mod prime numbers to higher powers of those prime numbers using methods described in the constructive proof of Hensel's lemma. It then uses this system of congruences to find a root of the polynomial mod the product of these powers of primes. This is also based on a constructive proof, this time of the Chinese Remainder Theorem.

`factor` will hopefully tie all of this together to find roots of polynomials over finite fields and rings. I just need to get on top of my studies before I finish it.
                                                                          
All of these programs can be built with the Makefile in the src directory:
`make <program-name>`

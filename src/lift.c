/**
 * @file    lift.c
 * @brief   Implementation of Hensel's lemma and the Chinese Remainder Theorem.
 *
 * Hensel's lemma is used to lift roots of polynomials mod p to higher powers of
 * p, and the remainder theorem is used to solve systems of congruences under
 * the condition that the divisors are pairwise coprime.
 *
 * @author  L. Foxcroft
 * @date    2022-04-01
 */

#include <stdlib.h>
#include <stdio.h>
#include "euclid.h"
#include "lift.h"

/* --- function prototypes ---------------------------------------------------*/

int evaluate(Polynomial *f, int x, int m);
int inverse(int a, int m);

/* --- lift interface --------------------------------------------------------*/

int is_simple_root(Polynomial *f, int root, int m)
{
	/* Evaluate f'(root) mod m */
	Polynomial *f_prime = get_formal_derivative(f, m);
	int val = evaluate(f_prime, root, m);
	free_polynomial(f_prime);
	return val != 0;
}

int hensel(int *power, Polynomial *f, int root, int p, int k)
{
	int new_root = root;
	int m = p;
	int f_x;

	/* Calculate [f'(x)]^-1 */
	Polynomial *f_prime = get_formal_derivative(f, m);
	int f_prime_x = evaluate(f_prime, root, m);
	int f_prime_x_inv = inverse(f_prime_x, m);
	
	for (int i = 1; i <= k; i++) {
		/* Calculate f(new_root) */
		f_x = evaluate(f, new_root, m);

		/* Calculate new root */
		new_root = new_root - f_x * f_prime_x_inv;
		new_root = mod(new_root, m);

		/* Update divisor */
		if (i != k) {
			m *= p;
		}
	}

	/* Store p^k */
	*power = m;

	return new_root;
}

int chinese_remainder(int *product, int num_congruences, int *remainders,
		int *divisors)
{
	/* Bezout coefficients, current 2 solutions and moduli */
	int s, t, a1, a2, n1, n2;
	a1 = remainders[0];
	a2 = remainders[1];
	n1 = divisors[0];
	n2 = divisors[1];

	/* Iterate over equations and build up solution */
	for (int i = 1; i < num_congruences; i++) {
		a2 = remainders[i];
		n2 = divisors[i];
		extended_gcd_z(&s, &t, n1, n2);
		a1 = a1*t*n2 + a2*s*n1;
		n1 = n1*n2;
		a1 = mod(a1, n1);
	}

	/* Store product of divisors */
	*product = n1;

	return a1;
}

/* --- utility functions -----------------------------------------------------*/

/** Evaluate polynomial f(x) mod m */
int evaluate(Polynomial *f, int x, int m)
{
	int val = 0, power = 1;
	for (int i = 0; i <= f->degree; i++) {
		val += f->coefficients[i] * power;
		power *= x;
	}
	return mod(val, m);
}

/** Find the inverse of a in Z_m (assumes it exists) */
int inverse(int a, int m)
{
	int s, t;
	extended_gcd_z(&s, &t, a, m);
	return s;
}


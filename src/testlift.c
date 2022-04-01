/**
 * @file    testlift.c
 * @brief   A driver program to test the implementation of Hensel's lemma and
 *          the remainder theorem.
 * @author  L. Foxcroft
 * @date    2022-04-01
 */

#include <stdlib.h>
#include <stdio.h>
#include "euclid.h"
#include "lift.h"

/* --- main routine ----------------------------------------------------------*/

int main()
{
	/**
	 * Sage is very useful for constructing test cases. Eg:
	 * x = PolynomialRing(GF(7), 'x').gen()
	 * f = x^4 + x^3 + 2*x^2 + x - 13
	 * f.roots()
	 * (https://doc.sagemath.org/html/en/constructions/polynomials.html)
	 * Integers(n) can also be used for rings instead of fields
	 */

	/* Read in polynomial, and simple roots mod prime numbers which should be
	 * lifted */
	printf("Polynomial f(x)\n");
	Polynomial *f = scan_polynomial();

	int n;
	printf("Number of roots:\n");
	scanf("%d", &n);

	int *remainders = malloc(sizeof(int) * n);
	int *divisors = malloc(sizeof(int) * n);
	int *exponents = malloc(sizeof(int) * n);
	printf("simple root, prime divisor, power to lift to\n");
	for (int i = 0; i < n; i++) {
		scanf("%d %d %d", remainders + i, divisors + i, exponents + i);
	}

	/* Print polynomial to make output easier to follow */
	printf("f(x) = ");
	print_polynomial(f);
	printf("\n");

	/* Check roots are simple, and lift them if they are */
	int simple_roots = TRUE;
	int helper, new_root;
	for (int i = 0; i < n; i++) {
		if (!is_simple_root(f, remainders[i], divisors[i])) {
			printf("%d is not a simple root of f(x) mod %d\n",
					remainders[i], divisors[i]);
			simple_roots = FALSE;
		} else {
			printf("f(%d) = 0 mod %d -> ", remainders[i], divisors[i]);
			new_root = hensel(&helper, f, remainders[i], divisors[i],
					exponents[i]);
			remainders[i] = new_root;
			divisors[i] = helper;
			printf("f(%d) = 0 mod %d\n", remainders[i], divisors[i]);
		}
	}

	/* If all the roots were simple, assume they were lifted successfully, and
	 * use the remainder theorem to find the root mod the product of the
	 * divisors */
	if (simple_roots) {
		int x = chinese_remainder(&helper, n, remainders, divisors);
		printf("f(%d) = 0 mod %d\n", x, helper);
	}

	return EXIT_SUCCESS;
}

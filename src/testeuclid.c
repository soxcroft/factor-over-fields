/**
 * @file    file.c
 * @brief   A driver program to test the implementation of some methods in 
 *          euclid.c
 */

#include <stdlib.h>
#include <stdio.h>
#include "euclid.h"

int main()
{
	/* TODO make this more interactive */

	/* get p for field Z_p */
	int p;
	printf("P for Z_p? ");
	scanf("%d", &p);

	/* Scan polynomials */
	Polynomial *poly1, *poly2;
	poly1 = scan_polynomial();
	poly2 = scan_polynomial();

	/* Print what program scanned, to make sure it checks out */
	printf("Working in Z_%d\n", p);
	print_polynomial(poly1);
	printf("\n");
	print_polynomial(poly2);
	printf("\n");

	/* Get formal derivatives */
	Polynomial *d1, *d2;
	d1 = get_formal_derivative(poly1, p);
	d2 = get_formal_derivative(poly2, p);
	printf("Derivative of p1: ");
	print_polynomial(d1);
	printf("\nDerivative of p2: ");
	print_polynomial(d2);
	printf("\n");

	/* Divide p1 by p2 and print the quotient and remainder after division */
	Polynomial *q, *r;
	long_div(&q, &r, poly1, poly2, p);
	print_polynomial(poly1);
	printf(" = (");
	print_polynomial(poly2);
	printf(")(");
	print_polynomial(q);
	printf(") + ");
	print_polynomial(r);
	printf("\n");

	/* Calculate the gcd of p1 and p2 */
	Polynomial *gcd = gcd_p(poly1, poly2, p);
	printf("gcd(p1, p2) = ");
	print_polynomial(gcd);
	printf("\n");

	/* free allocated memory */
	free_polynomial(poly1);
	free_polynomial(poly2);
	free_polynomial(d1);
	free_polynomial(d2);
	free_polynomial(q);
	free_polynomial(r);
	free_polynomial(gcd);

	return EXIT_SUCCESS;
}

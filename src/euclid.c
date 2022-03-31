/**
 * @file    euclid.c
 * @brief   Implementations of Euclid's division and gcd algorithm for Z_p and
 *          polynomials over Z_p.
 */

#include <stdlib.h>
#include <stdio.h>
#include "euclid.h"

/* --- function prototypes --------------------------------------------------*/

static int lc(Polynomial *p);
static int degree(Polynomial *p);
static int is_zero(Polynomial *p);

/* --- euclid interface -----------------------------------------------------*/

/* INTEGER FUNCTIONS */

int mod(int a, int p)
{
	return (a % p + p) % p;
}

void extended_gcd_z(int *s, int *t, int a, int m)
{
	int s0 = 1, s1 = 0;
	int t0 = 0, t1 = 1;
	int r0 = a, r1 = m;
	int q, helper;

	while (r1 != 0) {
		q = r0 / r1;

		helper = s0 - q*s1;
		s0 = s1;
		s1 = helper;

		helper = t0 - q*t1;
		t0 = t1;
		t1 = helper;
		
		helper = r0 - q*r1;
		r0 = r1;
		r1 = helper;
	}

	*s = s0;
	*t = t0;
}

/* POLYNOMIAL FUNCTIONS */

void free_polynomial(Polynomial *polynomial)
{
	free(polynomial->coefficients);
	free(polynomial);
}

Polynomial *init_polynomial(int degree)
{
	Polynomial *p = malloc(sizeof(Polynomial));
	p->degree = degree;
	p->coefficients = malloc(sizeof(int) * (degree + 1));
	for (int i = 0; i <= degree; i++) {
		p->coefficients[i] = 0;
	}
	return p;
}

Polynomial *copy_polynomial(Polynomial *poly)
{
	Polynomial *new_poly = malloc(sizeof(Polynomial));
	new_poly->degree = poly->degree;
	new_poly->coefficients = malloc(sizeof(int) * (new_poly->degree + 1));
	for (int i = 0; i <= poly->degree; i++) {
		new_poly->coefficients[i] = poly->coefficients[i];
	}
	return new_poly;
}

Polynomial *scan_polynomial()
{
	printf("Enter the degree of your polynomial:\n");
	Polynomial *polynomial = malloc(sizeof(Polynomial));
	scanf("%d", &polynomial->degree);

	printf("Enter the coefficients of your polynomial (from lowest order term to highest):\n");
	polynomial->coefficients = malloc(sizeof(int) * (polynomial->degree + 1));
	for (int i = 0; i <= polynomial->degree; i++) {
		scanf("%d", polynomial->coefficients + i);
	}

	return polynomial;
}


void print_polynomial(Polynomial *polynomial)
{
	int printed_first_term = FALSE;
	int *p = polynomial->coefficients;
	int n = polynomial->degree;

	/* print polynomial */
	for (int i = 0; i <= n; i++) {
		if (p[i] != 0) {
			if (!printed_first_term) {
				printf("%d", p[i]);
				printed_first_term = TRUE;
			} else {
				printf(" + %d", p[i]);
			}
			if (i != 0) {
				printf("*x^%d", i);
			}
		}
	}
}

Polynomial *get_formal_derivative(Polynomial *p, int m)
{
	Polynomial *derivative = malloc(sizeof(Polynomial));
	derivative->degree = p->degree; /* not -1 in case p is already constant */
	derivative->coefficients = malloc(sizeof(int) * (p->degree + 1));

	for (int i = 0; i < p->degree; i++) {
		derivative->coefficients[i] = mod((i + 1) * p->coefficients[i + 1], m);
	}
	/* highest powers coefficient falls away */
	derivative->coefficients[p->degree] = 0;

	return derivative;
}

void long_div(Polynomial **q, Polynomial **r, Polynomial *p1, Polynomial *p2, 
		int m)
{
	/* initialize quotient */
	*q = init_polynomial(p1->degree); /* has a max degree of deg(r)-deg(p2) */

	/* initialize remainder */
	*r = init_polynomial(p1->degree);
	for (int i = 0; i <= p1->degree; i++) {
		(*r)->coefficients[i] = p1->coefficients[i];
	}
	Polynomial *sb = init_polynomial(p1->degree + p2->degree);

	/* initialize helper variables */
	int d = degree(p2);
	int c = lc(p2);
	int deg_r = degree(*r);
	int mult_factor, s, t;

	while (deg_r >= d && !is_zero(*r)) {
		/* calculate polynomial s*p2 (stored in sb) */
		extended_gcd_z(&s, &t, mod(c, m), m);
		mult_factor = mod(lc(*r) * s, m);

		/* calculate s*b, = lc(r)/c * x^(deg(r)-d) * b */
		for (int i = 0; i <= sb->degree; i++) {
			sb->coefficients[i] = 0;
		}
		for (int i = 0; i <= p2->degree; i++) {
			sb->coefficients[i + (deg_r - d)] = p2->coefficients[i] * mult_factor;
		}

		/* q = q + s */
		(*q)->coefficients[deg_r - d] = mod((*q)->coefficients[deg_r - d] +
				mult_factor, m);

		/* r = r - sb */
		for (int i = 0; i <= (*r)->degree; i++) {
			(*r)->coefficients[i] = mod((*r)->coefficients[i] - sb->coefficients[i], m);
		}

		/* recalculate deg_r */
		deg_r = degree(*r);
	}

	/* free sb, q and r should be handled after outside of function */
	free_polynomial(sb);
}

Polynomial *gcd_p(Polynomial *p1, Polynomial *p2, int m)
{
	/* initialize remainders to p1 and p2 */
	Polynomial *r0 = init_polynomial(p1->degree);
	for (int i = 0; i <= p1->degree; i++) {
		r0->coefficients[i] = p1->coefficients[i];
	}
	Polynomial *r1 = init_polynomial(p2->degree);
	for (int i = 0; i <= p2->degree; i++) {
		r1->coefficients[i] = p2->coefficients[i];
	}
	Polynomial *helper, *q;

	while (!is_zero(r1)) {
		long_div(&q, &helper, r0, r1, m);
		free_polynomial(r0);
		free_polynomial(q);
		r0 = r1;
		r1 = helper;
	}

	free_polynomial(r1);

	return r0;
}

/* --- utility functions -----------------------------------------------------*/

/** Returns the leading coefficient of a polynomial */
int lc(Polynomial *p)
{
	for (int i = p->degree; i >= 0; i--) {
		if (p->coefficients[i] != 0) {
			return p->coefficients[i];
		}
	}
	return 0;
}

/** Returns the actual degree of a polynomial. */
int degree(Polynomial *p)
{
	for (int i = p->degree; i >= 0; i--) {
		if (p->coefficients[i] != 0) {
			return i;
		}
	}
	return 0;
}

/** Checks if a polynomial is the zero polynomial */
int is_zero(Polynomial *p)
{
	for (int i = 0; i <= p->degree; i++) {
		if (p->coefficients[i] != 0) {
			return FALSE;
		}
	}
	return TRUE;
}

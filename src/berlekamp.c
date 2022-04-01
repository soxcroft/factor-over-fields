/**
 * @file    berlekamp.c
 * @brief   Implementation of Berlekamp's algorithm for factorising polynomials
 *          over finite fields.
 */

#include <stdlib.h>
#include <stdio.h>
#include "euclid.h"
#include "berlekamp.h"

/* --- function prototypes ---------------------------------------------------*/

int is_constant(Polynomial *p);

/* --- berlekamp interface ---------------------------------------------------*/

int **get_berlekamp_matrix(Polynomial *p, int m)
{
	/* initialize helper polynomial and matrix */
	int degree = p->degree;
	Polynomial *helper = init_polynomial(m * degree);
	int **matrix;
	matrix = malloc(sizeof(int *) * degree);

	/* compute powers x^(2i) mod p */
	Polynomial *q, *r; /* quotient and remainder */
	for (int i = 0; i < degree; i++) {
		if (i != 0) {
			helper->coefficients[m*(i - 1)] = 0;
		}
		helper->coefficients[m*i] = 1;

		long_div(&q, &r, helper, p, m);
		matrix[i] = r->coefficients;

		free_polynomial(q);
		free(r);
	}

	free_polynomial(helper);

	return matrix;
}

void transpose(int ***A, int m, int n)
{
	/* initialize a new matrix which will replace A */
	int **B = malloc(sizeof(int *) * n);
	for (int i = 0; i < n; i++) {
		B[i] = malloc(sizeof(int) * m);
	}

	/* copy values from A into B */
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			B[j][i] = (*A)[i][j];
		}
	}

	/* free old matrix and assign B to *A */
	free_matrix(*A, m);

	*A = B;
}

void subtract_identity(int **A, int m, int n, int p)
{
	/* Check A is square */
	if (m != n) {
		printf("ERROR: Cannot subtract identity, A is not square");
		return;
	}

	/* Subtract I from A */
	for (int i = 0; i < m; i++) {
		A[i][i] = mod(A[i][i] - 1, p);
	}
}

void gauss_jordan(int **A, int m, int n, int p)
{
	/* Initialize and declare variables */
	int lead = 0;
	int i, j, s, t, factor, *helper;

	for (int r = 0; r < m && lead < n; r++) {
		/* Find row with pivot element in 'lead' column */
		i = r;
		while (A[i][lead] % p == 0) {
			i++;
			/* If we have exhausted rows, increment lead and start over */
			if (i == m) {
				i = r;
				lead++;
				if (lead == n) {
					return;
				}
			}
		}

		/* Swap rows i and r */
		if (i != r) {
			helper = A[i];
			A[i] = A[r];
			A[r] = helper;
		}

		/* Multiply row r by inverse of A[r][lead] */
		extended_gcd_z(&s, &t, p, A[r][lead]);
		for (i = lead; i < n; i++) {
			A[r][i] = mod(A[r][i] * t, p);
		}

		/* Make sure col lead only has an element in row r */
		for (i = 0; i < m; i++) {
			if (i != r && A[i][lead] % p != 0) {
				factor = A[i][lead];
				for (j = 0; j < n; j++) {
					A[i][j] = mod(A[i][j] - factor * A[r][j], p);
				}
			}
		}

		lead++;
	}
}

int **null_space(int *rank, int **R, int m, int n, int p)
{
	/* First store the pivot elements positions and count them to get the rank */
	*rank = 0;
	int col = 0, row = 0;
	int *pivot = malloc(sizeof(int) * m);
	for (int i = 0; i < m; i++) {
		pivot[i] = -1; /* corresponds to no pivot in row */
	}
	while (row < m && col < n) {
		if (R[row][col] != 0) {
			pivot[row] = col;
			(*rank)++;
			row++;
			col++;
		} else {
			col++;
		}
	}

	/* Nullity is cols - rank */
	int **kernel = malloc(sizeof(int *) * (n - *rank));
	for (int i = 0; i < n - *rank; i++) {
		kernel[i] = malloc(sizeof(int) * n);
		for (int j = 0; j < n; j++) {
			kernel[i][j] = 0;
		}
	}

	/* Now we find the free variables and update the kernel */
	int free_variables = 0;
	col = 0;
	row = -1;
	while (row < m && col < n) {
		if (row + 1 < m && R[row+1][col] != 0) {
			/* Pivot element, so go to next row */
			row++;
			col++;
		} else {
			/* Free variable */
			for (int i = 0; i < m; i++) {
				if (row >= 0 && pivot[i] != -1) {
					kernel[free_variables][pivot[i]] = mod(-R[i][col], p);
				}
			}
			kernel[free_variables][col] = 1;
			free_variables++;
			col++;
		}
	}

	free(pivot);

	return kernel;
}

Polynomial **kernel_to_arr(int **kernel, int m, int n)
{
	Polynomial **arr = malloc(sizeof(Polynomial *) * m);
	for (int i = 0; i < m; i++) {
		arr[i] = malloc(sizeof(Polynomial));
		arr[i]->degree = n - 1;
		arr[i]->coefficients = malloc(sizeof(int) * n);
		for (int j = 0; j < n; j++) {
			arr[i]->coefficients[j] = kernel[i][j];
		}
	}
	return arr;
}

Polynomial **factors(Polynomial *p, Polynomial **subalgebra, int nullity, int m)
{
	/* p has nullity distinct factors */
	Polynomial **facs = malloc(sizeof(Polynomial *) * nullity);
	int counter = 0;

	/* Find index of a non trivial polynomial in the subalgebra */
	int ip = -1; /* XXX change back to -1 */
	for (int i = 0; i < nullity && ip == -1; i++) {
		if (!is_constant(subalgebra[i])) {
			ip = i;
		}
	}

	/* NULL indicates there are no non-trivial factors */
	if (ip == -1) {
		free(facs);
		return NULL;
	}

	/* f(x) = product[s in F_q](gcd(p, g(x)-s), where F_q is field */
	Polynomial *factor;
	for (int s = 0; s < nullity; s++) {
		subalgebra[ip]->coefficients[0] -= s;
		factor = gcd_p(p, subalgebra[ip], m);
		subalgebra[ip]->coefficients[0] += s;
		facs[counter] = factor;
		counter++;
	}

	/* Return what we found and worry about discarding/reducing factors later */
	return facs;
}

void free_factors(Polynomial **factors, int nullity)
{
	for (int i = 0; i < nullity; i++) {
		free_polynomial(factors[i]);
	}
	free(factors);
}

void free_matrix(int **matrix, int m)
{
	for (int i = 0; i < m; i++) {
		free(matrix[i]);
	}
	if (matrix) {
		free(matrix);
	}
}

Polynomial **berlekamp(int *num_factors, Polynomial *poly, int m)
{
	/* Get Berlekamp subalgebra */
	int **matrix = get_berlekamp_matrix(poly, m);
	subtract_identity(matrix, poly->degree, poly->degree, m);
	transpose(&matrix, poly->degree, poly->degree);
	gauss_jordan(matrix, poly->degree, poly->degree, m);

	int **kernel, rank;
	kernel = null_space(&rank, matrix, poly->degree, poly->degree, m);

	*num_factors = poly->degree - rank;

	if (*num_factors == 0 || *num_factors == 1) {
		/* free memory allocated so far */
		free_matrix(matrix, poly->degree);
		free_matrix(kernel, *num_factors);
		/* polynomial is irreducible */
		*num_factors = 1;
		Polynomial **facs = malloc(sizeof(Polynomial *));
		*facs = copy_polynomial(poly);
		return facs;
	}

	Polynomial **subalgebra = kernel_to_arr(kernel, poly->degree - rank, poly->degree);

	/* Now, find factors of poly and recursively call berlekamp on them until we
	 * are left with irreducible polynomial factors */
	Polynomial **check = factors(poly, subalgebra, poly->degree - rank, m);
	if (!check) {
		/* No non-trivial factors */
		/* TODO not sure if this will ever be executed */
		*num_factors = 1;
		Polynomial **facs = malloc(sizeof(Polynomial *));
		*facs = copy_polynomial(poly);
		return facs;
	}

	Polynomial **facs = malloc(sizeof(Polynomial *) * (poly->degree - rank));
	Polynomial **reduced_facs; /* use to store factors of factors */
	int counter, helper;
	counter = 0;

	/* iterate over factors, discard trivial ones, reduce reducible ones */
	for (int i = 0; i < poly->degree - rank; i++) {
		if (!is_constant(check[i])) {
			reduced_facs = berlekamp(&helper, check[i], m);
			for (int j = 0; j < helper; j++) {
				facs[counter++] = reduced_facs[j];
			}
			free(reduced_facs);
		}
	}

	/* FIXME num_factors should already be correct ...  this is just to try fix
	 * poly6 */
	*num_factors = counter;

	/* Free allocated memory */
	free_matrix(matrix, poly->degree);
	free_matrix(kernel, poly->degree - rank);
	free_factors(subalgebra, poly->degree - rank);
	free_factors(check, poly->degree - rank);

	return facs;
}

/* --- utility functions -----------------------------------------------------*/

/** Return true if polynomial is a constant */
int is_constant(Polynomial *p)
{
	int trivial = TRUE;
	for (int i = 1; i <= p->degree; i++) {
		if (p->coefficients[i] != 0) {
			trivial = FALSE;
		}
	}
	return trivial;
}

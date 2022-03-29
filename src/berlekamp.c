/**
 * @file    berlekamp.c
 * @brief   Implementation of Berlekamp's algorithm for factorising polynomials
 *          over finite fields.
 */

#include <stdlib.h>
#include <stdio.h>
#include "berlekamp.h"

/* TODO put these in euclid.c (change yuns.c, organise methods in terms of
 * polynomials and integers) */

/* Mod, like % but always returns positive values */
int mod(int a, int p)
{
	return (a % p + p) % p;
}

/* From yuns.c */
void extended_gcd(int *s, int *t, int a, int m)
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

/* I/O METHODS */

void scan_matrix(int ***arr, int m, int n)
{
	/* initialize matrix */
	*arr = malloc(sizeof(int *) * m);
	for (int i = 0; i < m; i++) {
		(*arr)[i] = malloc(sizeof(int) * n);
	}

	/* scan for values */
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			scanf("%d", (*arr)[i] + j);
		}
	}
}

void print_matrix(int **arr, int m, int n)
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n - 1; j++) {
			printf("%d ", arr[i][j]);
		}
		printf("%d\n", arr[i][n-1]);
	}
}

/* LINEAR ALGEBRA */

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
	for (int i = 0; i < m; i++) {
		free((*A)[i]);
	}
	free(*A);
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
		extended_gcd(&s, &t, p, A[r][lead]);
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
				if (row >= 0 && pivot[row] != -1) {
					kernel[free_variables][pivot[i]] = mod(-R[i][col], p);
				}
			}
			kernel[free_variables][col] = 1;
			free_variables++;
			col++;
		}
	}

	return kernel;
}


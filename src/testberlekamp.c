/**
 * @file    testberlekamp.c
 * @brief   A driver program to test the implementation of Berlekamp's algorithm
 */

#include <stdlib.h>
#include <stdio.h>
#include "euclid.h"
#include "berlekamp.h"

/* --- function prototypes ---------------------------------------------------*/

void scan_matrix(int ***arr, int m, int n);
void print_matrix(int **arr, int m, int n);

/* --- main routine ----------------------------------------------------------*/

int main()
{
	int p;
	printf("P for Z_p? ");
	scanf("%d", &p);
	Polynomial *polynomial = scan_polynomial();

	printf("Working in Z_%d\n", p);
	print_polynomial(polynomial);
	printf("\n");

	int **matrix = get_berlekamp_matrix(polynomial, p);
	int m = polynomial->degree;
	printf("Berlekamp matrix:\n");
	print_matrix(matrix, m, m);

	subtract_identity(matrix, m, m, p);
	printf("Matrix - I\n");
	print_matrix(matrix, m, m);

	transpose(&matrix, m, m);
	printf("Transposed\n");
	print_matrix(matrix, m, m);

	gauss_jordan(matrix, m, m, p);
	printf("Row reduced\n");
	print_matrix(matrix, m, m);

	int **kernel, rank;
	kernel = null_space(&rank, matrix, m, m, p);
	printf("Kernel, rank %d\n", rank);
	print_matrix(kernel, m - rank, m);

	Polynomial **subalgebra = kernel_to_arr(kernel, m - rank, m);
	printf("Subalgebra\n");
	for (int i = 0; i < m - rank; i++) {
		print_polynomial(subalgebra[i]);
		printf("\n");
	}

	Polynomial **facs = factors(polynomial, subalgebra, m - rank, p);
	printf("Factors, nullity %d\n", m - rank);
	if (facs == NULL) {
		printf("no non-trivial factors\n");
	} else {
		for (int i = 0; i < m - rank; i++) {
			print_polynomial(facs[i]);
			printf("\n");
		}
	}

	int num_factors;
	facs = berlekamp(&num_factors, polynomial, p);
	/* TODO doesn't always create enough factors... */
	printf("Berlekamp, %d factors\n", num_factors);
	for (int i = 0; i < num_factors; i++) {
		print_polynomial(facs[i]);
		printf("\n");
	}

	return EXIT_SUCCESS;
}

/* --- functions -------------------------------------------------------------*/

/** Scans m x n integers into a matrix */
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

/** Prints an m x n matrix to the terminal */
void print_matrix(int **arr, int m, int n)
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n - 1; j++) {
			printf("%d ", arr[i][j]);
		}
		printf("%d\n", arr[i][n-1]);
	}
}


/**
 * @file    testberlekamp.c
 * @brief   A driver program to test the implementation of Berlekamp's algorithm
 */

#include <stdlib.h>
#include <stdio.h>
#include "berlekamp.h"

int main()
{
	int m, n, p, **matrix, **kernel, rank;

	scanf("%d %d %d", &m, &n, &p);
	printf("m n p\n");
	printf("%d %d %d\n", m, n, p);

	scan_matrix(&matrix, m, n);
	printf("Matrix\n");
	print_matrix(matrix, m, n);

	subtract_identity(matrix, m, n, p);
	printf("Matrix - I\n");
	print_matrix(matrix, m, n);

	transpose(&matrix, m, n);
	printf("Transposed\n");
	print_matrix(matrix, n, m);

	int helper = m;
	m = n;
	n = helper;

	gauss_jordan(matrix, m, n, p);
	printf("Row reduced\n");
	print_matrix(matrix, m, n);

	kernel = null_space(&rank, matrix, m, n, p);
	printf("Kernel, rank %d\n", rank);
	print_matrix(kernel, n - rank, m);

	return EXIT_SUCCESS;
}

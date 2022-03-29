/**
 * @file 	berlekamp.h
 * @brief 	Prototypes for methods used in Berlekamp's algorithm 
 */

#ifndef BERLEKAMP
#define BERLEKAMP

int mod(int a, int m); /* TODO move this to euclid stuff */
void extended_gcd(int *s, int *t, int a, int m); /* TODO also move */

/**
 * Allocates memory for and reads matrix into *arr.
 *
 * @param[in]   arr
 *     pointer to a double pointer which should store the scanned matrix
 * @param[in]   m
 *     the number of rows in the matrix
 * @param[in]   n
 *     the number of columns in the matrix
 */
void scan_matrix(int ***arr, int m, int n);

/**
 * Prints an m x n matrix of integers to the terminal.
 *
 * @param[in]   arr
 *     double pointer to a matrix
 * @param[in]   m
 *     the number of rows in the matrix
 * @param[in]   n
 *     the number of columns in the matrix
 */
void print_matrix(int **arr, int m, int n);

/**
 * Transpose an m x n matrix, in place.
 *
 * @param[in]   A
 *     pointer to a double pointer which should store transposed matrix
 * @param[in]   m
 *     the number of rows in the matrix
 * @param[in]   n
 *     the number of columns in the matrix
 */
void transpose(int ***A, int m, int n);

/**
 * Subtracts the identity matrix from the provided matrix if possible.
 *
 * @param[in]   A
 *     double pointer to a matrix
 * @param[in]   m
 *     the number of rows in the matrix
 * @param[in]   n
 *     the number of columns in the matrix
 * @param[in]   p
 *     the modulus we are working with (Z_p is a field)
 */
void subtract_identity(int **A, int m, int n, int p);

/**
 * Performs Gauss-Jordan elimination on the given matrix to get it in reduced
 * row echelon form. This is done in place, so the original matrix passed as
 * input is lost.
 *
 * @param[in]   A
 *     double pointer to a matrix
 * @param[in]   m
 *     the number of rows in the matrix
 * @param[in]   n
 *     the number of columns in the matrix
 * @param[in]   p
 *     the modulus we are working with (Z_p)
 */
void gauss_jordan(int **A, int m, int n, int p);

/**
 * Finds the (right) null space of a matrix given its reduced row echelon form.
 *
 * @param[in]   rank
 *     pointer to an integer which we should store the rank in for later use
 * @param[in]   R
 *     double pointer to a matrix in reduced row echelon form
 * @param[in]   m
 *     the number of rows in the matrix
 * @param[in]   n
 *     the number of columns in the matrix
 * @param[in]   p
 *     the modulus we are working with (Z_p is a field)
 * @return      a matrix whose rows contain the basis vectors for Rs null space
 */
int **null_space(int *rank, int **R, int m, int n, int p);

#endif

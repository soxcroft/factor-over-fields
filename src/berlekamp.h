/**
 * @file 	berlekamp.h
 * @brief 	Prototypes for methods used in Berlekamp's algorithm 
 */

#ifndef BERLEKAMP
#define BERLEKAMP

#include "euclid.h"

/**
 * Find Berlekamp matrix. Berlekamp subalgebra is kernel of matrix derived from
 * this.
 *
 * @param[in] p
 *     polynomial f(x) which we will use to calculate matrix
 * @param[in] m
 *     prime number so that we work over field Z_m
 * @return    Berlekamp matrix derived from p = f(x)
 */
int **get_berlekamp_matrix(Polynomial *p, int m);

/**
 * Transpose an m x n matrix, in place.
 *
 * @param[in] A
 *     pointer to a double pointer which should store transposed matrix
 * @param[in] m
 *     the number of rows in the matrix
 * @param[in] n
 *     the number of columns in the matrix
 */
void transpose(int ***A, int m, int n);

/**
 * Subtracts the identity matrix from the provided matrix if possible.
 *
 * @param[in] A
 *     double pointer to a matrix
 * @param[in] m
 *     the number of rows in the matrix
 * @param[in] n
 *     the number of columns in the matrix
 * @param[in] p
 *     the modulus we are working with (Z_p is a field)
 */
void subtract_identity(int **A, int m, int n, int p);

/**
 * Performs Gauss-Jordan elimination on the given matrix to get it in reduced
 * row echelon form. This is done in place, so the original matrix passed as
 * input is lost.
 *
 * @param[in] A
 *     double pointer to a matrix
 * @param[in] m
 *     the number of rows in the matrix
 * @param[in] n
 *     the number of columns in the matrix
 * @param[in] p
 *     the modulus we are working with (Z_p)
 */
void gauss_jordan(int **A, int m, int n, int p);

/**
 * Finds the (right) null space of a matrix given its reduced row echelon form.
 * The null space of B - I, where B is a Berlekamp matrix, is the the Berlekamp
 * subalgebra.
 *
 * @param[in] rank
 *     pointer to an integer which we should store the rank in for later use
 * @param[in] R
 *     double pointer to a matrix in reduced row echelon form
 * @param[in] m
 *     the number of rows in the matrix
 * @param[in] n
 *     the number of columns in the matrix
 * @param[in] p
 *     the modulus we are working with (Z_p is a field)
 * @return    a matrix whose rows contain the basis vectors for Rs null space
 */
int **null_space(int *rank, int **R, int m, int n, int p);

/**
 * Converts a matrix to an array of pointers to polynomials. Each row in the
 * matrix corresponds to the array of coefficients for the polynomials.
 *
 * @param[in] kernel
 *     a matrix of polynomial's coefficients
 * @param[in] m
 *     the number of rows in the matrix
 * @param[in] n
 *     the number of columns in the matrix
 * @return    an array of pointers to polynomials
 */
Polynomial **kernel_to_arr(int **kernel, int m, int n);

/**
 * Takes in a polynomial and the polynomials in its berlekamp subalgebra to find
 * (potentially trivial) factors of a polynomial.
 *
 * @param[in] p
 *     pointer to the polynomial to be factorised
 * @param[in] subalgebra
 *     array of pointers to polynomials in p's Berlekamp subalgebra
 * @param[in] nullity
 *     the number of polynomials in the subalgebra
 * @param[in] m
 *     prime number for field Z_m
 * @return    NULL if p is irreducible, else array of pointers to factors found
 */
Polynomial **factors(Polynomial *p, Polynomial **subalgebra, int nullity, int m);

/**
 * Berlekamp's algorithm. Takes in a polynomial defined over Z_m as input, finds
 * its square free factorization, and returns an array of its factors.
 *
 * @param[in] num_factors
 *     pointer to the number of factors found, written to in function
 * @param[in] poly
 *     pointer to the polynomial over Z_m to be factorised
 * @param[in] m
 *     prime number so that we can work over field Z_m
 * @return    an array of pointers to polynomial factors of poly
 */
Polynomial **berlekamp(int *num_factors, Polynomial *poly, int m);

#endif

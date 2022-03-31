/**
 * @file    euclid.h
 * @brief   Prototypes for methods involving euclid's division algorithms, used
 *          with integers and polynomials in this project.
 */

#ifndef EUCLID
#define EUCLID

#define FALSE 0
#define TRUE 1

typedef struct polynomial {
	int degree;
	int *coefficients;
} Polynomial;

/* INTEGERS */

/** 
 * Mod, like % but always returns positive values. 
 *
 * @param[in] a
 *     the dividend
 * @param[in] p
 *     the divisor
 * return     signed remainder after dividing a by p
 */
int mod(int a, int p);

/**
 * Extended Euclidean algorithm. Calculates the gcd(a, m) and stores the Bezout
 * coefficients in *s and *t. a and m should be integers.
 *
 * @param[in] s
 *     pointer to where the first Bezout coefficient should be stored
 * @param[in] t
 *     pointer to where the second Bezout coefficient should be stored
 * @param[in] a
 *     first integer in gcd
 * @param[in] b
 *     second integer in gcd
 */
void extended_gcd_z(int *s, int *t, int a, int m);

/* POLYNOMIALS */

/**
 * Free the memory allocated for a polynomial.
 *
 * @param[in] polynomial
 *     the polynomial to be freed
 */
void free_polynomial(Polynomial *polynomial);

/**
 * Prints a polynomial to the terminal, from lowest order coefficient to
 * highest.
 *
 * @param[in] polynomial
 *     pointer to the polynomial which should be printed
 */
void print_polynomial(Polynomial *polynomial);

/**
 * Allocates memory for and returns a polynomial. Initializes all its 
 * coefficients to zero.
 *
 * @param[in] degree
 *     the degree of the polynomial to be initialized
 */
Polynomial *init_polynomial(int degree);

/**
 * Allocates memory for and returns a pointer to a new polynomial, which is a
 * copy of the one passed as input.
 *
 * @param[in] poly
 *     the polynomial to be copied
 * @return    a pointer to a copy of poly
 */
Polynomial *copy_polynomial(Polynomial *poly);

/**
 * Scans for the degree of a polynomial, then allocates memory for a new
 * polynomial, scans in all its coefficients and finally returns a pointer to it.
 *
 * @return    the polynomial read in from the command line
 */
Polynomial *scan_polynomial(void);

/**
 * Gets the formal derivative. Calculated the same way that we are used to for
 * polynomials, just a different name because we can't use limits when working
 * with integers.
 *
 * @param[in] p
 *     the polynomial whose derivative should be calculated
 * @param[in] m
 *     prime number so that we can ensure coefficients are in Z_m
 * @return    the formal derivative of p
 */
Polynomial *get_formal_derivative(Polynomial *p, int m);

/** 
 * Euclidean division of polynomial 1 by polynomial 2 over a finite field (Z_m, 
 * where m is prime). Writes to q, the quotient, and r, the remainder.
 *
 * @param[in] q
 *     double pointer to a polynomial, where quotient should be written
 * @param[in] r
 *     double pointer to a polynomial, where remainder should be written
 * @param[in] p1
 *     pointer to the dividend
 * @param[in] p2
 *     pointer to the divisor
 * @param[in] m
 *     prime number so that we can work with field Z_m
 */
void long_div(Polynomial **q, Polynomial **r, Polynomial *p1, Polynomial *p2,
		int m);

/**
 * Euclid's algorithm for calculating the gcd of 2 polynomials.
 * 
 * @param[in] p1
 *     pointer to the first polynomial
 * @param[in] p2
 *     pointer to the second polynomial
 * @param[in] m
 *     prime number, so we work with field Z_m
 * @return    the gcd of polynomials p1 and p2
 */
Polynomial *gcd_p(Polynomial *p1, Polynomial *p2, int m);

#endif

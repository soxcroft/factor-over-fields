/**
 * @file    lift.c
 * @brief   Functions for implementing Hensel's lemma and the Chinese Remainder
 *          Theorem.
 * @author  L. Foxcroft
 * @date    2022-04-01
 */

#ifndef LIFT
#define LIFT

#include "euclid.h"

/** 
 * Checks if a root of a polynomial has multiplicity 1 by evaluating its
 * derivative at root. If f'(root) == 0, it is a simple root. Does not check if
 * root is actually a root of the polynomial though.
 *
 * @param[in] f
 *     pointer to the polynomial, f(x)
 * @param[in] root
 *     a root of the polynomial
 * @param[in] m
 *     the divisor, which specifies the ring Z_m that we are working over
 * @return    true if root is a simple root of f, else false
 */
int is_simple_root(Polynomial *f, int root, int m);


/** 
 * Lift a root of a polynomial mod p to a root mod p^k using Hensel's lemma. The
 * proof uses a constructive argument. Assumes that root is a simple root.
 *
 * @param[out] power
 *     pointer to where p^k should be written
 * @param[in]  f
 *     the polynomial
 * @param[in]  root
 *     a simple root of the polynomial
 * @param[in]  p
 *     specifies ring of integers, Z_p, that we are working over
 * @param[in]  k
 *     the power that p should be raised to
 * @return     a root of the polynomial mod p^k
 */
int hensel(int *power, Polynomial *f, int root, int p, int k);

/**
 * Takes in the remainders after the Euclidean division of some integer x by the
 * corresponding divisors, and returns the remainder of x divided by the product
 * of these divisors (mod the product). Assumes the divisors are coprime. This
 * algorithm is based on the constructive proof of the Chinese remainder
 * theorem.
 *
 * @param[out] product
 *     pointer to where the product of the divisors should be stored
 * @param[in] num_congruences
 *     the number of equations in the system of congruences
 * @param[in] remainder
 *     array of integers containing remainders after division
 * @param[in] divisors
 *     array of integers containing the divisors used in equations
 * @return    x mod the product of the divisors
 */
int chinese_remainder(int *product, int num_congruences, int *remainders,
		int *divisors);

#endif

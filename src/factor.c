/**
 * @file    factor.c
 * @brief   Finds the factors of a polynomial mod m.
 * @author  L. Foxcroft
 * @date    TODO
 */

#include <stdlib.h>
#include <stdio.h>
#include "euclid.h"
#include "berlekamp.h"
#include "lift.h"

/* --- type definitions ------------------------------------------------------*/

typedef struct factor_s {
	int p;
	int exp;
	struct factor_s *next;
} Factor;

/* --- function prototypes ---------------------------------------------------*/

Factor *init_factor(int p, int exp);
Factor *factorize(int *n, int x);
void free_factors(Factor *factors);

/* --- main routine ----------------------------------------------------------*/

int main()
{
	/* Read in polynomial and divisor */
	//Polynomial *p = scan_polynomial();
	int m;
	printf("Factor polynomial mod ...\n");
	scanf("%d", &m);

	/* Factorize divisor, and convert linked list into parallel array, 1
	 * containing the factors (divisors) and the other the exponents */
	int num_factors;
	Factor *factors = factorize(&num_factors, m);
	Factor *helper = factors;

	int *divisors = malloc(sizeof(int) * num_factors);
	int *exponents = malloc(sizeof(int) * num_factors);
	for (int i = 0; i < num_factors; i++) {
		divisors[i] = helper->p;
		exponents[i] = helper->exp;
		helper = helper->next;
	}

	free_factors(factors);

	/* TODO delete, for debugging */
	printf("%d factors\n", num_factors);
	for (int i = 0; i < num_factors; i++) {
		printf("%d^%d, ", divisors[i], exponents[i]);
	}
	printf("\n");

	/* TODO Call berlekamp with polynomial and each divisor to find roots */

	/* TODO Use Hensel's lemma to lift roots */

	/* TODO Use remainder theorem to solve systems of congruences with roots */

	free(exponents);
	free(divisors);

	return EXIT_SUCCESS;
}

/* --- utility functions -----------------------------------------------------*/

/** Allocates memory for and returns a factor which represents p^exp */
Factor *init_factor(int p, int exp)
{
	Factor *factor = malloc(sizeof(Factor));
	factor->p = p;
	factor->exp = exp;
	factor->next = NULL;
	return factor;
}

/** Use trial division to create a linked list which contains x's prime
 * factorization. Also count number of distinct prime factors and store in
 * variable pointed to by n. */
Factor *factorize(int *n, int x)
{
	/* Initialize head and tail of linked list, and number of factors */
	Factor *head = NULL, *tail = NULL;
	*n = 0;
	
	/* Count how many times 2 divides x, and update linked list if necessary */
	int p = 2, cnt = 0;
	while (x % 2 == 0) {
		cnt++;
		x /= 2;
	}
	if (cnt > 0) {
		head = init_factor(2, cnt);
		tail = head;
		(*n)++;
	}

	/* Now iterate over odd numbers and add them to list if they divide x */
	p = 3;
	cnt = 0;
	while (x > 1) {
		if (x % p == 0) {
			x /= p;
			cnt++;
		} else {
			if (cnt != 0) {
				if (tail) {
					tail->next = init_factor(p, cnt);
					tail = tail->next;
				} else {
					head = init_factor(p, cnt);
					tail = head;
				}
				(*n)++;
			}
			p += 2;
			cnt = 0;
		}
	}

	/* Add last factor */
	if (cnt != 0) {
		if (tail) {
			tail->next = init_factor(p, cnt);
			tail = tail->next;
		} else {
			head = init_factor(p, cnt);
			tail = head;
		}
		(*n)++;
	}

	return head;
}

/** Frees memory allocated to linked list of factors */
void free_factors(Factor *factors) {
	Factor *helper = factors;
	while (factors) {
		helper = factors;
		factors = factors->next;
		free(helper);
	}
}

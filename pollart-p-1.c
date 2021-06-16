/*****************************************************************************/
/* This program is free software; you can redistribute it and/or modify	     */
/* it under the terms of the GNU General Public License as published by	     */
/* the Free Software Foundation; either version 2 of the License, or	     */
/* (at your option) any later version.					     */
/* This program is distributed in the hope that it will be useful,	     */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of	     */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	     */
/* GNU General Public License for more details.				     */
/*									     */
/* You should have received a copy of the GNU General Public License	     */
/* along with this program; if not, write to the Free Software		     */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/*****************************************************************************/
/*	(c) 2003 by Paolo Ardoino  <paolo.ardoino@gmail.com>                 */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <gmp.h>

#define MAX_B 1000L  /* MAX b */
	
int main(int argc, char *argv[])
{
  float b = 0.;
  mpz_t N, a, GCD, tmp, k;
  struct timeval tm0, tm1;
  gmp_randstate_t state;
 
  printf("SPMI - Simple Pollard p-1 Method Implementaion.\n");
  printf("(c) 2003 by Paolo Ardoino  <paolo.ardoino@gmail.com>\n");
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <N> <b>\n", *argv);
    fprintf(stderr, "\t<N>: integer to factorize.\n");
    fprintf(stderr, "\t<b>: small integer for computation of k.\n");
    exit(EXIT_FAILURE);
  } else if (*(argv + 1) && *(argv + 2)){
    mpz_init_set_str(N, *(argv + 1), 10);
    if (mpz_cmp_ui(N, 1) <= 0) {
      fprintf(stderr, "Errot: Please insert N >= 2");
      exit(EXIT_FAILURE);
    }
    gmp_printf("N = %Zd\n", N);
    b = atof(*(argv + 2));
    if (b > MAX_B) {
      printf("Warning: b too large. Setting to %ld\n", MAX_B);
      b = (float)MAX_B;
    }
    printf("b = %.0f\n", b);
    mpz_init(tmp);
    gettimeofday(&tm0, NULL);

    /* Tries to compute m = N mod 2 */
    /* if m == 0 => 2|N [2 is a factor of N] */
    while (mpz_mod_ui(tmp, N, 2) == 0) {
      printf("Factor = 2\n");
      mpz_div_ui(N, N, 2);
    }

    /* Checks if N == 1 */
    if (mpz_cmp_ui(N, 1) == 0) {
      mpz_clear(tmp);
      mpz_clear(N);
      gettimeofday(&tm1, NULL);
      printf("Factorization has been completed in %ld seconds.\n",\
      tm1.tv_sec - tm0.tv_sec); 
      exit(EXIT_SUCCESS);
    }

    /* Checks if N is prime */
    /* Uses a probility primality test that has */
    /* probabity of failure == 0.25 ^ x [here x == 10] */
    if (mpz_probab_prime_p(N, 10) > 0) {
      gmp_printf("Factor = %Zd\n", N);
      mpz_clear(tmp);
      mpz_clear(N);
      gettimeofday(&tm1, NULL);
      printf("Factorization has been completed in %ld seconds.\n",\
      tm1.tv_sec - tm0.tv_sec); 
      exit(EXIT_SUCCESS);
    }
    if (mpz_perfect_power_p(N) != 0) {
      printf("N is a perfect root.\n");
      mpz_clear(tmp);
      mpz_clear(N);
      gettimeofday(&tm1, NULL);
      printf("Factorization has been completed in %ld seconds.\n",\
      tm1.tv_sec - tm0.tv_sec); 
      exit(EXIT_SUCCESS);
    }
    mpz_init(a);
    mpz_init(GCD);
    mpz_sub_ui(tmp, N, 1); /* tmp = N - 1 */
    mpz_init(k);
    mpz_fac_ui(k, b); /* k = b! */
    gmp_printf("k = %Zd\n", k);
    gmp_randinit_default(state);
    while (1) {
      mpz_sub_ui(tmp, N, 1);
      mpz_urandomm(a, state, tmp); /* 0 < a < N - 2 */
      if (mpz_cmp_ui(a, 1) <= 0)
        mpz_set_ui(a, 2);
      mpz_powm(tmp, a, k, N); /* computes a^k (mod(N)) */
      mpz_sub_ui(tmp, tmp, 1); /* a^k - 1 (mod(N)) */
      mpz_abs(tmp, tmp);
      mpz_gcd(GCD, tmp, N); /* GCD(a^k - 1 (mod(N)), N) */
      if (mpz_cmp_ui(GCD, 1) > 0) { /* GCD > 1 */
        if (mpz_probab_prime_p(GCD, 10) > 0) { /* GCD is prime */
	  gmp_printf("Factor = %Zd\n", GCD); /* GCD is a factor of N */
	  mpz_div(N, N, GCD);
	}
      }
      if (mpz_cmp_ui(N, 1) == 0) {
        mpz_clear(a);
        mpz_clear(GCD);
        mpz_clear(tmp);
        mpz_clear(N);
        mpz_clear(k);
        gettimeofday(&tm1, NULL);
        printf("Factorization has been completed in %ld seconds.\n",\
	tm1.tv_sec - tm0.tv_sec); 
        exit(EXIT_SUCCESS);
      }
      if (mpz_probab_prime_p(N, 10) > 0) {
        gmp_printf("Factor = %Zd\n", N);
        mpz_clear(a);
        mpz_clear(GCD);
        mpz_clear(tmp);
        mpz_clear(N);
        mpz_clear(k);
        gettimeofday(&tm1, NULL);
        printf("Factorization has been completed in %ld seconds.\n",\
	tm1.tv_sec - tm0.tv_sec); 
	exit(EXIT_SUCCESS);
      }
    }
  }
  return 0;
}


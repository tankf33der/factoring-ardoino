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

int main(int argc, char *argv[])
{
  unsigned long long int tmpui = 0;
  mpz_t N, sqrt, tmp, ctr;
  struct timeval tm0, tm1;
 
  if (argc != 2) {
    fprintf(stderr, "CFMI - Classical Factorization Method \
    Implementaion.\n");
    fprintf(stderr, "Usage: %s <N>\n", *argv);
    fprintf(stderr, "\t<N>: integer to factorize.\n");
    exit(EXIT_FAILURE);
  } else if (*(argv + 1)){
    mpz_init_set_str(N, *(argv + 1), 10);
    if (mpz_cmp_ui(N, 1) <= 0) {
      fprintf(stderr, "Errot: Please insert N >= 2");
      exit(EXIT_FAILURE);
    }
    gmp_printf("N = %Zd\n", N);
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
    mpz_init(sqrt);
    mpz_init_set_ui(ctr, 3); /* sets ctr = 3 */
    mpz_sqrt(sqrt, N);
    
    if (mpz_odd_p(sqrt) == 0)
      mpz_add_ui(sqrt, sqrt, 1);

    /* while ctr < sqrt(N) */
    while (mpz_cmp(ctr, sqrt) < 0) {
      while (1) {
        mpz_mod(tmp, N, ctr);
	if (mpz_cmp_ui(tmp, 0) == 0) {
          gmp_printf("Factor = %Zd\n", ctr);
	  mpz_div(N, N, ctr);
	  mpz_sqrt(sqrt, N);
         if (mpz_odd_p(sqrt) == 0)
           mpz_add_ui(sqrt, sqrt, 1);
	} else
	  break;
      }
      mpz_add_ui(ctr, ctr, 2);
      if (mpz_cmp_ui(N, 1) == 0) {
        mpz_clear(tmp);
	mpz_clear(N);
	mpz_clear(ctr);
	mpz_clear(sqrt);
        gettimeofday(&tm1, NULL);
        printf("Factorization has been completed in %ld seconds.\n",\
        tm1.tv_sec - tm0.tv_sec); 
	exit(EXIT_SUCCESS);
      }
      if (mpz_probab_prime_p(N, 10) > 0) {
        gmp_printf("Factor = %Zd\n", N);
        mpz_clear(tmp);
	mpz_clear(N);
	mpz_clear(ctr);
	mpz_clear(sqrt);
        gettimeofday(&tm1, NULL);
        printf("Factorization has been completed in %ld seconds.\n",\
        tm1.tv_sec - tm0.tv_sec); 
        exit(EXIT_SUCCESS);
      }
    }

  }
  return 0;
}

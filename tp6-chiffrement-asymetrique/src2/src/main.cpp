//
//  TP6_RSA
//

#include <stdio.h>
#include <iostream>
#include <gmp.h>

#define BITSTRENGTH 14              /* size of modulus (n) in bits */
#define PRIMESIZE (BITSTRENGTH / 2) /* size of the primes p and q  */

/* Declare global variables */

mpz_t d, e, n;
mpz_t M, c;

/* Main subroutine */
int main()
{
    /* Initialize the GMP integers */
    mpz_init(d);
    mpz_init(e);
    mpz_init(n);

    mpz_init(M);
    mpz_init(c);

    /* This function creates the keys. The basic algorithm is...
     *
     *  1. Generate two large distinct primes p and q randomly
     *  2. Calculate n = pq and x = (p-1)(q-1)
     *  3. Select a random integer e (1<e<x) such that gcd(e,x) = 1
     *  4. Calculate the unique d such that ed = 1(mod x)
     *  5. Public key pair : (e,n), Private key pair : (d,n)
     *
     */

    /*
     *  Step 1 : Get two large primes.
     */
    mpz_t p, q;
    mpz_init(p);
    mpz_init(q);

    mpz_init_set_str(p, "47", 0);
    mpz_init_set_str(q, "71", 0);
    char p_str[1000];
    char q_str[1000];
    mpz_get_str(p_str, 10, p);
    mpz_get_str(q_str, 10, q);

    std::cout << "Random Prime 'p' = " << p_str << std::endl;
    std::cout << "Random Prime 'q' = " << q_str << std::endl;

    /*
     *  Step 2 : Calculate n (=pq) ie the 1024 bit modulus
     *  and x (=(p-1)(q-1)).
     */
    char n_str[1000];
    mpz_t x;
    mpz_init(x);

    /* Calculate n... */
    mpz_mul(n, p, q);
    mpz_get_str(n_str, 10, n);
    std::cout << "\t n = " << n_str << std::endl;

    /* Calculate x... */
    mpz_t p_minus_1, q_minus_1;
    mpz_init(p_minus_1);
    mpz_init(q_minus_1);

    mpz_sub_ui(p_minus_1, p, (unsigned long int)1);
    mpz_sub_ui(q_minus_1, q, (unsigned long int)1);

    mpz_mul(x, p_minus_1, q_minus_1);
    char phi_str[1000];
    mpz_get_str(phi_str, 10, x);
    std::cout << "\t phi(n) = " << phi_str << std::endl;

    /*
     *  Step 3 : Get small odd integer e such that gcd(e,x) = 1.
     */
    mpz_init_set_str(e, "79", 0);
    char e_str[1000];
    mpz_get_str(e_str, 10, e);
    std::cout << "\t e = " << e_str << std::endl;

    /*
     *  Step 4 : Calculate unique d such that ed = 1(mod x)
     */
    mpz_init_set_str(d, "1019", 0);
    char d_str[1000];
    mpz_get_str(d_str, 10, d);
    std::cout << "\t d = " << d_str << std::endl
              << std::endl;

    /*
     *  Step 5 : Print the public and private key pairs...
     */
    std::cout << "Public Keys  (e,n): ( " << e_str << " , " << n_str << " )" << std::endl;
    std::cout << "Private Keys (d,n): ( " << d_str << " , " << n_str << " )" << std::endl;
    /*
     *  Encrypt
     */

    //TODO
    mpz_init_set_str(M, "11", 0);
    char c_str[1000];
    mpz_powm(c, M, e, n);
    mpz_get_str(c_str, 10, c);
    std::cout << "Ciphered message c: " << c_str << std::endl;

    /*
     * Decrypt
     */
    mpz_mul(x,p_minus_1,q_minus_1);
    mpz_invert(d,e,x);
    mpz_get_str(d_str, 10, d);

    std::cout << "d: " << d_str << std::endl;


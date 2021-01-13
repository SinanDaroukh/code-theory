//
//  TP6_RSA
//

#include <stdio.h>
#include <iostream>
#include <gmp.h>
#include <cmath>

#define BITSTRENGTH 14              /* size of modulus (n) in bits */
#define PRIMESIZE (BITSTRENGTH / 2) /* size of the primes p and q  */

/* Declare global variables */

mpz_t d, e, n;
mpz_t M, c;

void exponentiation_by_squaring(mpz_t &m, const mpz_t &g, const mpz_t &k, const mpz_t &p)
{
    mpz_t g_tmp, k_tmp, p_tmp;
    mpz_init(g_tmp);
    mpz_init(k_tmp);
    mpz_init(p_tmp);
    mpz_set(g_tmp, g);
    mpz_set(k_tmp, k);
    mpz_set(p_tmp, p);

    if (mpz_cmp_si(k_tmp, 0) < 0)
    { //si k<0
        mpz_t one;
        mpz_init(one);
        mpz_set_ui(one, 1);

        mpz_fdiv_q(g_tmp, one, g_tmp); //g=1/g
        mpz_mul_si(k_tmp, k_tmp, -1);  //k=k*(-1)
        mpz_clear(one);
    }
    if (mpz_cmp_si(k_tmp, 0) == 0)
    {
        mpz_set_ui(m, 1);
        return;
    }

    mpz_t y;
    mpz_init(y);
    mpz_set_ui(y, 1);

    // while k > 1
    while (mpz_cmp_si(k_tmp, 1) > 0)
    {
        // if k is even
        if (mpz_even_p(k_tmp) != 0)
        {
            // g = g * g
            mpz_mul(g_tmp, g_tmp, g_tmp);
            // g = g mod p
            mpz_mod(g_tmp, g_tmp, p_tmp);
            // k = k / 2
            mpz_fdiv_q_ui(k_tmp, k_tmp, 2);
        }
        else
        {
            // y = g * y
            mpz_mul(y, g_tmp, y);
            // g = g * g
            mpz_mul(g_tmp, g_tmp, g_tmp);
            // k = (k - 1) / 2
            mpz_sub_ui(k_tmp, k_tmp, 1);
            mpz_fdiv_q_ui(k_tmp, k_tmp, 2);
        }
    }

    // m = g * y
    // std::cout << "g: " << mpz_get_ui(g_tmp) << std::endl;
    // std::cout << "y: " << mpz_get_ui(y) << std::endl;

    mpz_mul(m, g_tmp, y);
    mpz_mod(m, m, p_tmp);
}

void rabin_miller_step_1_dif(mpz_t &s, mpz_t &t, const mpz_t &n)
{
    mpz_t local_n;
    char t_str[1000];
    char s_str[1000];
    mpz_init(local_n);

    mpz_set(t, n);
    mpz_sub_ui(t, t, 1);
    mpz_set_ui(s, 0);

    while (mpz_even_p(t))
    {
        mpz_fdiv_q_ui(t, t, 2);
        mpz_add_ui(s, s, 1);
        mpz_get_str(t_str, 10, t);
        mpz_get_str(s_str, 10, s);
        // std::cout << " t : " << t_str << std::endl;
        // std::cout << " s : " << s_str << std::endl;
    }
}

bool rabin_miller(int k, const mpz_t &n)
{
    // std::cout << "Rabin Miller" << std::endl;

    char n_str[1000];
    mpz_get_str(n_str, 10, n);
    // std::cout << "n: " << n_str << std::endl;

    if (mpz_get_si(n) <= 2 && mpz_odd_p(n))
        return false;

    gmp_randstate_t state;
    mpz_t t, s, a, x, r, tmp_oprd;

    char t_str[1000];
    char s_str[1000];

    mpz_init(t);
    mpz_init(s);
    mpz_init(a);
    mpz_init(x);
    mpz_init(r);
    mpz_init(tmp_oprd);
    rabin_miller_step_1_dif(s, t, n);

    mpz_get_str(t_str, 10, t);
    mpz_get_str(s_str, 10, s);
    // std::cout << "s: " << s_str << std::endl;
    // std::cout << "t: " << t_str << std::endl;
    mpz_t seed;
    mpz_init_set_str(seed, std::to_string(1000 + rand() % 100000).c_str(), 0);
    gmp_randinit_default(state);
    gmp_randseed(state, seed);

    for (int i = 0; i < k; ++i)
    {
        mpz_sub_ui(tmp_oprd, n, 2);
        mpz_urandomm(tmp_oprd, state, tmp_oprd);
        mpz_add_ui(a, tmp_oprd, 2);
        mpz_powm(x, a, t, n);
        mpz_sub_ui(tmp_oprd, n, 1);
        if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, tmp_oprd) == 0)
        {
            goto loop;
        }

        for (mpz_set_ui(r, 1); mpz_cmp(r, s) < 0; mpz_add_ui(r, r, 1))
        {
            // x = (x * x) % n;
            mpz_mul(x, x, x);
            mpz_mod(x, x, n);
            // std::cout << "x: " << x << std::endl;
            if (mpz_cmp_ui(x, 1) == 0)
                return false;
            mpz_sub_ui(tmp_oprd, n, 1);
            if (mpz_cmp(x, tmp_oprd) == 0)
                goto loop;
        }
        return false;

    loop:
        continue;
    }
    return true;
}

void nextprime(mpz_t &rop, const mpz_t op)
{
    mpz_set(rop, op);
    if (mpz_even_p(op) != 0)
    {
        mpz_add_ui(rop, op, 1);
    } // If op is even
    while (rabin_miller(100000, rop) != true)
    {
        mpz_add_ui(rop, rop, 2);
    }
}

void euclide_etendu(const mpz_t a, const mpz_t b, mpz_t &r, mpz_t &u, mpz_t &v)
{
    mpz_t r_prime, u_prime, v_prime, q, rs, us, vs, q_mult_r_prime, q_mult_u_prime, q_mult_v_prime;

    mpz_init_set(r, a);
    mpz_init_set_ui(u, 1);
    mpz_init(v);
    mpz_init_set(r_prime, b);
    mpz_init(u_prime);
    mpz_init_set_ui(v_prime, 1);
    mpz_init(rs);
    mpz_init(us);
    mpz_init(vs);
    mpz_init(q);
    mpz_init(q_mult_r_prime);
    mpz_init(q_mult_u_prime);
    mpz_init(q_mult_v_prime);

    while (mpz_cmp_ui(r_prime, 0) != 0)
    {
        mpz_fdiv_q(q, r, r_prime);

        mpz_set(rs, r);
        mpz_set(us, u);
        mpz_set(vs, v);

        mpz_set(r, r_prime);
        mpz_set(u, u_prime);
        mpz_set(v, v_prime);

        mpz_mul(q_mult_r_prime, q, r_prime);
        mpz_sub(r_prime, rs, q_mult_r_prime);

        mpz_mul(q_mult_u_prime, q, u_prime);
        mpz_sub(u_prime, us, q_mult_u_prime);

        mpz_mul(q_mult_v_prime, q, v_prime);
        mpz_sub(v_prime, vs, q_mult_v_prime);
    }

    mpz_clear(r_prime);
    mpz_clear(u_prime);
    mpz_clear(v_prime);
    mpz_clear(q);
    mpz_clear(q_mult_r_prime);
    mpz_clear(q_mult_u_prime);
    mpz_clear(q_mult_v_prime);
}

void invert(mpz_t &rop, const mpz_t op1, const mpz_t op2)
{
    mpz_t r, u, v;

    char r_str[1000], u_str[1000], v_str[1000];

    mpz_init(r);
    mpz_init(u);
    mpz_init(v);

    euclide_etendu(op1, op2, r, u, v);

    mpz_get_str(r_str, 10, r);
    mpz_get_str(u_str, 10, u);
    mpz_get_str(v_str, 10, v);

    std::cout << "r: " << r_str << std::endl;
    std::cout << "u: " << u_str << std::endl;
    std::cout << "v :" << v_str << std::endl;

    mpz_abs(rop, u);

    mpz_clear(r);
    mpz_clear(u);
    mpz_clear(v);
}

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
    mpz_t p, p_tmp, q, q_tmp, seed;
    gmp_randstate_t state;

    char p_str[1000];
    char p_tmp_str[1000];
    char q_str[1000];
    char q_tmp_str[1000];

    mpz_init(p);
    mpz_init(p_tmp);
    mpz_init(q);
    mpz_init(q_tmp);
    mpz_init(seed);

    // mpz_init_set_str(p, "47", 0);
    // mpz_init_set_str(q, "71", 0);

    // -- Generation and randomisation of a prime number p
    srand(time(NULL)); // -- Initialization of a seed based on time

    mpz_init_set_str(seed, std::to_string(1000 + rand() % 100000).c_str(), 0);

    gmp_randinit_default(state);
    gmp_randseed(state, seed);

    mpz_urandomm(p_tmp, state, seed);
    mpz_urandomm(q_tmp, state, seed);

    // mpz_nextprime(p, p_tmp);
    // mpz_nextprime(q, q_tmp);
    nextprime(p, p_tmp);
    nextprime(q, q_tmp);

    mpz_get_str(p_tmp_str, 10, p_tmp);
    mpz_get_str(p_str, 10, p);

    mpz_get_str(q_tmp_str, 10, q_tmp);
    mpz_get_str(q_str, 10, q);

    std::cout << "**************************************" << std::endl;
    std::cout << "Random 'p' tmp = " << p_tmp_str << std::endl;
    std::cout << "Random Prime 'p' = " << p_str << std::endl;
    std::cout << "Random 'q' tmp = " << q_tmp_str << std::endl;
    std::cout << "Random Prime 'q' = " << q_str << std::endl;
    std::cout << "**************************************" << std::endl;

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

    mpz_t e_tmp, pgcd;
    mpz_init(e_tmp);
    mpz_init(pgcd);

    char e_str[1000];

    srand(time(NULL));
    do
    {
        mpz_init_set_str(seed, std::to_string(rand()).c_str(), 0);

        gmp_randinit_default(state);
        gmp_randseed(state, seed);

        mpz_urandomm(e_tmp, state, seed);
        mpz_init_set_str(e, std::to_string(mpz_get_ui(e_tmp) % mpz_get_ui(x)).c_str(), 0);

        mpz_gcd(pgcd, e, x);
    } while (mpz_get_ui(pgcd) != 1);

    // mpz_init_set_str(e, "79", 0);
    mpz_get_str(e_str, 10, e);
    std::cout << "\t e = " << e_str << std::endl;

    /*
     *  Step 4 : Calculate unique d such that ed = 1(mod x)
     */
    // if (mpz_invert(d, e, x) == 0)
    //     std::cout << "Error while inversing e mod x" << std::endl;

    char x_str[1000];
    mpz_get_str(x_str, 10, x);
    std::cout << "\t x = " << x_str << std::endl;

    invert(d, e, x);

    // mpz_init_set_str(d, "1019", 0);
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
    mpz_t m;
    mpz_t c;

    mpz_init(m);
    mpz_init(c);

    char m_str[1000];
    char c_str[1000];

    std::cout << "Message à chiffrer " << std::endl;
    std::cin >> m_str;
    // TODO: check m < n
    mpz_set_str(m, m_str, 10);
    // mpz_powm(c, m, e, n);
    exponentiation_by_squaring(c, m, e, n);
    mpz_get_str(c_str, 10, c);
    std::cout << "Message chiffré: " << c_str << std::endl;

    /*
     *  Decrypt
     */
    mpz_t c2;
    mpz_t m2;

    mpz_init(m2);
    mpz_init(c2);

    char m2_str[1000];
    char c2_str[1000];

    std::cout << "Message à déchiffrer ? " << std::endl;
    std::cin >> c2_str;
    mpz_set_str(c2, c2_str, 10);
    // mpz_powm(m2, c2, d, n);
    exponentiation_by_squaring(m2, c2, d, n);
    mpz_get_str(m2_str, 10, m2);
    std::cout << "Message déchiffré: " << m2_str << std::endl;

    /* Clean up the GMP integers */
    mpz_clear(p_minus_1);
    mpz_clear(q_minus_1);
    mpz_clear(x);
    mpz_clear(p);
    mpz_clear(q);

    mpz_clear(d);
    mpz_clear(e);
    mpz_clear(n);

    mpz_clear(M);
    mpz_clear(c);
}

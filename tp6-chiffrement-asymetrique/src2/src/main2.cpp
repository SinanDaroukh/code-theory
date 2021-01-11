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

/*
Détermine u, v tels que a*u+b*v =PGCD(a,b)
Inspiré de https://fr.wikipedia.org/wiki/Algorithme_d%27Euclide_%C3%A9tendu
*/
void Euclide_Etendu(mpz_t u, mpz_t v, const mpz_t a, const mpz_t b)
{
    mpz_t r, r_bis, u_bis, v_bis, q;
    mpz_t rs, us, vs, mul1, mul2, mul3; //Variables de stockage intermédiaires

    //Initialisation
    mpz_init(r);
    mpz_init(r_bis);
    mpz_init(u_bis);
    mpz_init(v_bis);
    mpz_init(q);
    mpz_init(rs);
    mpz_init(us);
    mpz_init(vs);
    mpz_init(mul1);
    mpz_init(mul2);
    mpz_init(mul3);

    //Initialisation
    mpz_set(r, a);
    mpz_set(r_bis, b);
    mpz_set_ui(u, 1);
    mpz_set_ui(v, 0);
    mpz_set_ui(u_bis, 0);
    mpz_set_ui(v_bis, 1);

    while (mpz_cmp_ui(r_bis, 0) != 0)
    {
        mpz_fdiv_q(q, r, r_bis);
        mpz_set(rs, r);
        mpz_set(us, u);
        mpz_set(vs, v);
        mpz_set(r, r_bis);
        mpz_set(u, u_bis);
        mpz_set(v, v_bis);
        mpz_mul(mul1, q, r_bis);
        mpz_sub(r_bis, rs, mul1);
        mpz_mul(mul2, q, u_bis);
        mpz_sub(u_bis, us, mul2);
        mpz_mul(mul3, q, v_bis);
        mpz_sub(v_bis, vs, mul3);
    }

    // Nettoyage
    mpz_clear(r);
    mpz_clear(r_bis);
    mpz_clear(u_bis);
    mpz_clear(v_bis);
    mpz_clear(q);
    mpz_clear(rs);
    mpz_clear(us);
    mpz_clear(vs);
    mpz_clear(mul1);
    mpz_clear(mul2);
    mpz_clear(mul3);
}

/* Détermine l'entier d tel que ed = 1(mod x) */
int Inverse_modulaire(mpz_t d, const mpz_t e, const mpz_t x)
{
    mpz_t r, u, v;

    mpz_init(r);
    mpz_init(u);
    mpz_init(v);

    Euclide_Etendu(u, v, e, x);

    if (mpz_cmp_ui(u, 0) < 0) // Si Euclide_Etendu renvoie un négatif on rajoute x au résultat pour le rendre positif
    {
        mpz_add(u, u, x);
    }
    mpz_set(d, u);

    mpz_clear(r);
    mpz_clear(u);
    mpz_clear(v);

    return 1;
}

//fonction de test de primalite
//https://fr.wikipedia.org/wiki/Test_de_primalit%C3%A9_de_Miller-Rabin
bool rabinMiller(mpz_t k, mpz_t n)
{ //n => le nombre a tester , k => le nombre de repetition de test pour la detection du temoin de miller
    if (mpz_get_si(n) > 2)
    {

        mpz_t a;
        mpz_t x;
        mpz_t d;

        mpz_init(a); //temoin de miller
        mpz_init(x);
        mpz_init(d);

        mpz_t exp;
        mpz_init(exp);
        mpz_set_si(exp, 2);

        //t <-- n-1
        long t = mpz_get_si(n) - 1;
        //calcul de s et d tel que t=pow(2,s)*d
        unsigned s = 0;
        while (t % 2 == 0)
        {
            t = t / 2;
            s++;
        }
        mpz_set_si(d, t); //initialiser la valeur de d ==> maintenant nous avons s et d tel que n-1=pow(2,s)*d

        for (unsigned i = 0; i < mpz_get_ui(k); i++)
        { //tester k temoin de miller => si on trouve un seul temoin on sort
            //tirer a aleatoirement dans l'intervalle [2,n-1]
            mpz_set_si(a, (rand() % (mpz_get_si(n) - 1)) + 2);
            //calculer  x=pow(a,d) mod n
            mpz_powm(x, a, d, n);
            if ((mpz_get_si(x) == 1) || (mpz_get_si(x) == mpz_get_si(n) - 1))
            {
                // a n'est pas un temoin de miller, tirer un autre nombre a
                continue;
            }

            for (unsigned r = 1; r < s; r++)
            {
                mpz_powm(x, x, exp, n); //x=pow(x,2) mod n
                if (mpz_get_si(x) == 1)
                {
                    mpz_clear(a);
                    mpz_clear(d);
                    mpz_clear(x);
                    return false; // a est un temoin de miller donc n est composite
                }
                if (mpz_get_si(x) == mpz_get_si(n) - 1)
                { //si x=n-1
                    goto LOOP;
                }
            }
            mpz_clear(a);
            mpz_clear(d);
            mpz_clear(x);
            return false; //n est composite
        LOOP:
            continue;
        }
        mpz_clear(a);
        mpz_clear(d);
        mpz_clear(x);
        return true; //n est probablement premier
    }
    return false; //n est composite
}

//fonction de calcul de GCD
unsigned long GCD(mpz_t a, mpz_t b)
{
    if (mpz_cmp_ui(b, 0) == 0)
    {
        return mpz_get_ui(a);
    }
    else
    {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mod(tmp, a, b); //tmp=a mod b
        return GCD(b, tmp); //appel recursif
    }
}

//fonction d'exponentiation by squaring pour le calcul de m=pow(g,k) mod p
unsigned long EBS(mpz_t &m, mpz_t g, mpz_t k, mpz_t p)
{
    if (mpz_cmp_si(k, 0) < 0)
    { //si k<0
        mpz_t tmp;
        mpz_init(tmp);
        mpz_set_ui(tmp, 1);

        mpz_fdiv_q(g, tmp, g); //g=1/g
        mpz_mul_si(k, k, -1);  //k=k*(-1)
        mpz_clear(tmp);
    }
    if (mpz_cmp_si(k, 0) == 0)
    {
        mpz_set_si(m, 1); //m=1
        return mpz_get_ui(m);
    }
    mpz_t y;
    mpz_init(y);
    mpz_set_si(y, 1); //y=1

    mpz_t exp;
    mpz_init(y);
    mpz_set_si(exp, 2);

    while (mpz_cmp_si(k, 1) > 0)
    { //tant que  k>1
        if (mpz_even_p(k) != 0)
        {                                     //si k est pair
            mpz_powm(g, g, exp, p);           //g=g*g mod p
            mpz_set_si(k, mpz_get_si(k) / 2); //k=k/2
        }
        else
        {                                           //si k est impair
            mpz_mul(y, g, y);                       //y=g*y
            mpz_mul(g, g, g);                       //g=g*g
            mpz_set_si(k, (mpz_get_si(k) - 1) / 2); //(k-1)/2
        }
    }

    mpz_mul(m, g, y); //m=g*y
    mpz_mod(m, m, p); //m=m mod p
    mpz_clear(y);
    mpz_clear(exp);
    return mpz_get_ui(m);
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

    //******************************************Step 1 : Get two large primes.*******************************************
    mpz_t p, q;
    mpz_init(p);
    mpz_init(q);

    //***************************************************Utilisation de GMP

    /*
    mpz_t op;
    mpz_init(op);
    char op_str[1000];
   	char p_str[1000];
  	char q_str[1000];

	std::cout<<"Donnez une valeur de op1:"<<std::endl;
	std::cin>>op_str;
	mpz_set_str (op, op_str,10);
	mpz_nextprime(p, op);
   	mpz_get_str(p_str,10,p);

	std::cout<<"Donnez une valeur de op2:"<<std::endl;
	std::cin>>op_str;
	mpz_set_str (op, op_str,10);
	mpz_nextprime(q, op);
   	mpz_get_str(q_str,10,q);

    std::cout << "Random Prime 'p' = " << p_str <<  std::endl;
  	std::cout << "Random Prime 'q' = " << q_str <<  std::endl;
    */

    //*********************************************utilisation de nos propres fonctions
    //generate random numbers
    mpz_t taille;
    mpz_init(taille);
    mpz_set_ui(taille, 100000000000);

    int seed = time(NULL);
    srand(seed);

    //Génération de p
    //initialisation du nombre d'iteration k
    mpz_t k;
    mpz_init(k);
    mpz_set_ui(k, 10);

    while ((mpz_cmp_ui(p, 2) <= 0) || (mpz_even_p(p) != 0) || (rabinMiller(k, p) == false))
    { //TQ le nombre p est inferieur a 2, paire ou n'est pas premier
        mpz_set_si(p, rand() % mpz_get_ui(taille));
    }

    while ((mpz_cmp(p, q) == 0) || (mpz_cmp_ui(q, 2) <= 0) || (mpz_even_p(q) != 0) || (rabinMiller(k, q) == false))
    { //TQ le nombre q est égale à p, est inferieur a 2, paire ou n'est pas premier
        mpz_set_si(q, rand() % mpz_get_ui(taille));
    }

    unsigned long firstPrime = mpz_get_ui(p);
    unsigned long secondPrime = mpz_get_ui(q);

    std::cout << "Random Prime 'p' = " << firstPrime << std::endl;
    std::cout << "Random Prime 'q' = " << secondPrime << std::endl;

    //***********************************************Step 2 : Calculate n (=pq) ie the 1024 bit modulus*********************
    //  and x (=(p-1)(q-1)).

    char n_str[1000];
    mpz_t x;
    mpz_init(x);

    // Calculate n...
    mpz_mul(n, p, q);
    mpz_get_str(n_str, 10, n);
    std::cout << "\t n = " << n_str << std::endl;

    // Calculate x
    mpz_t p_minus_1, q_minus_1;
    mpz_init(p_minus_1);
    mpz_init(q_minus_1);

    mpz_sub_ui(p_minus_1, p, (unsigned long int)1);
    mpz_sub_ui(q_minus_1, q, (unsigned long int)1);

    mpz_mul(x, p_minus_1, q_minus_1);
    unsigned long phi = mpz_get_ui(x);
    char phi_str[1000];
    mpz_get_str(phi_str, 10, x);
    std::cout << "\t phi(n) = " << phi_str << std::endl;

    //****************************Step 3 : Get small odd integer e such that gcd(e,x) = 1.**********************************

    //***************************************************Utilisation de GMP
    /*
    mpz_init_set_str(e, "79", 0);
    char e_str[1000];
    mpz_get_str(e_str,10,e);
    std::cout << "\t e = " << e_str << std::endl;
    */

    //***************************utilisation de nos propres fonctions

    mpz_set_ui(taille, phi); //entier naturel strictement inferieur a phi(n)
    while (GCD(e, x) != 1)
    { //TQ e et x ne sont pas premiers entre eux
        mpz_set_si(e, rand() % mpz_get_ui(taille));
    }
    char e_str[1000];
    mpz_get_str(e_str, 10, e);
    std::cout << "\t e = " << e_str << std::endl;

    //*********************  Step 4 : Calculate unique d such that ed = 1(mod x)**************************

    //***************************************************Utilisation de GMP
    /*
    mpz_invert(d, e, x);
    char d_str[1000];
    mpz_get_str(d_str,10,d);
    std::cout << "\t d = " << d_str << std::endl << std::endl;
    */

    //***************************utilisation de nos propres fonctions

    Inverse_modulaire(d, e, x);
    char d_str[1000];
    mpz_get_str(d_str, 10, d);
    std::cout << "\t d = " << d_str << std::endl
              << std::endl;

    //******************** *  Step 5 : Print the public and private key pairs...******************************

    std::cout << "Public Keys  (e,n): ( " << e_str << " , " << n_str << " )" << std::endl;
    std::cout << "Private Keys (d,n): ( " << d_str << " , " << n_str << " )" << std::endl;

    //***************************************** Encrypt  **************************************************

    mpz_t m;
    mpz_t c;

    char m_str[1000];
    char c_str[1000];
    std::cout << "Message à chiffrer " << std::endl;
    std::cin >> m_str;
    mpz_set_str(m, m_str, 10);
    mpz_powm(c, m, e, n);
    mpz_get_str(c_str, 10, c);
    std::cout << "Message chiffré: " << c_str << std::endl;

    //*****************************************  Decrypt   ***********************************************

    mpz_t c2;
    mpz_t m2;
    char m2_str[1000];
    char c2_str[1000];
    std::cout << "Message à déchiffrer ? " << std::endl;
    std::cin >> c2_str;
    mpz_set_str(c2, c2_str, 10);

    mpz_powm(m2, c2, d, n);

    mpz_get_str(m2_str, 10, m2);
    std::cout << "Message déchiffré: " << m2_str << std::endl;

    // ****************************************** Clean up the GMP integers  *****************************
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

#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void E( mpz_t c, mpz_t m, mpz_t g, mpz_t n ) {
    mpz_t r;
    mpz_t n2;
    
    mpz_init( r );
    mpz_init( n2 );
    
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, time(NULL) );
    
    mpz_mul( n2, n, n );
    
    mpz_urandomm( r, rstate, n );
    mpz_powm( c, g, m, n2 );
    mpz_powm( r, r, n, n2 );
    mpz_mul( c, c, r );    
    mpz_mod( c, c, n2 );     // c1 = (g^m1)(r^n) mod n^2
}

void D( mpz_t m, mpz_t c, mpz_t lambda, mpz_t micro, mpz_t n ) {
    mpz_t r;
    mpz_t n2;
    mpz_t L;
    
    mpz_init( n2 );
    mpz_init( L );
    
    mpz_mul( n2, n, n );
    
    mpz_powm( L, c, lambda, n2 );
    mpz_sub_ui( L, L, 1 );
    mpz_div( L, L, n );
    
    mpz_mul( m, L, micro );
    mpz_mod( m, m, n );
    
}
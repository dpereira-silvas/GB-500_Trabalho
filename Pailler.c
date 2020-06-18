#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


int main() {
    int nbits = 2048;
    
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, time(NULL) );
    
    
    mpz_t p;
    mpz_t q;
    mpz_t n;
    mpz_t n2;
    mpz_t g;
    mpz_t lambda;
    mpz_t micro;
    
    
    mpz_init( p );
    mpz_init( q );
    mpz_init( n );
    mpz_init( n2 );
    mpz_init( g );
    mpz_init( lambda );
    mpz_init( micro );
    
    mpz_urandomb( p, rstate, nbits/2 );
    mpz_urandomb( q, rstate, nbits/2 );
    
    mpz_nextprime( p, p );
    mpz_nextprime( q, q );
        
    mpz_mul( n, p, q );    // n = p * q

    mpz_sub_ui( p, p, 1 ); // p = p - 1
    mpz_sub_ui( q, q, 1 ); // q = q - 1
    
    mpz_mul( lambda, p, q ); // lambda = (p-1)(q-1)
    
    
    mpz_invert( micro, lambda, n );  // micro = lambda^-1 mod n
    
    mpz_add_ui( g, n , 1 );  // g = n + 1 
    
    mpz_mul( n2, n, n );
    
    // publica : (n,g)
    // privada : (lambda, micro)

    FILE *arq;
    arq = fopen("./keys/pubkey.txt","w");
    gmp_fprintf(arq,"%ZX\n",n);
    gmp_fprintf(arq,"%ZX\n",g);

    fclose(arq);

    arq = fopen("./keys/privkey.txt","w");
    gmp_fprintf(arq,"%ZX\n",lambda);
    gmp_fprintf(arq,"%ZX\n",micro);


    fclose(arq);

    mpz_clear( n );
    mpz_clear( p );
    mpz_clear( q );
    mpz_clear( n2 );
    mpz_clear( g );
    mpz_clear( lambda );
    mpz_clear( micro );
    
    
    return 0;
}

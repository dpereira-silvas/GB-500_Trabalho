#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Pailler.h"
#include "utils.h"

int main(int argc, char **argv)
{
    FILE *arq;

    mpz_t n;
    mpz_t n2;
    mpz_t g;
    mpz_t lambda;
    mpz_t micro;
    
    mpz_init( n );
    mpz_init( n2 );
    mpz_init( g );
    mpz_init( lambda );
    mpz_init( micro );

    // Lendo as chaves publica e privadas dos arquivos
    arq = fopen("./keys/pubkey.txt","r");
    gmp_fscanf(arq,"%ZX\n",n);
    gmp_fscanf(arq,"%ZX\n",g);

    mpz_mul( n2, n, n );

    fclose(arq);

    arq = fopen("./keys/privkey.txt","r");

    gmp_fscanf(arq,"%ZX\n",lambda);
    gmp_fscanf(arq,"%ZX\n",micro);

    fclose(arq);
    // publica : (n,g)
    // privada : (lambda, micro)


    int dx;
    int dt;
    int n_nodes;
    int t_final;
    int mul;

    arq = fopen("./in/param.txt","r");
    
    fscanf(arq,"%d\n",&n_nodes);
    fscanf(arq,"%d\n",&dx);
    fscanf(arq,"%d\n",&dt);
    fscanf(arq,"%d\n",&t_final);
    fscanf(arq,"%d\n",&mul);

    fclose(arq);



    mpz_t *Uant;
    Uant = malloc(sizeof(mpz_t) * n_nodes);

    // Condições de contorno

    int cc;
    arq = fopen("./in/CC.txt","r");
    
    fscanf(arq,"%d\n",&cc);
    mpz_init_set_ui(Uant[0],cc);

    fscanf(arq,"%d\n",&cc);
    mpz_init_set_ui(Uant[n_nodes-1],cc);

    fclose(arq);

    arq = fopen("./in/IC.txt","w");

    mpz_t tmp;
    mpz_init(tmp);
    E( tmp, Uant[0], g, n );
    gmp_fprintf (arq, "%ZX\n",tmp);

    // Condição inicial
    for(int i = 1; i < (n_nodes-1); i++)
    {   
        //Exemplo 1
        mpz_init_set_str(Uant[i],"20000000",10); // 20.000000
        E( tmp, Uant[i], g, n );
        gmp_fprintf (arq, "%ZX\n",tmp);

    }
    E( tmp, Uant[n_nodes-1], g, n );
    gmp_fprintf (arq, "%ZX\n",tmp);

    clear_vector(Uant,n_nodes);


    mpz_clear(n);
    mpz_clear(n2);
    mpz_clear(g);
    mpz_clear(lambda);
    mpz_clear(micro);

    fclose(arq);


    return 0;
}
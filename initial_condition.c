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


void init_vector(mpz_t **v, int n)
{
    (*v) = malloc(sizeof(mpz_t) * n);

    // inicializando vetor
    for(int i = 0; i < n; i++)
    {
        mpz_init((*v)[i]); 
    }
    
}

void print(mpz_t *v, int n)
{
    for(int i = 0; i < (n-1); i++)
    {
        printf("%f ",(1/1000000.0)*mpz_get_ui(v[i] ));
    }
    printf("%f\n",(1/1000000.0)*mpz_get_ui(v[(n-1)]));
}


void write_vector(mpz_t *v, int n, FILE *f,int op)
{

    switch (op)
    {
        case 0:
            for(int i = 0; i < (n-1); i++)
            {
                fprintf (f, "%f\t",(1/1000000.0)*mpz_get_ui(v[i] ));
            }
            fprintf(f,"%f\n",(1/1000000.0)*mpz_get_ui(v[(n-1)]));  
            break;
        case 1:
            for(int i = 0; i < (n-1); i++)
            {
                fprintf (f, "%ld\n",mpz_get_ui(v[i]));
            }
            fprintf(f,"%ld\n",mpz_get_ui(v[(n-1)]));
            break;
    
    default:
        break;
    }
    
}


void clear_vector(mpz_t *v, int n)
{
    for(int i = 0; i < n; i++)
    {
        mpz_clear(v[i]);
    }
    free(v);

}

void copy_vector(mpz_t *u, mpz_t *v, int n)
{
    for(int i = 0; i < n; i++)
    {
        mpz_set(u[i], v[i]);
    }
}


int main(int argc, char **argv)
{
    if(argc < 2)
    {
        printf("./initial_condition <n_nodes>\n");
        return 0;
    }
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


    int dx = 5;
    int dt = 5;
    // Exemplo 1
    int n_nodes = atoi(argv[1]);
 
    // Exemplo 2
    // int n_nodes = 7;

    int t_final = 500;

    // numero de Fourier 1/mul
    int mul = 5;

    mpz_t *Uant;
    Uant = malloc(sizeof(mpz_t) * n_nodes);

    // Condições de contorno
    // Exemplo 1
    int cc;
    arq = fopen("./in/CC.txt","r");
    
    fscanf(arq,"%d\n",&cc);
    mpz_init_set_ui(Uant[0],cc);

    fscanf(arq,"%d\n",&cc);
    mpz_init_set_ui(Uant[n_nodes-1],cc);

    fclose(arq);


    
    // Exemplo 2
    // mpz_init_set_str(Uant[0],"20000000",10);
    // mpz_init_set_str(Uant[n_nodes-1],"50000000",10);



    arq = fopen("./in/IC.txt","w");

    mpz_t tmp;
    mpz_init(tmp);
    E( tmp, Uant[0], g, n );
    gmp_fprintf (arq, "%ZX\n",tmp);
    // gmp_fprintf (arq, "%ZX\n",Uant[0]);
    // Condição inicial
    for(int i = 1; i < (n_nodes-1); i++)
    {   
        //Exemplo 1
        mpz_init_set_str(Uant[i],"20000000",10); // 20.000000
        // gmp_fprintf (arq, "%ZX\n",Uant[i]);
        // Exemplo 2
        // mpz_init_set_ui(Uant[i],(60-2*(i*dx))*1000000  ); // 60 -2*x
        E( tmp, Uant[i], g, n );
        gmp_fprintf (arq, "%ZX\n",tmp);

    }
    E( tmp, Uant[n_nodes-1], g, n );
    gmp_fprintf (arq, "%ZX\n",tmp);
    // gmp_fprintf (arq, "%ZX\n",Uant[n_nodes-1]);

    clear_vector(Uant,n_nodes);


    mpz_clear(n);
    mpz_clear(n2);
    mpz_clear(g);
    mpz_clear(lambda);
    mpz_clear(micro);

    fclose(arq);


    return 0;
}
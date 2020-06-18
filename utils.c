#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>

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
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
                fprintf (f, "%ld\t",mpz_get_ui(v[i]));
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

// MEF homomorfico no cenario 2
double homomorphic_sc1() {
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
    // gmp_printf("%ZX\n",n);
    // gmp_printf("%ZX\n",g);
    // gmp_printf("%ZX\n",lambda);
    // gmp_printf("%ZX\n",micro);

    

    
    int dx = 5;
    int dt = 5;

    // Exemplo 1
    int n_nodes;

    // Exemplo 2
    // int n_nodes = 7;

    int t_final;

    // numero de Fourier 1/mul
    int mul;

    // Lendo os parametros do problema

    arq = fopen("./in/param.txt","r");
    
    fscanf(arq,"%d\n",&n_nodes);
    fscanf(arq,"%d\n",&dx);
    fscanf(arq,"%d\n",&dt);
    fscanf(arq,"%d\n",&t_final);
    fscanf(arq,"%d\n",&mul);

    fclose(arq);




    mpz_t *U;
    init_vector(&U, n_nodes);

    mpz_t *C;
    init_vector(&C, n_nodes);

    mpz_t *U_dec;
    init_vector(&U_dec, n_nodes);

    arq = fopen("./in/IC.txt","r");

    // Cifrando a condição inicial
    for(int i = 0; i < n_nodes; i++)
    {
        gmp_fscanf(arq,"%ZX\n",C[i]);
    }
    fclose(arq);

    // for(int i = 0; i < n_nodes; i++)
    // {
    //     D( U_dec[i], C[i], lambda, micro, n );

    // }

    // print(U_dec,n_nodes);

    arq = fopen("./out/out_sc1.txt","w");

    // write_vector(C, n_nodes, arq,0);

    // passando para o vetor U as condições de contorno criptografadas
    mpz_set(U[0],C[0]);
    mpz_set(U[n_nodes-1],C[n_nodes-1]);

    mpz_t tmp;
    mpz_init(tmp);

    int j = dt; // Iteração começa do instante de tempo 0+dt
    
    clock_t time_1, time_2;
	time_1 = clock();

    while(j <= t_final)
    {
        /////   SERVIDOR
        for(int i = 1; i < (n_nodes-1); i++)
        {
            mpz_powm_ui( tmp, C[i], 2, n2 );      // 
            mpz_invert(tmp, tmp, n2);            //  -2*Uant[i]

            mpz_powm_ui( U[i], C[i], mul, n2 );  // 5*Uant[i]

            mpz_mul( U[i], U[i], C[i+1] );       //  U[i] = U[i] + Uant[i+1]
            mpz_mod( U[i], U[i], n2 );

            mpz_mul( U[i], U[i], C[i-1] );       //  U[i] = U[i] + Uant[i-1]
            mpz_mod( U[i], U[i], n2 );

            mpz_mul( U[i], U[i], tmp );          //  U[i] = U[i] -  2*Uant[i]
            mpz_mod( U[i], U[i], n2 );
        }


        /////   CLIENTE
        // Decifra o U
        for(int i = 0; i < n_nodes; i++)
        {
            D( U_dec[i], U[i], lambda, micro, n );
            if((i != 0) && (i != (n_nodes-1)) )
                mpz_div_ui(U_dec[i],U_dec[i],mul); 
        }

        // print(U_dec,n_nodes);
        write_vector(U_dec, n_nodes, arq,0);

        // Cifra novamente e manda para o servidor
        for(int i = 0; i < n_nodes; i++)
        {
            E( C[i], U_dec[i], g, n );
        }
        j = j + dt;


    }
    time_2 = clock();

    mpz_clear(n);
    mpz_clear(n2);
    mpz_clear(g);
    mpz_clear(lambda);
    mpz_clear(micro);

    mpz_clear(tmp);
    clear_vector(U,n_nodes);
    clear_vector(C,n_nodes);
    clear_vector(U_dec,n_nodes);
    fclose(arq);

    return  ((time_2 - time_1)/ (double)CLOCKS_PER_SEC);
}

// MEF homomorfico no cenario 2
double homomorphic_sc2()
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

    // Exemplo 1
    int n_nodes;

    // Exemplo 2
    // int n_nodes = 7;

    int t_final;

    // numero de Fourier 1/mul
    int mul;

    // Lendo os parametros do problema

    arq = fopen("./in/param.txt","r");
    
    fscanf(arq,"%d\n",&n_nodes);
    fscanf(arq,"%d\n",&dx);
    fscanf(arq,"%d\n",&dt);
    fscanf(arq,"%d\n",&t_final);
    fscanf(arq,"%d\n",&mul);

    fclose(arq);


    mpz_t *U;
    init_vector(&U, n_nodes);

    mpz_t *C;
    init_vector(&C, n_nodes);

    mpz_t *U_dec;
    init_vector(&U_dec, n_nodes);

    arq = fopen("./in/IC.txt","r");

    // Cifrando a condição inicial
    for(int i = 0; i < n_nodes; i++)
    {
        gmp_fscanf(arq,"%ZX\n",C[i]);
    }
    fclose(arq);

    // for(int i = 0; i < n_nodes; i++)
    // {   
    //     D( U_dec[i], C[i], lambda, micro, n );

    // }
    // print(U_dec,n_nodes);
    // return 0;
    arq = fopen("./out/out_sc2.txt","w");

    // write_vector(C, n_nodes, arq,0);
    // clear_vector(Uant,n_nodes);

    // passando para o vetor U as condições de contorno criptografadas
    mpz_set(U[0],C[0]);
    mpz_set(U[n_nodes-1],C[n_nodes-1]);

    mpz_t tmp;
    mpz_init(tmp);
    mpz_t tmp1;
    mpz_init(tmp1);

    int j = dt; // Iteração começa do instante de tempo 0+dt
    
    int iter = 0;
    mpz_ui_pow_ui(tmp1,mul,iter);

    clock_t time_1, time_2;
	time_1 = clock();
    while(j <= t_final)
    {   
        mpz_powm( C[0]           , C[0],           tmp1, n2 );  
        mpz_powm( C[(n_nodes-1)] , C[(n_nodes-1)], tmp1, n2 );

        /////   SERVIDOR
        for(int i = 1; i < (n_nodes-1); i++)
        {
            mpz_powm_ui( tmp, C[i], 2, n2 );      // 
            mpz_invert(tmp, tmp, n2);            //  -2*Uant[i]

            mpz_powm_ui( U[i], C[i], mul, n2 );  // 5*Uant[i]

            mpz_mul( U[i], U[i], C[i+1] );       //  U[i] = U[i] + Uant[i+1]
            mpz_mod( U[i], U[i], n2 );

            mpz_mul( U[i], U[i], C[i-1] );       //  U[i] = U[i] + Uant[i-1]
            mpz_mod( U[i], U[i], n2 );

            mpz_mul( U[i], U[i], tmp );          //  U[i] = U[i] -  2*Uant[i]
            mpz_mod( U[i], U[i], n2 );
        }

        iter +=1;
        mpz_ui_pow_ui(tmp1,mul,iter);
        ///   CLIENTE
        // Decifra o U
        for(int i = 0; i < n_nodes; i++)
        {   
            D( U_dec[i], U[i], lambda, micro, n );
            if((i != 0) && (i != (n_nodes-1)) )
                    mpz_div(U_dec[i],U_dec[i],tmp1);
        }
        copy_vector(C,U,n_nodes);
        // write_vector(U, n_nodes, arq);
        // print(U_dec,n_nodes);
        write_vector(U_dec, n_nodes, arq,0);
        j = j + dt;
    }
        time_2 = clock();



    mpz_clear(n);
    mpz_clear(n2);
    mpz_clear(g);
    mpz_clear(lambda);
    mpz_clear(micro);

    mpz_clear(tmp);
    mpz_clear(tmp1);
    clear_vector(U,n_nodes);
    clear_vector(C,n_nodes);
    clear_vector(U_dec,n_nodes);
    fclose(arq);


    return  ((time_2 - time_1)/ (double)CLOCKS_PER_SEC);    
}


// MEF não homomorfico
double n_homomorphic() {
        

    FILE *arq;
    
    int dx = 5;
    int dt = 5;

    // Exemplo 1
    int n_nodes;

    // Exemplo 2
    // int n_nodes = 7;

    int t_final;

    // numero de Fourier 1/mul
    int mul;

    // Lendo os parametros do problema

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

    // Condição inicial
    for(int i = 1; i < (n_nodes-1); i++)
    {
        mpz_init_set_str(Uant[i],"20000000",10); // 20.000000
    }


    arq = fopen("./out/NH_out.txt","w");
    // print(Uant,n_nodes);
    write_vector(Uant, n_nodes, arq,0);

    mpz_t *U;
    U = malloc(sizeof(mpz_t) * n_nodes);

    // inicializando vetor
    for(int i = 0; i < n_nodes; i++)
    {
        mpz_init(U[i]); 
    }
    mpz_set(U[0],Uant[0]);
    mpz_set(U[n_nodes-1],Uant[n_nodes-1]);

    
    mpz_t tmp;
    mpz_init(tmp);

    int j = dt; // Iteração começa do instante de tempo 0+dt

    clock_t time_1, time_2;
	time_1 = clock();
    while(j <= t_final)
    {
        /////   SERVIDOR
        for(int i = 1; i < (n_nodes-1); i++)
        {
            
            mpz_mul_ui(tmp,Uant[i],2);          // tmp = 2*Uant[i]
            mpz_div_ui(tmp,tmp,mul);            // tmp = (0.2)*2*Uant[i]

            mpz_add(U[i],Uant[i+1],Uant[i-1]);  // U[i] = Uant[i+1]+Uant[i-1]
            mpz_div_ui(U[i],U[i],mul);          // U[i] = (0.2)*U[i]
            mpz_add(U[i],U[i],Uant[i]);         // U[i] = (0.2)*(Uant[i+1]+Uant[i-1]) + Uant[i]

            mpz_sub(U[i], U[i], tmp);           // U[i] = (0.2)*(Uant[i+1]+Uant[i-1]) + Uant[i] - (0.2)*2*Uant[i]
        }
        // print(U,n_nodes);
        copy_vector(Uant,U,n_nodes);
        
        write_vector(U, n_nodes, arq,0);

        j = j + dt;
    }

    time_2 = clock();

    clear_vector(U,n_nodes);
    fclose(arq);

    return  ((time_2 - time_1)/ (double)CLOCKS_PER_SEC);
}
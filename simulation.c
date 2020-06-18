#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "trabalho.h"
#include "Pailler.h"
#include "utils.h"



int main() {
    double total_time = 0;
    printf("Nao Homomorfico\n");
    for(int i=0; i < 1; i++){
        total_time = total_time + n_homomorphic();
    }
    printf("Total Time %f seconds\n",total_time/10.0);

    total_time = 0;
    printf("Homomorfico Cenario 1\n");
    for(int i=0; i < 1; i++){
        total_time = total_time + homomorphic_sc1();;
    }
    printf("Total Time %f seconds\n",total_time/10.0);


    total_time = 0;
    printf("Homomorfico  Cenario 2\n");
    for(int i=0; i < 1; i++){
        total_time = total_time + homomorphic_sc2();;
    }
    printf("Total Time %f seconds\n",total_time/10.0);

    return 0;
}

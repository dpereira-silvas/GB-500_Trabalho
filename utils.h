#ifndef UTILS_H
#define UTILS_H

void init_vector(mpz_t **v, int n);
void write_vector(mpz_t *v, int n, FILE *f,int op);
void clear_vector(mpz_t *v, int n);
void copy_vector(mpz_t *u, mpz_t *v, int n);

#endif /* UTILS_H */

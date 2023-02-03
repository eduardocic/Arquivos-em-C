#ifndef _MAT_H_
#define _MAT_H_


/* cria números randômicos */
void randu(int size, double minValue, double maxValue, int seed, matrix *M);
void randn(int size, int seed, matrix *M);

/* Demais funções */
matrix* linspace(double x1, double x2, int n);
int fatorial(int x); 
double min(double a, double b);
double max(double a, double b);



#endif

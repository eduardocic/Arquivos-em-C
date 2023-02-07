#ifndef _CALC_NUM_H_
#define _CALC_NUM_H_

/* Estruturas de dados para as funções */
typedef struct {
    matrix *temp;
    matrix *k1;
    matrix *k2;
    matrix *hh;
    matrix *meio;
} RK2_s;


typedef struct {
    matrix *temp;
    matrix *k1;
    matrix *k2;
    matrix *k3;
    matrix *k4;
    matrix *hh;
    matrix *meio;
    matrix *dois;
    matrix *umsexto;
} RK4_s;


/* Definindo ponteiro para função */
typedef void   (*ptrFuncao_MMM)(double t, matrix *x, matrix *u, matrix *xnew);
typedef double (*ptrFuncao_D  )(double x);

/* Métodos de cálculo de derivada */
void RK2(ptrFuncao_MMM function, double t, double h, matrix *x, matrix *u, matrix *xnew);
void RK4(ptrFuncao_MMM function, double t, double h, matrix *x, matrix *u, matrix *xnew);
void RK2_init(int n, RK2_s *ptr);
void RK4_init(int n, RK4_s *ptr);

/* Métodos de cálculo de Integral */
double integral(ptrFuncao_D function, double t0, double tf);

//void funcao(double t, matrix *x, matrix *u, matrix *y);
//void VanderPol(double t, matrix *x, matrix *u, matrix *y);

#endif

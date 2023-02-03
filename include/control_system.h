#ifndef _CONTROL_SYSTEM_H_
#define _CONTROL_SYSTEM_H_

/* ================================================================
 *
 *  					Estruturas diversas
 *
 * ============================================================== */ 
/* Estrutura para um sistema linear */
typedef struct{
	matrix *A;
	matrix *B;
	matrix *C;
	matrix *D;
} lin_sys;

/* Estrutura de dados com as variáveis que rodam na função principal de Kalman */
typedef struct {
    matrix *PHIk_t;             /* Transposta de PHIk, ou seja, (PHIk)'                */
    matrix *Mk;                 /* Variável derivada na EQ. 1                          */
    matrix *H_t;                /* Transposta de H, ou seja, (H)'                      */
    matrix *temp_nxp;           /* Variável de auxílio aos cálculos                    */
    matrix *temp_nxn;           /* Variável de auxílio aos cálculos                    */
    matrix *temp_pxn;           /* Variável de auxílio aos cálculos                    */
    matrix *soma;               /* Matriz soma na EQ. 2 que fará a inversão            */
    matrix *inv_soma;           /* Matriz inversa da soma                              */
    matrix *Inxn;               /* Ponteiro para uma matriz identidade de dimensão nxn */
} SKalman;




/* ==================================================================
 *
 *						PROTÓTIPOS DE FUNÇÕES
 *
 * ================================================================ */
 
/* ==================== */
/*   FUNÇÕES LINEARES   */
/* ==================== */
lin_sys ss(matrix *a, matrix *b, matrix *c, matrix *d);	        /* Função 'ss': State Space 			   		   */
void fun_lsim(double t, matrix *x, matrix *u, matrix *xnew);    /* Função 'fun_lsim': utilizada pela 'lsim' apenas */
matrix* lsim(lin_sys sys, matrix *u, matrix *t, matrix *x0);	/* Função 'lsim': linear simulation		           */

/* ==================== */
/*   FILTRO DE KALMAN   */
/* ==================== */
void kalman_init(int n, int p, SKalman *ptr);					/* Função 'kalman_init': chamada por 'kalman       */
void kalman(matrix *P_k_1, matrix *PHI_k, matrix *Q_k, 
            matrix *R_k, matrix *H, matrix *K_k);				/* Função 'kalman': função de kalman               */




#endif

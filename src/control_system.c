#include <stdio.h>
#include <stdlib.h>


#include "matrix.h"
#include "control_system.h"
#include "calc_num.h"
#include "matematica.h"


/* ======================================================== */
/* 				Função para soma de variáveis 				*/
/* ======================================================== */
void fun_lsim(double t, matrix *x, matrix *u, matrix *xnew){    

	double x1, x2, uk, x1d, x2d;
	x1 = get(x, 0, 0);
	x2 = get(x, 1, 0);
	uk = get(u, 0, 0);
	
	x1d = x2;
	x2d = -3.0 * x1 - 2.0 * x2 + 3.0 * uk;
	
	set( x1d, xnew, 0, 0);
	set( x2d, xnew, 1, 0);
}




/* ======================================================== */
/* 						Função 'lsim'	    				*/
/* ======================================================== */
matrix* lsim(lin_sys sys, matrix *u, matrix *t, matrix *x0){
	
	/* Variáveis para auxílio */
	int num_states  = sys.A->lines;		/* Número de estados 			   */
	int num_control = sys.B->columns;   /* Número de variáveis de controle */
	int num_output  = sys.C->lines; 	/* Número de variáveis de saída    */
	
	/* Criando variável auxiliar para cálculo */
	matrix *yout  = new_matrix(t->lines, num_output);      /* vetor de saída (completo)        */
	matrix *ytemp = new_matrix( num_output, 1);            /* vetor de saída intermediário     */
	matrix *x     = new_matrix( num_states, 1);		       /* vetor de estados em 'k'		   */
	matrix *xnew  = new_matrix( num_states, 1);            /* vetor de estados em 'k+1' 	   */
	matrix *uu    = new_matrix(num_control, 1);			   /* vetor de entradas em 'k'		   */
		
	/* Variáveis de inicialização de tempo e de intervalo de integração do sistema */
	double h     = (get(t, (t->lines)-1, 0) - 0.0)/t->lines;
	double time  = get(t, 0, 0);
	
	/* Inicializando o sistema */
	isequal(x, x0);
	for (int cont = 0 ; cont < (t->lines) ; cont++){
		
		/* Pegar o extrato de 'u' e salvar em 'uu' */
		for (int i = 0 ; i < num_control ; i++){
			set(get(u, cont, i), uu, i, 0);
		}		
		
		/* Executa o loop */				
		RK4(fun_lsim, time, h, x, uu, xnew);			
		
		/* Pega a saída */
		prod(sys.C, x, ytemp);
		for (int i = 0 ; i < num_output ; i++){
			set(get(ytemp, i, 0), yout, cont, i);
		}
		
		/* Atualiza o estado e o contador */
		isequal(x, xnew);
		time += h;
	}
	
	/* Libera o heap */
	free(ytemp);
	free(x);
	free(xnew);
	free(uu);
		
	return yout;
}




lin_sys ss(matrix *a, matrix *b, matrix *c, matrix *d){
/* ========================================================
 *  Esta função salva na estrutura do tipo 'lin_sys' os pa-
 *  râmetros de um sistema linear
 * ======================================================== */
	/* Cria uma estrutura */
	lin_sys e;
	
	/* Variáveis para auxílio */
	int num_states  = a->lines;		/* Número de estados 			   */
	int num_control = b->columns;   /* Número de variáveis de controle */
	int num_output  = c->lines; 	/* Número de variáveis de saída    */

	/* É necessário que se faça a inserção de um código aqui que possibi-
	   lite tratar excessões de dimensões de matrizes estranhas 	   */

	/* Alimenta ela com os parâmetros do sistema */
	e.A   = a;
	e.B   = b;
	e.C   = c;
	e.D   = d;
	
	/* Retorna a estrutura criada */
	return e;
}





void kalman(matrix *P_k_1, matrix *PHI_k, matrix *Q_k, matrix *R_k, matrix *H, matrix *K_k){
/* =========================================================================
 *
 *  REFERẼNCIA: Fundamentals of Kalman Filtering: A Practical Approach, 
 * ============ 3th Edition
 *              Paul Zarchan & Howard Musoff
 *              Página 131
 *
 *
 *  EQUAÇÕES A SE RESOLVER:
 * ========================
 *
 *  ---> Mk = PHI_k * P_k_1 * (PHI_k)' + Q_k             (EQ. 1)
 *  ---> K_k = Mk * (H)' * inv{(H * Mk * (H)' + R_k)}    (EQ. 2)
 *  ---> P_k = (I - K_k * H) * Mk                        (EQ. 3)
 *
 *
 *  ENTRADAS:
 * ==========
 *
 *    1. P_k_1 ..........: tem dimensão 'n x n';
 *    2. PHI_k ..........: tem dimensão 'n x n';
 *    3. Q_k ............: tem dimensão 'n x n';
 *    4. R_k ............: tem dimensão 'p x p';
 *    5. H ..............: tem dimensão 'p x n';
 *    6. K_k ............: tem dimensão 'n x p'; e
 *    7. ptr ............: ponteiro para estrutura do tipo 'SKalman'.
 *
 * Eduardo H. Santos
 * 25/01/2023
 * ======================================================================= */

    /* Necessária definição da estrutura de Kalman */
    /* =========================================== */
    static bool flag = false;
    static SKalman ptr = { .PHIk_t   = NULL, 
                           .Mk       = NULL,
                           .H_t      = NULL, 
                           .temp_nxp = NULL,
                           .temp_nxn = NULL, 
                           .temp_pxn = NULL, 
                           .soma     = NULL, 
                           .inv_soma = NULL, 
                           .Inxn     = NULL};

    /* Inicializa as variáveis static a serem utilizadas apenas dentro aqui */
    /* ==================================================================== */
    if (!flag){
        /* Variáveis auxiliares */    
        int n = P_k_1->lines;
        int p = H->lines;

        /* Inicialização dos parâmetros de Kalman */
        kalman_init(n, p, &ptr);

        /* Faz a operação de transposição de matrizes pertinentes */
        transpose(PHI_k, ptr.PHIk_t);
        transpose(    H, ptr.H_t   );
        flag = true;
    }

    /* (EQ. 1)    Mk = PHI_k * P_k_1 * (PHI_k)' + Q_k    */
    /* ================================================= */
    prod(       PHI_k,      P_k_1, ptr.temp_nxn);
    prod(ptr.temp_nxn, ptr.PHIk_t,       ptr.Mk);
    add(       ptr.Mk,        Q_k,       ptr.Mk);

    /* (EQ. 2)  Mk * (H)' * inv{(H * Mk * (H)' + R_k)}   */
    /* ================================================= */
    prod(            H,         ptr.Mk,  ptr.temp_pxn);
    prod( ptr.temp_pxn,        ptr.H_t,      ptr.soma);
    add(      ptr.soma,            R_k,      ptr.soma);
    inv(      ptr.soma,  ptr.inv_soma);
    prod(       ptr.Mk,        ptr.H_t,  ptr.temp_nxp);
    prod( ptr.temp_nxp,   ptr.inv_soma,           K_k);

    /* (EQ. 3)     P_k = (I - K_k * H) * Mk              */
    /* ================================================= */
    prod(          K_k,        H,         P_k_1);
    sub(      ptr.Inxn,    P_k_1,  ptr.temp_nxn); 
    prod( ptr.temp_nxn,   ptr.Mk,         P_k_1);
}





void kalman_init(int n, int p, SKalman *ptr){
/* =========================================================================
 * Inicialização de estrutura do tipo 'SKalman' com parâmetros estáticos
 * os quais são utilizados recursivamente pela função de 'kalman'
 * ======================================================================= */
 
	/* Aloca espaço para matrizes utilizadas no programa. */
	matrix *PHIk_t    = new_matrix(n, n);
	matrix *Mk        = new_matrix(n, n);
	matrix *H_t       = new_matrix(n, p);
	matrix *temp_nxp  = new_matrix(n, p);
	matrix *temp_nxn  = new_matrix(n, n);
	matrix *temp_pxn  = new_matrix(p, n);
	matrix *soma      = new_matrix(p, p);  
	matrix *inv_soma  = new_matrix(p, p); 
	matrix *Inxn      = eye(n);

	/* Faz o ponteiro apontar para as variáveis */
	ptr->PHIk_t    = PHIk_t;
	ptr->Mk        = Mk;
	ptr->H_t       = H_t;
	ptr->temp_nxp  = temp_nxp;
	ptr->temp_nxn  = temp_nxn;
	ptr->temp_pxn  = temp_pxn;
	ptr->soma      = soma;
	ptr->inv_soma  = inv_soma;
	ptr->Inxn      = Inxn;
}





//matrix* phi_k(matrix* PHI_t, double ts, int ordem){
//
//    /* Dimensão do sistema */
//    int n = PHI_t->lines;
//
//    /* Recussão */
//    if(ordem == 0){
//        return eye(n);    
//    } else {
//
//        /* Matrizes auxiliares e de saída*/
//        matrix *PHI_k_x_Ts  = new_matrix(n, n);
//        matrix *soma        = new_matrix(n, n);
//        matrix *temp1       = eye(n);
//        matrix *temp2       = eye(n);
//        matrix *PHI_k       = eye(n);
//        matrix *Ts          = new_constant(ts);
//        matrix *fat         = new_constant(1.0);
//
//        /* A ser utilizada recorrentemente (PHI * Ts)*/
//        prod(Ts, PHI_t, PHI_k_x_Ts);   
//        
//        /* Faz a recursão */
//        for (int i = 0 ; i < ordem ; i++){   
//            /* Calculo o fatorial */ 
//            set(1.0/fatorial(i+1), fat, 0, 0);     
//                
//            /* calculo o produto da matriz por ela mesma */
//            prod(temp1, PHI_k_x_Ts, temp2);
//            isequal(temp1, temp2);
//            prod(fat, temp2, temp2);
//            add(temp2, PHI_k, PHI_k);
//        }
//
//        /* Libera a memória do heap */
//        free(PHI_k_x_Ts);
//        free(soma);
//        free(temp1);
//        free(temp2);
//        free(Ts);
//        free(fat);
//
//        return PHI_k;
//    }
//}


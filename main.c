#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#include "matrix.h"
#include "calc_num.h"
#include "matematica.h"
#include "control_system.h"



int main(void){

	clock_t start_t, end_t;
	double total_t;

	start_t = clock();

	/* ================================== */	
	/*    Criando a estrutura de dados    */
	/* ================================== */
	matrix *A = new_matrix(2, 2);
	matrix *B = new_matrix(2, 1);
	matrix *C = new_matrix(1, 2);
	matrix *D = new_matrix(1, 1);


	/* ================================== */	
	/* Ajustando os parâmetros do sistema */
	/* ================================== */
	set( 0.0, A, 0, 0);
	set( 1.0, A, 0, 1);
	set(-3.0, A, 1, 0);
	set(-2.0, A, 1, 1);
	set( 0.0, B, 0, 0);
	set( 3.0, B, 1, 0);	
	set( 1.0, C, 0, 0);
	set( 0.0, C, 0, 1);
	set( 0.0, D, 0, 0);

	/* Inicialização do sistema linear */
	LINEAR_SYSTEM_s sys = ss(A, B, C, D);

	/* Tempo total de simulação 	   */	
	matrix *t = linspace(0.0, 8.0, 201);
	
	/* Vetor de entradas 			   */
	matrix *u = linspace(0.0, 0.0, 201);
	for(int i = 0 ; i < 201 ; i++){
		double MAX, MIN;
		MIN = min( (get(t, i, 0) - 1.0), 1.0);
		MAX = max(	    			0.0, MIN);
		set(MAX, u, i, 0);
	}
	
	/* Condição inicial 		*/
	matrix *x0 = new_matrix(2, 1);
	set(0.0, x0, 0, 0);
	set(0.0, x0, 1, 0);
	
	/* Saída 					*/
	matrix *y = lsim(sys, u, t, x0);
	
	
	
	FILE *fp;
    fp = fopen("lsim.txt", "w");
    fprintf(fp, "t    y    u\n");
	for (int i = 0 ; i < t->lines; i++){
		fprintf(fp, "%f  %f  %f\n", 
				get(t, i, 0), 
                get(y, i, 0),
                get(u, i, 0));
	}
	fclose(fp);
	
	end_t = clock();
	total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	printf("Tempo total tomado pela CPU: %f\n", total_t);
	printf("Saindo do programa ...\n");
	return 0;
}













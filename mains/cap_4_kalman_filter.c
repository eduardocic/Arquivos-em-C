#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "matrix.h"
#include "matematica.h"
#include "calc_num.h"
#include "control_system.h"


void funKalman(double t, matrix *x, matrix *u, matrix *y){
    /* ============================
       | y1 |   |  x2 |
       |    |   |     |
       | y2 | = |  x3 |
       |    |   |     |
       | y3 |   | 0.0 | 
       ============================ */
    double x2 = get(x, 1, 0);
    double x3 = get(x, 2, 0);
     
    // y1 = x2 
    set(x2, y, 0, 0);

    // y2 = x3
    set(x3, y, 1, 0);     

    // y3 = 0
    set(0.0, y, 2, 0);     
}



int main(){

    // =================================================================
    //
    //       Variáveis constantes e utilizadas para a simulação
    //
    // =================================================================
    // Gravidade
    double g = 32.2;

    // Tempo de simulação e parâmetros derivados
    double Ttotal = 30.0;                      // tempo total de amostragem
    double fs = 10.0;                          // frequencia de amostragem
    double h  = 1.0/fs;                        // passo de amostragem 
    double Ts = h;                             // tempo de amostragem
    int nT = (int)(fs*Ttotal + 1.0);           // intervalos total de simulação


    // =================================================================
    //
    //           Variáveis de estado, de controle e de saída
    //
    // =================================================================
    // x(k)
    matrix *x_k = new_matrix(3, 1);
    set(400000.0, x_k, 0, 0);
    set( -6000.0, x_k, 1, 0);
    set(      -g, x_k, 2, 0);      /* nova */
    
    // x(k+1)
    matrix *x_k1 = new_matrix(3, 1);

    // u
    matrix *u = new_constant(g);

    // y
    matrix *y_k = new_constant(0.0);


    // =================================================================
    //  
    //                       Simulação do ruído
    //
    // =================================================================
    // No caso, o radar mede a posição do objeto, com ruído. 
    //
    // --> O primeiro ponto é modelar o ruído:
    matrix *desvio_padrao = new_constant(1000.0);
    matrix *noise = new_matrix(nT, 1);
    randn(nT, 0, noise);
    prod(desvio_padrao, noise, noise);

    // =================================================================
    //  
    //         Aloca espaço para as matrizes do filtro de kalman
    //
    // ================================================================= 
    // PHI_k
    matrix *PHI_k = new_matrix(3, 3);
    set(    1.0, PHI_k, 0, 0);
    set(     Ts, PHI_k, 0, 1);
    set(Ts*Ts/2, PHI_k, 0, 2);
    set(    0.0, PHI_k, 1, 0);
    set(    1.0, PHI_k, 1, 1);
    set(     Ts, PHI_k, 1, 2);
    set(    0.0, PHI_k, 2, 0);
    set(    0.0, PHI_k, 2, 1);
    set(    1.0, PHI_k, 2, 2);

    // Q_k
    matrix *Q_k = new_matrix(3, 3);
    double phiS = 0.0;
    set(  pow(Ts, 5.0) * (phiS/20), Q_k, 0, 0);
    set(  pow(Ts, 4.0) *  (phiS/8), Q_k, 0, 1);
    set(  pow(Ts, 3.0) *  (phiS/6), Q_k, 0, 2);
    set(  pow(Ts, 4.0) *  (phiS/8), Q_k, 1, 0);
    set(  pow(Ts, 3.0) *  (phiS/3), Q_k, 1, 1);
    set(  pow(Ts, 2.0) *  (phiS/2), Q_k, 1, 2);
    set(  pow(Ts, 3.0) *  (phiS/6), Q_k, 2, 0);
    set(  pow(Ts, 2.0) *  (phiS/2), Q_k, 2, 1);
    set(          Ts *        phiS, Q_k, 2, 2);

    // H
    matrix *H = new_matrix(1, 3);
    set(1.0, H, 0, 0);

    // R_k
    matrix *R_k = new_matrix(1, 1);
    prod(desvio_padrao, desvio_padrao, R_k);

    // P_(k-1)
    matrix *P_k_1 = eye(3);
    matrix *c = new_constant(999999999999.0);
    prod(c, P_k_1, P_k_1);

    // K_k
    matrix *K_k = new_matrix(3, 1);



    // =================================================================
    //  
    //          Variáveis temporárias para auxílio no cálculo
    //
    // =================================================================
    int n = P_k_1->lines;
    int p = H->lines;
    matrix *temp_pxp  = new_matrix (p, p);
    matrix *temp_pxn  = new_matrix (p, n);
    matrix *temp_nxp  = new_matrix (n, p);
    matrix *res       = new_matrix (p, p);


    // Estimativa inicial
    // xhat(k) e xhat(k+1)
    matrix *xhat_k  = new_matrix(3, 1);
    matrix *xhat_k1 = new_matrix(3, 1);


    // =================================================================
    //  
    //                       Executa a simulação
    //
    // =================================================================
    int contNoise = 0;
    FILE *fp;
    fp = fopen("test2.txt", "w");
    fprintf(fp, "t    x_k1    x_k2    xh1    xh2    xh3    p11    p22    p33\n");
    for (double t = 0 ; t <= Ttotal ; t = t+h){

        // Executa o algoritmo de Kalman
        kalman(P_k_1, PHI_k, Q_k, R_k, H, K_k);

        // x(k+1) e xhat(k+1)
        RK4(funKalman, t, h,    x_k, u,    x_k1);
        RK4(funKalman, t, h, xhat_k, u, xhat_k1);

        // Adiciona o ruído à medida real
        contNoise++;
        set(get(x_k, 0, 0) + get(noise, contNoise, 0), y_k, 0, 0); 

        // Determinação do resíduo e produto com a matriz K_k
        prod(       H,    PHI_k, temp_pxn);
        prod(temp_pxn,   xhat_k, temp_pxp);
        sub(      y_k, temp_pxp,      res);        
        prod(     K_k,      res, temp_nxp);  

        // Corrige o valor de xhat_k(k)
        add(xhat_k1, temp_nxp, xhat_k1);

        // Salva o resultado
        fprintf(fp, "%f  %f  %f  %f  %f  %f  %f  %f  %f \n", 
                t, 
                get(x_k, 0, 0),
                get(x_k, 1, 0),
                get(xhat_k, 0, 0),
                get(xhat_k, 1, 0),
                get(xhat_k, 2, 0), 
                pow(get(P_k_1, 0, 0), 0.5),
                pow(get(P_k_1, 1, 1), 0.5),
                pow(get(P_k_1, 2, 2), 0.5));

        // Atualizo as variáveis para o próximo passo
        isequal(   x_k,    x_k1);    /* x(k)    = x(k+1)    */
        isequal(xhat_k, xhat_k1);    /* xhat(k) = xhat(k+1) */
    }  
    fclose(fp);


    /* Libera a memória de todo mundo */
    /* ============================== */
    free(x_k);
    free(x_k1);
    free(u);
    free(y_k);

    free(desvio_padrao);
    free(noise);

    free(PHI_k);
    free(Q_k);
    free(H);
    free(R_k);
    free(P_k_1);
    free(c);
    free(K_k);

    return 0;
}

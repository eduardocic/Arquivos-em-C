#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "calc_num.h"


void RK2_init(int n, RK2_struct *ptr){

    /* Aloca espaço para matrizes utilizadas no programa. */
    matrix *temp = new_matrix(n, 1);
    matrix *k1   = new_matrix(n, 1);
    matrix *k2   = new_matrix(n, 1);
    matrix *hh   = new_matrix(1, 1);
    matrix *meio = new_matrix(1, 1);

    /* Faz o ponteiro apontar para as variáveis */
    ptr->temp    = temp;
    ptr->k1      = k1;
    ptr->k2      = k2;
    ptr->hh      = hh;
    ptr->meio    = meio;
}



void RK2(ptrFuncao_MMM function, double t, double h, matrix *x, matrix *u, matrix *xnew){
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Procedimento que resolve o seguinte problema :
    _
   |  
   |    k1 = h * f[          x_n ,        t_n) ]        (EQ. 1) 
 --|    k2 = h * f[ (x_n + k1/2) , (t_n + h/2) ]        (EQ. 2) 
   | x_n1  = x_n + k2                                   (EQ. 3) 
   |_

em que:

   i) x_n1 .........: simbologia para 'x_{n+1}';
  ii) x_n ..........: simbologia para 'x_{n}'; e
 iii) h ............: simbologia para o passo de simulação.  


 PARAMETROS DA FUNÇÃO
======================

    1. function .........: é um ponteiro para função;
    2. t ................: variável de tempo corrente;
    3. h ................: passo de interação;
    4. x ................: vetor de estados de entrada;
    5. u ................: vetor de controle; e
    6. xnew .............: vetor de estados de saída.

Eduardo H. Santos
18/01/2023
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    /*  Necessária definição da estrutura de RK2   */
    /* =========================================== */
    static bool flag = false;
    static RK2_struct ptr = { .temp  = NULL, 
                              .k1    = NULL,
                              .k2    = NULL, 
                              .hh    = NULL,
                              .meio  = NULL};

    /* Inicializa as variáveis static a serem utilizadas apenas dentro aqui */
    /* ==================================================================== */
    if (!flag){

        // Inicialização dos parâmetros de RK2
        RK2_init(x->lines, &ptr);

        // Faz a inicialização das matrizes pertinentes.
        ptr.hh->data[0]   = h;
        ptr.meio->data[0] = 0.5;
        flag = true;
    }

    /*  (EQ. 1)          k1 = h * f[ x_n , t_n) ]            */
    /* ===================================================== */
    (*function)(t, x, u, xnew);
    prod(ptr.hh, xnew, ptr.k1);  

    /*  (EQ. 2)   k2 = h * f[ (x_n + k1/2) , (t_n + h/2) ]   */ 
    /* ===================================================== */        
    prod(ptr.meio, ptr.k1, ptr.temp);               
    add(x, ptr.temp, xnew);
    (*function)((t + 0.5*h), x, u, ptr.temp);
    prod(ptr.hh, ptr.temp, ptr.k2);
    
    /*  (EQ. 3)            x_n1  = x_n + k2                  */
    /* ===================================================== */ 
    add(x, ptr.k2, xnew);
}



void RK4_init(int n, RK4_struct *ptr){

    /* Aloca espaço para matrizes utilizadas no programa. */
    matrix *temp     = new_matrix(n, 1);
    matrix *k1       = new_matrix(n, 1);
    matrix *k2       = new_matrix(n, 1);
    matrix *k3       = new_matrix(n, 1);
    matrix *k4       = new_matrix(n, 1);
    matrix *hh       = new_matrix(1, 1);
    matrix *meio     = new_matrix(1, 1);
    matrix *dois     = new_matrix(1, 1);
    matrix *umsexto  = new_matrix(1, 1);

    /* Faz o ponteiro apontar para as variáveis */
    ptr->temp     = temp;
    ptr->k1       = k1;
    ptr->k2       = k2;
    ptr->k3       = k3;
    ptr->k4       = k4;
    ptr->hh       = hh;
    ptr->meio     = meio;
    ptr->dois     = dois;
    ptr->umsexto  = umsexto;
}



void RK4(ptrFuncao_MMM function, double t, double h, matrix *x, matrix *u, matrix *xnew){
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Procedimento que resolve o seguinte problema:
    _
   |  
   |    k1 = h * f[          x_n ,        t_n  ]        (EQ. 1) 
   |    k2 = h * f[ (x_n + k1/2) , (t_n + h/2) ]        (EQ. 2)
 --|    k3 = h * f[ (x_n + k2/2) , (t_n + h/2) ]        (EQ. 3) 
   |    k4 = h * f[   (x_n + k3) ,   (t_n + h) ]        (EQ. 4)  
   | x_n1  = x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4)        (EQ. 5)
   |_

Eduardo H. Santos
18/01/2023
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    /* Criando uma matriz temporária a ser utilizada recorrentemente */
    /* ============================================================= */
    static bool flag = false;
    static RK4_struct ptr = { .temp     = NULL, 
                              .k1       = NULL,
                              .k2       = NULL,
                              .k3       = NULL, 
                              .k4       = NULL,
                              .hh       = NULL,
                              .meio     = NULL,
                              .dois     = NULL,
                              .umsexto  = NULL};

    /* Inicializa as variáveis static a serem utilizadas apenas dentro aqui */
    /* ==================================================================== */
    if (!flag){
        // Inicialização dos parâmetros de RK2
        RK4_init(x->lines, &ptr);

        // Faz a inicialização das matrizes pertinentes.
        ptr.hh->data[0]      = h;
        ptr.meio->data[0]    = 0.5;
        ptr.dois->data[0]    = 2.0;
        ptr.umsexto->data[0] = 1.0/6;
        flag = true;
    }

    /*  (EQ. 1)           k1 = h * f[ x_n, t_n]              */
    /* ===================================================== */  
    (*function)(t, x, u, ptr.temp);
    prod(ptr.hh, ptr.temp, ptr.k1);

    /*  (EQ. 2)  k2 = h * f[ (x_n + k1/2) , (t_n + h/2) ]    */
    /* ===================================================== */  
    prod(ptr.meio, ptr.k1, ptr.temp);               
    add(x, ptr.temp, ptr.temp);
    (*function)((t + 0.5*h), ptr.temp, u, ptr.temp);
    prod(ptr.hh, ptr.temp, ptr.k2);

    /*  (EQ. 3)  k3 = h * f[ (x_n + k2/2) , (t_n + h/2) ]    */
    /* ===================================================== */         
    prod(ptr.meio, ptr.k2, ptr.temp);
    add(x, ptr.temp, ptr.temp);         
    (*function)((t + 0.5*h), ptr.temp, u, ptr.temp);
    prod(ptr.hh, ptr.temp, ptr.k3);

    /*  (EQ. 4)    k4 = h * f[ x_n + k3,   (t_n + h) ]       */
    /* ===================================================== */       
    add(x, ptr.k3, ptr.temp);
    (*function)((t + h), ptr.temp, u, ptr.temp);
    prod(ptr.hh, ptr.temp, ptr.k4);

    /*  (EQ. 5)  x_n1  = x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4) */
    /* ===================================================== */ 
    prod(ptr.dois, ptr.k2, ptr.k2);
    prod(ptr.dois, ptr.k3, ptr.k3);
    add(ptr.k1, ptr.k2, ptr.temp);
    add(ptr.temp, ptr.k3, ptr.temp);
    add(ptr.temp, ptr.k4, ptr.temp);
    prod(ptr.umsexto, ptr.temp, xnew);
    add(x, xnew, xnew); 
}




/* Integral genérica */
/* ================= */
double integral(ptrFuncao_D function, double t0, double tf){
/* ================================================================
 *
 *  --> Todas as informações contidas nesta função eu peguei do se-
 *      guinte link: 
 * 
 * https://www.ufrgs.br/reamat/CalculoNumerico/livro-sci/in.html
 *
 *  --> A partir desse site, com direcionamento para outros links,
 *      que eu li e programei os seguintes métodos; e
 *  --> A programação apresentada segui a apresentação dos tópicos
 *      apresentados no site supracitado.
 *
 * Eduardo H. Santos
 * 26/01/2023
 * ===============================================================*/    

    /* Total de separação entre 't0' e 'tf' */
    int N    = 10001;   
    double h = (tf - t0)/N;
    double soma = 0.0;

    /* ===============================================
     *   Soma de Riemann
     * 
     * -- A regra implementada é a do 'ponto médio'; 
     * -- Divisão do intervalo (tf - t0) em N partes;
     * -- No ponto médio das partes que calculamos o 
     *    valor da função.
     * -- A função implementada não passou por simpli-
     *    ficação.
     * -- Para maiores informações, consultar o link
     *    supramencionado na aba 'Soma Riemann'.
     * ============================================== */
    double xi;       /* x_i        */
    double xi1;      /* x_{i+1}    */
    double xm;       /* x médio    */
    double fx;       /* f(xm)      */
    for (int i = 1 ; i <= N ; i++){
        xi  = t0 + (i-1) * h;    
        xi1 = xi + h;
        xm  = (xi + xi1) / 2;
        fx = (*function)(xm);
        soma += (fx * h);
    }
    printf("Soma pelo método de Riemann: %f\n", soma);


    /* ===============================================
     *   Newton-Cotes: Regra do ponto médio simples
     * ============================================= */
    soma = 0.0;
    double x = (tf + t0)/2;
    soma = (tf - t0) * ((*function)(x));
    printf("Soma pelo método de Newton-Cotes (Ponto Médio simples): %f\n", soma);


    /* ===============================================
     *   Newton-Cotes: Regra do ponto médio composto
     * ============================================= */
    soma = 0.0;
    double pm;
    for(double x = t0 ; x < tf ; x += h){
       pm = (2.0*x + h)/2;
       soma += (h) * ((*function)(pm));
    }
    printf("Soma pelo método de Newton-Cotes (Ponto Composto): %f\n", soma);


    /* ===============================================
     *         Newton-Cotes: Regra do trapézio
     * ============================================= */
    soma = 0.0;
    double fa, fb;    
    for (double x = t0 ; x < tf ; x += h){
        fa = (*function)(x);
        fb = (*function)(x + h);
        soma += (fb + fa)*h/2;
    }
    printf("Soma pelo método de Newton-Cotes (Trapézio): %f\n", soma);


    /* ===============================================
     *     Newton-Cotes: Regra de Simpson Simples
     * ============================================= */
    soma = 0.0;
    xm   = (t0 + tf) / 2;    
    soma += (*function)(t0) * (1.0/3);
    soma += (*function)(xm) * (4.0/3);
    soma += (*function)(tf) * (1.0/3);
    soma *= (tf - t0)/2;
    printf("Soma pelo método de Newton-Cotes (Regra de Simpson Simples): %f\n", soma);


    /* ===============================================
     *     Newton-Cotes: Regra de Simpson Composto
     * ============================================= */
    soma = 0.0;
    double soma_temp;
    double a, b;
    
    for (double x = t0 ; x < tf ; x += h){
        a  = x;
        b  = x + h;
        xm = (a + b) / 2;    
        soma_temp = 0.0;
        soma_temp += (*function)(a)  * (1.0/3);
        soma_temp += (*function)(xm) * (4.0/3);
        soma_temp += (*function)(b)  * (1.0/3);
        soma += soma_temp *(b - a)/2;
    }
    printf("Soma pelo método de Newton-Cotes (Regra de Simpson Composta): %f\n", soma);

    return soma;
}






















//void funcao(double t, matrix *x, matrix *u, matrix *xnew){
//    /* Variável auxiliar temporária */
//    double temp;
//    double calc;
//
//    /* ==========================
//                 2000 - 2v(t)
//        v'(t) = -------------
//                  200 - t 
//       ========================== */
//    temp = (2000 - 2*get_num(x_in, 0, 0))/(200 - t);
//    set_num(temp, x_out, 0, 0);  
//
//    /* ==========================
//        x'(t) = x^2 + t^2       
//       ========================== */
//    calc = get(x, 0, 0);    
//    temp = calc*calc + t*t;
//    set(temp, x, 0, 0);     
//}
//






//void VanderPol(double t, matrix *x, matrix *u, matrix *xnew){
//    /* ===========================================
//        y[0] = x[1] 
//        y[1] = -u[0]*(x[0]*x[0] - 1)*x[1] - x[0];
//       =========================================== */
//    double x0 = get(x, 0, 0);
//    double x1 = get(x, 1, 0);
//    double u0 = get(u, 0, 0);
//     
//    // y[0] = x[1] 
//    set(x1, y, 0, 0);
//
//    // y[1] = -u[0]*(x[0]*x[0] - 1)*x[1] - x[0];
//    double temp;
//    temp = -u0*(x0*x0 - 1)*x1 - x0;
//    set(temp, y, 1, 0);     
//}




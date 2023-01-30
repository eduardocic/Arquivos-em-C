#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"


/* Get values */
/* ========== */
double get(matrix *A, 
           int x, 
           int y){

    return A->data[x*A->columns + y];
}


/*  Set values  */
/* ============ */
void set(double value, 
         matrix *A, 
         int x, 
         int y){

    A->data[ x*A->columns + y] = value;
}



/*  Print matriz  */
/* ============== */
void display(matrix *A){
    for (int i = 0 ; i < A->lines ; i++){
        printf("|  ");        
        for (int j = 0 ; j < A->columns ; j++){
            printf("%f  ", get(A, i, j));
        }
        printf("|\n");
    }
}


/*  Create a new matrix  */
/* ===================== */
matrix* new_matrix(int x, 
                   int y){

    /* Alocação de memória (FAM) */ 
    matrix *e = NULL;
    if ( !(e = (matrix *)malloc(sizeof(matrix) + x*y*sizeof(double))) ){
        printf("Problema em criar matriz -- dentro da função 'new_matrix' \n");
    }
    else {
        /* Inicializa a estrutura de dados */
        for (int i = 0 ; i < x ; i++){
            for (int j = 0 ; j < y ; j++){
                 set(0.0, e, i, j);
            }
        }

        e->lines   = x;
        e->columns = y;  
    }
    return e;
}


/*  Create a new constant  */
/* ======================= */
matrix* new_constant(double x){
    
    /* Alocação de memória (FAM) */ 
    matrix *e = NULL;
    if ( !(e = (matrix *)malloc(sizeof(matrix) + sizeof(double))) ){
        printf("Problema em criar matriz -- dentro da função 'new_constant' \n");
    }
    else {
        e->data[0] = x;
        e->lines   = 1;
        e->columns = 1;  
    }
    return e;
}




/* Create a identity matrix of nxn deminsion */
/* ========================================= */
matrix* eye(int x){

    /* Alocação de memória (FAM) */ 
    matrix *e = NULL;
    if ( !(e = (matrix *)malloc(sizeof(matrix) + x*x*sizeof(double))) ){
        printf("Problema em criar matriz -- dentro da função 'eye' \n");
    }
    else {
        /* Inicializa a estrutura de dados */
        int cont = 0;   
        for (int i = 0 ; i < x ; i++){
            for (int j = 0 ; j < x ; j++){
                if (i == j){
                    e->data[cont] = 1.0;
                } else {
                    e->data[cont] = 0.0;
                }
                cont++;
            }
        }
        e->lines   = x;
        e->columns = x;  
    }
    return e;
}


/***********************************************************
                Produto entre matrizes. 
************************************************************/
void prod(matrix *A, 
          matrix *B, 
          matrix *Output){
       
    if ( ((A->lines == 1) && (A->columns == 1)) ){
        int n = B->lines * B->columns;

        /* Faz o produto entre os termos */
        for (int i = 0 ; i < n ; i++){
            Output->data[i] = A->data[0] * B->data[i];
        }
        Output->lines   = B->lines;
        Output->columns = B->columns; 
    } 
    else if( (B->lines == 1) && (B->columns == 1) ){
        int n = A->lines * A->columns;

        /* Faz o produto entre os termos */
        for (int i = 0 ; i < n ; i++){
            Output->data[i] = B->data[0] * A->data[i];
        }
        Output->lines   = A->lines;
        Output->columns = A->columns; 
    } 
    else{ /* Caso de 'A' e 'B' serem matrizes */    

        /* Verifica se é possível a multiplicação entre as matrizes */
        if ( A->columns != B->lines ){
            printf("Impossível realizar o produto entre as matrizes. \n");
        }

        /* Realiza a multiplicação */
        int contador = 0;
        double sum;
        for (int indexA = 0 ; indexA < A->lines ; indexA++){
            for (int indexB = 0 ; indexB < B->columns ; indexB++){
                sum = 0.0;
                for (int cont = 0 ; cont < A->columns ; cont++){
                    sum = sum + (A->data[indexA * A->columns + cont]) * (B->data[indexB + cont * B->columns]);
                }   
                Output->data[contador] = sum;
                contador++;
            }
        }
        Output->lines   = A->lines;
        Output->columns = B->columns;
    }
}


/***********************************************************
             Soma de matriz 'A' com matriz 'B' 
 ***********************************************************/
void add(matrix *A, 
         matrix *B, 
         matrix *Output){

    if ( (A->lines != B->lines) || (A->columns != B->columns) )
        printf("Dimensão diferente entre as matrizes \n");

    for (int i = 0 ; i < (A->lines * A->columns) ; i++){
        Output->data[i] =  A->data[i] + B->data[i];
    }
}


/***********************************************************
            Subtração de matriz 'A' com matriz 'B' 
 ***********************************************************/
void sub(matrix *A, 
         matrix *B, 
         matrix *Output){

    if ( (A->lines != B->lines) || (A->columns != B->columns) )
        printf("Dimensão diferente entre as matrizes \n");

    for (int i = 0 ; i < (A->lines * A->columns) ; i++){
        Output->data[i] =  A->data[i] - B->data[i];
    }
}


/***********************************************************
                   Transposição de matriz 'A'
 ***********************************************************/
void transpose(matrix *A,  
               matrix *Output){

    for (int i = 0 ; i < A->lines ; i++){
        for (int j = 0 ; j < A->columns ; j++){
            set(get(A, i, j), Output, j, i);
        }
    }
}


/***********************************************************
               Troca de linhas de matriz 'A'
 ***********************************************************/
void switch_lines(matrix *A, 
                  int line_1, 
                  int line_2){
    
    /* Código de proteção */
    if ( (line_1 >= A->lines) || (line_2 >= A->lines) ){
        printf("Erro na dimensão das matrizes para troca de linhas "); 
    }

    double valor_1, valor_2;
    for (int i = 0 ; i < A->columns ; i++){
        valor_1 = get(A, line_1, i);
        valor_2 = get(A, line_2, i);

        set(valor_2, A, line_1, i);
        set(valor_1, A, line_2, i);
    }
} 


/***********************************************************
                  Determinante da matriz 'A'
 ***********************************************************/
double det(matrix *A){

    /* Inicio da variável do determinante */
    double determinante = 1.0;

    /* Matriz identidade */
    matrix *I = eye(A->lines);

    /* Faz uma matriz ser igual a outra */
    isequal(I, A);

    /* Diagonalizar a matriz para execução das próximas etapas */
    double v1;
    double v2;
    for (int i = 0 ; i < A->lines ; i++){

        /* Subtrai os números subsequentes */
        for (int k = i+1 ; k < A->lines ; k++){

            /* pivor de referência para não ter divisão por zero */
            bool flag = true;
            int  cont = 0;
            while( flag ){
                v1 = get(I, i, i);
                if (v1 == 0.0){
                    cont++;
                    switch_lines(I, i, i+cont);
                    determinante *= -1.0;
                } else{
                    flag = false;              
                }
            }
            v2 = get(I, k, i);
            for (int j = 0 ; j < A->columns ; j++){
                set( get(I, k, j) - (v2/v1) * get(I, i, j), I, k, j);
            }
        }
    }

    /* Calcula o determinante */
    for (int i = 0; i < A->lines ; i++){
        determinante *= get(I, i, i);
    }
    
    /* Libera o espaço para a matriz I*/
    free(I);

    /* Retorna o parâmetro do determinante */
    return determinante;
}


/***********************************************************
                  Inversa de matriz 'A'
 ***********************************************************/
void inv(matrix *A, 
         matrix *Ainv){

    double det = 1.0;

        /* Nova Matriz I */
        matrix *I = eye(A->lines);
        
        if ((A->lines == 1) && (A->data[0] != 0.0) ){
            set(1.0/A->data[0], Ainv, 0, 0);
        } else{
            /* Diagonalizar a matriz 'I' para execução das próximas etapas */
            double v1;
            double v2;
            for (int i = 0 ; i < A->lines ; i++){

                /* Subtrai os números subsequentes */
                for (int k = i+1 ; k < A->lines ; k++){

                    /* Pivor de referência para não ter divisão por zero */
                    bool flag = true;
                    int  cont = 0;
                    while( flag ){
                        v1 = get(I, i, i);
                        if (v1 == 0){
                            cont++;
                            switch_lines(   I, i, i+cont);
                            switch_lines(Ainv, i, i+cont);
                            det *= -1.0;
                        } else{
                            flag = false;              
                        }
                    }
                    v2 = get(I, k, i);
                    for (int j = 0 ; j < A->columns ; j++){
                        set( get(   I, k, j) - (v2/v1) * get(   I, i, j),    I, k, j);
                        set( get(Ainv, k, j) - (v2/v1) * get(Ainv, i, j), Ainv, k, j);
                    }
                }
            }


            /* Calcula o determinante */
            for (int i = 0; i < A->lines ; i++){
                det *= get(I, i, i);
            }

            /* Se o determinante não for nulo, existe inversa */
            if (det != 0.0){

                /* Normaliza as matrizes e subtrai os números subsequentes */
                double pivor;
                double valor_ref;
                for (int i = 0 ; i < A->lines ; i++){
                    pivor = get(I, i, i);              

                    /* Divide todo mundo pelo pivor */
                    for (int j = 0 ; j < A->columns ; j++){    
                        set(get(   I, i, j)/pivor,    I, i, j);
                        set(get(Ainv, i, j)/pivor, Ainv, i, j);
                    }
                }

                /* Faz os cálculos restantes para obtenção da matriz identidade de um lado */  
                for (int i = 0 ; i < (A->lines - 1); i++){
                    for (int k = (i + 1) ; k < (A->columns) ; k++){ 
                        valor_ref = get(I, i, k);
                        for (int j = 0 ; j < A->columns ; j++){
                            set( get(   I, i, j) - valor_ref * get(   I, k, j),    I, i, j);   
                            set( get(Ainv, i, j) - valor_ref * get(Ainv, k, j), Ainv, i, j);
                        }
                    }
                }
            } 
            else {
                printf("\n\n Não há matriz inversa!!!\n\n");
            }

            /* Libera o espaço para a matriz I*/
            free(I); 
        }
}


/***********************************************************
                Carrega a matriz B na matriz A
 ***********************************************************/
void isequal(matrix *A, 
             matrix *B){

    int n = B->lines * B->columns;
    for (int i = 0 ; i < n ; i++){
        A->data[i] = B->data[i];
    }
    A->lines   = B->lines;
    A->columns = B->columns;
}
































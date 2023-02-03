#ifndef _MATRIX_H_
#define _MATRIX_H_


/***********************************************************
    Macros a serem empregadas na minha biblioteca  
************************************************************/
#define sign(x)  (x > 0) ? (1) : (-1)


/***********************************************************
    Enumerate a ser empregado em várias partes do código  
************************************************************/
typedef enum { false = 0, 
               true  = 1} bool;

/***********************************************************
               Estrutura para matriz de dados  
************************************************************/
typedef struct {
  int lines;
  int columns;
  double data[];
} matrix;


/***********************************************************
                   Protótipos de funções
************************************************************/
double get(matrix *A, int i, int j);
void set(double num, matrix *M, int i, int j);
void display(matrix *A);

matrix* new_matrix(int x, int y);
matrix* new_constant(double x);
matrix* eye(int x);

void prod(matrix *A, matrix *B, matrix *Output);
void add( matrix *A, matrix *B, matrix *output);
void sub( matrix *A, matrix *B, matrix *output);
void transpose( matrix *A, matrix *Output);
void switch_lines(matrix *A, int line_1, int line_2);
double det(matrix *A);
void inv(matrix *A, matrix *Ainv);
void isequal(matrix *A, matrix *B);


#endif

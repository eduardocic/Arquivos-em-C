#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "matrix.h"

void save(const char* nameFile, const int num, ...)
{
	/* ==============================================================
	   Objetos do tipo 'va_list' (lista variádica) é necessário para 
	   guardar as informações necessárias para as macros:
	     
	     -- va_start;
	     -- va_copy;
	     -- va_arg; e
	     -- va_end.
	 ==============================================================  */
    va_list args;
    va_start(args, num);
 
 
	/* Criei o nome o arquivo */ 
 	FILE *fp; 
 	fp = fopen(nameFile, "w");
 	
 	for (int i = 0 ; i < num ; i++){
		if (va_arg(args, matrix *)) {
			printf("Imprime matriz\n\n");
		}
		if (va_arg(args, int)) {
			printf("Imprime inteiro\n\n");
		}
    }
    va_end(args);
}
 
int main(void)
{
	matrix *a = new_constant(5.0);
//    save("m.txt", 1, a); 
    save("m.txt", 1, 1);
}











//void save(const char* fmt, ...)
//{
//	/* ==============================================================
//	   Objetos do tipo 'va_list' (lista variádica) é necessátio para 
//	   guardar as informações necessárias para as macros:
//	     
//	     -- va_start;
//	     -- va_copy;
//	     -- va_arg; e
//	     -- va_end.
//	 ==============================================================  */
//    va_list args;
//    
//    
//    va_start(args, fmt);
// 
//    while (*fmt != '\0') {
//    
//        if (*fmt == 'd') {
//            int i = va_arg(args, int);
//            printf("%d\n", i);
//        } else if (*fmt == 'c') {
//            // A 'char' variable will be promoted to 'int'
//            // A character literal in C is already 'int' by itself
//            int c = va_arg(args, int);
//            printf("%c\n", c);
//        } else if (*fmt == 'f') {
//            double d = va_arg(args, double);
//            printf("%f\n", d);
//        }
//        ++fmt;
//    }
// 
//    va_end(args);
//}















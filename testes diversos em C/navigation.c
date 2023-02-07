#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "matrix.h"
#include "navigation.h"



void rotZ(double x, matrix *R){
/* ================================================ 
 * Equação 2.51, pág. 35 
 *
 * O que esta equação faz é uma rotação no entorno
 * do eixo Z.
 * 
 *  VARIÁVEIS
 * -----------  
 *
 *  1. x ............: ângulo de rotação (rad); e
 *  2. R ............: matrix resultado da rotação
 *					   dada (matrix 3x3).
 * =============================================== */
	/* Pega os parâmetros do sistema */
	double cosx, sinx;
	cosx = cos(x);
	sinx = sin(x);

	/* Linha 1, 2 e 3 de W */
	set(  cosx, R, 0, 0);
	set(  sinx, R, 0, 1);
	set(   0.0, R, 0, 2);
	set( -sinx, R, 1, 0);
	set(  cosx, R, 1, 1);
	set(   0.0, R, 1, 2);
	set(   0.0, R, 2, 0);
	set(   0.0, R, 2, 1);
	set(   1.0, R, 2, 2);
}





void rotX(double x, matrix *R){
/* ================================================ 
 * Equação 2.53, pág. 36 
 *
 * O que esta equação faz é uma rotação no entorno
 * do eixo X.
 * 
 *  VARIÁVEIS
 * -----------
 *  
 *  1. x ............: ângulo de rotação (rad); e
 *  2. R ............: matrix resultado da rotação
 *					   dada (matrix 3x3).
 * =============================================== */
	/* Pega os parâmetros do sistema */
	double cosx, sinx;
	cosx = cos(x);
	sinx = sin(x);

	/* Linha 1, 2 e 3 de W */
	set(   1.0, R, 0, 0);
	set(   0.0, R, 0, 1);
	set(   0.0, R, 0, 2);
	set(   0.0, R, 1, 0);
	set(  cosx, R, 1, 1);
	set(  sinx, R, 1, 2);
	set(   0.0, R, 2, 0);
	set( -sinx, R, 2, 1);
	set(  cosx, R, 2, 2);
}





void rotY(double x, matrix *R){
/* ================================================ 
 * Equação 2.55, pág. 36 
 *
 * O que esta equação faz é uma rotação no entorno
 * do eixo Y.
 * 
 *  VARIÁVEIS
 * -----------
 *  
 *  1. x ............: ângulo de rotação (rad); e
 *  2. R ............: matrix resultado da rotação
 *					   dada (matrix 3x3).
 * =============================================== */	
	/* Pega os parâmetros do sistema */
	double cosx, sinx;
	cosx = cos(x);
	sinx = sin(x);

	/* Linhx 1, 2 e 3 de W */
	set(  cosx, R, 0, 0);
	set(   0.0, R, 0, 1);
	set( -sinx, R, 0, 2);
	set(   0.0, R, 1, 0);
	set(   1.0, R, 1, 1);
	set(   0.0, R, 1, 2);
	set(  sinx, R, 2, 0);
	set(   0.0, R, 2, 1);
	set(  cosx, R, 2, 2);
}






void skew_symmetric(matrix *w, matrix *W){
	/* Separa os parâmetros de w */
	double wx, wy, wz;
	wx = get(w, 0, 0);
	wy = get(w, 1, 0);
	wz = get(w, 2, 0);
	
	/* Linha 1, 2 e 3 de W */
	set(  0.0, W, 0, 0);
	set(  -wz, W, 0, 1);
	set(   wy, W, 0, 2);
	set(   wz, W, 1, 0);
	set(  0.0, W, 1, 1);
	set(  -wx, W, 1, 2);
	set(  -wy, W, 2, 0);
	set(   wx, W, 2, 1);
	set(  0.0, W, 2, 2);
}






double Rn(double LAT){
/* =======================================
 * 	  Raio normal (Eq. 2.109, pág. 48)	
 *
 *  VARIÁVEL
 * ----------
 *
 *  1. LAT .........: latitude (rad).	
 * ====================================== */
	double e2   = WGS84_e * WGS84_e;
	double sin2 = sin(LAT) * sin(LAT);	
	
	double den  = sqrt(1.0 - e2 * sin2);
	double res  = WGS84_a/den;
	
	return res; 	
}






double Rm(double LAT){
/* =======================================
 *   Raio meridional (Eq. 2.110, pág. 48)	
 *
 *  VARIÁVEL
 * ----------
 *
 *  1. LAT .........: latitude (rad).
 * ====================================== */
	double e2   = WGS84_e * WGS84_e;
	double sin2 = sin(LAT) * sin(LAT);	
	
	double num  = WGS84_a * (1.0 - e2);
	double den  = pow((1.0 - e2 * sin2), 1.5);
	double res  = num/den;
	
	return res; 	
}





void lla2ecef(double LAT, double LONG, double H, matrix *geod){
/* ========================================================
 * Esta função faz a conversão das coordenadas Geodésica do 
 * ECEF (representada pelas informações LATITUDE, LOGITUDE e
 * ALTITUDE) para os coordenadas Retangulares (x, y e z), 
 * aos moldes de (Eq. 2.111), na pág. 50.
 *
 *  (*) Ler Referência a partir da pág. 49 até 52.	
 *
 *  VARIÁVEIS
 * -----------
 *
 *  1. LAT .........: latitude (rad);
 *  2. LONG ........: longitude (rad);
 *  3. H ...........: altitude (m); e
 *  4. geod ........: vetor coluna com a conversão das 
 *					  coordenadas ECEF LLA para coorde-
 * 					  nadas ECEF retangulares. Todas as
 * 					  variáveis do geod estão expressas
 *					  em metros (m).
 * ====================================================== */
 
 	double rn = Rn(LAT);
 	double e2 = WGS84_e * WGS84_e;

	set(  (rn + H) * cos(LAT) * cos(LONG), geod, 0, 0);
	set(  (rn + H) * cos(LAT) * sin(LONG), geod, 1, 0);
	set( (rn * (1.0 - e2) + H) * sin(LAT), geod, 2, 0);
}




void ecef2lla(double x, double y, double z, matrix *LLA){
/* ========================================================
 * Esta função faz a conversão das coordenadas Retangulares 
 * do ECEF (x, y, z) para os coordenadas geodésicas (LATI-
 * TUDE, LOGITUDE e ALTITUDE aos moldes do procedimento
 * de iteração numérica disposto na seção 2.5.4 do livro
 * (a partir da página 50).
 *
 *  (*) Ler Referência a partir da pág. 49 até 52.	
 *
 *  VARIÁVEIS
 * -----------
 *
 *  1. x ........: coordenada retangular x (m);
 *  2. y ........: coordenada retangular y (m);
 *  3. z ........: coordenada retangular z (m); e
 *  4. LLA ......: vetor coluna com a conversão das coor-
 *				   denadas ECEF retangular para ECEF LLA.
 * 				   As variáveis de latitude e longitude
 *				   estão expressas em radiano (rad) e a 
 *	 			   altura em metros (m).
 * ====================================================== */
	/* Captura dos parâmetros do WBS-84 */
	double e = WGS84_e;
	double a = WGS84_a;

	/* Inicialização */
	double LONG = atan2(y, x);
	double h    = 0.0;
	double rn   = a;
	double LAT  = atan(z / ((Rn(0.0) + h) * (1.0 - e*e) ));	
		
	int cont = 0;
	while ( cont < 10 ){
		/* Iteratividade */
		rn   = a / sqrt(1.0 - e*e*sin(LAT)*sin(LAT));
		h    = sqrt(x*x + y*y)/cos(LAT) - rn;
	    LAT  = z / (sqrt(x*x + y*y) * (1.0 - (e*e)*rn/(rn+h)));
			
		/* Verifica o resultado */
		/* printf("Cont: %d\n", cont);	*/
		/* printf("LAT, LONG, h: %f, %f, %f \n", LAT, LONG, h); */
		cont++;
	}
	/* printf("Tolerância alcançada! \n"); */
	
	set(LAT , LLA, 0, 0);
	set(LONG, LLA, 1, 0);
	set(h   , LLA, 2, 0);
	display(LLA);
}






void llt2ecef(double LAT, double LONG, matrix *Rel){
/* =======================================================
 * Rotação que leva um corpo do Local Level Frame (LLT)
 * para o Earth Center Earth Fixed (ECEF). Essa trans-
 * formação busca alinhar o 'l-frame' com o 'e-frame', 
 * em que o l-frame é a referência para o LLT e o 
 * 'e-frame' para o ECEF.
 * 
 * --> Equação (2.70), página 40
 * ======================================================= */

	double cLAT   = cos(LAT);
	double cLONG  = sin(LONG);
	double sLAT   = sin(LAT);
	double sLONG  = sin(LONG);
	
	
	set(        -sLONG, Rel, 0, 0);
	set( -sLAT * cLONG, Rel, 0, 1);
	set(  cLAT * cLONG, Rel, 0, 2);
	set(  		 cLONG, Rel, 1, 0);
	set( -sLAT * sLONG, Rel, 1, 1);
	set(  cLAT * cLONG, Rel, 1, 2);
	set(  		   0.0, Rel, 2, 0);
	set(          cLAT, Rel, 2, 1);
	set(  		  sLAT, Rel, 2, 2);
}




void body2llf(double p, double r, double y, matrix *Rlb){
/* =======================================================
 * Rotação que leva as coordenadas do corpo (body) para o 
 * Local Level Frame. Essa transformação busca alinhar o 
 * 'b-frame' com o 'e-frame', 
 * em que o b-frame é a referência para o BODY (corpo) e o 
 * 'l-frame' para o LLT.
 * 
 * --> Equação (2.82), pág. 42
 *
 * As entradas dizem respeito ao 'pitch', 'roll' e 'yaw', 
 * todas elas expressas em radiano (rad), em que:
 *
 *  --> pitch é simbolizado por 'theta';
 *  --> roll é simbolizado por 'phi'; e
 *  --> yaw é simbolizado por 'psi'.
 * ======================================================= */
	double cosp   = cos(p);
	double cosr   = cos(r);
	double cosy   = cos(y);
	double sinp   = sin(p);
	double sinr   = sin(r);
	double siny   = sin(y);
	
	set( (cosp*cosr - siny*sinp*sinr), Rlb, 0, 0);
	set(     			   -siny*cosp, Rlb, 0, 1);
	set( (cosy*sinr + siny*sinp*cosr), Rlb, 0, 2);
	set( (siny*cosr + cosy*sinp*sinr), Rlb, 1, 0);
	set(				  (cosy*cosp), Rlb, 1, 1);
	set( (siny*sinr - cosy*sinp*cosr), Rlb, 1, 2);
	set(    	   		   -cosp*sinr, Rlb, 2, 0);
	set(          			     sinp, Rlb, 2, 1);
	set(  		            cosp*cosr, Rlb, 2, 2);
}





void eci2ecef(double t, matrix *Rei){
/* =======================================================
 * Rotação que leva as coordenadas ECI (Earth Centered 
 * Inertial) para o ECEF (Earth Centered Earth Fixed). Ou
 * seja, faz a conversão de coordenadas entre o 'i-frame'
 * em que o 'e-frame'
 * 
 * --> Equação (2.65), pág. 39.
 * ======================================================= */


	set(  cos(WGS84_omega_e * t), Rei, 0, 0);
	set(  sin(WGS84_omega_e * t), Rei, 0, 1);
	set( 				     0.0, Rei, 0, 2);
	set( -sin(WGS84_omega_e * t), Rei, 1, 0);
	set(  cos(WGS84_omega_e * t), Rei, 1, 1);
	set( 					 0.0, Rei, 1, 2);
	set(    	   		     0.0, Rei, 2, 0);
	set(          			 0.0, Rei, 2, 1);
	set(  		             1.0, Rei, 2, 2);
	
}


void dot_Rmk (matrix *w, matrix *Rmk, matrix *dot_Rmk){
/* ==========================================================
 *
 * 		Derivada temporal da matriz de transformação
 *     ----------------------------------------------
 *
 * Existe a explicação desse assunto em:
 * 
 * (1) Fundamentals of Inertial Navigation, Satellite-based
 *     Positioning and their Integration (pág. 43 a 45); e
 * (2) Applied mathematics in integrated navigation systems
 *     (pág. 25 a 28) --- MELHOR EXPLICADO
 *
 * A simbologia que irei utilizar aqui segue as nomenclatura 
 * da primeira referência (pág. 45 da Ref. 1).
 *
 * --> O que essa matriz oferece é a taxa de rotação de um
 *     determinado frame (k-frame) em relação aos eixos
 *	   x, y e z de um outro frame (m-frame).
 *
 *            +----------------------------+
 *      	  | dot_Rmk = Rmk * Omega_k_mk |
 *            +----------------------------+
 *
 * ======================================================== */
 
	/* Criação da matrix do tipo antissimétrica */
	matrix *Omega_k_mk = new_matrix(3, 3);
	skew_symmetric(w, Omega_k_mk);

	/* Calcula o produto entre as matrizes */
	prod(Rmk,  Omega_k_mk, dot_Rmk);
	
	/* Libera a memória do heap */
	free(Omega_k_mk);
}

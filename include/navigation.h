#ifndef _NAVIGATION_H_
#define _NAVIGATION_H_

/* ==================================================================
 * 
 *  REFERÊNCIA
 * ------------
 *  
 * --> 'Fundamentals of Inertial Navigation, Satellite-based Positio-
 *     ning and their Integration', 
 * --> Autores: Aboelmagd Noureldin, Tashfeen B. Karamat e Jacques 
 *	   Georgy
 *
 *
 *  INFORMAÇÔES IMPORTANTES
 * -------------------------	   
 *
 * --> vetores são simbolizados por um vetor COLUNA;
 * -->
 *
 * =============================================================== */
#define  PI  	       3.141592653589793
#define  RAD2DEG(x)	   (x * 57.295779513082320)
#define  DEG2RAD(x)	   (x *  0.017453292519943)    

/* ============================= */
/*      Constantes diversas      */
/*								 */
/* --> pag 47 do livro.		     */
/* ============================= */ 
#define WGS84_omega_e  0.00007292115		/* velocidade de rotação da terra (rad/s)  */ 
#define WGS84_a        6378137.0		    /* semieixo maior - raio equatorial (em m) */
#define WGS84_GM       398600441800000      /* constante gravitacional (em m^3/s^2)    */
#define WGS84_f 	   0.00335281       	/* valor de 'f'							   */
#define WGS84_b        6356752.3142         /* semieixo menor (em m)				   */
#define WGS84_e        0.08181919   	    /* excentricidade 					       */

/* ============================ */
/*   Derivações para o WGS-84   */
/* ============================ */
double Rn(double Lat);					 	                     /* raio normal (em m)					   */
double Rm(double Lat);						                     /* raio meridional (em m)				   */
void lla2ecef(double Lat, double Long, double H, matrix *geod);  /* Conversão lat-long-H para coordenada 
																    retangular no ECEF                     */
void ecef2lla(double x, double y, double z, matrix *LLA);		 /* Conversão coordenada retangular x-y-z 
																    em coordenada geodésica no ECEF        */



/* ============================ */
/*     Álgebra linear básica    */
/* ============================ */ 
void rotZ(double ang, matrix *R);			/* rotação no entorno do eixo Z 		   */
void rotX(double ang, matrix *R);			/* rotação no entorno do eixo X 		   */
void rotY(double ang, matrix *R);			/* rotação no entorno do eixo Y 		   */
void skew_symmetric(matrix *w, matrix *W);  /* matrix antissimétrica 				   */


/* ============================ */
/*  	Rotação entre frames    */
/* ============================ */ 
void llt2ecef(double Lat, double Long, matrix *R_e_l);
void body2llf(double p, double r, double y, matrix *R_l_b);
void eci2ecef(double t, matrix *Rei);
/* Falta estudar direitinho sobre o frame de Wander */


/* ================================= */
/*    Taxa de Rotação entre frames   */
/* ================================= */ 
void dot_Rmk (matrix *w, matrix *Rmk, matrix *dot_Rmk);

#endif

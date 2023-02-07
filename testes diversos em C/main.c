#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "navigation.h"

int main(){

	/* teste da primeira rotina */
	matrix *lla  = new_matrix(3, 1);
	double x = 4510731;
	double y = 4510731;
	double z = 0;
	
	ecef2lla(x, y, z, lla);
	printf("lat: %f, long: %f, h: %f\n\n", get(lla, 0, 0)*180/PI, get(lla, 1, 0)*180/PI, get(lla, 2, 0));
	
	
	/* teste da primeira rotina */
	matrix *geod = new_matrix(3, 1);
	double LAT  = DEG2RAD(45.0);
	double LONG = DEG2RAD(90.0);
	double H    = 2000;
	
	lla2ecef(LAT, LONG, H, geod);
	printf("x: %f, y: %f, z: %f\n\n", get(geod, 0, 0), get(geod, 1, 0), get(geod, 2, 0));
	
	return 0;
}

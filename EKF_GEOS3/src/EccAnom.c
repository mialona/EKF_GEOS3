/** @file EccAnom.c
 *  @brief Eccentric anomaly for elliptic orbits.
 *
 *  This driver contains the code for the 
 *  computation of the eccentric anomaly for elliptic orbits.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double EccAnom(double M, double e) {
	int maxit = 15;
	int i = 1;

	// Starting value
	M = fmod(M, 2.0*M_PI);
	if(M < 0)
		M = M + 2.0*M_PI;

	double E;
	if(e<0.8) {
		E = M; 
	}
	else {
		E = M_PI;
	}

	double f = E - e*sin(E) - M;
	E = E - f / (1.0 - e*cos(E));

	// Iteration
	double eps = 2.220446049250313e-14;
	while(fabs(f) > eps) {
		f = E - e*sin(E) - M;
		E = E - f / (1.0 - e*cos(E));
		i = i+1;
		if(i == maxit) {
			printf("Convergence problems in EccAnom\n");
			exit(EXIT_FAILURE);
		}
	}
	
	return E;
}


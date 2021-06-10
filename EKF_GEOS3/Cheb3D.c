/** @file Cheb3D.c
 *  @brief Chebyshev approximation of 3-dimensional vectors.
 *
 *  This driver contains the code for the calculation of the
 *  Chebyshev approximation of 3-dimensional vectors.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "m_utils.h"

#include <stdio.h>
#include <stdlib.h>


double *Cheb3D(double t, int N, double Ta, double Tb, double *Cx, double *Cy, double *Cz) {
	// Check validity
	if((t<Ta) || (Tb<t)) {
		printf("ERROR: Time out of range in Cheb3D::Value\n");
        exit(EXIT_FAILURE);
	}

	// Clenshaw algorithm
	double tau = (2.0*t-Ta-Tb)/(Tb-Ta);  

	double *f1 = v_create(3);
	double *f2 = v_create(3);
	double *old_f1 = v_create(3);

	for(int i = N-1; i>0; i--) {
		old_f1[0] = f1[0];
		old_f1[1] = f1[1];
		old_f1[2] = f1[2];
		
		f1[0] = 2.0*tau*f1[0]-f2[0]+Cx[i];
		f1[1] = 2.0*tau*f1[1]-f2[1]+Cy[i];
		f1[2] = 2.0*tau*f1[2]-f2[2]+Cz[i];
		
		f2[0] = old_f1[0];
		f2[1] = old_f1[1];
		f2[2] = old_f1[2];
	}

	double *ChebApp = v_create(3);
	ChebApp[0] = tau*f1[0]-f2[0]+Cx[0];
	ChebApp[1] = tau*f1[1]-f2[1]+Cy[0];
	ChebApp[2] = tau*f1[2]-f2[2]+Cz[0];
	
	return ChebApp;
}


/** @file G_AccelHarmonic.c
 *  @brief Gradient of the Earth's harmonic gravity field.
 *
 *  This driver contains the code for the computation of the
 *  gradient of the Earth's harmonic gravity field.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"
#include "../includes/AccelHarmonic.h"

#include <stdio.h>
#include <math.h>


double **G_AccelHarmonic(double *r, double **U, int n_max, int m_max) {
	double d = 1.0;   // Position increment [m]

	double **G = m_zeros(3,3);
	double *dr = v_create(3);
	double *da;

	// Gradient
	for(int i=0; i<3; i++) {
		// Set offset in i-th component of the position vector
		dr[0] = 0.0;
		dr[1] = 0.0;
		dr[2] = 0.0;
		dr[i] = d;
		// Acceleration difference
		da = v_sum(AccelHarmonic(v_sum(r,3,v_mul_scalar(dr,3,1/2.0),3),U,n_max,m_max),3,v_mul_scalar(AccelHarmonic(v_sum(r,3,v_mul_scalar(dr,3,-1/2.0),3),U,n_max,m_max),3,-1),3);
		// Derivative with respect to i-th axis
		G = v_assign(v_mul_scalar(da,3,1/d),3,G,3,3,i);
	}
	
	v_free(dr,3);
	v_free(da,3);
	
	return G;
}


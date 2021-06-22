/** @file AccelPointMass.c
 *  @brief Perturbational acceleration due to a point mass.
 *
 *  This driver contains the code for the computation of the
 *  perturbational acceleration due to a point mass.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"

#include <stdio.h>
#include <math.h>


double *AccelPointMass(double *r, double *s, double GM) {
	// Relative position vector of satellite w.r.t. point mass 
	double *d = v_sum(r,3,v_mul_scalar(s,3,-1),3);

	// Acceleration 
	double *d_aux = v_mul_scalar(d,3,-GM/(pow(v_norm(d,3),3.0)));
	double *v_aux = v_mul_scalar(s,3,-GM/(pow(v_norm(s,3),3.0)));
	return v_sum(d_aux,3,v_aux,3);
}

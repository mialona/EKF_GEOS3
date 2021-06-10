/** @file MeasUpdate.c
 *  @brief Meas updater.
 *
 *  This driver contains the code for the 
 *  meas updater.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug m = 1 and n = 6 have been taken.
 */

#include "../includes/m_utils.h"

#include <stdio.h>


void MeasUpdate(double z, double g, double s, double *G, double **K, double **x, double ***P) {
	double **Inv_W = m_zeros(1,1);

	Inv_W[0][0] = s*s;    // Inverse weight (measurement covariance)
	
	// Kalman gain
	double aux = v_dot(v_dot_m(G,6,*P,6,6),6,G,6);
	Inv_W[0][0] += aux;
	*K = v_mul_scalar(m_dot_v(*P,6,6,G,6),6,m_inv(Inv_W,1)[0][0]);

	// State update
	*x = v_sum(*x,6,v_mul_scalar(*K,6,(z-g)),6);
	
	// Covariance update
	double **eye = m_eye(6);
	for(int i=0; i<6; i++) {
		for(int j=0; j<6; j++) {
			eye[i][j] -= (*K)[i]*G[j];
		}
	}
	*P = m_dot(eye,6,6,*P,6,6);
}


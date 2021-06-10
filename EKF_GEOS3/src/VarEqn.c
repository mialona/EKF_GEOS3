/** @file VarEqn.c
 *  @brief Variational equations.
 *
 *  This driver contains the code for the computation of the
 *  variational equations.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/global.h"
#include "../includes/const.h"
#include "../includes/IERS.h"
#include "../includes/timediff.h"
#include "../includes/PrecMatrix.h"
#include "../includes/NutMatrix.h"
#include "../includes/m_utils.h"
#include "../includes/PoleMatrix.h"
#include "../includes/GHAMatrix.h"
#include "../includes/AccelHarmonic.h"
#include "../includes/G_AccelHarmonic.h"

#include <stdio.h>
#include <math.h>


void VarEqn(double x, double *yPhi, double **yPhip) {
	extern Param AuxParam;
	
	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	IERS(AuxParam.Mjd_UTC,'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);
	
	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
	timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);
	
	double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400;
	
	// Transformation matrix
	double **P = PrecMatrix((MJD_J2000),AuxParam.Mjd_TT + x/86400.0);
	double **N = NutMatrix(AuxParam.Mjd_TT + x/86400.0);
	double **T = m_dot(N,3,3,P,3,3);
	double **E = m_dot(m_dot(PoleMatrix(x_pole,y_pole),3,3,GHAMatrix(Mjd_UT1),3,3),3,3,T,3,3);
	
	// State vector components
	double *r = v_create(3);
	double *v = v_create(3);
	for(int i=0; i<3; i++) {
		r[i] = yPhi[i];
		v[i] = yPhi[i+3];
	}
	double **Phi = m_zeros(6,6);

	// State transition matrix
	for(int j=0; j<6; j++) {
		for(int i=0; i<6; i++) {
			Phi[i][j] = yPhi[6*(j+1)+i];
		}
	}
	
	// Acceleration and gradient
	double *a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
	double **G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
	
	// Time derivative of state transition matrix
	*yPhip = v_create(42);
	double **dfdy = m_zeros(6,6);

	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			dfdy[i][j] = 0.0;                 // dv/dr(i,j)
			dfdy[i+3][j] = G[i][j];           // da/dr(i,j)
			if(i==j) {
				dfdy[i][j+3] = 1;
			}
			else {
				dfdy[i][j+3] = 0;             // dv/dv(i,j)
			}
			dfdy[i+3][j+3] = 0.0;             // da/dv(i,j)
		}
	}

	double **Phip = m_dot(dfdy,6,6,Phi,6,6);
	// Derivative of combined state vector and state transition matrix
	for(int i=0; i<3; i++) {
		(*yPhip)[i]   = v[i];                 // dr/dt(i)
		(*yPhip)[i+3] = a[i];                 // dv/dt(i)
	}

	for(int i=0; i<6; i++) {
		for(int j=0; j<6; j++) {
			(*yPhip)[6*(j+1)+i] = Phip[i][j];     // dPhi/dt(i,j)
		}
	}
}


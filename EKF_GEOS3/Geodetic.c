/** @file Geodetic.c
 *  @brief Geodetic coordinates.
 *
 *  This driver contains the code for the 
 *  Geodetic coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "const.h"
#include "m_utils.h"

#include <stdio.h>
#include <math.h>


void Geodetic(double *r, double *lon, double *lat, double *h) {
	double R_equ = R_Earth;
	double f     = f_Earth;

	double eps = 2.22044604925031e-16;
	double epsRequ = (eps)*R_equ;        // Convergence criterion
	double e2      = f*(2.0-f);        // Square of eccentricity

	double X = r[0];                   // Cartesian coordinates
	double Y = r[1];
	double Z = r[2];
	double rho2 = X*X + Y*Y;           // Square of distance from z-axis

	// Check validity of input data
	if(v_norm(r,3)==0.0) {
		*lon = 0.0;
		*lat = 0.0;
		*h   = -R_Earth;
	}

	// Iteration 
	double dZ = e2*Z;
	
	double ZdZ, Nh, SinPhi, N, dZ_new;
	while(1) {
		ZdZ    =  Z + dZ;
		Nh     =  sqrt(rho2 + ZdZ*ZdZ); 
		SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
		N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
		dZ_new =  N*e2*SinPhi;
		if(fabs(dZ-dZ_new) < epsRequ) {
			break;
		}
		dZ = dZ_new;
	}

	// Longitude, latitude, altitude
	*lon = atan2(Y, X);
	*lat = atan2(ZdZ, sqrt(rho2));
	*h   = Nh - N;
}

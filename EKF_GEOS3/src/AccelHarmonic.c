/** @file AccelHarmonic.c
 *  @brief Acceleration due to the harmonic gravity field
 *  of the  central body.
 *
 *  This driver contains the code for the computation of the
 *  acceleration due to the harmonic gravity field of the 
 *  central body.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/global.h"
#include "../includes/Legendre.h"
#include "../includes/m_utils.h"

#include <stdio.h>
#include <math.h>


double *AccelHarmonic(double *r, double **E, int n_max, int m_max) {
	extern double **Cnm, **Snm;
	
	double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
	double gm    = 398600.4415e9; // [m^3/s^2]; GGM03S
	
	// Body-fixed position 
	double *r_bf = m_dot_v(E,3,3,r,3);

	// Auxiliary quantities
	double d = v_norm(r_bf,3);                     // distance
	double latgc = asin(r_bf[2]/d);
	double lon = atan2(r_bf[1],r_bf[0]);
	
	double **pnm, **dpnm;
	Legendre(n_max,m_max,latgc,&pnm,&dpnm);
	
	double dUdr = 0;
	double dUdlatgc = 0;
	double dUdlon = 0;
	double q3 = 0, q2 = q3, q1 = q2;
	double b1,b2,b3;
	for(int n=0; n<=n_max; n++) {
		b1 = (-gm/(d*d))*pow((r_ref/d),n)*(n+1);
		b2 =  (gm/d)*pow((r_ref/d),n);
		b3 =  (gm/d)*pow((r_ref/d),n);
		for(int m=0; m<=m_max; m++) {
			q1 = q1 + pnm[n+1][m+1]*(Cnm[n+1][m+1]*cos(m*lon)+Snm[n+1][m+1]*sin(m*lon));
			q2 = q2 + dpnm[n+1][m+1]*(Cnm[n+1][m+1]*cos(m*lon)+Snm[n+1][m+1]*sin(m*lon));
			q3 = q3 + m*pnm[n+1][m+1]*(Snm[n+1][m+1]*cos(m*lon)-Cnm[n+1][m+1]*sin(m*lon));
		}
		dUdr     = dUdr     + q1*b1;
		dUdlatgc = dUdlatgc + q2*b2;
		dUdlon   = dUdlon   + q3*b3;
		q3 = 0; q2 = q3; q1 = q2;
	}
	
	// Body-fixed acceleration
	double r2xy = r_bf[0]*r_bf[0]+r_bf[1]*r_bf[1];

	double ax = (1/d*dUdr-r_bf[2]/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf[0]-(1/r2xy*dUdlon)*r_bf[1];
	double ay = (1/d*dUdr-r_bf[2]/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf[1]+(1/r2xy*dUdlon)*r_bf[0];
	double az =  1/d*dUdr*r_bf[2]+sqrt(r2xy)/(d*d)*dUdlatgc;

	double *a_bf = v_create(3);
	a_bf[0] = ax; a_bf[1] = ay; a_bf[2] = az;

	// Inertial acceleration
	return m_dot_v(m_trans(E,3),3,3,a_bf,3);
}

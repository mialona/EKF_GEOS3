/** @file elements.c
 *  @brief Osculating Keplerian elements from the satellite state
 *  vector for elliptic orbits.
 *
 *  This driver contains the code for the computation of the 
 *  Osculating Keplerian elements from the satellite state
 *  vector for elliptic orbits.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "const.h"
#include "m_utils.h"

#include <stdio.h>
#include <math.h>


void elements(double *y, double *p, double *a, double *e, double *i,
			  double *Omega, double *omega, double *M) {
	double *r = v_create(3);                           // Position
	r[0] = y[0];
	r[1] = y[1];
	r[2] = y[2];
	double *v = v_create(3);                           // Velocity
	v[0] = y[3];
	v[1] = y[4];
	v[2] = y[5];

	double *h = v_cross(r,3,v,3);                          // Areal velocity
	double magh = v_norm(h,3);
	*p = magh*magh/(GM_Earth);

	*Omega = atan2(h[0],-h[1]);                        // Long. ascend. node 
	*Omega = fmod(*Omega,pi2);
	if(*Omega < 0)
		*Omega = *Omega + pi2;
	*i     = atan2(sqrt(h[0]*h[0]+h[1]*h[1]), h[2]);   // Inclination        
	double u = atan2(r[2]*magh, -r[0]*h[1]+r[1]*h[0]); // Arg. of latitude   

	double R = v_norm(r,3);                            // Distance           

	*a = 1/(2.0/R-v_dot(v,3,v,3)/(GM_Earth));            // Semi-major axis

	double eCosE = 1-R/(*a);                           // e*cos(E)           
	double eSinE = v_dot(r,3,v,3)/sqrt((GM_Earth)*(*a)); // e*sin(E)

	double e2 = eCosE*eCosE + eSinE*eSinE;
	*e = sqrt(e2);                                     // Eccentricity 
	double E = atan2(eSinE,eCosE);                     // Eccentric anomaly  

	*M = fmod(E-eSinE,pi2);                            // Mean anomaly
	if(*M < 0)
		*M = *M + pi2;

	double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);   // True anomaly

	*omega = fmod(u-nu,pi2);
	if(*omega < 0)
		*omega = *omega + pi2;
}


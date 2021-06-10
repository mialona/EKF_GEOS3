/** @file AzElPa.c
 *  @brief Azimuth, elevation and partials from local tangent coordinates.
 *
 *  This driver contains the code for the calculation of
 *  azimuth, elevation and partials from local tangent coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "const.h"
#include "m_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void AzElPa(double *s, double *Az, double *El, double **dAds, double **dEds) {
	double rho = sqrt(s[0]*s[0]+s[1]*s[1]);
	
	// Angles
	*Az = atan2(s[0],s[1]);
	
	if(*Az < 0.0) {
		*Az = *Az+pi2;
	}
	
	*El = atan(s[2]/rho);
	
	// Partials
	*dAds = v_create(3);
	(*dAds)[0] = s[1]/(rho*rho);
	(*dAds)[1] = -s[0]/(rho*rho);
	(*dAds)[2] = 0.0;
	
	*dEds = v_create(3);
	(*dEds)[0] = (-s[0]*s[2]/rho)/v_dot(s,3,s,3);
	(*dEds)[1] = (-s[1]*s[2]/rho)/v_dot(s,3,s,3);
	(*dEds)[2] = (rho)/v_dot(s,3,s,3);
}

/** @file hgibbs.c
 *  @brief Herrick-gibbs approximation for orbit determination.
 *
 *  This driver contains the code for the 
 *  Herrick-gibbs approximation for orbit determination.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/const.h"
#include "../includes/m_utils.h"
#include "../includes/unit.h"
#include "../includes/angl.h"

#include <stdio.h>
#include <math.h>


void hgibbs(double *r1, double *r2, double *r3, double Mjd1, double Mjd2, double Mjd3, double **v2, double *theta, double *theta1, double *copa, char **error) {
	*error = "          ok";
	*theta= 0.0;
	*theta1= 0.0;
	
	double magr1 = v_norm(r1,3);
	double magr2 = v_norm(r2,3);
	double magr3 = v_norm(r3,3);
	*v2 = v_create(3);
	
	double tolangle = 0.01745329251994;
	double dt21 = (Mjd2-Mjd1)*86400.0;
	double dt31 = (Mjd3-Mjd1)*86400.0;
	double dt32 = (Mjd3-Mjd2)*86400.0;
	
	double *p = v_cross(r2,3,r3,3);
	double *pn = unit(p);
	double *r1n = unit(r1);
	*copa = asin(v_dot(pn,3,r1n,3));

	if(fabs(v_dot(r1n,3,pn,3)) > 0.017452406) {
		*error = "not coplanar";
	}
	
	*theta  = angl(r1,r2);
	*theta1 = angl(r2,r3);
	
	if ((*theta > tolangle) | (*theta1 > tolangle)) {
		*error = "   angl > 1ø";
	}
	
	double gm = GM_Earth;
	double term1 = -dt32*(1.0/(dt21*dt31) + gm/(12.0*magr1*magr1*magr1));
	double term2 = (dt32-dt21)*(1.0/(dt21*dt32) + gm/(12.0*magr2*magr2*magr2));
	double term3 = dt21*(1.0/(dt32*dt31) + gm/(12.0*magr3*magr3*magr3));

	*v2 =  v_sum(v_mul_scalar(r1,3,term1),3,v_sum(v_mul_scalar(r2,3,term2),3,v_mul_scalar(r3,3,term3),3),3);
}


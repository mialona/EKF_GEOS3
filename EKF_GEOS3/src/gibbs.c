/** @file gibbs.c
 *  @brief Gibbs method of orbit determination.
 *
 *  This driver contains the code for the 
 *  Gibbs method of orbit determination.
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


void gibbs(double *r1, double *r2, double *r3, double **v2, double *theta, double *theta1, double *copa, char **error) {
	double small= 0.00000001;
	*theta= 0.0;
	*error = "          ok";
	*theta1= 0.0;

	double magr1 = v_norm(r1,3);
	double magr2 = v_norm(r2,3);
	double magr3 = v_norm(r3,3);
	*v2 = v_create(3);

	double *p = v_cross(r2,3,r3,3);
	double *q = v_cross(r3,3,r1,3);
	double *w = v_cross(r1,3,r2,3);
	double *pn = unit(p);
	double *r1n = unit(r1);
	*copa = asin(v_dot(pn,3,r1n,3));

	if(fabs(v_dot(r1n,3,pn,3)) > 0.017452406) {
		*error = "not coplanar";
	}

	double *d = v_sum(p,3,v_sum(q,3,w,3),3);
	double magd = v_norm(d,3);
	double *n = v_sum(v_mul_scalar(p,3,magr1),3,v_sum(v_mul_scalar(q,3,magr2),3,v_mul_scalar(w,3,magr3),3),3);
	double magn = v_norm(n,3);
	double *nn = unit(n);
	double *dn = unit(d);

	// -------------------------------------------------------------
	// determine if  the orbit is possible. both d and n must be in
	// the same direction, and non-zero.
	// -------------------------------------------------------------
	if ((fabs(magd)<small) || (fabs(magn)<small) || (v_dot(nn,3,dn,3) < small)) {
		*error = "  impossible";
	}
	else {
		*theta  = angl(r1,r2);
		*theta1 = angl(r2,r3);

		// ----------- perform gibbs method to find v2 -----------
		double r1mr2 = magr1-magr2;
		double r3mr1 = magr3-magr1;
		double r2mr3 = magr2-magr3;
		double *s = v_sum(v_mul_scalar(r3,3,r1mr2),3,v_sum(v_mul_scalar(r2,3,r3mr1),3,v_mul_scalar(r1,3,r2mr3),3),3);
		double *b  = v_cross(d,3,r2,3);
		double gm = GM_Earth;
		double l  = sqrt(gm/(magd*magn));
		double tover2 = l/magr2;
		*v2 =  v_sum(v_mul_scalar(b,3,tover2),3,v_mul_scalar(s,3,l),3);
	}
}


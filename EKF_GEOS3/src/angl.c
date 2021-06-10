/** @file angl.c
 *  @brief Angle between the two vectors.
 *
 *  This driver contains the code for the 
 *  angle between the two vectors (-pi to pi).
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"

#include <stdio.h>
#include <math.h>


double angl(double *vec1, double *vec2) {
	double small     = 0.00000001;
	double undefined = 999999.1;

	double magv1 = v_norm(vec1,3);
	double magv2 = v_norm(vec2,3);

	if(magv1*magv2 > small^2) {
		double temp = v_dot(vec1,3,vec2,3)/(magv1*magv2);
		if(fabs(temp) > 1.0) {
			temp = temp/fabs(temp);
		}
		return acos(temp);
	}
	else {
		return undefined;
	}
}


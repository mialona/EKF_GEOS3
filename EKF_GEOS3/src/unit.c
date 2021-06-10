/** @file unit.c
 *  @brief Unit vector.
 *
 *  This driver contains the code for the calculation of
 *  a unit vector given the original vector.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"

#include <stdio.h>
#include <math.h>


double *unit(double *vec) {
	double small = 0.000001,
		   *outvec = v_create(3),
		   magv = v_norm(vec,3);

	if(magv > small) {
		for(int i=0; i<3; i++) {
			outvec[i]= vec[i]/magv;
		}
	}
	else {
		for(int i=0; i<3; i++) {
			outvec[i]= 0.0;
		}
	}
	
	return outvec;
}

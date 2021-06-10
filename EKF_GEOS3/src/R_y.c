/** @file R_y.c
 *  @brief Rotation about y axis utilities code driver.
 *
 *  This driver contains the code for the rotation about y axis.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double **R_y(double angle) {
	double C, S, **rotmat;
	
	C = cos(angle);
	S = sin(angle);
	rotmat = m_create(3,3);
	
	rotmat[0][0] = C;  rotmat[0][1] =    0.0;  rotmat[0][2] = -1.0*S;
	rotmat[1][0] = 0.0;  rotmat[1][1] =      1.0;  rotmat[1][2] =   0.0;
	rotmat[2][0] = S;  rotmat[2][1] = 0.0;  rotmat[2][2] =   C;
	
	return rotmat;
}

/** @file PoleMatrix.c
 *  @brief Pseudo Earth-fixed to Earth-fixed coordinates for a given date.
 *
 *  This driver contains the code for the transformation from pseudo
 *  Earth-fixed to Earth-fixed coordinates for a given date.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "m_utils.h"
#include "R_x.h"
#include "R_y.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double **PoleMatrix(double xp, double yp) {
	return m_dot(R_y(-xp),3,3,R_x(-yp),3,3);
}
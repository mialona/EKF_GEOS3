/** @file GHAMatrix.c
 *  @brief Transformation from true equator and equinox to Earth 
 *  equator and Greenwich meridian system.
 *
 *  This driver contains the code for the 
 *  transformation from true equator and equinox to Earth 
 *  equator and Greenwich meridian system.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/R_z.h"
#include "../includes/gast.h"

#include <stdio.h>


double **GHAMatrix(double Mjd_UT1) {
	return R_z(gast(Mjd_UT1));
}


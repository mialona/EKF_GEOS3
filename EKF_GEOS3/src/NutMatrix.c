/** @file NutMatrix.c
 *  @brief Transformation from mean to true equator and equinox code driver.
 *
 *  This driver contains the code for the transformation
 *  from mean to true equator and equinox.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/MeanObliquity.h"
#include "../includes/NutAngles.h"
#include "../includes/m_utils.h"
#include "../includes/R_x.h"
#include "../includes/R_z.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>


double **NutMatrix(double Mjd_TT) {
	// Mean obliquity of the ecliptic
	double eps = MeanObliquity(Mjd_TT);

	// Nutation in longitude and obliquity
	double dpsi, deps;
	NutAngles(Mjd_TT, &dpsi, &deps);

	// Transformation from mean to true equator and equinox
	 return m_dot(m_dot(R_x(-eps-deps),3,3,R_z(-dpsi),3,3),3,3,R_x(eps),3,3);
}

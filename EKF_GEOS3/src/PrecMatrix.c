/** @file PrecMatrix.c
 *  @brief Precession transformation of equatorial coordinates.
 *
 *  This driver contains the code for the precession
 *  transformation of equatorial coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/const.h"
#include "../includes/m_utils.h"
#include "../includes/R_y.h"
#include "../includes/R_z.h"

#include <stdio.h>
#include <stdlib.h>


double **PrecMatrix(double Mjd_1, double Mjd_2) {
	double T  = (Mjd_1-(MJD_J2000))/36525.0;
	double dT = (Mjd_2-Mjd_1)/36525.0;

	// Precession angles
	double zeta = ((2306.2181+(1.39656-0.000139*T)*T)+
				  ((0.30188-0.000344*T)+0.017998*dT)*dT)*dT/(Arcs);
	double z = zeta + ((0.79280+0.000411*T)+0.000205*dT)*dT*dT/(Arcs);
	double theta = ((2004.3109-(0.85330+0.000217*T)*T)-
				  ((0.42665+0.000217*T)+0.041833*dT)*dT)*dT/(Arcs);

	// Precession matrix
	return m_dot(m_dot(R_z(-z),3,3,R_y(theta),3,3),3,3,R_z(-zeta),3,3);
}


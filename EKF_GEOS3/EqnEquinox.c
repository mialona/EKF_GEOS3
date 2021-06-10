/** @file EqnEquinox.c
 *  @brief Equation of the equinoxes.
 *
 *  This driver contains the code for the 
 *  computation of the equation of the equinoxes.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "NutAngles.h"
#include "MeanObliquity.h"

#include <stdio.h>
#include <math.h>


double EqnEquinox(double Mjd_TT) {
	double dpsi, deps;
	
	// Nutation in longitude and obliquity
	NutAngles(Mjd_TT, &dpsi, &deps);
	
	// Equation of the equinoxes
	return dpsi * cos(MeanObliquity(Mjd_TT));
}
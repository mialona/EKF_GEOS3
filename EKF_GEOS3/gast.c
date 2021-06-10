/** @file gast.c
 *  @brief Greenwich Apparent Sidereal Time.
 *
 *  This driver contains the code for the 
 *  Greenwich Apparent Sidereal Time.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "gmst.h"
#include "EqnEquinox.h"

#include <stdio.h>
#include <math.h>


double gast(double Mjd_UT1) {
	double gstime  = fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2.0*M_PI);
	
	if(gstime < 0)
		gstime = gstime + 2.0*M_PI;
	
	return gstime;
}


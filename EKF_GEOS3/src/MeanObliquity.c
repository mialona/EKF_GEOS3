/** @file MeanObliquity.c
 *  @brief Mean obliquity of the ecliptic code driver.
 *
 *  This driver contains the code for the computation of
 *  the mean obliquity of the ecliptic.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"
#include "../includes/const.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double MeanObliquity(double Mjd_TT) {
	double T, MOblq;
	
	T = (Mjd_TT-MJD_J2000)/36525.0;

	MOblq = Rad *( 84381.448/3600.0-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0);
	
	return MOblq;
}

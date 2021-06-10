/** @file Mjday.c
 *  @brief Modified julian date code driver.
 *
 *  This driver contains the code for the calculation of the
 *  Modified julian date.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug Removed pre-initiated hours, min, sec.
 */

#include "../includes/m_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Mjday(int yr, int mon, int day, int hr, int min, int sec) {
	double jd;

	jd = 367.0 * yr
		- floor( (7.0 * (yr + floor( (mon + 9.0) / 12.0) ) ) * 0.25 )
		+ floor( 275.0 * mon / 9.0 )
		+ day + 1721013.5
		+ ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;
	
	return jd-2400000.5;
}

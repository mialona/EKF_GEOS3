/** @file timediff.h
 *  @brief Function prototypes for the time differences.
 *
 *  This header file contains the prototypes for the 
 *  time differences [s].
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include <stdio.h>


/** @brief Time differences [s].
 *
 *  @param [in] UT1_UTC Time.
 *  @param [in] TAI_UTC Time.
 *  @param [out] UT1-TAI Time difference [s].
 *  @param [out] UTC_GPS Time difference [s].
 *  @param [out] UT1-GPS Time difference [s].
 *  @param [out] TT-UTC Time difference [s].
 *  @param [out] GPS-UTC Time difference [s].
 */
void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI, double *UTC_GPS,
			  double *UT1_GPS, double *TT_UTC, double *GPS_UTC);


#endif
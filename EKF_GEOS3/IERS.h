/** @file IERS.h
 *  @brief Function prototypes for the 
 *  IERS time and polar motion data.
 *
 *  This header file contains the prototypes for the 
 *  management of IERS time and polar motion data.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _IERS_
#define _IERS_

#include "const.h"
#include "global.h"
#include "m_utils.h"

#include <stdio.h>
#include <math.h>


/** @brief IERS time and polar motion data.
 *
 *  @param [in] Mjd_UTC Modified Julian Date UTC.
 *  @param [in] interp
 *  @param [out] x_pole Pole coordinate [rad]
 *  @param [out] y_pole Pole coordinate [rad]
 *  @param [out] UT1_UTC UT1-UTC time difference [s]
 *  @param [out] LOD Length of day [s]
 *  @param [out] dpsi
 *  @param [out] deps
 *  @param [out] dx_pole Pole coordinate [rad]
 *  @param [out] dy_pole Pole coordinate [rad]
 *  @param [out] TAI_UTC TAI-UTC time difference [s]
 */
void IERS(double Mjd_UTC, char interp,
			double *x_pole, double *y_pole, double *UT1_UTC,
			double *LOD, double *dpsi, double *deps, double *dx_pole,
			double *dy_pole, double *TAI_UTC);


#endif
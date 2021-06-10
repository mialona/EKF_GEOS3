/** @file JPL_Eph_DE430.h
 *  @brief Function prototypes for the computation of the
 *  sun, moon, and nine major planets' equatorial position.
 *
 *  This header file contains the prototypes for the computation of the
 *  sun, moon, and nine major planets' equatorial position using JPL Ephemerides.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _JPL_
#define _JPL_

#include "m_utils.h"
#include "Cheb3D.h"

#include <stdio.h>
#include <math.h>


/** @brief Sun, moon, and nine major planets' equatorial position using JPL Ephemerides.
 *
 *  @param [in] Mjd_TDB Modified julian date of TDB.
 *  @param [out] r_Mercury Vector.
 *  @param [out] r_Venus Vector.
 *  @param [out] r_Earth Vector.
 *  @param [out] r_Mars Vector.
 *  @param [out] r_Jupiter Vector.
 *  @param [out] r_Saturn Vector.
 *  @param [out] r_Uranus Vector.
 *  @param [out] r_Neptune Vector.
 *  @param [out] r_Pluto Vector.
 *  @param [out] r_Moon Vector.
 *  @param [out] r_Sun Vector.
 */
void JPL_Eph_DE430(double Mjd_TDB, double **r_Mercury, double **r_Venus,
				   double **r_Earth, double **r_Mars, double **r_Jupiter,
				   double **r_Saturn, double **r_Uranus, double **r_Neptune,
				   double **r_Pluto, double **r_Moon, double **r_Sun);


#endif
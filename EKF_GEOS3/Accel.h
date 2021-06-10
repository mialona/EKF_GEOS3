/** @file Accel.h
 *  @brief Function prototypes for the computation of the
 *  acceleration of an Earth orbiting satellite.
 *
 *  This header file contains the prototypes for the computation of the
 *  acceleration of an Earth orbiting satellite due to 
 *  the Earth's harmonic gravity field, 
 *  the gravitational perturbations of the Sun and Moon
 *  the solar radiation pressure and the atmospheric drag.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _ACCEL_
#define _ACCEL_

#include "global.h"
#include "const.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "m_utils.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Mjday_TDB.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"

#include <stdio.h>
#include <math.h>


/** @brief Acceleration of an Earth orbiting satellite.
 *
 *  @param [in] x Terrestrial Time (Modified Julian Date).
 *  @param [in] Y Satellite state vector in the ICRF/EME2000 system.
 *  @param [out] dY Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system.
 */
void Accel(double x, double *Y, double **dY);


#endif

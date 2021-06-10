/** @file AzElPa.h
 *  @brief Function prototypes for the calculation of
 *  azimuth, elevation and partials from local tangent coordinates.
 *
 *  This header file contains the prototypes for the calculation of
 *  azimuth, elevation and partials from local tangent coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _AZELPA_
#define _AZELPA_

#include "const.h"
#include "m_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/** @brief Calculation of azimuth, elevation and partials from local
 *  tangent coordinates.
 *
 *  @param [in] s Topocentric local tangent coordinates (East-North-Zenith frame).
 *  @param [out] A Azimuth [rad].
 *  @param [out] E Elevation [rad].
 *  @param [out] dAds Partials of azimuth w.r.t. s.
 *  @param [out] dEds Partials of elevation w.r.t. s.
 */
void AzElPa(double *s, double *Az, double *El, double **dAds, double **dEds);


#endif
/** @file Geodetic.h
 *  @brief Function prototypes for the 
 *  Geodetic coordinates.
 *
 *  This header file contains the prototypes for the 
 *  Geodetic coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _GEODETIC_
#define _GEODETIC_


/** @brief Geodetic coordinates.
 *
 *  @param [in] r Position vector [m].
 *  @param [out] lon Longitude [rad].
 *  @param [out] lat Latitude [rad].
 *  @param [out] h Altitude [m].
 */
void Geodetic(double *r, double *lon, double *lat, double *h);


#endif

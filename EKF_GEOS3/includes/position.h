/** @file position.h
 *  @brief Function prototypes for the
 *  position vector from geodetic coordinates.
 *
 *  This header file contains the prototypes for the position vector
 *  from geodetic coordinates (longitude, latitude and altitude).
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _POSITION_
#define _POSITION_


/** @brief Position vector (r[m]) from geodetic coordinates..
 *
 *  @param [in] lon Longitude [rad].
 *  @param [in] lat Latitude [rad].
 *  @param [in] h Altitude [m].
 *  @return Vector result.
 */
double *position(double lon, double lat, double h);


#endif

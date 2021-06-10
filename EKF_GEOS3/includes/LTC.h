/** @file LTC.h
 *  @brief Transformation from Greenwich meridian system to 
 *  local tangent coordinates.
 *
 *  This header file contains the prototypes for the 
 *  transformation from Greenwich meridian system to 
 *  local tangent coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _LTC_
#define _LTC_


/** @brief Transformation from Greenwich meridian system to 
 *  local tangent coordinates.
 *
 *  @param [in] lon Geodetic East longitude [rad].
 *  @param [in] lat Geodetic latitude [rad].
 *  @return Rotation matrix from the Earth equator and Greenwich meridian
 *  to the local tangent (East-North-Zenith) coordinate system.
 */
double **LTC(double lon, double lat);


#endif

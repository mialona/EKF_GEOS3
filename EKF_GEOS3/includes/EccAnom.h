/** @file EccAnom.h
 *  @brief Function prototypes for the 
 *  computation of the eccentric anomaly for elliptic orbits.
 *
 *  This header file contains the prototypes for the 
 *  computation of the eccentric anomaly for elliptic orbits.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _ECCANOM_
#define _ECCANOM_


/** @brief Eccentric anomaly for elliptic orbits.
 *
 *  @param [in] M Mean anomaly in [rad].
 *  @param [in] e Eccentricity of the orbit [0,1].
 *  @return Eccentric anomaly in [rad].
 */
double EccAnom(double M, double e);


#endif

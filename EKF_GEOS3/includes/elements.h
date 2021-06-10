/** @file elements.h
 *  @brief Function prototypes for the computation of the 
 *  Osculating Keplerian elements.
 *
 *  This header file contains the prototypes for the computation of the 
 *  Osculating Keplerian elements from the satellite state
 *  vector for elliptic orbits.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _ELEMENTS_
#define _ELEMENTS_


/** @brief Osculating Keplerian elements from the satellite state
 *  vector for elliptic orbits.
 *  
 *  The function cannot be used with state vectors describing a circular
 *  or non-inclined orbit.
 *
 *  @param [in] y State vector (x,y,z,vx,vy,vz).
 *  @param [out] p Semilatus rectum [m].
 *  @param [out] a Semimajor axis.
 *  @param [out] e Eccentricity.
 *  @param [out] i Inclination [rad].
 *  @param [out] Omega Longitude of the ascending node [rad].
 *  @param [out] omega Argument of pericenter [rad].
 *  @param [out] M Mean anomaly [rad].
 */
void elements(double *y, double *p, double *a, double *e, double *i,
			  double *Omega, double *omega, double *M);


#endif

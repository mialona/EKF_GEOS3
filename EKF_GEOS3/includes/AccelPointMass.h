/** @file AccelPointMass.h
 *  @brief Function prototypes for the computation of the
 *  perturbational acceleration due to a point mass.
 *
 *  This header file contains the prototypes for the 
 *  perturbational acceleration due to a point mass.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _ACCELPOINTMASS_
#define _ACCELPOINTMASS_

#include "m_utils.h"

#include <stdio.h>


/** @brief Perturbational acceleration due to a point mass.
 *
 *  @param [in] r Satellite position vector.
 *  @param [in] s Point mass position vector.
 *  @param [in] GM Gravitational coefficient of point mass.
 *  @return Acceleration (a=d^2r/dt^2).
 */
double *AccelPointMass(double *r, double *s, double GM);


#endif
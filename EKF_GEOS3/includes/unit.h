/** @file unit.h
 *  @brief Function prototypes for the calculation of
 *  a unit vector given the original vector.
 *
 *  This header file contains the prototypes for the calculation of
 *  a unit vector given the original vector.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _UNIT_
#define _UNIT_


/** @brief Calculation of a unit vector given the original vector.
 *  If a zero vector is input, the vector is set to zero.
 *
 *  @param [in] vec Vector.
 *  @return Unit vector.
 */
double *unit(double *vec);


#endif

/** @file Cheb3D.h
 *  @brief Function prototypes for the calculation of the
 *  Chebyshev approximation of 3-dimensional vectors.
 *
 *  This header file contains the prototypes for the calculation of the
 *  Chebyshev approximation of 3-dimensional vectors.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _CHEB3D_
#define _CHEB3D_

#include "m_utils.h"

#include <stdio.h>
#include <stdlib.h>


/** @brief Chebyshev approximation of 3-dimensional vectors.
 *
 *  @param [in] t Initial value.
 *  @param [in] N Number of coefficients.
 *  @param [in] Ta Begin interval.
 *  @param [in] Tb End interval.
 *  @param [in] Cx Coefficients of Chebyshev polyomial (x-coordinate).
 *  @param [in] Cy Coefficients of Chebyshev polyomial (y-coordinate).
 *  @param [in] Cz Coefficients of Chebyshev polyomial (z-coordinate).
 *  @return Chebyshev approximation.
 */
double *Cheb3D(double t, int N, double Ta, double Tb, double *Cx, double *Cy, double *Cz);


#endif
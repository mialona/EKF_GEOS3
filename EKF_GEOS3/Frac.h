/** @file Frac.h
 *  @brief Function prototypes for the fractional
 *  part of a number (y=x-[x]).
 *
 *  This header file contains the prototypes for the fractional
 *  part of a number (y=x-[x]).
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _FRAC_
#define _FRAC_

#include <stdio.h>
#include <math.h>


/** @brief Fractional part of a number.
 *
 *  @param [in] x Number.
 *  @return y=x-[x].
 */
double Frac(double x);


#endif
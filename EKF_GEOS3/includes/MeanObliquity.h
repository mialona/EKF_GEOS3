/** @file MeanObliquity.h
 *  @brief Function prototypes for the mean obliquity of the ecliptic.
 *
 *  This header file contains the prototypes for the computation of
 *  the mean obliquity of the ecliptic.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _MEANOBLIQUITY_
#define _MEANOBLIQUITY_

#include "const.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/** @brief Computation of the mean obliquity of the ecliptic.
 *
 *  @param [in] Mjd_TT Modified Julian Date (Terrestrial Time).
 *  @return Mean obliquity of the ecliptic [rad].
 */
double MeanObliquity(double Mjd_TT);


#endif
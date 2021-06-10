/** @file NutAngles.h
 *  @brief Function prototypes for the calculation
 *  of the nutation in longitude and obliquity.
 *
 *  This header file contains the prototypes for the calculation
 *  of the nutation in longitude and obliquity.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _NUTANGLES_
#define _NUTANGLES_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/** @brief Calculation of the nutation in longitude and obliquity.
 *
 *  @param [in] year Year.
 *  @param [out] dpsi Nutation Angle.
 *  @param [out] deps Nutation Angle.
 */
void NutAngles(double Mjd_TT, double *dpsi, double *deps);


#endif
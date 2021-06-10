/** @file gast.h
 *  @brief Function prototypes for the 
 *  Greenwich Apparent Sidereal Time.
 *
 *  This header file contains the prototypes for the 
 *  Greenwich Apparent Sidereal Time.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _GAST_
#define _GAST_

#include "gmst.h"
#include "EqnEquinox.h"

#include <stdio.h>
#include <math.h>


/** @brief Greenwich Apparent Sidereal Time.
 *
 *  @param [in] Mjd_UT1 Modified Julian Date UT1.
 *  @return GAST in [rad].
 */
double gast(double Mjd_UT1);


#endif
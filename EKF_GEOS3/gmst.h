/** @file gmst.h
 *  @brief Function prototypes for the Greenwich
 *  Mean Sidereal Time.
 *
 *  This header file contains the prototypes for the Greenwich
 *  Mean Sidereal Time.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _GMST_
#define _GMST_

#include "const.h"
#include "m_utils.h"
#include "R_y.h"
#include "R_z.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/** @brief Greenwich Mean Sidereal Time.
 *
 *  @param [in] Mjd_UT1 Modified Julian Date UT1.
 *  @return GMST in [rad].
 */
double gmst(double Mjd_UT1);


#endif
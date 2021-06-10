/** @file R_x.h
 *  @brief Function prototypes for the rotation about x axis.
 *
 *  This header file contains the prototypes for the 
 *  rotation about x axis code driver.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _RX_
#define _RX_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/** @brief Vector rotation about x axis.
 *
 *  @param [in] c angle of rotation [rad].
 *  @return Vector result.
 */
double **R_x(double angle);


#endif
/** @file R_z.h
 *  @brief Function prototypes for the rotation about z axis.
 *
 *  This header file contains the prototypes for the 
 *  rotation about z axis code driver.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _RZ_
#define _RZ_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/** @brief Vector rotation about z axis.
 *
 *  @param [in] c angle of rotation [rad].
 *  @return Vector result.
 */
double **R_z(double angle);


#endif
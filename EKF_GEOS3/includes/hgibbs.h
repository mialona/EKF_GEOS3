/** @file hgibbs.h
 *  @brief Function prototypes for the 
 *  Herrick-gibbs approximation for orbit determination.
 *
 *  This header file contains the prototypes for the 
 *  Herrick-gibbs approximation for orbit determination.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _HGIBBS_
#define _HGIBBS_

#include "const.h"
#include "m_utils.h"
#include "unit.h"
#include "angl.h"

#include <stdio.h>
#include <math.h>


/** @brief Herrick-gibbs approximation for orbit determination.
 *  
 *  Finds the middle velocity vector for the 3 given
 *  position vectors.
 *
 *  @param [in] r1 ijk position vector 1.
 *  @param [in] r2 ijk position vector 2.
 *  @param [in] r3 ijk position vector 3.
 *  @param [in] Mjd1 Julian date of 1st sighting.
 *  @param [in] Mjd2 Julian date of 2st sighting.
 *  @param [in] Mjd3 Julian date of 3st sighting.
 *  @param [out] v2 ijk velocity vector for r2.
 *  @param [out] theta Angl between vectors.
 *  @param [out] theta1.
 *  @param [out] copa.
 *  @param [out] error Flag indicating success.
 */
void hgibbs(double *r1, double *r2, double *r3, double Mjd1, double Mjd2, double Mjd3, double **v2, double *theta, double *theta1, double *copa, char **error);


#endif
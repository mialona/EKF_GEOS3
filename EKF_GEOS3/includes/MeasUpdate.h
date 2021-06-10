/** @file MeasUpdate.h
 *  @brief Function prototypes for the 
 *  meas updater.
 *
 *  This header file contains the prototypes for the 
 *  meas updater.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "m_utils.h"

#include <stdio.h>


/** @brief Meas updater.
 *
 *  @param [in] z Scalar.
 *  @param [in] g Scalar.
 *  @param [in] s Scalar.
 *  @param [in] G Vector (6 components).
 *  @param [out] K Vector (6 components).
 *  @param [in,out] x Vector (6 components).
 *  @param [in,out] P Matrix 6x6.
 *  @return GAST in [rad].
 */
void MeasUpdate(double z, double g, double s, double *G, double **K, double **x, double ***P);


#endif
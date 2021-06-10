/** @file G_AccelHarmonic.h
 *  @brief Function prototypes for the computation of the
 *  gradient of the Earth's harmonic gravity field.
 *
 *  This header file contains the prototypes for the computation of the
 *  gradient of the Earth's harmonic gravity field.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _G_ACCELHARMONIC_
#define _G_ACCELHARMONIC_


/** @brief Gradient of the Earth's harmonic gravity field.
 *
 *  @param [in] r Satellite position vector in the true-of-date system.
 *  @param [in] U Transformation matrix to body-fixed system.
 *  @param [in] n_max Gravity model degree.
 *  @param [in] m_max Gravity model order.
 *  @return Gradient (G=da/dr) in the true-of-date system.
 */
double **G_AccelHarmonic(double *r, double **U, int n_max, int m_max);


#endif

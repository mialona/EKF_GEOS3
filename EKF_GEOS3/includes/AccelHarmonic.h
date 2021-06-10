/** @file AccelHarmonic.h
 *  @brief Function prototypes for the computation of the
 *  acceleration due to the harmonic gravity field.
 *
 *  This header file contains the prototypes for the computation of the
 *  acceleration due to the harmonic gravity field of the 
 *  central body.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_


/** @brief Acceleration due to the harmonic gravity field of the 
 *  central body.
 *
 *  @param [in] r Satellite position vector in the inertial system.
 *  @param [in] E Transformation matrix to body-fixed system.
 *  @param [in] n_max Maximum degree.
 *  @param [in] m_max Maximum order (m_max<=n_max; m_max=0 for zonals, only).
 *  @return Acceleration (a=d^2r/dt^2).
 */
double *AccelHarmonic(double *r, double **E, int n_max, int m_max);


#endif

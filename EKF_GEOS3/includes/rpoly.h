/** @file rpoly.h
 *  @brief Function prototypes for the computation of the
 *  roots of a polynomial with real coefficients.
 *
 *  This header file contains the prototypes for the computation of the roots
 *  of a polynomial with real coefficients.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _RPOLY_
#define _RPOLY_


/** @brief Acceleration of an Earth orbiting satellite.
 *
 *  @param [in] p Coefficients of the polynomial.
 *  @param [in] degree Polynomial degree.
 *  @param [out] zeror Real coefficients of the roots.
 *  @param [out] zeroi Imaginary coefficients of the roots.
 *  @return Number of roots.
 */
int real_poly_roots(double *p, int degree, double **zeror, double **zeroi);


#endif

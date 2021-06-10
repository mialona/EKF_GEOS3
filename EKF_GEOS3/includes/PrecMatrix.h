/** @file PrecMatrix.h
 *  @brief Function prototypes for the precession
 *  transformation of equatorial coordinates.
 *
 *  This header file contains the prototypes for the precession
 *  transformation of equatorial coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _PRECMATRIX_
#define _PRECMATRIX_


/** @brief Precession transformation of equatorial coordinates.
 *
 *  @param [in] Mjd_1 Epoch given (Modified Julian Date TT).
 *  @param [in] MjD_2 Epoch to precess to (Modified Julian Date TT).
 *  @return Precession transformation matrix.
 */
double **PrecMatrix(double Mjd_1, double Mjd_2);


#endif

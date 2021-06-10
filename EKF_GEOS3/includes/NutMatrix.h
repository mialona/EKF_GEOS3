/** @file NutMatrix.h
 *  @brief Function prototypes for the transformation
 *  from mean to true equator and equinox.
 *
 *  This header file contains the prototypes for the transformation
 *  from mean to true equator and equinox.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _NUTMATRIX_
#define _NUTMATRIX_


/** @brief Transformation from mean to true equator and equinox.
 *
 *  @param [in] Mjd_TT Modified Julian Date (Terrestrial Time).
 *  @return Nutation matrix.
 */
double **NutMatrix(double Mjd_TT);


#endif

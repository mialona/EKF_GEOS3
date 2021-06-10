/** @file PoleMatrix.h
 *  @brief Function prototypes for the transformation from pseudo
 *  Earth-fixed to Earth-fixed coordinates for a given date.
 *
 *  This header file contains the prototypes for the transformation from pseudo
 *  Earth-fixed to Earth-fixed coordinates for a given date.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _POLEMATRIX_
#define _POLEMATRIX_


/** @brief Pseudo Earth-fixed to Earth-fixed coordinates for a given date.
 *
 *  @param [in] xp Pole coordinte.
 *  @param [in] yp Pole coordinte.
 *  @return Pole Matrix.
 */
double **PoleMatrix(double xp, double yp);


#endif

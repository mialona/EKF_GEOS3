/** @file Legendre.h
 *  @brief Function prototypes for the 
 *  calculation of Legendre coefficients.
 *
 *  This header file contains the prototypes for the 
 *  calculation of Legendre coefficients.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _LEGENDRE_
#define _LEGENDRE_


/** @brief Legendre coefficients.
 *
 *  @param [in] n Rows.
 *  @param [in] m Columns.
 *  @param [in] fi Angle [rad].
 *  @param [out] pnm Result matrix.
 *  @param [out] dpnm Result matrix.
 */
void Legendre(int n, int m, double fi, double ***pnm, double ***dpnm);


#endif

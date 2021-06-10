/** @file GHAMatrix.h
 *  @brief Function prototypes for the 
 *  transformation from true equator and equinox to Earth 
 *  equator and Greenwich meridian system.
 *
 *  This header file contains the prototypes for the 
 *  transformation from true equator and equinox to Earth 
 *  equator and Greenwich meridian system.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _GHAMATRIX_
#define _GHAMATRIX_


/** @brief Transformation from true equator and equinox to Earth 
 *  equator and Greenwich meridian system.
 *
 *  @param [in] Mjd_UT1 Modified Julian Date UT1.
 *  @return Greenwich Hour Angle matrix.
 */
double **GHAMatrix(double Mjd_UT1);


#endif

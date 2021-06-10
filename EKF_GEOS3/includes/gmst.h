/** @file gmst.h
 *  @brief Function prototypes for the Greenwich
 *  Mean Sidereal Time.
 *
 *  This header file contains the prototypes for the Greenwich
 *  Mean Sidereal Time.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _GMST_
#define _GMST_


/** @brief Greenwich Mean Sidereal Time.
 *
 *  @param [in] Mjd_UT1 Modified Julian Date UT1.
 *  @return GMST in [rad].
 */
double gmst(double Mjd_UT1);


#endif

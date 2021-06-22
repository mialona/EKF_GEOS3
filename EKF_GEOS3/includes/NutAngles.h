/** @file NutAngles.h
 *  @brief Function prototypes for the calculation
 *  of the nutation in longitude and obliquity.
 *
 *  This header file contains the prototypes for the calculation
 *  of the nutation in longitude and obliquity.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _NUTANGLES_
#define _NUTANGLES_


/** @brief Calculation of the nutation in longitude and obliquity.
 *
 *  @param [in] Mjd_TT Modified Julian Date (Terrestrial Time).
 *  @param [out] dpsi Nutation Angle.
 *  @param [out] deps Nutation Angle.
 */
void NutAngles(double Mjd_TT, double *dpsi, double *deps);


#endif

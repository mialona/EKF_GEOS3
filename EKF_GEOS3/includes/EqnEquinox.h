/** @file EqnEquinox.h
 *  @brief Function prototypes for the 
 *  computation of the equation of the equinoxes.
 *
 *  This header file contains the prototypes for the 
 *  computation of the equation of the equinoxes.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _EQNEQUINOX_
#define _EQNEQUINOX_


/** @brief Equation of the equinoxes.
 *  
 *  The equation of the equinoxes dpsi*cos(eps) is the right ascension of
 *  the mean equinox referred to the true equator and equinox and is equal
 *  to the difference between apparent and mean sidereal time.
 *
 *  @param [in] Mjd_TT Modified Julian Date (Terrestrial Time).
 *  @return Equation of the equinoxes.
 */
double EqnEquinox(double Mjd_TT);


#endif

/** @file const.h
 *  @brief Definition of astronomical and mathematical constants.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _CONST_
#define _CONST_

#include <math.h>

// Mathematical constants
         
#define pi2         2*M_PI                // 2pi
#define Rad         M_PI/180              // Radians per degree
#define Deg         180/M_PI              // Degrees per radian
#define Arcs        3600*180/M_PI         // Arcseconds per radian

// General
#define MJD_J2000   51544.5             // Modified Julian Date of J2000
#define T_B1950     -0.500002108        // Epoch B1950
#define c_light     299792458.000000000 // Speed of light  [m/s]; DE430
#define AU          149597870700.000000 // Astronomical unit [m]; DE430

// Physical parameters of the Earth, Sun and Moon

// Equatorial radius and flattening
#define R_Earth     6378.1363e3      // Earth's radius [m]; DE430
#define f_Earth     1/298.257223563  // Flattening; WGS-84
#define R_Sun       696000e3         // Sun's radius [m]; DE430
#define R_Moon      1738e3           // Moon's radius [m]; DE430

// Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
#define omega_Earth   15.04106717866910/3600*Rad   // [rad/s]; WGS-84

// Gravitational coefficients
#define GM_Earth      398600.435436e9                  // [m^3/s^2]; DE430
#define GM_Sun        132712440041.939400e9            // [m^3/s^2]; DE430
#define GM_Moon       GM_Earth/81.30056907419062       // [m^3/s^2]; DE430
#define GM_Mercury    22031.780000e9                   // [m^3/s^2]; DE430
#define GM_Venus      324858.592000e9                 // [m^3/s^2]; DE430
#define GM_Mars       42828.375214e9                   // [m^3/s^2]; DE430
#define GM_Jupiter    126712764.800000e9               // [m^3/s^2]; DE430
#define GM_Saturn     37940585.200000e9                // [m^3/s^2]; DE430
#define GM_Uranus     5794548.600000e9                 // [m^3/s^2]; DE430
#define GM_Neptune    6836527.100580e9                 // [m^3/s^2]; DE430
#define GM_Pluto      977.0000000000009e9              // [m^3/s^2]; DE430

// Solar radiation pressure at 1 AU
#define P_Sol         1367/c_light // [N/m^2] (~1367 W/m^2); IERS 96

#endif

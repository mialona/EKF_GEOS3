/** @file position.c
 *  @brief Position vector from geodetic coordinates code driver.
 *
 *  This driver contains the code for the position vector
 *  from geodetic coordinates (longitude, latitude and altitude).
 *  latitude [rad], altitude [m]).
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/const.h"
#include "../includes/m_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *position(double lon, double lat, double h) {
    double R_equ, f, e2, CosLat, SinLat, N, *r;
    
    r = v_create(3);

    R_equ = R_Earth;
    f     = f_Earth;

    e2     = f*(2.0 - f);    // Square of eccentricity
    CosLat = cos(lat);     // (Co)sine of geodetic latitude
    SinLat = sin(lat);

    // Position vector
    N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

    r[0] =  (         N+h)*CosLat*cos(lon);
    r[1] =  (         N+h)*CosLat*sin(lon);
    r[2] =  ((1.0-e2)*N+h)*SinLat;
    
    return r;
}


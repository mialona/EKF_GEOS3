/** @file Accel.c
 *  @brief Acceleration of an Earth orbiting satellite.
 *
 *  This driver contains the code for the computation of thecomputation of the
 *  acceleration of an Earth orbiting satellite due to 
 *  the Earth's harmonic gravity field, 
 *  the gravitational perturbations of the Sun and Moon
 *  the solar radiation pressure and the atmospheric drag.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/global.h"
#include "../includes/const.h"
#include "../includes/IERS.h"
#include "../includes/timediff.h"
#include "../includes/PrecMatrix.h"
#include "../includes/NutMatrix.h"
#include "../includes/m_utils.h"
#include "../includes/PoleMatrix.h"
#include "../includes/GHAMatrix.h"
#include "../includes/Mjday_TDB.h"
#include "../includes/JPL_Eph_DE430.h"
#include "../includes/AccelHarmonic.h"
#include "../includes/AccelPointMass.h"

#include <stdio.h>
#include <math.h>


void Accel(double x, double *Y, double **dY) {
	extern Param AuxParam;

	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	IERS(AuxParam.Mjd_UTC + x/86400.0,'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);
	
	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
	timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);
	
	double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400.0 + UT1_UTC/86400.0;
	double Mjd_TT = AuxParam.Mjd_UTC + x/86400.0 + TT_UTC/86400.0;
	
	double **P = PrecMatrix((MJD_J2000),Mjd_TT);
	double **N = NutMatrix(Mjd_TT);
	double **T = m_dot(N,3,3,P,3,3);
	double **E = m_dot(m_dot(PoleMatrix(x_pole,y_pole),3,3,GHAMatrix(Mjd_UT1),3,3),3,3,T,3,3);

	double Mjday = Mjday_TDB(Mjd_TT);
	double *r_Mercury, *r_Venus, *r_Earth, *r_Mars, *r_Jupiter, *r_Saturn, *r_Uranus, *r_Neptune, *r_Pluto, *r_Moon, *r_Sun;
	JPL_Eph_DE430(Mjday,&r_Mercury,&r_Venus,&r_Earth,&r_Mars,&r_Jupiter,&r_Saturn,&r_Uranus,&r_Neptune,&r_Pluto,&r_Moon,&r_Sun);
	
	// Acceleration due to harmonic gravity field
	double *Y_aux = v_create(3);
	Y_aux[0] = Y[0]; Y_aux[1] = Y[1]; Y_aux[2] = Y[2];
	double *a = AccelHarmonic(Y_aux, E, AuxParam.n, AuxParam.m);

	// Luni-solar perturbations
	if(AuxParam.sun) {
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Sun,(GM_Sun)),3);
	}

	if(AuxParam.moon) {
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Moon,(GM_Moon)),3);
	}

	// Planetary perturbations
	if(AuxParam.planets) {
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Mercury,(GM_Mercury)),3);
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Venus,(GM_Venus)),3);
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Mars,(GM_Mars)),3);
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Jupiter,(GM_Jupiter)),3);
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Saturn,(GM_Saturn)),3);
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Uranus,(GM_Uranus)),3);
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Neptune,(GM_Neptune)),3);
		a = v_sum(a,3,AccelPointMass(Y_aux,r_Pluto,(GM_Pluto)),3);
	}

	*dY = v_create(6);
	(*dY)[0] = Y[3]; (*dY)[1] = Y[4]; (*dY)[2] = Y[5]; (*dY)[3] = a[0]; (*dY)[4] = a[1]; (*dY)[5] = a[2];
}


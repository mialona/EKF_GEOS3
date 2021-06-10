/** @file EKF_GEOS3.c
 *  @brief Main driver.
 *
 *  This driver contains the code for the Initial Orbit Determination
 *  using Gauss and Extended Kalman Filter methods.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "includes/global.h"
#include "includes/const.h"
#include "includes/m_utils.h"

#include "includes/position.h"
#include "includes/Mjday.h"
#include "includes/ode.h"
#include "includes/Accel.h"
#include "includes/LTC.h"
#include "includes/gmst.h"
#include "includes/R_z.h"
#include "includes/TimeUpdate.h"
#include "includes/AzElPa.h"
#include "includes/MeasUpdate.h"
#include "includes/IERS.h"
#include "includes/timediff.h"
#include "includes/VarEqn.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main()
{
	DE430Coeff(2285,1020);
	GGM03S(181);
	eop19620101(21413);
	GEOS3(46);
	
	double sigma_range = 92.5;    // [m]
	double sigma_az = 0.0224*Rad; // [rad]
	double sigma_el = 0.0139*Rad; // [rad]

	// Kaena Point station
	double lat = Rad*21.5748;     // [rad]
	double lon = Rad*(-158.2706); // [rad]
	double alt = 300.20;          // [m]

	double *Rs = position(lon, lat, alt);

	// // // Mjd1 = obs(1,1);
	// // // Mjd2 = obs(9,1);
	// // // Mjd3 = obs(18,1);
	// // // 
	// // // [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
	// // //                   Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
	// // // 
	// // // Y0_apr = [r2;v2];

	extern int n_eqn;
	n_eqn = 6;
	
	double *Y = v_create(n_eqn);
	Y[0] = 6221397.62857869;
	Y[1] = 2867713.77965741;
	Y[2] = 3006155.9850995;
	Y[3] = 4645.0472516175;
	Y[4] = -2752.21591588182;
	Y[5] = -7507.99940986939;

	double Mjd0 = Mjday(1995,1,29,2,38,0);

	extern double **obs;
	extern int fobs;
	double Mjd_UTC = obs[8][0];

	extern Param AuxParam;
	AuxParam.Mjd_UTC = Mjd_UTC;
	AuxParam.n       = 20;
	AuxParam.m       = 20;
	AuxParam.sun     = 1;
	AuxParam.moon    = 1;
	AuxParam.planets = 1;

	int iflag = 1;
	int iwork[5];
	double t = 0.0, relerr = 1e-13, abserr = 1e-6;
	double *work = v_create(100 + 21 * n_eqn);
	ode(Accel, n_eqn, Y, &t, -(obs[8][0]-Mjd0)*86400.0, relerr, abserr, &iflag, work, iwork);

	double **P = m_zeros(6,6);
	  
	for(int i=0; i<3; i++) {
		P[i][i]=1e8;
	}
	for(int i=3; i<6; i++) {
		P[i][i]=1e3;
	}

	double **LT = LTC(lon,lat);

	double *yPhi = v_create(42);
	double **Phi = m_zeros(6,6);

	// Measurement loop
	t = 0.0;
	
	double t_old;
	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
	double Mjd_TT, Mjd_UT1, t_aux;
	double theta, **U, *r = v_create(3), *s;
	double Azim, Elev, *dAds, *dEds, *aux, *K, *dAdY, *dEdY, Dist, *dDds, *dDdY;
	for(int i=0; i<fobs; i++) {
		// Previous step
		t_old = t;

		// Time increment and propagation
		Mjd_UTC = obs[i][0];                       // Modified Julian Date
		t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
		
		IERS(Mjd_UTC,'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);
		timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);
		
		Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
		Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
		AuxParam.Mjd_UTC = Mjd_UTC;
		AuxParam.Mjd_TT = Mjd_TT;
		
		for(int ii=0; ii<6; ii++) {
			yPhi[ii] = Y[ii];
			for(int jj=0; jj<6; jj++) { 
				if(ii==jj) {
					yPhi[6*(jj+1)+ii] = 1.0;
				}
				else {
					yPhi[6*(jj+1)+ii] = 0.0;
				}
			}
		}

		n_eqn = 42;
		t_aux = 0.0;
		iflag = 1;
		work = v_create(100 + 21 * n_eqn);
		ode(VarEqn, n_eqn, yPhi, &t_aux, t-t_old, relerr, abserr, &iflag, work, iwork);

		// Extract state transition matrices
		for(int j=0; j<6; j++) {
			for (int ii = 0; ii < 6; ii++){
				Phi[ii][j] = yPhi[6*(j+1)+ii];
			}
		}

		n_eqn = 6;
		t_aux = 0;
		iflag = 1;
		work = v_create(100 + 21 * n_eqn);
		ode(Accel, n_eqn, Y, &t_aux, t-t_old, relerr, abserr, &iflag, work, iwork);


		// Topocentric coordinates
		theta = gmst(Mjd_UT1);                    // Earth rotation
		U = R_z(theta);
		for(int j=0; j<3; j++) {
			r[j] = Y[j];
		}
		s = m_dot_v(LT,3,3,v_sum(m_dot_v(U,3,3,r,3),3,v_mul_scalar(Rs,3,-1.0),3),3); // Topocentric position [m]

		// Time update
		P = TimeUpdate(P, Phi, m_zeros(6,6));

		// Azimuth and partials
		AzElPa(s, &Azim, &Elev, &dAds, &dEds);     // Azimuth, Elevation
		aux = v_dot_m(dAds,3,m_dot(LT,3,3,U,3,3),3,3);
		dAdY = v_create(6);
		for(int j=0; j<3; j++) {
			dAdY[j] = aux[j];
			dAdY[j+3] = 0.0;
		}

		// Measurement update
		MeasUpdate(obs[i][1],Azim,sigma_az,dAdY,&K,&Y,&P);

		// Elevation and partials
		for(int j=0; j<3; j++) {
			r[j] = Y[j];
		}
		s = m_dot_v(LT,3,3,v_sum(m_dot_v(U,3,3,r,3),3,v_mul_scalar(Rs,3,-1.0),3),3); // Topocentric position [m]
		AzElPa(s, &Azim, &Elev, &dAds, &dEds);     // Azimuth, Elevation
		aux = v_dot_m(dEds,3,m_dot(LT,3,3,U,3,3),3,3);
		dEdY = v_create(6);
		for(int j=0; j<3; j++) {
			dEdY[j] = aux[j];
			dEdY[j+3] = 0.0;
		}

		// Measurement update
		MeasUpdate(obs[i][2],Elev,sigma_el,dEdY,&K,&Y,&P);

		// Range and partials
		for(int j=0; j<3; j++) {
			r[j] = Y[j];
		}
		s = m_dot_v(LT,3,3,v_sum(m_dot_v(U,3,3,r,3),3,v_mul_scalar(Rs,3,-1.0),3),3); // Topocentric position [m]
		Dist = v_norm(s,3);
		dDds = v_mul_scalar(s,3,1/Dist);         // Range
		aux = v_dot_m(dDds,3,m_dot(LT,3,3,U,3,3),3,3);
		dDdY = v_create(6);
		for(int j=0; j<3; j++) {
			dDdY[j] = aux[j];
			dDdY[j+3] = 0.0;
		}
		// Measurement update
		MeasUpdate(obs[i][3],Dist,sigma_range,dDdY,&K,&Y,&P);
	}

	IERS(obs[45][0],'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);
	timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	AuxParam.Mjd_UTC = Mjd_UTC;
	AuxParam.Mjd_TT = Mjd_TT;

	t_aux = 0.0;
	iflag = 1;
	work = v_create(100 + 21 * n_eqn);
	ode(Accel, n_eqn, Y, &t_aux, -(obs[45][0]-obs[0][0])*86400.0, relerr, abserr, &iflag, work, iwork);

	double *Y_true = v_create(n_eqn);
	Y_true[0] = 5753.173e3; Y_true[1] = 2673.361e3; Y_true[2] = 3440.304e3;
	Y_true[3] = 4.324207e3; Y_true[4] = -1.924299e3; Y_true[5] = -5.728216e3;

	printf("\nError of Position Estimation\n");
	printf("dX	%10.1lf [m]\n",Y[0]-Y_true[0]);
	printf("dY	%10.1lf [m]\n",Y[1]-Y_true[1]);
	printf("dZ	%10.1lf [m]\n",Y[2]-Y_true[2]);
	printf("\nError of Velocity Estimation\n");
	printf("dVx	%10.1lf [m/s]\n",Y[3]-Y_true[3]);
	printf("dVy	%10.1lf [m/s]\n",Y[4]-Y_true[4]);
	printf("dVz	%10.1lf [m/s]\n",Y[5]-Y_true[5]);
	
    return 0;
}

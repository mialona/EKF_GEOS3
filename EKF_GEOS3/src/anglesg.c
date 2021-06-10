/** @file anglesg.c
 *  @brief Orbit determination using three optical sightings.
 *
 *  This driver contains the code for the 
 *  orbit determination using three optical sightings.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/const.h"
#include "../includes/m_utils.h"
#include "../includes/Geodetic.h"
#include "../includes/LTC.h"
#include "../includes/IERS.h"
#include "../includes/timediff.h"
#include "../includes/PrecMatrix.h"
#include "../includes/NutMatrix.h"
#include "../includes/PoleMatrix.h"
#include "../includes/GHAMatrix.h"
#include "../includes/rpoly.h"
#include "../includes/gibbs.h"
#include "../includes/hgibbs.h"
#include "../includes/elements.h"
#include "../includes/angl.h"

#include <stdio.h>
#include <math.h>
#include <string.h>


void anglesg(double az1, double az2, double az3, double el1,
		 	 double el2, double el3, double Mjd1, double Mjd2,
			 double Mjd3, double *Rs1, double *Rs2, double *Rs3,
			 double **r, double **v) {

	double *L1 = v_create(3);
	L1[0] = cos(el1)*sin(az1);
	L1[1] = cos(el1)*cos(az1);
	L1[2] = sin(el1);
	double *L2 = v_create(3);
	L2[0] = cos(el2)*sin(az2);
	L2[1] = cos(el2)*cos(az2);
	L2[2] = sin(el2);
	double *L3 = v_create(3);
	L3[0] = cos(el3)*sin(az3);
	L3[1] = cos(el3)*cos(az3);
	L3[2] = sin(el3);

	double lon1, lat1, h1;
	Geodetic(Rs1, &lon1, &lat1, &h1);
	//double lon2, lat2, h2;
	//Geodetic(Rs2, &lon2, &lat2, &h2);
	//double lon3, lat3, h3;
	//Geodetic(Rs3, &lon3, &lat3, &h3);

	double **M1 = LTC(lon1, lat1);
	//double **M2 = LTC(lon2, lat2);
	//double **M3 = LTC(lon3, lat3);

	// body-fixed system
	double *Lb1 = m_dot_v(m_trans(M1,3),3,3,L1,3);
	double *Lb2 = m_dot_v(m_trans(M1,3),3,3,L2,3);
	double *Lb3 = m_dot_v(m_trans(M1,3),3,3,L3,3);

	// mean of date system (J2000)
	double Mjd_UTC = Mjd1;

	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	IERS(Mjd_UTC,'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);

	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
	timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);

	double Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

	double **P = PrecMatrix((MJD_J2000),Mjd_TT);
	double **N = NutMatrix(Mjd_TT);
	double **T = m_dot(N,3,3,P,3,3);
	double **E = m_dot(m_dot(PoleMatrix(x_pole,y_pole),3,3,GHAMatrix(Mjd_UT1),3,3),3,3,T,3,3);

	double *Lm1 = m_dot_v(m_trans(E,3),3,3,Lb1,3);
	Rs1 = m_dot_v(m_trans(E,3),3,3,Rs1,3);

	Mjd_UTC = Mjd2;
	IERS(Mjd_UTC,'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);
	timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

	P = PrecMatrix((MJD_J2000),Mjd_TT);
	N = NutMatrix(Mjd_TT);
	T = m_dot(N,3,3,P,3,3);
	E = m_dot(m_dot(PoleMatrix(x_pole,y_pole),3,3,GHAMatrix(Mjd_UT1),3,3),3,3,T,3,3);

	double *Lm2 = m_dot_v(m_trans(E,3),3,3,Lb2,3);
	Rs2 = m_dot_v(m_trans(E,3),3,3,Rs2,3);

	Mjd_UTC = Mjd3;
	IERS(Mjd_UTC,'l',&x_pole,&y_pole,&UT1_UTC,&LOD,&dpsi,&deps,&dx_pole,&dy_pole,&TAI_UTC);
		timediff(UT1_UTC,TAI_UTC,&UT1_TAI,&UTC_GPS,&UT1_GPS,&TT_UTC,&GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

	P = PrecMatrix((MJD_J2000),Mjd_TT);
	N = NutMatrix(Mjd_TT);
	T = m_dot(N,3,3,P,3,3);
	E = m_dot(m_dot(PoleMatrix(x_pole,y_pole),3,3,GHAMatrix(Mjd_UT1),3,3),3,3,T,3,3);

	double *Lm3 = m_dot_v(m_trans(E,3),3,3,Lb3,3);
	Rs3 = m_dot_v(m_trans(E,3),3,3,Rs3,3);

	// geocentric inertial position
	double tau1 = (Mjd1-Mjd2)*86400.0;
	double tau3 = (Mjd3-Mjd2)*86400.0;

	double a1 = tau3/(tau3-tau1);
	double a3 =-tau1/(tau3-tau1);

	double b1 = tau3/(6*(tau3-tau1))*((tau3-tau1)*(tau3-tau1)-tau3*tau3);
	double b3 =-tau1/(6*(tau3-tau1))*((tau3-tau1)*(tau3-tau1)-tau1*tau1);

	double **D1 = m_create(3,3);
	D1[0] = Lm1;
	D1[1] = Lm2;
	D1[2] = Lm3;
	D1 = m_trans(D1,3);
	double **D2 = m_create(3,3);
	D2[0] = Rs1;
	D2[1] = Rs2;
	D2[2] = Rs2;
	D2 = m_trans(D2,3);
	double **D = m_dot(m_inv(D1,3),3,3,D2,3,3);

	double d1s = (D[1][0])*a1-(D[1][1])+(D[1][2])*a3;
	double d2s = D[1][0]*b1+D[1][2]*b3;

	double Ccye = 2*v_dot(Lm2,3,Rs2,3);

	double *poly = v_create(9);
	poly[0]=  1.0;  // R2^8... polynomial
	poly[1]=  0.0;
	poly[2]=  -(d1s*d1s + d1s*Ccye + (v_norm(Rs2,3))*(v_norm(Rs2,3)));
	poly[3]=  0.0;
	poly[4]=  0.0;
	poly[5]=  -(GM_Earth)*(d2s*Ccye + 2*d1s*d2s);
	poly[6]=  0.0;
	poly[7]=  0.0;
	poly[8]=  -(GM_Earth)*(GM_Earth)*d2s*d2s;

	double *zeror, *zeroi;
	int nr = real_poly_roots(poly, 8, &zeror, &zeroi);

	double bigr2= -99999990.0;

	for(int j=0; j<nr; j++) {
	    if((zeror[j]>bigr2) && (zeroi[j]==0.0)) {
	        bigr2 = zeror[j];
	    }
	}

	double u = (GM_Earth)/(bigr2*bigr2*bigr2);

	double *C = v_create(3);
	C[0] = a1+b1*u;
	C[1] = -1;
	C[2] = a3+b3*u;

	double *temp = m_dot_v(D,3,3,v_mul_scalar(C,3,-1.0),3);
	double rho1 = temp[0]/(a1+b1*u);
	double rho2 = -temp[1];
	double rho3 = temp[2]/(a3+b3*u);

	double rhoold1 = rho1;
	double rhoold2 = rho2;
	double rhoold3 = rho3;

	rho2 = 99999999.9;
	double ll = 0;

	double *r1, *r2, *r3, magr1, magr2, magr3, *v2, theta, theta1, copa, *y, p, a, e, i, Omega, omega, M;
	double rdot, udot, tausqr, f1, g1, f3, g3, H1, H2, H3, G1, G2, G3, *Daux, tempAux;
	char *error;
	while((fabs(rhoold2-rho2) > 1e-12) && (ll <= 0 )) {
	    ll = ll + 1;
	    rho2 = rhoold2;

	    r1 = v_sum(Rs1,3,v_mul_scalar(Lm1,3,rho1),3);
	    r2 = v_sum(Rs2,3,v_mul_scalar(Lm2,3,rho2),3);
	    r3 = v_sum(Rs3,3,v_mul_scalar(Lm3,3,rho3),3);

	    magr1 = v_norm(r1,3);
	    magr2 = v_norm(r2,3);
	    magr3 = v_norm(r3,3);

	    gibbs(r1,r2,r3,&v2,&theta,&theta1,&copa,&error);

	    if ((strcmp(error, "          ok")) && (copa < M_PI/180.0)) {
	        hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3,&v2,&theta,&theta1,&copa,&error);
	    }

	    y = v_create(6);
	    y[0] = r2[0]; y[1] = r2[1]; y[2] = r2[2];
	    y[3] = r2[0]; y[4] = r2[1]; y[5] = r2[2];
	    elements(y, &p, &a, &e, &i, &Omega, &omega, &M);

	    if(ll <= 8) {
	        u = (GM_Earth)/(magr2*magr2*magr2);
	        rdot = v_dot(r2,3,v2,3)/magr2;
	        udot= (-3*(GM_Earth)*rdot)/(magr2*magr2*magr2*magr2);

	        tausqr = tau1*tau1;
	        f1 =  1 - 0.5*u*tausqr -(1/6.0)*udot*tausqr*tau1  - (1/24.0) * u*u*tausqr*tausqr - (1/30.0)*u*udot*tausqr*tausqr*tau1;
	        g1 = tau1 - (1/6.0)*u*tau1*tausqr - (1/12.0) * udot*tausqr*tausqr - (1/120.0)*u*u*tausqr*tausqr*tau1 - (1/120.0)*u*udot*tausqr*tausqr*tausqr;
	        tausqr = tau3*tau3;
	        f3 = 1 - 0.5*u*tausqr -(1/6.0)*udot*tausqr*tau3 - (1/24.0) * u*u*tausqr*tausqr - (1/30.0)*u*udot*tausqr*tausqr*tau3;
	        g3 = tau3 - (1/6.0)*u*tau3*tausqr - (1/12.0) * udot*tausqr*tausqr - (1/120.0)*u*u*tausqr*tausqr*tau3 - (1/120.0)*u*udot*tausqr*tausqr*tausqr;
	    }
	    else {

	        theta  = angl(r1,r2);
	        theta1 = angl(r2,r3);

	        f1= 1 - ( (magr1*(1 - cos(theta)) / p ) );
	        g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );
	        f3= 1 - ( (magr3*(1 - cos(theta1)) / p ) );
	        g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );
	    }

	    C[0] = g3/(f1*g3-f3*g1);
	    C[1] = -1;
	    C[2] =-g1/(f1*g3-f3*g1);

	    H1 = (GM_Earth)*tau3/12.0;
	    H3 =-(GM_Earth)*tau1/12.0;
	    H2 = H1-H3;

	    G1 = -tau3/(tau1*(tau3-tau1));
	    G3 = -tau1/(tau3*(tau3-tau1));
	    G2 = G1-G3;

	    Daux = v_create(3);
	    Daux[0] = G1+H1/(magr1*magr1*magr1);
	    Daux[1] = G2+H2/(magr2*magr2*magr2);
	    Daux[2] = G3+H3/(magr3*magr3*magr3);

	    tempAux = -v_dot(Daux,3,C,3);
	    rhoold1 = tempAux/(a1+b1*u);
	    rhoold2 = -tempAux;
	    rhoold3 = tempAux/(a3+b3*u);

	    r1 = v_sum(Rs1,3,v_mul_scalar(Lm1,3,rhoold1),3);
		r2 = v_sum(Rs2,3,v_mul_scalar(Lm2,3,rhoold2),3);
		r3 = v_sum(Rs3,3,v_mul_scalar(Lm3,3,rhoold3),3);
	}
	
	r1 = v_sum(Rs1,3,v_mul_scalar(Lm1,3,rho1),3);
	r2 = v_sum(Rs2,3,v_mul_scalar(Lm2,3,rho2),3);
	r3 = v_sum(Rs3,3,v_mul_scalar(Lm3,3,rho3),3);

	*r = r2;
	*v = v2;
}


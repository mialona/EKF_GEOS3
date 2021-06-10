/** @file JPL_Eph_DE430.c
 *  @brief Sun, moon, and nine major planets' equatorial position using JPL Ephemerides.
 *
 *  This driver contains the code for the computation of the
 *  sun, moon, and nine major planets' equatorial position using JPL Ephemerides.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/global.h"
#include "../includes/m_utils.h"
#include "../includes/Cheb3D.h"

#include <stdio.h>
#include <math.h>


void JPL_Eph_DE430(double Mjd_TDB, double **r_Mercury, double **r_Venus,
				   double **r_Earth, double **r_Mars, double **r_Jupiter,
				   double **r_Saturn, double **r_Uranus, double **r_Neptune,
				   double **r_Pluto, double **r_Moon, double **r_Sun) {
	extern double **PC;
	extern int fPC, cPC;
	
	double JD = Mjd_TDB + 2400000.5;
	
	int i = 0, j = 0;
	for(j=0; j<fPC; j++) {
		if(PC[j][0]<=JD && JD<=PC[j][1]) {
			i = j++;
			break;
		}
	}
	
	double *PCtemp = v_create(fPC);
	for(j=0; j<cPC; j++) {
		PCtemp[j] = PC[i][j];
	}

	double t1 = PCtemp[0]-2400000.5; // MJD at start of interval

	double dt = Mjd_TDB - t1;
	
	int temp[4];
	for(j=0; j<4; j++) {
		temp[j] = 231 + 13*j;
	}
	
	double *Cx_Earth = v_create(26);
	double *Cy_Earth = v_create(26);
	double *Cz_Earth = v_create(26);
	for(j=0; j<13; j++) {
		Cx_Earth[j] = PCtemp[temp[0] + j -1];
		Cy_Earth[j] = PCtemp[temp[1] + j -1];
		Cz_Earth[j] = PCtemp[temp[2] + j -1];
	}
	
	for(j=0; j<4; j++) {
		temp[j] += 39;
	}
	
	double *Cx = v_create(13);
	double *Cy = v_create(13);
	double *Cz = v_create(13);
	for(j=0; j<13; j++) {
		Cx[j] = PCtemp[temp[0] + j -1];
		Cy[j] = PCtemp[temp[1] + j -1];
		Cz[j] = PCtemp[temp[2] + j -1];
	}
	
	for(j=0; j<13; j++) {
		Cx_Earth[j + 13] = Cx[j];
		Cy_Earth[j + 13] = Cy[j];
		Cz_Earth[j + 13] = Cz[j];
	}
	
	double Mjd0;
	if((0<=dt) && (dt<=16)) {
		j = 0;
		Mjd0 = t1;
	}
	else if((16<dt) && (dt<=32)) {
		j = 1;
		Mjd0 = t1+16*j;
	}
	
	*r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, &Cx_Earth[13*j],
						  &Cy_Earth[13*j], &Cz_Earth[13*j]);
	
	*r_Earth = v_mul_scalar(*r_Earth,3,1e3);

	v_free(Cx_Earth,26);
	v_free(Cy_Earth,26);
	v_free(Cz_Earth,26);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 441 + 13*j;
	}
	
	double *Cx_Moon = v_create(104);
	double *Cy_Moon = v_create(104);
	double *Cz_Moon = v_create(104);
	for(j=0; j<13; j++) {
		Cx_Moon[j] = PCtemp[temp[0] + j -1];
		Cy_Moon[j] = PCtemp[temp[1] + j -1];
		Cz_Moon[j] = PCtemp[temp[2] + j -1];
	}
	
	for(i=0; i<7; i++) {
		for(j=0; j<4; j++) {
			temp[j] += 39;
		}
		for(j=0; j<13; j++) {
			Cx[i] = PCtemp[temp[0] + j -1];
			Cy[i] = PCtemp[temp[1] + j -1];
			Cz[i] = PCtemp[temp[2] + j -1];
		}
		for(j=0; j<13; j++) {
			Cx_Moon[13*(i+1)+j] = Cx[j];
			Cy_Moon[13*(i+1)+j] = Cy[j];
			Cz_Moon[13*(i+1)+j] = Cz[j];
		}
	}
	
	if((0<=dt) && (dt<=4)) {
		j = 0;
		Mjd0 = t1;
	}
	else if((4<dt) && (dt<=8)) {
		j = 1;
		Mjd0 = t1+4*j;
	}
	else if((8<dt) && (dt<=12)) {
		j = 2;
		Mjd0 = t1+4*j;
	}
	else if((12<dt) && (dt<=16)) {
		j = 3;
		Mjd0 = t1+4*j;
	}
	else if((16<dt) && (dt<=20)) {
		j = 4;
		Mjd0 = t1+4*j;
	}
	else if((20<dt) && (dt<=24)) {
		j = 5;
		Mjd0 = t1+4*j;
	}
	else if((24<dt) && (dt<=28)) {
		j = 6;
		Mjd0 = t1+4*j;
	}
	else if((28<dt) && (dt<=32)) {
		j = 7;
		Mjd0 = t1+4*j;
	}
	
	*r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, &Cx_Moon[13*j],
						  &Cy_Moon[13*j], &Cz_Moon[13*j]);

	*r_Moon = v_mul_scalar(*r_Moon,3,1e3);
	
	v_free(Cx_Moon,104);
	v_free(Cy_Moon,104);
	v_free(Cz_Moon,104);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 753 + 11*j;
	}
	
	double *Cx_Sun = v_create(22);
	double *Cy_Sun = v_create(22);
	double *Cz_Sun = v_create(22);
	for(j=0; j<11; j++) {
		Cx_Sun[j] = PCtemp[temp[0] + j -1];
		Cy_Sun[j] = PCtemp[temp[1] + j -1];
		Cz_Sun[j] = PCtemp[temp[2] + j -1];
	}
	
	for(j=0; j<4; j++) {
		temp[j] += 33;
	}
	
	for(j=0; j<11; j++) {
		Cx[j] = PCtemp[temp[0] + j -1];
		Cy[j] = PCtemp[temp[1] + j -1];
		Cz[j] = PCtemp[temp[2] + j -1];
	}
	
	for(j=0; j<11; j++) {
		Cx_Sun[j + 11] = Cx[j];
		Cx_Sun[j + 11] = Cy[j];
		Cx_Sun[j + 11] = Cz[j];
	}
	
	if((0<=dt) && (dt<=16)) {
		j = 0;
		Mjd0 = t1;
	}
	else if((16<dt) && (dt<=32)) {
		j = 1;
		Mjd0 = t1+16*j;
	}
	
	*r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, &Cx_Sun[11*j],
						  &Cy_Sun[11*j], &Cz_Sun[11*j]);

	*r_Sun = v_mul_scalar(*r_Sun,3,1e3);
	
	v_free(Cx_Sun,22);
	v_free(Cy_Sun,22);
	v_free(Cz_Sun,22);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 3 + 14*j;
	}
	
	double *Cx_Mercury = v_create(56);
	double *Cy_Mercury = v_create(56);
	double *Cz_Mercury = v_create(56);
	for(j=0; j<14; j++) {
		Cx_Mercury[j] = PCtemp[temp[0] + j -1];
		Cy_Mercury[j] = PCtemp[temp[1] + j -1];
		Cz_Mercury[j] = PCtemp[temp[2] + j -1];
	}
	
	for(i=0; i<3; i++) {
		for(j=0; j<4; j++) {
			temp[j] += 42;
		}
		for(j=0; j<14; j++) {
			Cx[i] = PCtemp[temp[0] + j -1];
			Cy[i] = PCtemp[temp[1] + j -1];
			Cz[i] = PCtemp[temp[2] + j -1];
		}
		for(j=0; j<14; j++) {
			Cx_Mercury[14+j] = Cx[j];
			Cy_Mercury[14+j] = Cy[j];
			Cz_Mercury[14+j] = Cz[j];
		}
	}
	
	if((0<=dt) && (dt<=8)) {
		j = 0;
		Mjd0 = t1;
	}
	else if((8<dt) && (dt<=16)) {
		j = 1;
		Mjd0 = t1+8*j;
	}
	else if((16<dt) && (dt<=24)) {
		j = 2;
		Mjd0 = t1+8*j;
	}
	else if((24<dt) && (dt<=32)) {
		j = 3;
		Mjd0 = t1+8*j;
	}
	
	*r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, &Cx_Mercury[14*j],
						  &Cy_Mercury[14*j], &Cz_Mercury[14*j]);
	
	*r_Mercury = v_mul_scalar(*r_Mercury,3,1e3);

	v_free(Cx_Mercury,56);
	v_free(Cy_Mercury,56);
	v_free(Cz_Mercury,56);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 171 + 10*j;
	}
	
	double *Cx_Venus = v_create(20);
	double *Cy_Venus = v_create(20);
	double *Cz_Venus = v_create(20);
	for(j=0; j<10; j++) {
		Cx_Venus[j] = PCtemp[temp[0] + j -1];
		Cy_Venus[j] = PCtemp[temp[1] + j -1];
		Cz_Venus[j] = PCtemp[temp[2] + j -1];
	}
	
	for(j=0; j<3; j++) {
		temp[j] += 30;
	}
	
	for(j=0; j<10; j++) {
		Cx[j] = PCtemp[temp[0] + j -1];
		Cy[j] = PCtemp[temp[1] + j -1];
		Cz[j] = PCtemp[temp[2] + j -1];
	}
	
	for(j=0; j<10; j++) {
		Cx_Venus[j + 10] = Cx[j];
		Cy_Venus[j + 10] = Cy[j];
		Cz_Venus[j + 10] = Cz[j];
	}
	
	if((0<=dt) && (dt<=16)) {
		j = 0;
		Mjd0 = t1;
	}
	else if((16<dt) && (dt<=32)) {
		j = 1;
		Mjd0 = t1+16*j;
	}
	
	*r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, &Cx_Venus[10*j],
						  &Cy_Venus[10*j], &Cz_Venus[10*j]);

	*r_Venus = v_mul_scalar(*r_Venus,3,1e3);
	
	v_free(Cx_Venus,56);
	v_free(Cy_Venus,56);
	v_free(Cz_Venus,56);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 309 + 11*j;
	}
	
	double *Cx_Mars = v_create(11);
	double *Cy_Mars = v_create(11);
	double *Cz_Mars = v_create(11);
	for(j=0; j<11; j++) {
		Cx_Mars[j] = PCtemp[temp[0] + j -1];
		Cy_Mars[j] = PCtemp[temp[1] + j -1];
		Cz_Mars[j] = PCtemp[temp[2] + j -1];
	}
	
	j = 0;
	Mjd0 = t1;
	
	*r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, &Cx_Mars[11*j],
						  &Cy_Mars[11*j], &Cz_Mars[11*j]);

	*r_Mars = v_mul_scalar(*r_Mars,3,1e3);
	
	v_free(Cx_Mars,11);
	v_free(Cy_Mars,11);
	v_free(Cz_Mars,11);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 342 + 8*j;
	}
	
	double *Cx_Jupiter = v_create(8);
	double *Cy_Jupiter = v_create(8);
	double *Cz_Jupiter = v_create(8);
	for(j=0; j<8; j++) {
		Cx_Jupiter[j] = PCtemp[temp[0] + j -1];
		Cy_Jupiter[j] = PCtemp[temp[1] + j -1];
		Cz_Jupiter[j] = PCtemp[temp[2] + j -1];
	}
	
	j = 0;
	Mjd0 = t1;
	
	*r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, &Cx_Jupiter[8*j],
						  &Cy_Jupiter[8*j], &Cz_Jupiter[8*j]);

	*r_Jupiter = v_mul_scalar(*r_Jupiter,3,1e3);
	
	v_free(Cx_Jupiter,8);
	v_free(Cy_Jupiter,8);
	v_free(Cz_Jupiter,8);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 366 + 7*j;
	}
	
	double *Cx_Saturn = v_create(7);
	double *Cy_Saturn = v_create(7);
	double *Cz_Saturn = v_create(7);
	for(j=0; j<7; j++) {
		Cx_Saturn[j] = PCtemp[temp[0] + j -1];
		Cy_Saturn[j] = PCtemp[temp[1] + j -1];
		Cz_Saturn[j] = PCtemp[temp[2] + j -1];
	}
	
	j = 0;
	Mjd0 = t1;
	
	*r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, &Cx_Saturn[7*j],
						  &Cy_Saturn[7*j], &Cz_Saturn[7*j]);

	*r_Saturn = v_mul_scalar(*r_Saturn,3,1e3);
	
	v_free(Cx_Saturn,7);
	v_free(Cy_Saturn,7);
	v_free(Cz_Saturn,7);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 387 + 6*j;
	}
	
	double *Cx_Uranus = v_create(6);
	double *Cy_Uranus = v_create(6);
	double *Cz_Uranus = v_create(6);
	for(j=0; j<6; j++) {
		Cx_Uranus[j] = PCtemp[temp[0] + j -1];
		Cy_Uranus[j] = PCtemp[temp[1] + j -1];
		Cz_Uranus[j] = PCtemp[temp[2] + j -1];
	}
	
	j = 0;
	Mjd0 = t1;
	
	*r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, &Cx_Uranus[6*j],
						  &Cy_Uranus[6*j], &Cz_Uranus[6*j]);

	*r_Uranus = v_mul_scalar(*r_Uranus,3,1e3);
	
	v_free(Cx_Uranus,6);
	v_free(Cy_Uranus,6);
	v_free(Cz_Uranus,6);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 405 + 6*j;
	}
	
	double *Cx_Neptune = v_create(6);
	double *Cy_Neptune = v_create(6);
	double *Cz_Neptune = v_create(6);
	for(j=0; j<6; j++) {
		Cx_Neptune[j] = PCtemp[temp[0] + j -1];
		Cy_Neptune[j] = PCtemp[temp[1] + j -1];
		Cz_Neptune[j] = PCtemp[temp[2] + j -1];
	}
	
	j = 0;
	Mjd0 = t1;
	
	*r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, &Cx_Neptune[6*j],
						  &Cy_Neptune[6*j], &Cz_Neptune[6*j]);

	*r_Neptune = v_mul_scalar(*r_Neptune,3,1e3);
	
	v_free(Cx_Neptune,6);
	v_free(Cy_Neptune,6);
	v_free(Cz_Neptune,6);
	
	
	//--------------------------------------------------------------------------
	
	
	for(j=0; j<4; j++) {
		temp[j] = 423 + 6*j;
	}
	
	double *Cx_Pluto = v_create(6);
	double *Cy_Pluto = v_create(6);
	double *Cz_Pluto = v_create(6);
	for(j=0; j<6; j++) {
		Cx_Pluto[j] = PCtemp[temp[0] + j -1];
		Cy_Pluto[j] = PCtemp[temp[1] + j -1];
		Cz_Pluto[j] = PCtemp[temp[2] + j -1];
	}
	
	j = 0;
	Mjd0 = t1;
	
	*r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, &Cx_Pluto[6*j],
						  &Cy_Pluto[6*j], &Cz_Pluto[6*j]);

	*r_Pluto = v_mul_scalar(*r_Pluto,3,1e3);
	
	v_free(Cx_Pluto,6);
	v_free(Cy_Pluto,6);
	v_free(Cz_Pluto,6);
	
	
	//--------------------------------------------------------------------------
	
	
//	for(j=0; j<2; j++) {
//		temp[j] = 819 + 10*j;
//	}
//
//	double *Cx_Nutations = v_create(40);
//	double *Cy_Nutations = v_create(40);
//	for(j=0; j<10; j++) {
//		Cx_Nutations[j] = PCtemp[temp[0] + j -1];
//		Cy_Nutations[j] = PCtemp[temp[1] + j -1];
//	}
//
//	for(i=0; i<3; i++) {
//		for(j=0; j<2; j++) {
//			temp[j] += 20;
//		}
//		Cx[i] = PCtemp[temp[0] + i -1];
//		Cy[i] = PCtemp[temp[1] + i -1];
//		for(j=0; j<10; j++) {
//			Cx_Nutations[10*(i+1)+j] = Cx[j];
//			Cy_Nutations[10*(i+1)+j] = Cy[j];
//		}
//	}
//
//	if((0<=dt) && (dt<=8)) {
//		j = 0;
//		Mjd0 = t1;
//	}
//	else if((8<dt) && (dt<=16)) {
//		j = 1;
//		Mjd0 = t1+8*j;
//	}
//	else if((16<dt) && (dt<=24)) {
//		j = 2;
//		Mjd0 = t1+8*j;
//	}
//	else if((24<dt) && (dt<=32)) {
//		j = 3;
//		Mjd0 = t1+8*j;
//	}
//
//	double *Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, &(Cx_Nutations[10*j]),
//						  &(Cx_Nutations[10*j]), v_create(10));
//
//	v_free(Cx_Nutations,10);
//	v_free(Cy_Nutations,10);
//
//
//	//--------------------------------------------------------------------------
//
//
//	for(j=0; j<3; j++) {
//		temp[j] = 899 + 10*j;
//	}
//
//	double *Cx_Librations = v_create(40);
//	double *Cy_Librations = v_create(40);
//	double *Cz_Librations = v_create(40);
//	for(j=0; j<10; j++) {
//		Cx_Librations[j] = PCtemp[temp[0] + j -1];
//		Cy_Librations[j] = PCtemp[temp[1] + j -1];
//		Cz_Librations[j] = PCtemp[temp[2] + j -1];
//	}
//
//	for(i=0; i<3; i++) {
//		for(j=0; j<3; j++) {
//			temp[j] += 30;
//		}
//		Cx[i] = PCtemp[temp[0] + i -1];
//		Cy[i] = PCtemp[temp[1] + i -1];
//		Cz[i] = PCtemp[temp[2] + i -1];
//		for(j=0; j<10; j++) {
//			Cx_Librations[10*(i+1)+j] = Cx[j];
//			Cy_Librations[10*(i+1)+j] = Cy[j];
//			Cz_Librations[10*(i+1)+j] = Cz[j];
//		}
//	}
//
//	if((0<=dt) && (dt<=8)) {
//		j = 0;
//		Mjd0 = t1;
//	}
//	else if((8<dt) && (dt<=16)) {
//		j = 1;
//		Mjd0 = t1+8*j;
//	}
//	else if((16<dt) && (dt<=24)) {
//		j = 2;
//		Mjd0 = t1+8*j;
//	}
//	else if((24<dt) && (dt<=32)) {
//		j = 3;
//		Mjd0 = t1+8*j;
//	}
//
//	double *Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, &(Cx_Librations[10*j]),
//						  &(Cx_Librations[10*j]), &(Cz_Librations[10*j]));
//
//	v_free(Cx_Librations,10);
//	v_free(Cy_Librations,10);
//	v_free(Cz_Librations,10);
	
	
	//--------------------------------------------------------------------------
	
	
	double EMRAT = 81.30056907419062; // DE430
	double EMRAT1 = 1/(1+EMRAT);
	*r_Earth = v_sum(*r_Earth,3,v_mul_scalar(*r_Moon,3,-EMRAT1),3);

	double *aux = v_mul_scalar(*r_Earth,3,-1.0);
	*r_Mercury = v_sum(aux,3,*r_Mercury,3);
	*r_Venus = v_sum(aux,3,*r_Venus,3);
	*r_Mars = v_sum(aux,3,*r_Mars,3);
	*r_Jupiter = v_sum(aux,3,*r_Jupiter,3);
	*r_Saturn = v_sum(aux,3,*r_Saturn,3);
	*r_Uranus = v_sum(aux,3,*r_Uranus,3);
	*r_Neptune = v_sum(aux,3,*r_Neptune,3);
	*r_Pluto = v_sum(aux,3,*r_Pluto,3);
	*r_Sun = v_sum(aux,3,*r_Sun,3);
	v_free(aux,3);
}


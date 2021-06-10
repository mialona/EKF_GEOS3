/** @file IERS.c
 *  @brief IERS time and polar motion data.
 *
 *  This driver contains the code for the 
 *  management of IERS time and polar motion data.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug Removed pre-initiated interp.
 */

#include "../includes/const.h"
#include "../includes/global.h"
#include "../includes/m_utils.h"

#include <stdio.h>
#include <math.h>


void IERS(double Mjd_UTC, char interp,
			double *x_pole, double *y_pole, double *UT1_UTC,
			double *LOD, double *dpsi, double *deps, double *dx_pole,
			double *dy_pole, double *TAI_UTC) {
	extern double** eopdata;

	double mjd = floor(Mjd_UTC);
	int i = 0;
	for(int j=0; j<ceopdata; j++) {
		if(mjd==eopdata[3][j]) {
			i = j++;
			break;
		}
	}
	if(i == 0) {
		printf("IERS: Not find mjd\n");
        exit(EXIT_FAILURE);
	}

	if(interp =='l') {
		// linear interpolation
		double *preeop = v_extract(eopdata,feopdata,ceopdata,i);
		double *nexteop = v_extract(eopdata,feopdata,ceopdata,i+1);
		double mfme = 1440.0*(Mjd_UTC-floor(Mjd_UTC));
		double fixf = mfme/1440.0;

		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		*x_pole  = (preeop[4]+(nexteop[4]-preeop[4])*fixf)/(Arcs); // Pole coordinate [rad]
		*y_pole  = (preeop[5]+(nexteop[5]-preeop[5])*fixf)/(Arcs); // Pole coordinate [rad]
		*UT1_UTC = (preeop[6]+(nexteop[6]-preeop[6])*fixf);
		*LOD     = (preeop[7]+(nexteop[7]-preeop[7])*fixf);
		*dpsi    = (preeop[8]+(nexteop[8]-preeop[8])*fixf)/(Arcs);
		*deps    = (preeop[9]+(nexteop[9]-preeop[9])*fixf)/(Arcs);
		*dx_pole = (preeop[10]+(nexteop[10]-preeop[10])*fixf)/(Arcs); // Pole coordinate [rad]
		*dy_pole = (preeop[11]+(nexteop[11]-preeop[11])*fixf)/(Arcs); // Pole coordinate [rad]
		*TAI_UTC = (preeop[12]);
	}
	else if (interp =='n') {
		double *auxeop = v_extract(eopdata,feopdata,ceopdata,i);
		
		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		*x_pole  = auxeop[4]/(Arcs);  // Pole coordinate [rad]
		*y_pole  = auxeop[5]/(Arcs);  // Pole coordinate [rad]
		*UT1_UTC = auxeop[6];             // UT1-UTC time difference [s]
		*LOD     = auxeop[7];             // Length of day [s]
		*dpsi    = auxeop[8]/(Arcs);
		*deps    = auxeop[9]/(Arcs);
		*dx_pole = auxeop[10]/(Arcs); // Pole coordinate [rad]
		*dy_pole = auxeop[11]/(Arcs); // Pole coordinate [rad]
		*TAI_UTC = auxeop[12];            // TAI-UTC time difference [s]
	}
}


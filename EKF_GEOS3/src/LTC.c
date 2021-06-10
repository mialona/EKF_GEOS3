/** @file LTC.c
 *  @brief Transformation from Greenwich meridian system to 
 *  local tangent coordinates.
 *
 *  This driver contains the code for the 
 *  transformation from Greenwich meridian system to 
 *  local tangent coordinates.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"
#include "../includes/R_y.h"
#include "../includes/R_z.h"

#include <stdio.h>


double **LTC(double lon, double lat) {
	double **M = m_dot(R_y(-lat),3,3,R_z(lon),3,3);

	double Aux;
	for(int j=0; j<3; j++) {
		Aux = M[0][j];
		M[0][j] = M[1][j];
		M[1][j] = M[2][j];
		M[2][j] = Aux;
	}
	
	return M;
}


/** @file Legendre.c
 *  @brief Legendre coefficients.
 *
 *  This driver contains the code for the 
 *  calculation of Legendre coefficients.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/m_utils.h"

#include <stdio.h>
#include <math.h>


void Legendre(int n, int m, double fi, double ***pnm, double ***dpnm) {
	
	*pnm = m_zeros(n+2,m+2);
	*dpnm = m_zeros(n+2,m+2);

	(*pnm)[1][1]=1.0;
	(*dpnm)[1][1]=0.0;
	(*pnm)[2][2]=sqrt(3.0)*cos(fi);
	(*dpnm)[2][2]=-sqrt(3.0)*sin(fi);

	// diagonal coefficients
	for(int i=2; i<=n; i++) {
		(*pnm)[i+1][i+1] = sqrt((2.0*i+1.0)/(2.0*i))*cos(fi)*(*pnm)[i][i];
	}
	for(int i=2; i<=n; i++) {
		(*dpnm)[i+1][i+1] = sqrt((2.0*i+1.0)/(2.0*i))*((cos(fi)*(*dpnm)[i][i])-(sin(fi)*(*pnm)[i][i]));
	}

	// horizontal first step coefficients
	for(int i=1; i<=n; i++) {
		(*pnm)[i+1][i] = sqrt(2.0*i+1.0)*sin(fi)*(*pnm)[i][i];
	}
	for(int i=1; i<=n; i++) {
		(*dpnm)[i+1][i] = sqrt(2.0*i+1.0)*((cos(fi)*(*pnm)[i][i])+(sin(fi)*(*dpnm)[i][i]));
	}

	// horizontal second step coefficients
	int j = 0;
	int k = 2;
	while(1) {
		for(int i=k; i<=n; i++) {
			(*pnm)[i+1][j+1]=sqrt((2.0*i+1.0)/((1.0*i-j)*(1.0*i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*(*pnm)[i][j+1])
				-(sqrt(((1.0*i+j-1.0)*(1.0*i-j-1.0))/(2.0*i-3.0))*(*pnm)[i-1][j+1]));
		}
		j = j+1;
		k = k+1;
		if(j>m) {
			break;
		}
	}

	j = 0;
	k = 2;
	while(1) {
		for(int i=k; i<=n; i++) {
			(*dpnm)[i+1][j+1] = sqrt((2.0*i+1.0)/((1.0*i-1.0*j)*(1.0*i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*(*dpnm)[i][j+1])+
								(sqrt(2.0*i-1.0)*cos(fi)*(*pnm)[i][j+1])-(sqrt(((1.0 * i + j-1.0)*(1.0*i-j-1.0))/(2.0*i-3.0))*(*dpnm)[i-1][j+1]));
		}
		j = j+1;
		k = k+1;
		if(j>m) {
			break;
		}
	}
}

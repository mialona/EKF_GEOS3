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
	
	*pnm = m_zeros(n+1,m+1);
	*dpnm = m_zeros(n+1,m+1);

	(*pnm)[0][0]=1;
	(*dpnm)[0][0]=0;
	(*pnm)[1][1]=sqrt(3)*cos(fi);
	(*dpnm)[1][1]=-sqrt(3)*sin(fi);
	
	// diagonal coefficients
	for(int i=2; i<=n; i++) {
		(*pnm)[i][i]= sqrt((2.0*i+1)/(2*i))*cos(fi)*(*pnm)[i-1][i-1];
	}
	for(int i=2; i<=n; i++) {
		(*dpnm)[i][i]= sqrt((2.0*i+1)/(2*i))*((cos(fi)*(*dpnm)[i-1][i-1])-(sin(fi)*(*pnm)[i-1][i-1]));
	}
	
	// horizontal first step coefficients
	for(int i=1; i<=n; i++) {
		(*pnm)[i][i-1]= sqrt(2.0*i+1)*sin(fi)*(*pnm)[i-1][i-1];
	}
	for(int i=1; i<=n; i++) {
		(*dpnm)[i][i-1]= sqrt(2.0*i+1)*((cos(fi)*(*pnm)[i-1][i-1])+(sin(fi)*(*dpnm)[i-1][i-1]));
	}
	
	// horizontal second step coefficients
	int j = 0;
	int k = 2;
	while(1) {
		for(int i=k; i<=n; i++) {
			(*pnm)[i][j]=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*(*pnm)[i-1][j])
				-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))*(*pnm)[i-2][j]));
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
			(*dpnm)[i][j]=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*(*dpnm)[i-1][j])
				 +(sqrt(2.0*i-1)*cos(fi)*(*pnm)[i-1][j])-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))*(*dpnm)[i-2][j]));
		}
		j = j+1;
		k = k+1;
		if(j>m) {
			break;
		}
	}
}

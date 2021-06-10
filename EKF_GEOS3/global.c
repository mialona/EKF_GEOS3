/** @file globals.c
 *  @brief Global variables and reading data files.
 *
 *  This driver contains the code for the 
 *  reading data files and global variables.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "m_utils.h"
#include "Mjday.h"
#include "const.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void DE430Coeff(int f, int c) {
	extern double **PC;
	extern int fPC, cPC;
	PC = m_create(f,c);
	fPC = f; //2285
	cPC = c; //1020
	
	FILE *fp = fopen("data/DE430Coeff.txt","r");
	if(fp == NULL) {
		printf("Fail open DE430Coeff.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<f; i++) {
		for(int j=0; j<c; j++) {
			fscanf(fp,"%lf",&PC[i][j]);
		}
	}
	
	fclose(fp);
}

void GGM03S(int n) {
	extern double **Cnm, **Snm;
	extern int fCnm, cCnm, fSnm, cSnm;
	Cnm = m_zeros(n,n);
	Snm = m_zeros(n,n);
	fCnm = n; //182
	cCnm = n; //182
	fSnm = n; //182
	cSnm = n; //182
	
	double aux1, aux2;
	int f, c;
	
	FILE *fp = fopen("data/GGM03S.txt","r");
	if(fp == NULL) {
		printf("Fail open GGM03S.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<n; i++) {
		for(int j=0; j<=i; j++) {
			fscanf(fp,"%d%d%lf%lf%lf%lf",&f,&c,&Cnm[i+1][j+1],&Snm[i+1][j+1],&aux1,&aux2);
		}
	}
	
	fclose(fp);
}

void eop19620101(int c) {
	extern double **eopdata;
	extern int feopdata, ceopdata;
	eopdata = m_create(13,c);
	feopdata = 13;
	ceopdata = c; //21413
	
	FILE *fp = fopen("data/eop19620101.txt","r");
	if(fp == NULL) {
		printf("Fail open eop19620101.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for(int j=0; j<c; j++) {
		fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&(eopdata[0][j])
			   ,&(eopdata[1][j]),&(eopdata[2][j]),&(eopdata[3][j]),&(eopdata[4][j])
			   ,&(eopdata[5][j]),&(eopdata[6][j]),&(eopdata[7][j]),&(eopdata[8][j])
			   ,&(eopdata[9][j]),&(eopdata[10][j]),&(eopdata[11][j]),&(eopdata[12][j]));
	}
	
	fclose(fp);
}

void GEOS3(int f) {
	extern double **obs;
	extern int fobs, cobs;
	obs = m_create(f,4);
	fobs = f; //46
	cobs = 4;
	
	FILE *fp = fopen("data/GEOS3.txt","r");
	if(fp == NULL) {
		printf("Fail open GEOS3.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	int Y, MO, D, H, MI, S;
	double AZ, EL, DIST;
	char line[55], y[5], mo[3], d[3], h[3], mi[3], s[7], az[9], el[9], dist[10];
	for(int i=0; i<f; i++) {
		fgets(line,sizeof(line)+2,fp);
		
		strncpy(y,&(line[0]),4);
		y[4] = '\0';
		Y = atoi(y);
		
		strncpy(mo,&(line[5]),2);
		mo[2] = '\0';
		MO = atoi(mo);
		
		strncpy(d,&(line[8]),2);
		d[2] = '\0';
		D = atoi(d);
		
		strncpy(h,&(line[12]),2);
		h[2] = '\0';
		H = atoi(h);
		
		strncpy(mi,&(line[15]),2);
		mi[2] = '\0';
		MI = atoi(mi);
		
		strncpy(s,&(line[18]),6);
		s[6] = '\0';
		S = atoi(s);
		
		strncpy(az,&(line[25]),8);
		az[8] = '\0';
		AZ = atof(az);
		
		strncpy(el,&(line[35]),8);
		el[8] = '\0';
		EL = atof(el);
		
		strncpy(dist,&(line[44]),9);
		dist[9] = '\0';
		DIST = atof(dist);
		
		obs[i][0] = Mjday(Y,MO,D,H,MI,S);
		obs[i][1] = Rad*AZ;
		obs[i][2] = Rad*EL;
		obs[i][3] = 1e3*DIST;
	}
	
	fclose(fp);
}

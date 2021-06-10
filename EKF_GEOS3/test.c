/** @file test.c
 *  @brief Tester driver.
 *
 *  This driver contains the code for the testing.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "includes/global.h"
#include "includes/m_utils.h"
#include "includes/R_x.h"
#include "includes/R_y.h"
#include "includes/R_z.h"
#include "includes/position.h"
#include "includes/Mjday_TDB.h"
#include "includes/MeanObliquity.h"
#include "includes/Mjday.h"
#include "includes/NutAngles.h"
#include "includes/NutMatrix.h"
#include "includes/AzElPa.h"
#include "includes/PoleMatrix.h"
#include "includes/PrecMatrix.h"
#include "includes/Frac.h"
#include "includes/gmst.h"
#include "includes/timediff.h"
#include "includes/unit.h"
#include "includes/sign_.h"
#include "includes/EqnEquinox.h"
#include "includes/gast.h"
#include "includes/Geodetic.h"
#include "includes/Legendre.h"
#include "includes/LTC.h"
#include "includes/GHAMatrix.h"
#include "includes/EccAnom.h"
#include "includes/Cheb3D.h"
#include "includes/elements.h"
#include "includes/IERS.h"
#include "includes/angl.h"
#include "includes/gibbs.h"
#include "includes/hgibbs.h"
#include "includes/TimeUpdate.h"
#include "includes/MeasUpdate.h"
#include "includes/AccelPointMass.h"
#include "includes/AccelHarmonic.h"
#include "includes/JPL_Eph_DE430.h"
#include "includes/Accel.h"
#include "includes/G_AccelHarmonic.h"
#include "includes/VarEqn.h"
#include "includes/ode.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <unistd.h>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

/** @brief Comparator of two vectors of c components.
 *
 *  @param [in] v First vector.
 *  @param [in] w Second vector.
 *  @param [in] c Number of components.
 *  @return 0=error, 1=pass.
 */
int equals_vector(double *v, double *w, int c, double p) {
	int equal = 1;
	
    for(int i = 0; i < c; i++)
		if(fabs(v[i]-w[i]) > p){
			printf("%2.20lf %2.20lf\n",v[i],w[i]);
			equal = 0;}
	
	return equal;
}

/** @brief Comparator of two matrices of f rows and c columns.
 *
 *  @param [in] A First matrix.
 *  @param [in] B Second matrix.
 *  @param [in] f Number of rows.
 *  @param [in] c Number of columns.
 *  @return 0=error, 1=pass.
 */
int equals_matrix(double **A, double **B, int f, int c, double p) {
	int equal = 1;
	
    for(int i = 0; i < f; i++)
		for(int j = 0; j < c; j++)
			if(fabs(A[i][j]-B[i][j]) > p){
				printf("%2.20lf %2.20lf\n",A[i][j],B[i][j]);
				equal = 0;}
	
	return equal;
}


/** @brief Unit test for function v_sum.
 *
 *  @return 0=error, 1=pass.
 */
int v_sum_01() {
    int c = 3;
	
	double *v = v_create(c);
	v[0] = 1; v[1] = 2; v[2] = 3;
	
	double *w = v_create(c);
	w[0] = 5; w[1] = 1; w[2] = 4;
	
	double *x = v_create(c);
	x[0] = 6; x[1] = 3; x[2] = 7;
    
	double *r = v_sum(v,c,w,c);
    _assert(equals_vector(r,x,3,1e-10));
    
	v_free(v,c);
	v_free(w,c);
	v_free(x,c);
	v_free(r,c);
    
    return 0;
}

/** @brief Unit test for function v_mul_scalar.
 *
 *  @return 0=error, 1=pass.
 */
int v_mul_scalar_01() {
    double s = 2.0;
    int c = 3;
	
	double *v = v_create(c);
	v[0] = 1; v[1] = 2; v[2] = 3;
	
	double *x = v_create(c);
	x[0] = 2; x[1] = 4; x[2] = 6;
    
    double *r = v_mul_scalar(v,c,s);
    _assert(equals_vector(r,x,3,1e-10));
	
	v_free(v,c);
	v_free(x,c);
	v_free(r,c);
    
    return 0;
}

/** @brief Unit test for function v_norm.
 *
 *  @return 0=error, 1=pass.
 */
int v_norm_01() {
    int c = 5;
	
	double *v = v_create(c);
	v[0] = 2; v[1] = -2; v[2] = 3; v[3] = 2; v[4] = -2;
    
	double r = v_norm(v,c);
    _assert(r == 5);
    
	v_free(v,c);
	
    return 0;
}

/** @brief Unit test for function v_dot.
 *
 *  @return 0=error, 1=pass.
 */
int v_dot_01() {
    int c = 3;
	
	double *v = v_create(c);
	v[0] = 2; v[1] = -1; v[2] = 3;
	
	double *w = v_create(c);
	w[0] = 1; w[1] = 5; w[2] = 4;
    
	double r = v_dot(v,c,w,c);
    _assert(r == 9);
    
	v_free(v,c);
	v_free(w,c);
	
    return 0;
}

/** @brief Unit test for function v_dot.
 *
 *  @return 0=error, 1=pass.
 */
int v_cross_01() {
    int c = 3;
	
	double *v = v_create(c);
	v[0] = 2; v[1] = -1; v[2] = 3;
	
	double *w = v_create(c);
	w[0] = 1; w[1] = 5; w[2] = 4;
	
	double *sol = v_create(c);
	sol[0] = -19; sol[1] = -5; sol[2] = 11;
    
	double *r = v_cross(v,c,w,c);
    _assert(equals_vector(sol,r,c,1e-10));
	
	v_free(v,c);
	v_free(w,c);
	v_free(sol,c);
	
    return 0;
}


/** @brief Unit test for function v_assign.
 *
 *  @return 0=error, 1=pass.
 */
int v_extract_01() {
    int f = 3;
    int c = 4;
	int k = 2;
	
	double *v = v_create(f);
	v[0] = 8; v[1] = 0; v[2] = 3;
	
	double **A = m_create(f,c);
	A[0][0] = 0; A[0][1] = 2; A[0][2] = 8; A[0][3] = 0;
	A[1][0] = 1; A[1][1] = -1; A[1][2] = 0; A[1][3] = 0;
	A[2][0] = 0; A[2][1] = 1; A[2][2] = 3; A[2][3] = 5;
    
	double *r = v_extract(A,f,c,k);
    _assert(equals_vector(v,r,f,1e-10));
    
	v_free(v,f);
	m_free(A,f,c);
	v_free(r,f);
	
    return 0;
}

/** @brief Unit test for function v_assign.
 *
 *  @return 0=error, 1=pass.
 */
int v_assign_01() {
    int f = 3;
    int c = 4;
	int k = 2;
	
	double *v = v_create(f);
	v[0] = 2; v[1] = -1; v[2] = 3;
	
	double **A = m_create(f,c);
	A[0][0] = 0; A[0][1] = 2; A[0][2] = 8; A[0][3] = 0;
	A[1][0] = 1; A[1][1] = -1; A[1][2] = 0; A[1][3] = 0;
	A[2][0] = 0; A[2][1] = 1; A[2][2] = 0; A[2][3] = 5;
	
	double **C = m_create(f,c);
	C[0][0] = 0; C[0][1] = 2; C[0][2] = 2; C[0][3] = 0;
	C[1][0] = 1; C[1][1] = -1; C[1][2] = -1; C[1][3] = 0;
	C[2][0] = 0; C[2][1] = 1; C[2][2] = 3; C[2][3] = 5;
    
	double **R = v_assign(v,f,A,f,c,k);
    _assert(equals_matrix(C,R,f,c,1e-10));
    
	v_free(v,f);
	m_free(A,f,c);
	m_free(C,f,c);
	m_free(R,f,c);
	
    return 0;
}

/** @brief Unit test for function v_dot_m.
 *
 *  @return 0=error, 1=pass.
 */
int v_dot_m_01() {
    int f = 3;
    int c = 4;
	
	double *v = v_create(f);
	v[0] = 2; v[1] = -1; v[2] = 3;
	
	double **A = m_create(f,c);
	A[0][0] = 0; A[0][1] = 2; A[0][2] = 8; A[0][3] = 0;
	A[1][0] = 1; A[1][1] = -1; A[1][2] = 0; A[1][3] = 0;
	A[2][0] = 0; A[2][1] = 1; A[2][2] = 0; A[2][3] = 5;
	
	double *w = v_create(c);
	w[0] = -1; w[1] = 8; w[2] = 16; w[3] = 15;
    
	double *r = v_dot_m(v,f,A,f,c);
    _assert(equals_vector(w,r,c,1e-10));
    
	v_free(v,f);
	m_free(A,f,c);
	v_free(w,c);
	v_free(r,c);
	
    return 0;
}

/** @brief Unit test for function m_dot_v.
 *
 *  @return 0=error, 1=pass.
 */
int m_dot_v_01() {
    int f = 3;
    int c = 4;
	
	double *v = v_create(c);
	v[0] = 2; v[1] = -1; v[2] = 3; v[3] = 3;
	
	double **A = m_create(f,c);
	A[0][0] = 0; A[0][1] = 2; A[0][2] = 8; A[0][3] = 0;
	A[1][0] = 1; A[1][1] = -1; A[1][2] = 0; A[1][3] = 0;
	A[2][0] = 0; A[2][1] = 1; A[2][2] = 0; A[2][3] = 5;
	
	double *w = v_create(f);
	w[0] = 22; w[1] = 3; w[2] = 14;
    
	double *r = m_dot_v(A,f,c,v,c);
    _assert(equals_vector(w,r,f,1e-10));
    
	v_free(v,c);
	m_free(A,f,c);
	v_free(w,f);
	v_free(r,f);
	
    return 0;
}



/** @brief Unit test for function m_zeros.
 *
 *  @return 0=error, 1=pass.
 */
int zeros_01() {
    int f = 3;
    int c = 4;
	
	double **A = m_create(f,c);
	A[0][0] = 0; A[0][1] = 0; A[0][2] = 0; A[0][3] = 0;
	A[1][0] = 0; A[1][1] = 0; A[1][2] = 0; A[1][3] = 0;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 0; A[2][3] = 0;
    
	double **R = m_zeros(f,c);
    _assert(equals_matrix(A,R,f,c,1e-10));
    
	m_free(A,f,c);
	m_free(R,f,c);
	
    return 0;
}

/** @brief Unit test for function m_eye.
 *
 *  @return 0=error, 1=pass.
 */
int eye_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 1; A[0][1] = 0; A[0][2] = 0;
	A[1][0] = 0; A[1][1] = 1; A[1][2] = 0;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 1;
    
	double **R = m_eye(n);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
    
    return 0;
}

/** @brief Unit test for function m_sum.
 *
 *  @return 0=error, 1=pass.
 */
int sum_01() {
    int f = 3;
    int c = 4;
	
	double **A = m_create(f,c);
	A[0][0] = 0; A[0][1] = 2; A[0][2] = 8; A[0][3] = 0;
	A[1][0] = 1; A[1][1] = -1; A[1][2] = 0; A[1][3] = 0;
	A[2][0] = 0; A[2][1] = 1; A[2][2] = 0; A[2][3] = 5;
	
	double **B = m_create(f,c);
	B[0][0] = 2; B[0][1] = 0; B[0][2] = 0; B[0][3] = 0;
	B[1][0] = 7; B[1][1] = -2; B[1][2] = 1; B[1][3] = 0;
	B[2][0] = 0; B[2][1] = -3; B[2][2] = 0; B[2][3] = 2;
	
	double **C = m_create(f,c);
	C[0][0] = 2; C[0][1] = 2; C[0][2] = 8; C[0][3] = 0;
	C[1][0] = 8; C[1][1] = -3; C[1][2] = 1; C[1][3] = 0;
	C[2][0] = 0; C[2][1] = -2; C[2][2] = 0; C[2][3] = 7;
    
	double **R = m_sum(A,f,c,B,f,c);
    _assert(equals_matrix(R,C,f,c,1e-10));
    
	m_free(A,f,c);
	m_free(B,f,c);
	m_free(C,f,c);
	m_free(R,f,c);
    
    return 0;
}

/** @brief Unit test for function m_dot.
 *
 *  @return 0=error, 1=pass.
 */
int m_dot_01() {
    int f = 3;
    int c = 4;
	
	double **A = m_create(f,c);
	
	A[0][0] = 0; A[0][1] = 2; A[0][2] = 8; A[0][3] = 0;
	A[1][0] = 1; A[1][1] = -1; A[1][2] = 0; A[1][3] = 0;
	A[2][0] = 0; A[2][1] = 1; A[2][2] = 0; A[2][3] = 5;
	
	double **B = m_create(c,f);
	
	B[0][0] = 2; B[0][1] = 0; B[0][2] = 0;
	B[1][0] = 7; B[1][1] = -2; B[1][2] = 1;
	B[2][0] = 0; B[2][1] = -3; B[2][2] = 0;
	B[3][0] = 1; B[3][1] = 0; B[3][2] = 0;
	
	double **C = m_create(f,f);
	
	C[0][0] = 14; C[0][1] = -28; C[0][2] = 2;
	C[1][0] = -5; C[1][1] = 2; C[1][2] = -1;
	C[2][0] = 12; C[2][1] = -2; C[2][2] = 1;
    
	double **R = m_dot(A,f,c,B,c,f);
    _assert(equals_matrix(R,C,f,f,1e-10));
    
	m_free(A,f,c);
	m_free(B,c,f);
	m_free(C,f,f);
	m_free(R,f,f);
    
    return 0;
}

/** @brief Unit test for function m_inv.
 *
 *  @return 0=error, 1=pass.
 */
int inv_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 1; A[0][1] = -1; A[0][2] = 0;
	A[1][0] = 0; A[1][1] = 1; A[1][2] = 0;
	A[2][0] = 2; A[2][1] = 0; A[2][2] = 1;
	
	double **C = m_create(n,n);
	C[0][0] = 1; C[0][1] = 1; C[0][2] = 0;
	C[1][0] = 0; C[1][1] = 1; C[1][2] = 0;
	C[2][0] = -2; C[2][1] = -2; C[2][2] = 1;
    
	double **R = m_inv(A,n);
    _assert(equals_matrix(R,C,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(C,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function m_trans.
 *
 *  @return 0=error, 1=pass.
 */
int trans_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 0; A[0][1] = 2; A[0][2] = 8;
	A[1][0] = 1; A[1][1] = -1; A[1][2] = 0;
	A[2][0] = 0; A[2][1] = 1; A[2][2] = 0;
	
	double **C = m_create(n,n);
	C[0][0] = 0; C[0][1] = 1; C[0][2] = 0;
	C[1][0] = 2; C[1][1] = -1; C[1][2] = 1;
	C[2][0] = 8; C[2][1] = 0; C[2][2] = 0;
    
	double **R = m_trans(A,n);
    _assert(equals_matrix(R,C,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(C,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function R_x.
 *
 *  @return 0=error, 1=pass.
 */
int R_x_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 1; A[0][1] = 0; A[0][2] = 0;
	A[1][0] = 0; A[1][1] = 0.54030230586814; A[1][2] = 0.841470984807897;
	A[2][0] = 0; A[2][1] = -0.841470984807897; A[2][2] = 0.54030230586814;
    
	double **R = R_x(1.0);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function R_y.
 *
 *  @return 0=error, 1=pass.
 */
int R_y_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 0.54030230586814; A[0][1] = 0; A[0][2] = -0.841470984807897;
	A[1][0] = 0; A[1][1] = 1; A[1][2] = 0;
	A[2][0] = 0.841470984807897; A[2][1] = 0; A[2][2] = 0.54030230586814;
    
	double **R = R_y(1.0);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function R_z.
 *
 *  @return 0=error, 1=pass.
 */
int R_z_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 0.54030230586814; A[0][1] = 0.841470984807897; A[0][2] = 0;
	A[1][0] = -0.841470984807897; A[1][1] = 0.54030230586814; A[1][2] = 0;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 1;
    
	double **R = R_z(1.0);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function position.
 *
 *  @return 0=error, 1=pass.
 */
int position_01() {
    double lon = -2.76234307910694;
    double lat = 0.376551295459273;
    double alt = 300.2;
	
    double *sol = v_create(3);
    sol[0] = -5512567.8400360728;
    sol[1] = -2196994.4466693182;
    sol[2] = 2330804.9661468901;
    
    double *R = position(lon, lat, alt);
    _assert(equals_vector(R,sol,3,1e-10));

    v_free(R,3);
    
    return 0;
}

/** @brief Unit test for function Mjday_TDB.
 *
 *  @return 0=error, 1=pass.
 */
int Mjday_TDB_01() {
    double sol = 49746.1097222304;
    
	double R = Mjday_TDB(49746.1097222222);
    _assert(sol - R < DBL_EPSILON);
	
    return 0;
}

/** @brief Unit test for function MeanObliquity.
 *
 *  @return 0=error, 1=pass.
 */
int MeanObliquity_01() {
    double sol = 0.409103979370901;
    
	double R = MeanObliquity(49746.1097222222);
    _assert(sol - R < DBL_EPSILON);
	
    return 0;
}

/** @brief Unit test for function Mjday.
 *
 *  @return 0=error, 1=pass.
 */
int Mjday_01() {
    double sol = 59333;
    
	double R = Mjday(2021,4,29,0,0,0);
    _assert(sol - R < DBL_EPSILON);
	
    return 0;
}

/** @brief Unit test for function NutAngles.
 *
 *  @return 0=error, 1=pass.
 */
int NutAngles_01() {
    double dpsi_sol = 6.23063736216799e-05;
    double deps_sol = -3.51110708894389e-05;
    
	double dpsi, deps;
	NutAngles(49746.1097222222, &dpsi, &deps);
    _assert(fabs(dpsi_sol - dpsi) < DBL_EPSILON && fabs(deps_sol - deps) < DBL_EPSILON);
	
    return 0;
}

/** @brief Unit test for function NutMatrix.
 *
 *  @return 0=error, 1=pass.
 */
int NutMatrix_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 0.9999999980589579; A[0][1] = -5.7164703144451515e-05; A[0][2] = -2.4784690905225753e-05;
	A[1][0] = 5.7165573326254723e-05; A[1][1] = 0.99999999774967985; A[1][2] = 3.511036246666723e-05;
	A[2][0] = 2.4782683776004558e-05; A[2][1] = -3.511177922960697e-05; A[2][2] = 0.99999999907649073;
    
	double **R = NutMatrix(49746.1097222222);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function AzElPa.
 *
 *  @return 0=error, 1=pass.
 */
int AzElPa_01() {
    int n = 3;
	
	double *s = v_create(n), Az, El, *dAds, *dEds, Az_sol = 1.05892995381517,
			    El_sol = 0.286534142298292, *dAds_sol = v_create(n), *dEds_sol = v_create(n);
	s[0] = 2159055.44810213; s[1] = 1212982.41102083; s[2] = 729669.130080208;
	dAds_sol[0] = 1.97784562210396e-07; dAds_sol[1] = -3.52047839037887e-07; dAds_sol[2] = 0.0;
	dEds_sol[0] = -9.54424040207182e-08; dEds_sol[1] = -5.36206503841483e-08; dEds_sol[2] = 3.71546961476313e-07;
    
	AzElPa(s, &Az, &El, &dAds, &dEds);
    _assert(fabs(Az_sol - Az) < 1e-10 &&
    		fabs(El_sol - El) < 1e-10 &&
			equals_vector(dAds_sol,dAds,n,1e-10) &&
			equals_vector(dEds_sol,dEds,n,1e-10));
    
	v_free(s,n);
	v_free(dAds,n);
	v_free(dEds,n);
	v_free(dAds_sol,n);
	v_free(dEds_sol,n);
	
    return 0;
}

/** @brief Unit test for function PoleMatrix.
 *
 *  @return 0=error, 1=pass.
 */
int PoleMatrix_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 0.955336489125606; A[0][1] = 0.141679934247038; A[0][2] = 0.259343380052231;
	A[1][0] = 0.0; A[1][1] = 0.877582561890373; A[1][2] = -0.479425538604203;
	A[2][0] = -0.29552020666134; A[2][1] = 0.458012710847292; A[2][2] = 0.838386643594204;
    
	double **R = PoleMatrix(0.3,0.5);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function PrecMatrix.
 *
 *  @return 0=error, 1=pass.
 */
int PrecMatrix_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 0.999999997546639; A[0][1] = 6.42297715142866e-05; A[0][2] = 2.79510023652174e-05;
	A[1][0] = -6.42297715142866e-05; A[1][1] = 0.999999997937268; A[1][2] = -8.97642805260184e-10;
	A[2][0] = -2.79510023652174e-05; A[2][1] = -8.97643692455259e-10; A[2][2] = 0.999999999609371;
    
	double **R = PrecMatrix(526,421);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function Frac.
 *
 *  @return 0=error, 1=pass.
 */
int Frac_01() {
    double sol = 0.7;
    
	double R = Frac(-2.3);
    _assert(fabs(sol - R) < 1e-15);
	
    return 0;
}

/** @brief Unit test for function gmst.
 *
 *  @return 0=error, 1=pass.
 */
int gmst_01() {
    double sol = 0.554091424936741;
    
	double R = gmst(523.5);
    _assert(fabs(sol - R) < 1e-15);
	
    return 0;
}

/** @brief Unit test for function timediff.
 *
 *  @return 0=error, 1=pass.
 */
int timediff_01() {
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC,
		   UT1_TAI_sol = -28.6742353018935,
		   UTC_GPS_sol = -10,
		   UT1_GPS_sol = -9.67423530189348,
		   TT_UTC_sol  = 61.184,
		   GPS_UTC_sol = 10;
    
	timediff(0.325764698106523,29, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    _assert(fabs(UT1_TAI_sol - UT1_TAI) < 1e-10 &&
			fabs(UTC_GPS_sol - UTC_GPS) < 1e-10 &&
			fabs(UT1_GPS_sol - UT1_GPS) < 1e-10 &&
			fabs(TT_UTC_sol - TT_UTC) < 1e-10 &&
			fabs(GPS_UTC_sol - GPS_UTC) < 1e-10);
	
    return 0;
}

/** @brief Unit test for function unit.
 *
 *  @return 0=error, 1=pass.
 */
int unit_01() {
    int n = 3;
	
	double *v = v_create(n);
	v[0] = 5; v[1] = 4; v[2] = 3;
	double *v_sol = v_create(n);
	v_sol[0] = 0.707106781186547; v_sol[1] = 0.565685424949238; v_sol[2] = 0.424264068711929;
    
	double *r = unit(v);
    _assert(equals_vector(v_sol,r,n,1e-10));
    
	v_free(v,n);
	v_free(v_sol,n);
	v_free(r,n);
	
    return 0;
}

/** @brief Unit test for function sign_.
 *
 *  @return 0=error, 1=pass.
 */
int sign_01() {
    double a = 5, b =-3, sol = -5;
    
	double r = sign_(a,b);
    _assert(fabs(sol - r) < 1e-15);
	
    return 0;
}

/** @brief Unit test for function EqnEquinox.
 *
 *  @return 0=error, 1=pass.
 */
int EqnEquinox_01() {
    double a = 5.6, sol = 2.65688073537765e-05;
    
	double r = EqnEquinox(a);
    _assert(fabs(sol - r) < 1e-15);
	
    return 0;
}

/** @brief Unit test for function gast.
 *
 *  @return 0=error, 1=pass.
 */
int gast_01() {
    double a = 5.6, sol = 4.83948153246421;
    
	double r = gast(a);
    _assert(fabs(sol - r) < 1e-10);
	
    return 0;
}

/** @brief Unit test for function Geodetic.
 *
 *  @return 0=error, 1=pass.
 */
int Geodetic_01() {
    int n = 3;
	
	double *v = v_create(n);
	v[0] = 5; v[1] = 4; v[2] = 3;
	double lon_sol = 0.674740942223553,
		   lat_sol = 1.57064687580039,
		   h_sol = -6356748.61611367;
    
	double lon, lat, h;
	Geodetic(v, &lon, &lat, &h);
    _assert(fabs(lon_sol - lon) < 1e-10 &&
			fabs(lat_sol - lat) < 1e-8 &&
			fabs(h_sol - h) < 1e-2);
    
	v_free(v,n);
	
    return 0;
}

/** @brief Unit test for function Legendre.
 *
 *  @return 0=error, 1=pass.
 */
int Legendre_01() {
    int n = 3;
	
	double **pnm_sol = m_create(n,n);
	pnm_sol[0][0] = 1.0; pnm_sol[0][1] = 0.0; pnm_sol[0][2] = 0.0;
	pnm_sol[1][0] = 0.830389391308554; pnm_sol[1][1] = 1.52001758503058; pnm_sol[1][2] = 0.0;
	pnm_sol[2][0] = -0.347097518865836; pnm_sol[2][1] = 1.62950155523887; pnm_sol[2][2] = 1.49139129468805;
	
	double **dpnm_sol = m_create(n,n);
	dpnm_sol[0][0] = 0.0; dpnm_sol[0][1] = 0.0; dpnm_sol[0][2] = 0.0;
	dpnm_sol[1][0] = 1.52001758503058; dpnm_sol[1][1] = -0.830389391308554; dpnm_sol[1][2] = 0.0;
	dpnm_sol[2][0] = 2.82237948468622; dpnm_sol[2][1] = 2.09258183254477; dpnm_sol[2][2] = -1.62950155523887;
    
	double **pnm, **dpnm;
	Legendre(2,2,0.5,&pnm,&dpnm);
    _assert(equals_matrix(pnm_sol,pnm,n,n,1e-10) &&
    		equals_matrix(dpnm_sol,dpnm,n,n,1e-10));
    
	m_free(pnm_sol,n,n);
	m_free(dpnm_sol,n,n);
	m_free(pnm,n,n);
	m_free(dpnm,n,n);
	
    return 0;
}

/** @brief Unit test for function LTC.
 *
 *  @return 0=error, 1=pass.
 */
int LTC_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = 0.370223471399199; A[0][1] = -0.928942722252092; A[0][2] = 0.0;
	A[1][0] = 0.341586711932422; A[1][1] = 0.136136938528208; A[1][2] = 0.929938305587722;
	A[2][0] = -0.863859421119156; A[2][1] = -0.344284987681776; A[2][2] = 0.367715580035218;
    
	double **R = LTC(-2.76234307910694,0.376551295459273);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function GHAMatrix.
 *
 *  @return 0=error, 1=pass.
 */
int GHAMatrix_01() {
    int n = 3;
	
	double **A = m_create(n,n);
	A[0][0] = -0.555497751197422; A[0][1] = -0.831518038538315; A[0][2] = 0.0;
	A[1][0] = 0.831518038538315; A[1][1] = -0.555497751197422; A[1][2] = 0.0;
	A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 1.0;
    
	double **R = GHAMatrix(0.5);
    _assert(equals_matrix(A,R,n,n,1e-10));
    
	m_free(A,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function EccAnom.
 *
 *  @return 0=error, 1=pass.
 */
int EccAnom_01() {
    double M = 0.5, e = 1, sol = 1.49730038909589;
    
	double r = EccAnom(M, e);
    _assert(fabs(sol - r) < 1e-10);
	
    return 0;
}

/** @brief Unit test for function Cheb3D.
 *
 *  @return 0=error, 1=pass.
 */
int Cheb3D_01() {
    int n = 3;
	
	double *sol = v_create(n);
	sol[0] = 3.26; sol[1] = 4.66; sol[2] = 1.1;
	
	double *Cx = v_create(n);
	Cx[0] = 5.0; Cx[1] = 4.0; Cx[2] = 3.0;
	double *Cy = v_create(n);
	Cy[0] = 6.0; Cy[1] = 7.0; Cy[2] = 8.0;
	double *Cz = v_create(n);
	Cz[0] = 2.0; Cz[1] = 1.0; Cz[2] = 0.0;
    
	double *R = Cheb3D(0.5,n,0,10,Cx,Cy,Cz);
    _assert(equals_vector(sol,R,n,1e-10));
    
	v_free(sol,n);
	v_free(Cx,n);
	v_free(Cy,n);
	v_free(Cz,n);
	v_free(R,n);
	
    return 0;
}

/** @brief Unit test for function elements.
 *
 *  @return 0=error, 1=pass.
 */
int elements_01() {
    int n = 6;
	
	double *y = v_create(n);
	y[0] = 1; y[1] = 2; y[2] = 3; y[3] = 4; y[4] = 5; y[5] = 6;
	
	double p_sol = 1.35474011564823e-13,
		   a_sol = 1.87082869338765,
		   e_sol = 0.999999999999964,
		   i_sol = 1.99133066207886,
		   Omega_sol = 3.6052402625906,
		   omega_sol = 5.21086941752228,
		   M_sol = 3.14159030993265;
    
	double p, a, e, i, Omega, omega, M;
	elements(y, &p, &a, &e, &i, &Omega, &omega, &M);
    _assert(fabs(p_sol - p) < 1e-10 &&
			fabs(a_sol - a) < 1e-10 &&
			fabs(e_sol - e) < 1e-10 &&
			fabs(i_sol - i) < 1e-10 &&
			fabs(Omega_sol - Omega) < 1e-10 &&
			fabs(omega_sol - omega) < 1e-10 &&
			fabs(M_sol - M) < 1e-10);
    
	v_free(y,n);
	
    return 0;
}

/** @brief Unit test for function IERS.
 *
 *  @return 0=error, 1=pass.
 */
int IERS_01() {
    double Mjd_UTC = 49746.1101504629;
	char interp = 'l';
	
	double x_pole_sol = -5.59386183152189e-07,
		   y_pole_sol = 2.33554438440373e-06,
		   UT1_UTC_sol = 0.325764698106523,
		   LOD_sol = 0.00272668635763815,
		   dpsi_sol = -1.16881960640421e-07,
		   deps_sol = -2.47881680412219e-08,
		   dx_pole_sol = -8.41764150670523e-10,
		   dy_pole_sol = -1.56618880121342e-09,
		   TAI_UTC_sol = 29;
    
	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	IERS(Mjd_UTC, interp, &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi,
		 &deps, &dx_pole, &dy_pole, &TAI_UTC);
	_assert(fabs(x_pole_sol - x_pole) < 1e-10 &&
			fabs(y_pole_sol - y_pole) < 1e-10 &&
			fabs(UT1_UTC_sol - UT1_UTC) < 1e-10 &&
			fabs(LOD_sol - LOD) < 1e-10 &&
			fabs(dpsi_sol - dpsi) < 1e-10 &&
			fabs(deps_sol - deps) < 1e-10 &&
			fabs(dx_pole_sol - dx_pole) < 1e-10 &&
			fabs(dy_pole_sol - dy_pole) < 1e-10 &&
			fabs(TAI_UTC_sol - TAI_UTC) < 1e-10);
	
    return 0;
}

/** @brief Unit test for function angl.
 *
 *  @return 0=error, 1=pass.
 */
int angl_01() {
    int n = 3;
	
	double *v = v_create(n);
	v[0] = 5.0; v[1] = 4.0; v[2] = 3.0;
	double *w = v_create(n);
	w[0] = 2.0; w[1] = 1.0; w[2] = 0.0;
    
	double R = angl(v,w), R_sol = 0.483361282211952;
    _assert(fabs(R_sol - R) < 1e-10);
    
	v_free(v,n);
	v_free(w,n);
	
    return 0;
}

/** @brief Unit test for function gibbs.
 *
 *  @return 0=error, 1=pass.
 */
int gibbs_01() {
    int n = 3;
	
	double *v = v_create(n);
	v[0] = 5.0; v[1] = -4.0; v[2] = 3.0;
	double *w = v_create(n);
	w[0] = -2.0; w[1] = 1.0; w[2] = 0.0;
	double *x = v_create(n);
	x[0] = -1.0; x[1] = -2.0; x[2] = 3.0;
    
	double *v2, theta, theta1, copa;
	char *error;
	gibbs(v,w,x,&v2,&theta,&theta1,&copa,&error);
	
	double *v2_sol = v_create(n);
	v2_sol[0] = -5585279.86512667; v2_sol[1] = -9245958.76630678; v2_sol[2] = 12496183.5615176;
	double theta_sol = 2.65823137137784,
		   theta1_sol = 1.5707963267949,
		   copa_sol = 0.101593180557725;
		   
	_assert(equals_vector(v2_sol,v2,n,1e-7) &&
			fabs(theta_sol - theta) < 1e-10 &&
			fabs(theta1_sol - theta1) < 1e-10 &&
			fabs(copa_sol - copa) < 1e-10);
	
    
	v_free(v,n);
	v_free(w,n);
	v_free(x,n);
	v_free(v2,n);
	v_free(v2_sol,n);
	
    return 0;
}

/** @brief Unit test for function hgibbs.
 *
 *  @return 0=error, 1=pass.
 */
int hgibbs_01() {
    int n = 3;
	
	double *v = v_create(n);
	v[0] = 5.0; v[1] = -4.0; v[2] = 3.0;
	double *w = v_create(n);
	w[0] = -2.0; w[1] = 1.0; w[2] = 0.0;
	double *x = v_create(n);
	x[0] = -1.0; x[1] = -2.0; x[2] = 3.0;
    
	double Mjd1 = 55,
		   Mjd2 = 44,
		   Mjd3 = 33,
		   *v2, theta, theta1, copa;
	char *error;
	hgibbs(v,w,x,Mjd1,Mjd2,Mjd3,&v2,&theta,&theta1,&copa,&error);
	
	double *v2_sol = v_create(n);
	v2_sol[0] = 1049113223860432400.0; v2_sol[1] = 848151707750756740.0; v2_sol[2] = -1540100720191158300.0;
	double theta_sol = 2.65823137137784,
		   theta1_sol = 1.5707963267949,
		   copa_sol = 0.101593180557725;
	
	_assert(equals_vector(v2_sol,v2,n,1e-1) &&
			fabs(theta_sol - theta) < 1e-10 &&
			fabs(theta1_sol - theta1) < 1e-10 &&
			fabs(copa_sol - copa) < 1e-10);
	
    
	v_free(v,n);
	v_free(w,n);
	v_free(x,n);
	v_free(v2,n);
	v_free(v2_sol,n);
	
    return 0;
}

/** @brief Unit test for function TimeUpdate.
 *
 *  @return 0=error, 1=pass.
 */
int TimeUpdate_01() {
    int n = 6;
	
	double **P = m_create(n,n);
	P[0][0] = 100000000; P[0][1] = 0; P[0][2] = 0; P[0][3] = 0; P[0][4] = 0; P[0][5] = 0;
	P[1][0] = 0; P[1][1] = 100000000; P[1][2] = 0; P[1][3] = 0; P[1][4] = 0; P[1][5] = 0;
	P[2][0] = 0; P[2][1] = 0; P[2][2] = 100000000; P[2][3] = 0; P[2][4] = 0; P[2][5] = 0;
	P[3][0] = 0; P[3][1] = 0; P[3][2] = 0; P[3][3] = 1000; P[3][4] = 0; P[3][5] = 0;
	P[4][0] = 0; P[4][1] = 0; P[4][2] = 0; P[4][3] = 0; P[4][4] = 1000; P[4][5] = 0;
	P[5][0] = 0; P[5][1] = 0; P[5][2] = 0; P[5][3] = 0; P[5][4] = 0; P[5][5] = 1000;
	double **Phi = m_create(n,n);
	Phi[0][0] = 1.00041922218367; Phi[0][1] = 0.000599210758843754; Phi[0][2] = 0.000737344012173523; Phi[0][3] = 37.0053480735614; Phi[0][4] = 0.00742082744985318; Phi[0][5] = 0.00907132659757866;
	Phi[1][0] = 0.000599205935189867; Phi[1][1] = 0.999703770692677; Phi[1][2] = 0.000418689926148395; Phi[1][3] = 0.00742079763278345; Phi[1][4] = 36.9963058593961; Phi[1][5] = 0.00509713273998443;
	Phi[2][0] = 0.000737334362385235; Phi[2][1] = 0.000418687815717319; Phi[2][2] = 0.999877413279163; Phi[2][3] = 0.0090712669637158; Phi[2][4] = 0.00509711970560663; Phi[2][5] = 36.9983505101528;
	Phi[3][0] = 2.34511597329531e-05; Phi[3][1] = 3.25246538718763e-05; Phi[3][2] = 3.97560066025033e-05; Phi[3][3] = 1.00044810751498; Phi[3][4] = 0.000604066118447109; Phi[3][5] = 0.000733457411125333;
	Phi[4][0] = 3.25240005603744e-05; Phi[4][1] = -1.61854966340218e-05; Phi[4][2] = 2.23408858389155e-05; Phi[4][3] = 0.000604061271608192; Phi[4][4] = 0.999697158529404; Phi[4][5] = 0.000407832040042641;
	Phi[5][0] = 3.97546995847966e-05; Phi[5][1] = 2.23405999335799e-05; Phi[5][2] = -7.22167527414373e-06; Phi[5][3] = 0.000733447714890626; Phi[5][4] = 0.000407829919202528; Phi[5][5] = 0.999855141838431;
	double **Qdt = m_zeros(n,n);
    
	
	double **R_sol = m_create(n,n);
	R_sol[0][0] = 101453348.207834; R_sol[0][1] = 120429.109556826; R_sol[0][2] = 148186.1448513; R_sol[0][3] = 39372.9209797587; R_sol[0][4] = 3284.21675106589; R_sol[0][5] = 4014.15727751921;
	R_sol[1][0] = 120429.109556826; R_sol[1][1] = 101309543.076907; R_sol[1][2] = 84141.6477758924; R_sol[1][3] = 3284.34773933075; R_sol[1][4] = 35369.9224513583; R_sol[1][5] = 2255.66799781441;
	R_sol[2][0] = 148186.1448513; R_sol[2][1] = 84141.6477758924; R_sol[2][2] = 101344434.103469; R_sol[2][3] = 4014.41933186659; R_sol[2][4] = 2255.72532205054; R_sol[2][5] = 36274.7873542153;
	R_sol[3][0] = 39372.9209797587; R_sol[3][1] = 3284.34773933075; R_sol[3][2] = 4014.41933186659; R_sol[3][3] = 1001.21615369228; R_sol[3][4] = 1.320962491756; R_sol[3][5] = 1.6045548104278;
	R_sol[4][0] = 3284.21675106589; R_sol[4][1] = 35369.9224513583; R_sol[4][2] = 2255.72532205054; R_sol[4][3] = 1.320962491756; R_sol[4][4] = 999.576829598137; R_sol[4][5] = 0.892927375360559;
	R_sol[5][0] = 4014.15727751921; R_sol[5][1] = 2255.66799781441; R_sol[5][2] = 36274.7873542153; R_sol[5][3] = 1.6045548104278; R_sol[5][4] = 0.892927375360559; R_sol[5][5] = 999.924178045366;
	double **R = TimeUpdate(P,Phi,Qdt);
    _assert(equals_matrix(R_sol,R,n,n,1e-5));
    
	m_free(P,n,n);
	m_free(Phi,n,n);
	m_free(Qdt,n,n);
	m_free(R_sol,n,n);
	m_free(R,n,n);
	
    return 0;
}

/** @brief Unit test for function MeasUpdate.
 *
 *  @return 0=error, 1=pass.
 */
int MeasUpdate_01() {
    int n = 6;
	
	double z = 1.0559084894933,
		   g = 1.05892995381517,
		   s = 0.00039095375244673;
	double *x = v_create(n);
	x[0] = 5738566.57769186; x[1] = 3123975.34092959; x[2] = 3727114.48156055; x[3] = 5199.63329181068; x[4] = -2474.43881044643; x[5] = -7195.16752553801;
	double *G = v_create(n);
	G[0] = 9.59123748603008e-08; G[1] = 2.16050345227537e-07; G[2] = -3.27382770920712e-07; G[3] = 0; G[4] = 0; G[5] = 0;
	double **P = m_create(n,n);
	P[0][0] = 101453348.207834; P[0][1] = 120429.109556826; P[0][2] = 148186.1448513; P[0][3] = 39372.9209797587; P[0][4] = 3284.21675106589; P[0][5] = 4014.15727751921;
	P[1][0] = 120429.109556826; P[1][1] = 101309543.076907; P[1][2] = 84141.6477758924; P[1][3] = 3284.34773933075; P[1][4] = 35369.9224513583; P[1][5] = 2255.66799781441;
	P[2][0] = 148186.1448513; P[2][1] = 84141.6477758924; P[2][2] = 101344434.103469; P[2][3] = 4014.41933186659; P[2][4] = 2255.72532205054; P[2][5] = 36274.7873542153;
	P[3][0] = 39372.9209797587; P[3][1] = 3284.34773933075; P[3][2] = 4014.41933186659; P[3][3] = 1001.21615369228; P[3][4] = 1.320962491756; P[3][5] = 1.6045548104278;
	P[4][0] = 3284.21675106589; P[4][1] = 35369.9224513583; P[4][2] = 2255.72532205054; P[4][3] = 1.320962491756; P[4][4] = 999.576829598137; P[4][5] = 0.892927375360559;
	P[5][0] = 4014.15727751921; P[5][1] = 2255.66799781441; P[5][2] = 36274.7873542153; P[5][3] = 1.6045548104278; P[5][4] = 0.892927375360559; P[5][5] = 999.924178045366;
	
	double *K;
	MeasUpdate(z,g,s,G,&K,&x,&P);
	
	double *K_sol = v_create(n);
	K_sol[0] = 582691.20646825; K_sol[1] = 1312775.3084207; K_sol[2] = -1989454.89979193; K_sol[3] = 190.367307502019; K_sol[4] = 433.242659522805; K_sol[5] = -660.433799143302;
	double *x_sol = v_create(n);
	x_sol[0] = 5736805.99700085; x_sol[1] = 3120008.83717257; x_sol[2] = 3733125.54856023; x_sol[3] = 5199.05810378301; x_sol[4] = -2475.74783768488; x_sol[5] = -7193.17204837694;
	double **P_sol = m_create(n,n);
	P_sol[0][0] = 95796502.3074957; P_sol[0][1] = -12624173.0726456; P_sol[0][2] = 19462086.3185101; P_sol[0][3] = 37524.8091307257; P_sol[0][4] = -921.762222304672; P_sol[0][5] = 10425.7388968351;
	P_sol[1][0] = -12624173.0726456; P_sol[1][1] = 72596566.3325921; P_sol[1][2] = 43597431.3213228; P_sol[1][3] = -879.359513633289; P_sol[1][4] = 25894.053787673; P_sol[1][5] = 16700.6535138616;
	P_sol[2][0] = 19462086.3185101; P_sol[2][1] = 43597431.3213228; P_sol[2][2] = 35401902.4365326; P_sol[2][3] = 10324.3398053662; P_sol[2][4] = 16615.9994846053; P_sol[2][5] = 14384.0288763645;
	P_sol[3][0] = 37524.8091307257; P_sol[3][1] = -879.35951363329; P_sol[3][2] = 10324.3398053662; P_sol[3][3] = 1000.61236892128; P_sol[3][4] = -0.0531459273905025; P_sol[3][5] = 3.6992415263927;
	P_sol[4][0] = -921.762222304671; P_sol[4][1] = 25894.053787673; P_sol[4][2] = 16615.9994846053; P_sol[4][3] = -0.0531459273905022; P_sol[4][4] = 996.449599435587; P_sol[4][5] = 5.66006757185727;
	P_sol[5][0] = 10425.7388968351; P_sol[5][1] = 16700.6535138616; P_sol[5][2] = 14384.0288763645; P_sol[5][3] = 3.6992415263927; P_sol[5][4] = 5.66006757185727; P_sol[5][5] = 992.657163955644;
	
	_assert(equals_vector(K_sol,K,n,1e-5) &&
			equals_vector(x_sol,x,n,1e-5) &&
			equals_matrix(P_sol,P,n,n,1e-5));
	
    v_free(x,n);
    v_free(G,n);
    v_free(K,n);
	m_free(P,n,n);
    v_free(K_sol,n);
    v_free(x_sol,n);
	m_free(P_sol,n,n);
	
    return 0;
}

/** @brief Unit test for function AccelPointMass.
 *
 *  @return 0=error, 1=pass.
 */
int AccelPointMass_01() {
    int n = 3;
	
	double *r = v_create(n);
	r[0] = 6221397.62857869; r[1] = 2867713.77965741; r[2] = 3006155.9850995;
	double *s = v_create(n);
	s[0] = 92298251728.4766; s[1] = -105375196079.054; s[2] = -45686367226.3533;
    
	double GM = 1.32712440041939e+20;
	double *a = AccelPointMass(r,s,GM);
	
	double *a_sol = v_create(n);
	a_sol[0] = -1.8685505934417e-07; a_sol[1] = -2.00332995883182e-07; a_sol[2] = -1.59993120755489e-07;
	_assert(equals_vector(a_sol,a,n,1e-17));
	
    
	v_free(r,n);
	v_free(s,n);
	v_free(a,n);
	v_free(a_sol,n);
	
    return 0;
}

/** @brief Unit test for function AccelHarmonic.
 *
 *  @return 0=error, 1=pass.
 */
int AccelHarmonic_01() {
    int n = 3;
	
	double *r = v_create(n);
	r[0] = 6221397.62857869; r[1] = 2867713.77965741; r[2] = 3006155.9850995;
	
	double **E = m_create(n,n);
	E[0][0] = -0.978185453896254; E[0][1] = 0.20773306636226; E[0][2] = -0.000436950239569363;
	E[1][0] = -0.207733028352522; E[1][1] = -0.978185550768511; E[1][2] = -0.000131145697267082;
	E[2][0] = -0.000454661708585098; E[2][1] = -3.75158169026289e-05; E[2][2] = 0.999999895937642;
    
	double n_max = 20,
		   m_max = 20;
	double *a = AccelHarmonic(r,E,n_max,m_max);
	
	double *a_sol = v_create(n);
	a_sol[0] = -5.92414856522537; a_sol[1] = -2.73076679296887; a_sol[2] = -2.86933544780686;
	
	_assert(equals_vector(a_sol,a,n,1e-10));
    
	
	v_free(r,n);
	m_free(E,n,n);
	v_free(a,n);
	v_free(a_sol,n);
	
    return 0;
}

/** @brief Unit test for function JPL_Eph_DE430.
 *
 *  @return 0=error, 1=pass.
 */
int JPL_Eph_DE430_01() {
    int n = 3;
	
	double Mjd_TDB = 49746.1119928785;
	
	double *r_Mercury, *r_Venus, *r_Earth, *r_Mars, *r_Jupiter, *r_Saturn,
		   *r_Uranus, *r_Neptune, *r_Pluto, *r_Moon, *r_Sun;
	JPL_Eph_DE430(Mjd_TDB, &r_Mercury, &r_Venus,
				   &r_Earth, &r_Mars, &r_Jupiter,
				   &r_Saturn, &r_Uranus, &r_Neptune,
				   &r_Pluto, &r_Moon, &r_Sun);
	
	double *r_Moon_sol = v_create(n);
	r_Moon_sol[0] = 89383372.3127192; r_Moon_sol[1] = -336603832.117946; r_Moon_sol[2] = -114648787.751995;
	double *r_Earth_sol = v_create(n);
	r_Earth_sol[0] = -92470961229.1923; r_Earth_sol[1] = 106394918389.493; r_Earth_sol[2] = 46130139909.4037;
	double *r_Mercury_sol = v_create(n);
	r_Mercury_sol[0] = 83775495895.6957; r_Mercury_sol[1] = -65291124913.4462; r_Mercury_sol[2] = -23391312101.2408;
	double *r_Venus_sol = v_create(n);
	r_Venus_sol[0] = -15229665573.9533; r_Venus_sol[1] = -110134992637.563; r_Venus_sol[2] = -41021803625.6266;
	double *r_Mars_sol = v_create(n);
	r_Mars_sol[0] = -88278413008.2432; r_Mars_sol[1] = 46964769778.2984; r_Mars_sol[2] = 29071026502.6713;
	double *r_Jupiter_sol = v_create(n);
	r_Jupiter_sol[0] = -298385936466.094; r_Jupiter_sol[1] = -754498258910.729; r_Jupiter_sol[2] = -314410518568.228;
	double *r_Saturn_sol = v_create(n);
	r_Saturn_sol[0] = 1482033999505.12; r_Saturn_sol[1] = -453872894236.397; r_Saturn_sol[2] = -249402247811.211;
	double *r_Uranus_sol = v_create(n);
	r_Uranus_sol[0] = 1412367984017.93; r_Uranus_sol[1] = -2511355045786.9; r_Uranus_sol[2] = -1118108651902.6;
	double *r_Neptune_sol = v_create(n);
	r_Neptune_sol[0] = 1871250770052.33; r_Neptune_sol[1] = -3928976313605.4; r_Neptune_sol[2] = -1655020476718.54;
	double *r_Pluto_sol = v_create(n);
	r_Pluto_sol[0] = -2171414794259.54; r_Pluto_sol[1] = -3915433128334.98; r_Pluto_sol[2] = -552716250355.845;
	double *r_Sun_sol = v_create(n);
	r_Sun_sol[0] = 92298251728.4766; r_Sun_sol[1] = -105375196079.054; r_Sun_sol[2] = -45686367226.3533;
	_assert(equals_vector(r_Moon_sol,r_Moon,n,1e-1) &&
			equals_vector(r_Earth_sol,r_Earth,n,1e-1) &&
			equals_vector(r_Mercury_sol,r_Mercury,n,1e-1) &&
			equals_vector(r_Venus_sol,r_Venus,n,1e-1) &&
			equals_vector(r_Mars_sol,r_Mars,n,1e-1) &&
			equals_vector(r_Jupiter_sol,r_Jupiter,n,1e-1) &&
			equals_vector(r_Saturn_sol,r_Saturn,n,1e-1) &&
			equals_vector(r_Uranus_sol,r_Uranus,n,1e-1) &&
			equals_vector(r_Neptune_sol,r_Neptune,n,1e-1) &&
			equals_vector(r_Pluto_sol,r_Pluto,n,1e-1) &&
			equals_vector(r_Sun_sol,r_Sun,n,1e-1));
    
	
	v_free(r_Moon,n);
	v_free(r_Earth,n);
	v_free(r_Mercury,n);
	v_free(r_Venus,n);
	v_free(r_Mars,n);
	v_free(r_Jupiter,n);
	v_free(r_Saturn,n);
	v_free(r_Uranus,n);
	v_free(r_Neptune,n);
	v_free(r_Pluto,n);
	v_free(r_Sun,n);
	v_free(r_Moon_sol,n);
	v_free(r_Earth_sol,n);
	v_free(r_Mercury_sol,n);
	v_free(r_Venus_sol,n);
	v_free(r_Mars_sol,n);
	v_free(r_Jupiter_sol,n);
	v_free(r_Saturn_sol,n);
	v_free(r_Uranus_sol,n);
	v_free(r_Neptune_sol,n);
	v_free(r_Pluto_sol,n);
	v_free(r_Sun_sol,n);
	
    return 0;
}

/** @brief Unit test for function Accel.
 *
 *  @return 0=error, 1=pass.
 */
int Accel_01() {
    int n = 6;
	
	extern Param AuxParam;
	AuxParam.Mjd_UTC = 49746.1112847221;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
	
	double Mjd_TT = 0.0;
	double *Y = v_create(n);
	Y[0] = 6221397.62857869; Y[1] = 2867713.77965741; Y[2] = 3006155.9850995;
	Y[3] = 4645.0472516175; Y[4] = -2752.21591588182; Y[5] = -7507.99940986939;
	
	double *dY;
	Accel(Mjd_TT,Y,&dY);
	
	double *dY_sol = v_create(n);
	dY_sol[0] = 4645.0472516175; dY_sol[1] = -2752.21591588182; dY_sol[2] = -7507.99940986939;
	dY_sol[3] = -5.92414951314002; dY_sol[4] = -2.73076669788114; dY_sol[5] = -2.86933570556258;
	_assert(equals_vector(dY_sol,dY,n,1e-10));
    
	
	v_free(Y,n);
	v_free(dY,n);
	v_free(dY_sol,n);
	
    return 0;
}

/** @brief Unit test for function G_AccelHarmonic.
 *
 *  @return 0=error, 1=pass.
 */
int G_AccelHarmonic_01() {
    int n = 3;
	
	double *r = v_create(n);
	r[0] = 5542555.93722869; r[1] = 3213514.86734919; r[2] = 3990892.97587674;
	
	double **U = m_create(n,n);
	U[0][0] = -0.976675972331716; U[0][1] = 0.214718082511189; U[0][2] = -0.000436019054674645;
	U[1][0] = -0.214718043811152; U[1][1] = -0.976676068937815; U[1][2] = -0.000134261271504216;
	U[2][0] = -0.000454677699074514; U[2][1] = -3.750859940872e-05; U[2][2] = 0.999999895930642;
    
	double n_max = 20,
		   m_max = 20;
	double **G = G_AccelHarmonic(r,U,n_max,m_max);
	
	
	double **G_sol = m_create(n,n);
	G_sol[0][0] = 5.70032034907797e-07; G_sol[0][1] = 8.67651590574781e-07; G_sol[0][2] = 1.08169354007259e-06;
	G_sol[1][0] = 8.67651592351137e-07; G_sol[1][1] = -4.23359107770693e-07; G_sol[1][2] = 6.27183704970946e-07;
	G_sol[2][0] = 1.08169353918441e-06; G_sol[2][1] = 6.27183701418232e-07; G_sol[2][2] = -1.46672925360747e-07;
	_assert(equals_matrix(G_sol,G,n,n,1e-12));
    
	
	v_free(r,n);
	m_free(U,n,n);
	m_free(G,n,n);
	m_free(G_sol,n,n);
	
    return 0;
}

/** @brief Unit test for function VarEqn.
 *
 *  @return 0=error, 1=pass.
 */
int VarEqn_01() {
    int n = 42;
	
	extern Param AuxParam;
	AuxParam.Mjd_UTC = 49746.1101504629;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
	AuxParam.Mjd_TT = 49746.1108586111;
	
	double x = 0.0;
	double *yPhi = v_create(n);
	yPhi[0] = 5542555.93722869; yPhi[1] = 3213514.86734919; yPhi[2] = 3990892.97587674;
	yPhi[3] = 5394.06842166295; yPhi[4] = -2365.21337882319; yPhi[5] = -7061.84554200204;
	yPhi[6] = 1.0; yPhi[13] = 1.0; yPhi[20] = 1.0; yPhi[27] = 1.0; yPhi[34] = 1.0; yPhi[41] = 1.0;
	
	double *yPhip;
	VarEqn(x,yPhi,&yPhip);
	
	double *yPhip_sol = v_create(n);
	yPhip_sol[0] = 5394.06842166295; yPhip_sol[1] = -2365.21337882319; yPhip_sol[2] = -7061.84554200204;
	yPhip_sol[3] = -5.13483678540858; yPhip_sol[4] = -2.97717622353621; yPhip_sol[5] = -3.70591776714193;
	yPhip_sol[9] = 5.70032034907797e-07; yPhip_sol[10] = 8.67651592351137e-07; yPhip_sol[11] = 1.08169353918441e-06;
	yPhip_sol[15] = 8.67651590574781e-07; yPhip_sol[16] = -4.23359107770693e-07; yPhip_sol[17] = 6.27183701418232e-07;
	yPhip_sol[21] = 1.08169354007259e-06; yPhip_sol[22] = 6.27183704970946e-07; yPhip_sol[23] = -1.46672925360747e-07;
	yPhip_sol[24] = 1.0; yPhip_sol[31] = 1.0; yPhip_sol[38] = 1.0;
	_assert(equals_vector(yPhip_sol,yPhip,n,1e-10));
    
	
	v_free(yPhi,n);
	v_free(yPhip,n);
	v_free(yPhip_sol,n);
	
    return 0;
}

/** @brief Unit test for function ode.
 *
 *  @return 0=error, 1=pass.
 */
int ode_01() {
    int n_eqn = 6, iflag = 1, iwork[5];
	double t = 0.0, relerr = 1e-13, abserr = 1e-6, obsMjd0 = -134.999991953373;
	double *work = v_create(100 + 21 * n_eqn);
	
	double *Y = v_create(n_eqn);
	Y[0] = 6221397.62857869; Y[1] = 2867713.77965741; Y[2] = 3006155.9850995;
	Y[3] = 4645.0472516175; Y[4] = -2752.21591588182; Y[5] = -7507.99940986939;
	
	ode(Accel, n_eqn, Y, &t, obsMjd0, relerr, abserr, &iflag, work, iwork);
	
	double *Y_sol = v_create(n_eqn);
	Y_sol[0] = 5542555.93722869; Y_sol[1] = 3213514.8673492; Y_sol[2] = 3990892.97587674;
	Y_sol[3] = 5394.06842166295; Y_sol[4] = -2365.21337882319; Y_sol[5] = -7061.84554200204;
	_assert(equals_vector(Y_sol,Y,n_eqn,1e-1));


	v_free(work,100 + 21 * n_eqn);
	v_free(Y,n_eqn);
	v_free(Y_sol,n_eqn);
	
    return 0;
}

/** @brief Unit test caller.
 *
 *  @return 0=error, 1=pass.
 */
int all_tests() {
	_verify(v_sum_01);
	_verify(v_mul_scalar_01);
    _verify(v_norm_01);
    _verify(v_dot_01);
    _verify(v_cross_01);
	
    _verify(v_extract_01);
    _verify(v_assign_01);
    _verify(v_dot_m_01);
    _verify(m_dot_v_01);
	
    _verify(zeros_01);
    _verify(eye_01);
    _verify(sum_01);
    _verify(m_dot_01);
    _verify(inv_01);
    _verify(trans_01);
	
	_verify(position_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(Mjday_TDB_01);
    _verify(MeanObliquity_01);
    _verify(Mjday_01);
    _verify(NutAngles_01);
    _verify(NutMatrix_01);
	
    _verify(AzElPa_01);
    _verify(PoleMatrix_01);
    _verify(PrecMatrix_01);
    _verify(Frac_01);
    _verify(gmst_01);
    _verify(timediff_01);
    _verify(unit_01);
    _verify(sign_01);
    _verify(EqnEquinox_01);
    _verify(gast_01);
    _verify(Geodetic_01);
    _verify(Legendre_01);
    _verify(LTC_01);
    _verify(GHAMatrix_01);
    _verify(EccAnom_01);
    _verify(Cheb3D_01);
    _verify(elements_01);
    _verify(IERS_01);
    _verify(angl_01);
    _verify(gibbs_01);
    _verify(hgibbs_01);
    _verify(TimeUpdate_01);
    _verify(MeasUpdate_01);
    _verify(AccelPointMass_01);
	_verify(AccelHarmonic_01);
	_verify(JPL_Eph_DE430_01);
	_verify(Accel_01);
	_verify(G_AccelHarmonic_01);
	_verify(VarEqn_01);
	
	_verify(ode_01);

    return 0;
}


//int main()
//{
//	DE430Coeff(2285,1020);
//	GGM03S(182);
//	eop19620101(21413);
//	GEOS3(46);
//
//    int result = all_tests();
//
//    if (result == 0)
//        printf("PASSED\n");
//
//    printf("Tests run: %d\n", tests_run);
//
//    return result != 0;
//}

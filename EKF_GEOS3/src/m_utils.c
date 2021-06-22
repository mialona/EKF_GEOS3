/** @file m_utils.c
 *  @brief Matrix utilities code driver.
 *
 *  This driver contains the code for the matrix 
 *  utilities.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *v_create(int c) {
	double *v = (double*)calloc(c, sizeof(double));
	
    if (v == NULL) {
		printf("vector create: error\n");
        exit(EXIT_FAILURE);
	}
	
	return v;
}

void v_free(double *v, int c) {
	free(v);
}

void v_show(double *v, int c) {
	for (int i = 0; i < c; i++) {
        printf("%2.20f ", v[i]);
        printf("\n");
    }
}

double *v_sum(double *v, int cv, double *w, int cw) {
	double *r = v_create(cv);
	
    if (cv != cw) {
        exit(EXIT_FAILURE);
    }
	
	for(int i = 0; i < cv; i++)
		r[i] = v[i]+w[i];
	
	return r;
}

double *v_mul_scalar(double *v, int c, double s) {
	double *r = v_create(c);
	
	for(int i = 0; i < c; i++)
		r[i] = v[i]*s;
	
	return r;
}

double v_norm(double *v, int c) {
	double r = 0.0;
	
	for (int i = 0; i < c; i++) {
        r += v[i]*v[i];
    }
	
	return (sqrt(r));
}

double v_dot(double *v, int cv, double *w, int cw) {
	double r = 0.0;
	
    if (cv != cw) {
		printf("vector dot: error\n");
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < cv; i++) {
        r += v[i]*w[i];
	}
	
	return r;
}

double *v_cross(double *v, int cv, double *w, int cw) {
	double *r = v_create(cv);
	
    if (cv != 3 || cw != 3) {
		printf("vector cross: error\n");
        exit(EXIT_FAILURE);
	}
	
	r[0] = v[1]*w[2] - v[2]*w[1];
	r[1] = v[2]*w[0] - v[0]*w[2];
	r[2] = v[0]*w[1] - v[1]*w[0];
	
	return r;
}


double **m_create(int f, int c) {
	double **L = (double **) malloc(f*sizeof(double *));
	for(int i = 0; i < f; i++)
		L[i] = (double *) malloc(c*sizeof(double));
	
    if (L == NULL) {
		printf("matrix create: error\n");
        exit(EXIT_FAILURE);
	}
	
	return L;
}

void m_free(double **L, int f, int c) {
	for(int i = 0; i < f; i++)
		free(L[i]);
	
    free(L);
}

void m_show(double **L, int f, int c) {
    for (int i = 0; i < f; i++) {
        for (int j = 0; j < c; j++)
            printf("%2.30f ", L[i][j]);
        printf("\n");
    }
}

double *v_extract(double **A, int fA, int cA, int k) {
	double *v = v_create(fA);
	
    if (k < 0 || k >= cA) {
		printf("vector extract: error\n");
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < fA; i++) {
		v[i] = A[i][k];
	}
	
	return v;
}

double **v_assign(double *v, int cv, double **A, int fA, int cA, int k) {
	double **L = m_create(fA,cA);
	
    if (cv != fA || k < 0 || k >= cA) {
		printf("vector assign: error\n");
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < fA; i++) {
		for(int j = 0; j < cA; j++) {
			if(j!=k) {
				L[i][j] = A[i][j];
			} else {
				L[i][j] = v[i];
			}
		}
	}
	
	return L;
}

double *v_dot_m(double *v, int cv, double **A, int fA, int cA) {
	double *x = v_create(cA);
	
    if (cv != fA) {
		printf("vector dot matrix: error\n");
        exit(EXIT_FAILURE);
	}
	
	for(int j = 0; j < cA; j++) {
		for(int i = 0; i < fA; i++) {
			x[j] += v[i]*A[i][j];
		}
	}
	
	return x;
}

double *m_dot_v(double **A, int fA, int cA, double *v, int cv) {
	double *x = v_create(fA);
	
    if (cv != cA){
		printf("matrix dot vector: error\n");
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < fA; i++) {
		for(int j = 0; j < cA; j++) {
			x[i] += v[j]*A[i][j];
		}
	}
	
	return x;
}

double **m_zeros(int f, int c) {
	double **L = m_create(f,c);
	
	for(int i = 0; i < f; i++) {
        for(int j = 0; j < c; j++) {
			L[i][j] = 0;
		}
	}
	
	return L;
}

double **m_eye(int n) {
	double **L = m_create(n,n);
	
	for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
			if(i == j) {
				L[i][j] = 1;
			} else {
				L[i][j] = 0;
			}
		}
	}
	
	return L;
}

double **m_sum(double **A, int fA, int cA, double **B, int fB, int cB) {
	double **L = m_zeros(fA,cA);
	
    if (fA != fB || cA != cB) {
		printf("matrix sum: error\n");
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < fA; i++) {
        for(int j = 0; j < cA; j++) {
			L[i][j] = A[i][j] + B[i][j];
		}
	}
	
	return L;
}

double **m_dot(double **A, int fA, int cA, double **B, int fB, int cB) {
	double **L = m_zeros(fA,cA);
	
    if (cA != fB) {
		printf("matrix dot: error\n");
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < fA; i++) {
        for(int j = 0; j < cB; j++) {
			for(int k = 0; k < cA; k++)
				L[i][j] += A[i][k]*B[k][j];
		}
	}
	
	return L;
}

double **m_inv(double **A, int n) {
	double **mat, ratio, **result;
    int i, j, k;
    
    mat = m_zeros(n, 2*n);
    result = m_zeros(n, n);
    
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            mat[i][j] = A[i][j];
    
    for(i = 0; i < n; i++){
         for(j = 0; j < n; j++){
              if(i == j){
                  mat[i][j+n] = 1.0;
              } else {
                  mat[i][j+n] = 0.0;
              }
         }
    }
    
    /* Applying Gauss Jordan Elimination */

    for(i = 0; i < n; i++){
         if(mat[i][i] == 0.0){
             printf("inv: error\n");
             exit(EXIT_FAILURE);
         }
         for(j = 0; j < n; j++){
              if(i != j){
                   ratio = mat[j][i]/mat[i][i];
                   for(k = 0; k < 2*n; k++){
                       mat[j][k] = mat[j][k] - ratio*mat[i][k];
                   }
              }
         }
    }

    /* Row Operation to Make Principal Diagonal to 1 */

    for(i = 0; i < n; i++){
         for(j = n; j < 2*n;j++){
             mat[i][j] = mat[i][j]/mat[i][i];
         }
    }

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = mat[i][j+n];
 
    m_free(mat, n, 2*n);
    
    return result;
}

double **m_trans(double **A, int n) {
	double **L = m_create(n, n);
	
	for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
			L[j][i] = A[i][j];
		}
	}
	
	return L;
}

/** @file m_utils.h
 *  @brief Function prototypes for the matrix utilities.
 *
 *  This header file contains the prototypes for the matrix 
 *  utilities code driver.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _MATUTILS_
#define _MATUTILS_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/** @brief Creating a vector of c components.
 *
 *  @param [in] c Number of components.
 *  @return Vector.
 */
double *v_create(int c);

/** @brief Release a vector of c components.
 *
 *  @param [in] v Vector.
 *  @param [in] c Number of components.
 */
void v_free(double *v, int c);

/** @brief Console printing of a vector of c components.
 *
 *  @param [in] v Vector.
 *  @param [in] c Number of components.
 */
void v_show(double *v, int c);

/** @brief Sum of two vector with the same number of components.
 *
 *  @param [in] v First vector.
 *  @param [in] cv Number of components in the first vector.
 *  @param [in] w Second vector.
 *  @param [in] cw Number of components in the second vector.
 *  @return Value of sum.
 */
double *v_sum(double *v, int cv, double *w, int cw);

/** @brief Product of a vector of c components by s scalar.
 *
 *  @param [in] v Vector.
 *  @param [in] cv Number of components.
 *  @param [in] s Scalar.
 *  @return Value of product.
 */
double *v_mul_scalar(double *v, int c, double s);

/** @brief Norm of a vector of c components.
 *
 *  @param [in] v Vector.
 *  @param [in] c Number of components.
 *  @return Value of norm.
 */
double v_norm(double *v, int c);

/** @brief Dot product of two vector with the same number of components.
 *
 *  @param [in] v First vector.
 *  @param [in] cv Number of components in the first vector.
 *  @param [in] w Second vector.
 *  @param [in] cw Number of components in the second vector.
 *  @return Value of product.
 */
double v_dot(double *v, int cv, double *w, int cw);

/** @brief Cross product of two vector with the same number of components.
 *
 *  @param [in] v First vector.
 *  @param [in] cv Number of components in the first vector.
 *  @param [in] w Second vector.
 *  @param [in] cw Number of components in the second vector.
 *  @return Value of product.
 */
double *v_cross(double *v, int cv, double *w, int cw);



/** @brief Extract a column of a matrix to a vector.
 *
 *  @param [in] A Matrix.
 *  @param [in] fA Number of rows in the matrix.
 *  @param [in] cA Number of columns in the matrix.
 *  @param [in] k Column to be extracted.
 *  @return Result vector.
 */
double *v_extract(double **A, int fA, int cA, int k);

/** @brief Assign a vector in a column of a matrix with the same number of components
 *  in the vector as the number of rows in the matrix.
 *
 *  @param [in] v Vector.
 *  @param [in] cv Number of components in the vector.
 *  @param [in] A Matrix.
 *  @param [in] fA Number of rows in the matrix.
 *  @param [in] cA Number of columns in the matrix.
 *  @param [in] k Column to be assigned.
 *  @return Result matrix.
 */
double **v_assign(double *v, int cv, double **A, int fA, int cA, int k);

/** @brief Vector dot product per matrix with the same number of components
 *  in the vector as the number of rows in the matrix.
 *
 *  @param [in] v Vector.
 *  @param [in] cv Number of components in the vector.
 *  @param [in] A Matrix.
 *  @param [in] fA Number of rows in the matrix.
 *  @param [in] cA Number of columns in the matrix.
 *  @return Result vector.
 */
double *v_dot_m(double *v, int cv, double **A, int fA, int cA);

/** @brief Matrix dot product per vector (column) with the same number of components
 *  in the vector as the number of columns in the matrix.
 *
 *  @param [in] A Matrix.
 *  @param [in] fA Number of rows in the matrix.
 *  @param [in] cA Number of columns in the matrix.
 *  @param [in] v Vector.
 *  @param [in] cv Number of components in the vector.
 *  @return Result vector (column).
 */
double *m_dot_v(double **A, int fA, int cA, double *v, int cv);



/** @brief Creating a matrix of f rows and c columns.
 *
 *  @param [in] f Number of rows.
 *  @param [in] c Number of columns.
 *  @return Matrix.
 */
double **m_create(int f, int c);

/** @brief Release a matrix of f rows and c columns.
 *
 *  @param [in] L Matrix.
 *  @param [in] f Number of rows.
 *  @param [in] c Number of columns.
 */
void m_free(double **L, int f, int c);

/** @brief Console printing of a matrix of f rows and c columns.
 *
 *  @param [in] L Matrix.
 *  @param [in] f Number of rows.
 *  @param [in] c Number of columns.
 */
void m_show(double **L, int f, int c);

/** @brief Creating a matrix of zeros of f rows and c columns.
 *
 *  @param [in] f Number of rows.
 *  @param [in] c Number of columns.
 *  @return Matrix.
 */
double **m_zeros(int f, int c);

/** @brief Creation of an identity matrix of order n.
 *
 *  @param [in] n Matrix order.
 *  @return Matrix.
 */
double **m_eye(int n);

/** @brief Sum of two matrices of f rows and c columns.
 *
 *  @param [in] A First matrix.
 *  @param [in] fA Number of rows in the first matrix.
 *  @param [in] cA Number of columns in the first matrix.
 *  @param [in] B Second matrix.
 *  @param [in] fB Number of rows in the second matrix.
 *  @param [in] cB Number of columns in the second matrix.
 *  @return Result matrix.
 */
double **m_sum(double **A, int fA, int cA, double **B, int fB, int cB);

/** @brief Dot product of two matrices with the same number of columns
 *  in the first matrix as the number of rows in the second.
 *
 *  @param [in] A First matrix.
 *  @param [in] fA Number of rows in the first matrix.
 *  @param [in] cA Number of columns in the first matrix.
 *  @param [in] B Second matrix.
 *  @param [in] fB Number of rows in the second matrix ().
 *  @param [in] cB Number of columns in the second matrix.
 *  @return Result matrix.
 */
double **m_dot(double **A, int fA, int cA, double **B, int fB, int cB);

/** @brief Inverse of a matrix of order n.
 *  Using Gauss Jordan Method
 *
 *  @param [in] L Matrix.
 *  @param [in] n Matrix order.
 *  @return Result matrix.
 */
double **m_inv(double **A, int n);

/** @brief Transpose of a matrix of order n.
 *
 *  @param [in] L Matrix.
 *  @param [in] n Matrix order.
 *  @return Result matrix.
 */
double **m_trans(double **A, int n);



#endif
/** @file VarEqn.h
 *  @brief Function prototypes for the computation of the
 *  variational equations.
 *
 *  This header file contains the prototypes for the computation of the
 *  variational equations, i.e. the derivative of the state vector
 *  and the state transition matrix.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _VAREQN_
#define _VAREQN_


/** @brief Variational equations.
 *
 *  @param [in] x Time since epoch in [s].
 *  @param [in] yPhi (6+36)-dim vector comprising the state vector (y) and
 *  the state transition matrix (Phi) in column wise storage order.
 *  @param [out] yPhip Derivative of yPhi.
 */
void VarEqn(double x, double *yPhi, double **yPhip);


#endif

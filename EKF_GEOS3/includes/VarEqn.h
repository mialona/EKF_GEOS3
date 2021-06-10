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

#include "global.h"
#include "const.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "m_utils.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"

#include <stdio.h>
#include <math.h>


/** @brief Variational equations.
 *
 *  @param [in] x Time since epoch in [s].
 *  @param [in] yPhi (6+36)-dim vector comprising the state vector (y) and
 *  the state transition matrix (Phi) in column wise storage order.
 *  @param [out] yPhip Derivative of yPhi.
 */
void VarEqn(double x, double *yPhi, double **yPhip);


#endif
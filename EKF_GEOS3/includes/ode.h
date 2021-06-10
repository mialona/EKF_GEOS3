/** @file ode.h
 *  @brief Function prototypes for the computation of the
 *  numerical integration methods.
 *
 *  This header file contains the prototypes for the numerical
 *  integration methods for ordinaray differential equations.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _ODE_
#define _ODE_


/** @brief Interface to an ordinary differential equation solver.
 *
 *  @param [in] f User-supplied function which accepts input
 *  values t and y, evaluates the right hand sides of the ODE,
 *  and stores the result in yp.
 *  @param [in] neqn Number of equations.
 *  @param [in,out] y Current vector solution.
 *  @param [in.out] t Current value of the independent variable.
 *  @param [in] tout Desired value of t on output.
 *  @param [in] relerr Relative error tolerances.
 *  @param [in] abserr Absolute error tolerances.
 *  @param [in,out] iflag Indicates the status of integration.
    On input, is normally 1.
 *  @param [in,out] work Workspace..
 *  @param [in,out] iwork Workspace..
 */
void ode ( void f ( double t, double *y, double **yp ), int neqn, double *y, 
  double *t, double tout, double relerr, double abserr, int *iflag, 
  double *work, int *iwork );


#endif

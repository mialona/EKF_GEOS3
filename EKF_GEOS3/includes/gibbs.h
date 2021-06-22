/** @file gibbs.h
 *  @brief Function prototypes for the 
 *  Gibbs method of orbit determination.
 *
 *  This header file contains the prototypes for the 
 *  Gibbs method of orbit determination. 
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _GIBBS_
#define _GIBBS_


/** @brief Gibbs method of orbit determination.
 *  
 *  This method determines the velocity at the middle point of the 3 given
 *  position vectors.
 *
 *  @param [in] r1 ijk position vector 1.
 *  @param [in] r2 ijk position vector 2.
 *  @param [in] r3 ijk position vector 3.
 *  @param [out] v2 ijk velocity vector for r2.
 *  @param [out] theta Angle between vectors.
 *  @param [out] theta1.
 *  @param [out] copa.
 *  @param [out] error Flag indicating success.
 */
void gibbs(double *r1, double *r2, double *r3, double **v2, double *theta, double *theta1, double *copa, char **error);


#endif

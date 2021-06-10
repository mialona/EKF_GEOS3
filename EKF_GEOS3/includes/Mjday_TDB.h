/** @file Mjday_TDB.h
 *  @brief Function prototypes for the Modified Julian Date
 *  for barycentric dynamical time.
 *
 *  This header file contains the prototypes for the
 *  calculation of the Modified Julian Date
 *  for barycentric dynamical time.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _MJDAY_TDB_
#define _MJDAY_TDB_


/** @brief Calculation of the Modified Julian Date
 *  for barycentric dynamical time.
 *  
 *  Reference:
 *  Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill;
 *  New York; 3rd edition(2007).
 *  
 *  @param [in] Mjd_TT Modified julian date (TT).
 *  @return Modified julian date (TDB).
 */
double Mjday_TDB(double Mjd_TT);


#endif

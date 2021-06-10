/** @file globals.h
 *  @brief Function prototypes for the 
 *  reading data files and global variables.
 *
 *  This header file contains the prototypes for the 
 *  reading data files and global variables.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _GLOBALES_
#define _GLOBALES_


typedef struct {
	double Mjd_UTC;
	int n, m, sun, moon, planets;
	double Mjd_TT;
} Param;


double **PC, **Cnm, **Snm, **eopdata, **obs;
int fPC, cPC, fCnm, cCnm, fSnm, cSnm, feopdata, ceopdata, fobs, cobs;
int n_eqn;
Param AuxParam;


/** @brief Read the GGM03S.txt file and store it in the matrix PC.
 *  
 *  @param [in] f Number of rows.
 *  @param [in] c Number of columns.
 */
void DE430Coeff(int f, int c);

/** @brief Read the GGM03S.txt file and store it in the Cnm and Snm
 *  matrices.
 *  
 *  @param [in] n Matrix order.
 */
void GGM03S(int n);

/** @brief Read the Earth orientation parameters in the eop19620101.txt file
 *  and store it in the matrix eopdata.
 *  
 *  @param [in] c Number of columns.
 */
void eop19620101(int c);

/** @brief Read the observations in the GEOS3.txt file and store
 *  it in the matrix obs.
 *  
 *  @param [in] f Number of rows.
 */
void GEOS3(int f);


#endif

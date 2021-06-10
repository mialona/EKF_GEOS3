/** @file anglesg.h
 *  @brief Function prototypes for the 
 *  orbit determination using three optical sightings.
 *
 *  This header file contains the prototypes for the 
 *  orbit determination using three optical sightings.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _ANGLESG_
#define _ANGLESG_


/** @brief Orbit determination using three optical sightings.
 *
 *  @param [in] az1 Azimuth at t1 [rad].
 *  @param [in] az2 Azimuth at t2 [rad].
 *  @param [in] az3 Azimuth at t3 [rad].
 *  @param [in] el1 Elevation at t1 [rad].
 *  @param [in] el2 Elevation at t2 [rad].
 *  @param [in] el3 Elevation at t3 [rad].
 *  @param [in] Mjd1 Modified julian date of t1.
 *  @param [in] Mjd2 Modified julian date of t2.
 *  @param [in] Mjd3 Modified julian date of t3.
 *  @param [in] Rs1 ijk site1 position vector [m].
 *  @param [in] Rs2 ijk site2 position vector [m].
 *  @param [in] Rs3 ijk site3 position vector [m].
 *  @param [out] r ijk position vector at t2 [m].
 *  @param [out] v ijk velocity vector at t2 [m/s].
 */
void anglesg(double az1, double az2, double az3, double el1,
		 	 double el2, double el3, double Mjd1, double Mjd2,
			 double Mjd3, double *Rs1, double *Rs2, double *Rs3,
			 double **r, double **v);


#endif

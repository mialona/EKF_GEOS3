/** @file Mjday.h
 *  @brief Function prototypes for the
 *  Modified julian date.
 *
 *  This header file contains the prototypes for the
 *  calculation of the Modified julian date.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _MJDAY_
#define _MJDAY_


/** @brief Calculation of the Modified julian date.
 *
 *  @param [in] year Year.
 *  @param [in] month Month.
 *  @param [in] day Day.
 *  @param [in] hr Universal time hour.
 *  @param [in] min Universal time min.
 *  @param [in] sec Universal time sec.
 *  @return Modified julian date.
 */
double Mjday(int yr, int mon, int day, int hr, int min, int sec);


#endif

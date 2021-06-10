/** @file TimeUpdate.h
 *  @brief Function prototypes for the 
 *  time updater.
 *
 *  This header file contains the prototypes for the 
 *  time updater.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No known bugs.
 */

#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "m_utils.h"

#include <stdio.h>


/** @brief Time updater.
 *
 *  @param [in] P Matrix 6x6.
 *  @param [in] Phi Matrix 6x6.
 *  @param [in] Qdt Matrix 6x6.
 *  @return Matrix 6x6.
 */
double **TimeUpdate(double **P, double **Phi, double **Qdt);


#endif
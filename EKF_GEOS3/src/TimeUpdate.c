/** @file TimeUpdate.c
 *  @brief Time updater.
 *
 *  This driver contains the code for the 
 *  time updater.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug Removed pre-initiated Qdt.
 */

#include "../includes/m_utils.h"

#include <stdio.h>


double **TimeUpdate(double **P, double **Phi, double **Qdt) {
	return m_sum(m_dot(Phi,6,6,m_dot(P,6,6,m_trans(Phi,6),6,6),6,6),6,6,Qdt,6,6);
}


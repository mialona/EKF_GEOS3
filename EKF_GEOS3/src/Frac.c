/** @file Frac.c
 *  @brief Fractional part of a number.
 *
 *  This driver contains the code for the fractional
 *  part of a number (y=x-[x]).
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include "../includes/Frac.h"

#include <stdio.h>
#include <math.h>


double Frac(double x) {
	return x-floor(x);
}

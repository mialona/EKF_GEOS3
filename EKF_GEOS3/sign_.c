/** @file sign_.c
 *  @brief Absolute value of a with sign of b.
 *
 *  This driver contains the code for the 
 *  absolute value of a with sign of b.
 *
 *  @author Miguel Alonso Angulo.
 *  @bug No know bugs.
 */

#include <math.h>


double sign_(double a, double b) {
	if(b >= 0.0) {
		return fabs(a);
	}
	else {
		return -fabs(a);
	}
}


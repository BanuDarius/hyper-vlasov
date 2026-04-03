#ifndef H_MATH_TOOLS_H
#define H_MATH_TOOLS_H

#include <math.h>

static inline double dot(const double *a, const double *b) {
	double x = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return x;
}

static inline double magnitude(const double *a) {
	double x = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	return x;
}

static inline void mult_vec(double *a, double *b, double c) {
	for(int i = 0; i < 3; i++)
		a[i] = b[i] * c;
}

#endif
#ifndef H_MATH_TOOLS_H
#define H_MATH_TOOLS_H

static inline double dot(const double *a, const double *b) {
	double x = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return x;
}

#endif
#ifndef FIT_ALGORITHM_H
#define FIT_ALGORITHM_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "sim_structs.h"

struct fit_data {
	struct test_particles *part;
	struct skyrme skm;
	int type, start, total;
};

void set_fit_function(struct fit_data *fit, struct test_particles *part, struct skyrme skm, int type, int start, int total);
int woods_saxon_f(const gsl_vector *x, void *p, gsl_vector *f);
int woods_saxon_df(const gsl_vector *x, void *p, gsl_matrix *J);
void minim_woods_saxon(struct test_particles *part, struct woods_saxon *ws, struct skyrme skm);


#endif
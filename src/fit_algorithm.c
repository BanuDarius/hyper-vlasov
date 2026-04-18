#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"
#include "fit_algorithm.h"

void set_fit_function(struct fit_data *fit, struct test_particles *part, struct skyrme skm, int type, int start, int total) {
	fit->skm = skm;
	fit->type = type;
	fit->part = part;
	fit->start = start;
	fit->total = total;
}

int woods_saxon_f(const gsl_vector *x, void *p, gsl_vector *f) {
	struct fit_data *fit = (struct fit_data*)p;
	double V0 = gsl_vector_get(x, 0), R12 = gsl_vector_get(x, 1), a = gsl_vector_get(x, 2);
	
	for(int i = 0; i < fit->total; i++) {
		int idx = fit->start + i;
		double r_vec[3];
		copy_particle_pos_to_vector(r_vec, *fit->part, idx);
		
		double r = magnitude(r_vec);
		double density_p = fit->part->density_p[idx];
		double density_n = fit->part->density_n[idx];
		double v_skyrme = skyrme_potential(fit->skm, density_p, density_n, fit->type);
		double v_woods_saxon = V0 / (1.0 + exp((r - R12) / a));
		
		gsl_vector_set(f, i, v_skyrme - v_woods_saxon);
	}
	return GSL_SUCCESS;
}

int woods_saxon_df(const gsl_vector *x, void *p, gsl_matrix *J) {
	struct fit_data *fit = (struct fit_data*)p;
	double V0 = gsl_vector_get(x, 0), R12 = gsl_vector_get(x, 1), a = gsl_vector_get(x, 2);
	
	for(int i = 0; i < fit->total; i++) {
		int idx = fit->start + i;
		double r_vec[3];
		copy_particle_pos_to_vector(r_vec, *fit->part, idx);
		
		double r = magnitude(r_vec);
		double exp_v = exp((r - R12) / a);
		double exp_v_squared = (1.0 + exp_v) * (1.0 + exp_v);
		
		double dV0 = -1.0 / (1.0 + exp_v);
		double dR12 = -(V0 * exp_v) / (a * exp_v_squared);
		double da = -(V0 * (r - R12) * exp_v) / (a * a * exp_v_squared);
		
		gsl_matrix_set(J, i, 0, dV0);
		gsl_matrix_set(J, i, 1, dR12);
		gsl_matrix_set(J, i, 2, da);
	}
	return GSL_SUCCESS;
}

void minim_woods_saxon(struct test_particles *part, struct woods_saxon *ws, struct skyrme skm) {
	const gsl_multifit_nlinear_type *T_MAGIC = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
	
	for(int type = 0; type <= 1; type++) {
		struct fit_data fit;
		int part_type = (type == 0) ? PROTONS : NEUTRONS;
		int start = (type == 0) ? 0 : part->protons;
		int total = (type == 0) ? part->protons : part->neutrons;
		set_fit_function(&fit, part, skm, part_type, start, total);
		
		gsl_multifit_nlinear_fdf magic_solver;
		magic_solver.f = woods_saxon_f;
		magic_solver.df = woods_saxon_df;
		magic_solver.fvv = NULL;
		magic_solver.n = total;
		magic_solver.p = 3;
		magic_solver.params = &fit;
		
		gsl_multifit_nlinear_workspace *magic_workspace = gsl_multifit_nlinear_alloc(T_MAGIC, &fdf_params, total, 3);
		
		gsl_vector *x = gsl_vector_alloc(3);
		gsl_vector_set(x, 0, ws[type].V0);
		gsl_vector_set(x, 1, ws[type].R12);
		gsl_vector_set(x, 2, ws[type].a);
		
		gsl_multifit_nlinear_init(x, &magic_solver, magic_workspace);
		
		int status, info;
		status = gsl_multifit_nlinear_driver(100, 1e-4, 1e-4, 1e-4, NULL, NULL, &info, magic_workspace);
		
		if (status != GSL_SUCCESS)
			printf("GSL Error %s\n", gsl_strerror(status));
		
		gsl_vector *fit_params = gsl_multifit_nlinear_position(magic_workspace);
		ws[type].V0 = gsl_vector_get(fit_params, 0);
		ws[type].R12 = gsl_vector_get(fit_params, 1);
		ws[type].a = gsl_vector_get(fit_params, 2);
		
		gsl_multifit_nlinear_free(magic_workspace);
		gsl_vector_free(x);
	}
}
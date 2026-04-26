/* MIT License

Copyright (c) 2026 Banu Darius-Matei

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */

#ifndef FIT_ALGORITHM_H
#define FIT_ALGORITHM_H

#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "sim_structs.hpp"

template <typename T>
struct FittingData {
	Skyrme<T> skm;
	int type, start, total;
	const TestParticles<T> *part;
};

template <typename T>
void set_fit_function(struct FittingData<T> *fit, const TestParticles<T> *part, const Skyrme<T> &skm, int type, int start, int total) {
	fit->skm = skm;
	fit->type = type;
	fit->part = part;
	fit->start = start;
	fit->total = total;
}

template <typename T>
int woods_saxon_f(const gsl_vector *x, void *p, gsl_vector *f) {
	struct FittingData<T> *fit = (FittingData<T>*)p;
	T V0 = T(gsl_vector_get(x, 0)), R12 = T(gsl_vector_get(x, 1)), a = T(gsl_vector_get(x, 2));
	
	for(int i = 0; i < fit->total; i++) {
		int idx = fit->start + i;
		T r_vec[3];
		copy_particle_pos_to_vector(r_vec, *fit->part, idx);
		
		T r = magnitude(r_vec);
		T density_p = fit->part->density_p[idx];
		T density_n = fit->part->density_n[idx];
		T v_skyrme = skyrme_potential(fit->skm, density_p, density_n, fit->type);
		T v_woods_saxon = V0 / (T(1.0) + std::exp((r - R12) / a));
		
		gsl_vector_set(f, i, v_skyrme - v_woods_saxon);
	}
	return GSL_SUCCESS;
}

template <typename T>
int woods_saxon_df(const gsl_vector *x, void *p, gsl_matrix *J) {
	struct FittingData<T> *fit = (struct FittingData<T>*)p;
	T V0 = T(gsl_vector_get(x, 0)), R12 = T(gsl_vector_get(x, 1)), a = T(gsl_vector_get(x, 2));
	
	for(int i = 0; i < fit->total; i++) {
		int idx = fit->start + i;
		T r_vec[3];
		copy_particle_pos_to_vector(r_vec, *fit->part, idx);
		
		T r = magnitude(r_vec);
		T exp_v = std::exp((r - R12) / a);
		T exp_v_squared = (T(1.0) + exp_v) * (T(1.0) + exp_v);
		
		T dV0 = -T(1.0) / (T(1.0) + exp_v);
		T dR12 = -(V0 * exp_v) / (a * exp_v_squared);
		T da = -(V0 * (r - R12) * exp_v) / (a * a * exp_v_squared);
		
		gsl_matrix_set(J, i, 0, dV0);
		gsl_matrix_set(J, i, 1, dR12);
		gsl_matrix_set(J, i, 2, da);
	}
	return GSL_SUCCESS;
}

template <typename T>
void minim_woods_saxon(TestParticles<T> *part, WoodsSaxon<T> *ws, Skyrme<T> skm) {
	const gsl_multifit_nlinear_type *T_MAGIC = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
	
	for(int type = 0; type <= 1; type++) {
		struct FittingData<T> fit;
		int part_type = (type == 0) ? PROTONS : NEUTRONS;
		int start = (type == 0) ? 0 : part->protons;
		int total = (type == 0) ? part->protons : part->neutrons;
		
		set_fit_function(&fit, part, skm, part_type, start, total);
		gsl_multifit_nlinear_fdf magic_solver;
		magic_solver.p = 3;
		magic_solver.n = total;
		magic_solver.fvv = NULL;
		magic_solver.params = &fit;
		magic_solver.f = woods_saxon_f<T>;
		magic_solver.df = woods_saxon_df<T>;
		
		gsl_multifit_nlinear_workspace *magic_workspace = gsl_multifit_nlinear_alloc(T_MAGIC, &fdf_params, total, 3);
		
		gsl_vector *x = gsl_vector_alloc(3);
		gsl_vector_set(x, 0, ws[type].V0);
		gsl_vector_set(x, 1, ws[type].R12);
		gsl_vector_set(x, 2, ws[type].a);
		
		gsl_multifit_nlinear_init(x, &magic_solver, magic_workspace);
		
		int status, info;
		status = gsl_multifit_nlinear_driver(100, 1e-4, 1e-4, 1e-4, NULL, NULL, &info, magic_workspace);
		
		if (status != GSL_SUCCESS)
			std::printf("GSL Error %s\n", gsl_strerror(status));
		
		gsl_vector *fit_params = gsl_multifit_nlinear_position(magic_workspace);
		ws[type].V0 = gsl_vector_get(fit_params, 0);
		ws[type].R12 = gsl_vector_get(fit_params, 1);
		ws[type].a = gsl_vector_get(fit_params, 2);
		
		gsl_multifit_nlinear_free(magic_workspace);
		gsl_vector_free(x);
	}
}

#endif
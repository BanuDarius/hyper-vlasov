// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#include <cstdio>

#include "fit_algorithm.hpp"
#include "tools.hpp"
#include "math_functions.hpp"
#include "physics_formulas.hpp"
#include "particle_in_cell.hpp"

template <std::floating_point T>
void set_fit_function(FittingData<T> *fit, const TestParticles<T> &part, const ScalarField<T> &density, const World<T> &world, const Skyrme<T> &skm, int type, int start, int total) {
	fit->skm = &skm;
	fit->type = type;
	fit->part = &part;
	fit->start = start;
	fit->total = total;
	fit->world = &world;
	fit->density = &density;
}

template <std::floating_point T>
int woods_saxon_f(const gsl_vector *x, void *p, gsl_vector *f) {
	FittingData<T> *fit = static_cast<FittingData<T>*>(p);
	T V0 = T(gsl_vector_get(x, 0)), R12 = T(gsl_vector_get(x, 1)), a = T(gsl_vector_get(x, 2));
	
	for(int i = 0; i < fit->total; i++) {
		int idx = fit->start + i;
		std::array<T, 3> r_vec;
		particle_pos_to_vector(r_vec, *fit->part, idx);
		
		T r = magnitude(r_vec);
		T density_p = scatter_scalar_field_cic(*(fit->density), r_vec, *(fit->world), is_proton);;
		T density_n = scatter_scalar_field_cic(*(fit->density), r_vec, *(fit->world), is_neutron);
		T v_skyrme = skyrme_potential(*(fit->skm), density_p, density_n, fit->type);
		T v_woods_saxon = V0 / (T(1.0) + std::exp((r - R12) / a));
		
		gsl_vector_set(f, i, v_skyrme - v_woods_saxon);
	}
	return GSL_SUCCESS;
}

template <std::floating_point T>
int woods_saxon_df(const gsl_vector *x, void *p, gsl_matrix *j) {
	FittingData<T> *fit = static_cast<FittingData<T>*>(p);
	T V0 = T(gsl_vector_get(x, 0)), R12 = T(gsl_vector_get(x, 1)), a = T(gsl_vector_get(x, 2));
	
	for(int i = 0; i < fit->total; i++) {
		int idx = fit->start + i;
		std::array<T, 3> r_vec;
		particle_pos_to_vector(r_vec, *(fit->part), idx);
		
		T r = magnitude(r_vec);
		T exp_v = std::exp((r - R12) / a);
		T exp_v_squared = (T(1.0) + exp_v) * (T(1.0) + exp_v);
		
		T dV0 = -T(1.0) / (T(1.0) + exp_v);
		T dR12 = -(V0 * exp_v) / (a * exp_v_squared);
		T da = -(V0 * (r - R12) * exp_v) / (a * a * exp_v_squared);
		
		gsl_matrix_set(j, i, 0, dV0);
		gsl_matrix_set(j, i, 1, dR12);
		gsl_matrix_set(j, i, 2, da);
	}
	return GSL_SUCCESS;
}

template <std::floating_point T>
void minim_woods_saxon(WoodsSaxon<T> &ws, TestParticles<T> &part, const ScalarField<T> &density, const Skyrme<T> &skm, const World<T> &world) {
	const gsl_multifit_nlinear_type *T_MAGIC = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
	
	for(int type = 0; type <= 1; type++) {
		FittingData<T> fit;
		int part_type = (type == 0) ? is_proton : is_neutron;
		int start = (type == 0) ? 0 : part.protons;
		int total = (type == 0) ? part.protons : part.neutrons;
		
		set_fit_function(&fit, part, density, world, skm, part_type, start, total);
		gsl_multifit_nlinear_fdf magic_solver;
		magic_solver.p = 3;
		magic_solver.n = total;
		magic_solver.fvv = NULL;
		magic_solver.params = &fit;
		magic_solver.f = woods_saxon_f<T>;
		magic_solver.df = woods_saxon_df<T>;
		
		gsl_multifit_nlinear_workspace *magic_workspace = gsl_multifit_nlinear_alloc(T_MAGIC, &fdf_params, total, 3);
		
		gsl_vector *x = gsl_vector_alloc(3);
		if(part_type == is_proton) {
			gsl_vector_set(x, 0, ws.V0_p);
			gsl_vector_set(x, 1, ws.R12_p);
			gsl_vector_set(x, 2, ws.a_p);
		} else {
			gsl_vector_set(x, 0, ws.V0_n);
			gsl_vector_set(x, 1, ws.R12_n);
			gsl_vector_set(x, 2, ws.a_n);
		}
		
		gsl_multifit_nlinear_init(x, &magic_solver, magic_workspace);
		
		int status, info;
		status = gsl_multifit_nlinear_driver(100, 1e-4, 1e-4, 1e-4, NULL, NULL, &info, magic_workspace);
		
		if (status != GSL_SUCCESS)
			std::fprintf(stderr, "GSL Error %s\n", gsl_strerror(status));
		
		gsl_vector *fit_params = gsl_multifit_nlinear_position(magic_workspace);
		if(part_type == is_proton) {
			ws.V0_p = gsl_vector_get(fit_params, 0);
			ws.R12_p = gsl_vector_get(fit_params, 1);
			ws.a_p = gsl_vector_get(fit_params, 2);
		} else {
			ws.V0_n = gsl_vector_get(fit_params, 0);
			ws.R12_n = gsl_vector_get(fit_params, 1);
			ws.a_n = gsl_vector_get(fit_params, 2);
		}
		
		gsl_multifit_nlinear_free(magic_workspace);
		gsl_vector_free(x);
	}
}

template void set_fit_function<double>(FittingData<double> *fit, const TestParticles<double> &part, const ScalarField<double> &density, const World<double> &world, const Skyrme<double> &skm, int type, int start, int total);
template int woods_saxon_f<double>(const gsl_vector *x, void *p, gsl_vector *f);
template int woods_saxon_df<double>(const gsl_vector *x, void *p, gsl_matrix *j);
template void minim_woods_saxon<double>(WoodsSaxon<double> &ws, TestParticles<double> &part, const ScalarField<double> &density, const Skyrme<double> &skm, const World<double> &world);

template void set_fit_function<float>(FittingData<float> *fit, const TestParticles<float> &part, const ScalarField<float> &density, const World<float> &world, const Skyrme<float> &skm, int type, int start, int total);
template int woods_saxon_f<float>(const gsl_vector *x, void *p, gsl_vector *f);
template int woods_saxon_df<float>(const gsl_vector *x, void *p, gsl_matrix *j);
template void minim_woods_saxon<float>(WoodsSaxon<float> &ws, TestParticles<float> &part, const ScalarField<float> &density, const Skyrme<float> &skm, const World<float> &world);
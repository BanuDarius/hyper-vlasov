// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#ifndef FIT_ALGORITHM_H
#define FIT_ALGORITHM_H

#include <concepts>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "sim_structs.hpp"

template <std::floating_point T>
struct FittingData {
	const Skyrme<T> *skm;
	const World<T> *world;
	int type, start, total;
	const TestParticles<T> *part;
	const ScalarField<T> *density;
};

template <std::floating_point T> void set_fit_function(FittingData<T> *fit, const TestParticles<T> &part, const ScalarField<T> &density, const World<T> &world, const Skyrme<T> &skm, int type, int start, int total);
template <std::floating_point T> int woods_saxon_f(const gsl_vector *x, void *p, gsl_vector *f);
template <std::floating_point T> int woods_saxon_df(const gsl_vector *x, void *p, gsl_matrix *j);
template <std::floating_point T> void minim_woods_saxon(WoodsSaxon<T> &ws, TestParticles<T> &part, const ScalarField<T> &density, const Skyrme<T> &skm, const World<T> &world);

#endif
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
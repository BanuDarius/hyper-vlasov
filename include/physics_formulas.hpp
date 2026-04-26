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

#ifndef PHYSICS_FORMULAS_H
#define PHYSICS_FORMULAS_H

#include <cmath>

#include "sim_structs.hpp"

template <typename T>
T nuclear_radius(int a) {
	T radius = T(1.5) * std::pow(T(a), T(1.0) / T(3.0));
	return radius;
}

template <typename T>
int max_particles(T r_max, T k_max, int part_per_nucleon) {
	T t = r_max * k_max, ct = T(2.0) * M_PI;
	T phase_space_volume = T(16.0 / 9.0) * M_PI * M_PI * (t * t * t);
	int max = part_per_nucleon * (int)std::floor(phase_space_volume / (ct * ct * ct) + T(0.5));
	return max;
}

template <typename T>
T kinetic_energy() {
	T hc2 = T(H_BAR_C) * T(H_BAR_C);
	T e_kin = hc2 / (T(2.0) * T(MC2));
	return e_kin;
}

template <typename T>
T fluctuation_energy(T sigma_k) {
	T e_fluc = T(3.0) * kinetic_energy<T>() * sigma_k * sigma_k;
	return e_fluc;
}

template <typename T>
T calc_sigma(T fwhm) {
	T t = T(2.0) * std::sqrt(T(2.0) * std::log(T(2.0)));
	T sigma = fwhm / t;
	return sigma;
}

template <typename T>
static inline T woods_saxon_potential(WoodsSaxon<T> ws, T r) {
	T v = ws.V0 / (T(1.0) + std::exp((r - ws.R12) / ws.a));
	return v;
}

template <typename T>
static inline T skyrme_potential(const Skyrme<T> &skm, T rho_p, T rho_n, int type) {
	T tau = (type == PROTONS) ? T(-1.0) : T(+1.0);
	T rho = rho_p + rho_n;
	T t = rho / T(RHO_0);
	T v = skm.A * t + skm.B * std::pow(t, T(skm.gamma)) + T(2.0) * tau * skm.C * ((rho_n - rho_p) / T(RHO_0));
	return v;
}

template <typename T>
static inline T coulomb_potential(const WoodsSaxon<T> &ws, T z, T r) {
	T R12 = ws.R12, v;
	if(r <= ws.R12)
		v = T(1.44) * (z - T(1.0)) / R12 * (T(1.5) - T(0.5) * (r / R12) * (r / R12));
	else
		v = T(1.44) * (z - T(1.0)) / r;
	return v;
}

#endif
// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#ifndef PHYSICS_FORMULAS_H
#define PHYSICS_FORMULAS_H

#include <cmath>
#include <concepts>

#include "sim_structs.hpp"

template <std::floating_point T>
T nuclear_radius(int a) noexcept {
	T radius = T(1.5) * std::pow(T(a), T(1.0) / T(3.0));
	return radius;
}

template <std::floating_point T>
int max_particles(T r_max, T k_max, int part_per_nucleon) noexcept {
	T t = r_max * k_max, ct = T(2.0) * pi<T>;
	T phase_space_volume = T(16.0) / T(9.0) * pi<T> * pi<T> * (t * t * t);
	int max = part_per_nucleon * static_cast<int>(std::floor(phase_space_volume / (ct * ct * ct) + T(0.5)));
	return max;
}

template <std::floating_point T>
T kinetic_energy() noexcept {
	T hc2 = h_bar_c<T> * h_bar_c<T>;
	T e_kin = hc2 / (T(2.0) * mc2<T>);
	return e_kin;
}

template <std::floating_point T>
T fluctuation_energy(T sigma_k) noexcept {
	T e_fluc = T(3.0) * kinetic_energy<T>() * sigma_k * sigma_k;
	return e_fluc;
}

template <std::floating_point T>
T calc_sigma(T fwhm) noexcept {
	T t = T(2.0) * std::sqrt(T(2.0) * std::log(T(2.0)));
	T sigma = fwhm / t;
	return sigma;
}

template <std::floating_point T>
inline T woods_saxon_potential(const WoodsSaxon<T> &ws, T r, int type) noexcept {
	if(type == is_proton)
		return ws.V0_p / (T(1.0) + std::exp((r - ws.R12_p) / ws.a_p));
	else
		return ws.V0_n / (T(1.0) + std::exp((r - ws.R12_n) / ws.a_n));
}

template <std::floating_point T>
inline T skyrme_potential(const Skyrme<T> &skm, T rho_p, T rho_n, int type) noexcept {
	T tau = (type == is_proton) ? T(-1.0) : T(+1.0);
	T rho = rho_p + rho_n;
	T t = rho / rho_0<T>;
	return skm.A * t + skm.B * std::pow(t, skm.gamma) + T(2.0) * tau * skm.C * ((rho_n - rho_p) / rho_0<T>);
}

template <std::floating_point T>
inline T coulomb_potential(const WoodsSaxon<T> &ws, T z, T r) noexcept {
	if(r <= ws.R12_p)
		return T(1.44) * (z - T(1.0)) / ws.R12_p * (T(1.5) - T(0.5) * (r / ws.R12_p) * (r / ws.R12_p));
	else
		return T(1.44) * (z - T(1.0)) / r;
}

#endif
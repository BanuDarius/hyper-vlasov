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

#ifndef TOOLS_H
#define TOOLS_H

#include <omp.h>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>

template <typename T>
void compute_particle_densities(TestParticles<T> &part, const Parameters<T> &param) {
	int part_per_nucleon = param.part_per_nucleon, protons = part.protons, neutrons = part.neutrons, total = protons + neutrons;
	T sigma_r = param.sigma_r, exp_term = T(1.0) / (T(2.0) * sigma_r * sigma_r);
	T cutoff_squared = T(16.0) * sigma_r * sigma_r;
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		std::array<T, 3> r_i, r_j, diff;
		T fact, dist_squared, density_p = T(0.0), density_n = T(0.0);
		copy_particle_pos_to_vector(r_i, part, i);
		
		for(int j = 0; j < total; j++) {
			copy_particle_pos_to_vector(r_j, part, j);
			
			diff = r_i - r_j;
			dist_squared = dot(diff, diff);
			if(dist_squared > cutoff_squared)
				continue;
			fact = exp(-dist_squared * exp_term);
			if(j < protons)
				density_p += fact;
			else
				density_n += fact;
		}
		part.density_p[i] = density_p;
		part.density_n[i] = density_n;
	}
	T term = (T(1.0) / part_per_nucleon) * (T(1.0) / std::pow(T(2.0) * pi<T> * sigma_r * sigma_r, T(1.5)));
	#pragma omp parallel for simd
	for(int i = 0; i < total; i++) {
		part.density_p[i] *= term;
		part.density_n[i] *= term;
	}
}

template <typename T>
T compute_energy(TestParticles<T> &part, const WoodsSaxon<T> &ws, T sigma_k, int z, int i) {
	std::array<T, 3> r_vec, k_vec;
	copy_particle_pos_to_vector(r_vec, part, i);
	copy_particle_vel_to_vector(k_vec, part, i);
	
	T r = magnitude(r_vec);
	T k = magnitude(k_vec);
	
	int type = (i < part.protons) ? is_proton : is_neutron;
	T energy = woods_saxon_potential(ws, r, type);
	energy += (k * k) * kinetic_energy<T>();
	energy += fluctuation_energy(sigma_k);
	
	if(i < part.protons)
		energy += coulomb_potential<T>(ws, z, r);
	return energy;
}

template <typename T>
void compute_particle_energies(TestParticles<T> &part, const WoodsSaxon<T> &ws, const Parameters<T> &param) {
	T sigma_k = param.sigma_k, z = param.z;
	#pragma omp parallel for
	for(int i = 0; i < part.protons + part.neutrons; i++)
		part.energy[i] = compute_energy(part, ws, sigma_k, z, i);
}

template <typename T>
void generate_random_particles(TestParticles<T> &part, T r_max) {
	int total = part.protons + part.neutrons, i = 0;
	while(i < total) {
		std::array<T, 3> r_new;
		random_vec(r_new, r_max);
		if(dot(r_new, r_new) < r_max * r_max) {
			copy_vector_to_particle_pos(part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < total) {
		std::array<T, 3> k_new;
		random_vec(k_new, k_max<T>);
		if(dot(k_new, k_new) < k_max<T> * k_max<T>) {
			copy_vector_to_particle_vel(part, k_new, i);
			i++;
		}
	}
}

template <typename T>
void generate_checking_particles(TestParticles<T> &part, const WoodsSaxon<T> &ws, const Parameters<T> &param, const Fermi<T> &fermi_levels) {
	T r_max = param.r_max, sigma_k = param.sigma_k, z = param.z, epsilon;
	int total = part.protons + part.neutrons, i = 0;
	while(i < total) {
		if(i < part.protons)
			epsilon = fermi_levels.epsilon_p;
		else
			epsilon = fermi_levels.epsilon_n;
		
		std::array<T, 3> r_new, k_new;
		random_vec(r_new, r_max);
		random_vec(k_new, k_max<T>);
		
		copy_vector_to_particle_pos(part, r_new, i);
		copy_vector_to_particle_vel(part, k_new, i);
		T energy = compute_energy(part, ws, sigma_k, z, i);
		if(energy < epsilon) {
			r_new = r_new * T(-1.0);
			k_new = k_new * T(-1.0);
			
			copy_vector_to_particle_pos(part, r_new, i + 1);
			copy_vector_to_particle_vel(part, k_new, i + 1);
			i+=2;
		}
	}
}

template <typename T>
void compute_coulomb_boundaries(ScalarField<T> &coulomb, const TestParticles<T> &part, const World<T> &world, int z) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	std::array<T, 3> cm_protons = center_of_mass(part, world, is_proton);
	#pragma omp parallel for collapse(3)
	for(int i = 0; i < nx; i++) {
		for(int j = 0; j < ny; j++) {
			for(int k = 0; k < nz; k++) {
				if(i == 0 || j == 0 || k == 0 || i == nx - 1 || j == ny - 1 || k == nz - 1) {
					int idx = grid_idx(i, j, k, nx, ny, nz);
					std::array<T, 3> r_vec;
					world_pos_to_vector(r_vec, world, idx);
					
					r_vec = r_vec - cm_protons;
					T r = magnitude(r_vec);
					coulomb.v[idx] = T(1.44) * z / r;
				}
			}
		}
	}
}

template <typename T>
T mean_squared_radius(const TestParticles<T> &part, const World<T> &world, int type) {
	int start, end, part_num = 0;
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	if(type == is_proton) { start = 0; end = part.protons; }
	else if(type == is_neutron) { start = part.protons; end = part.protons + part.neutrons; }
	
	T r_sqr = T(0.0);
	#pragma omp parallel for reduction(+:r_sqr, part_num)
	for(int i = start; i < end; i++) {
		std::array<T, 3> r_vec;
		copy_particle_pos_to_vector(r_vec, part, i);
		
		if(r_vec[0] < -d_max_x || r_vec[0] > +d_max_x
		|| r_vec[1] < -d_max_y || r_vec[1] > +d_max_y
		|| r_vec[2] < -d_max_z || r_vec[2] > +d_max_z)
			continue;
		
		T r2 = dot(r_vec, r_vec);
		r_sqr += r2;
		part_num++;
	}
	return r_sqr / static_cast<T>(part_num);
}

template <typename T>
std::array<T, 3> center_of_mass(const TestParticles<T> &part, const World<T> &world, int type) {
	int start, end, part_num = 0;
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	if(type == is_proton) { start = 0; end = part.protons; }
	else if(type == is_neutron) { start = part.protons; end = part.protons + part.neutrons; }
	
	T center_x = T(0.0), center_y = T(0.0), center_z = T(0.0);
	#pragma omp parallel for reduction(+:center_x, center_y, center_z, part_num)
	for(int i = start; i < end; i++) {
		std::array<T, 3> r_vec;
		copy_particle_pos_to_vector(r_vec, part, i);
		
		if(r_vec[0] < -d_max_x || r_vec[0] > +d_max_x
		|| r_vec[1] < -d_max_y || r_vec[1] > +d_max_y
		|| r_vec[2] < -d_max_z || r_vec[2] > +d_max_z)
			continue;
		
		center_x += r_vec[0];
		center_y += r_vec[1];
		center_z += r_vec[2];
		part_num++;
	}
	center_x /= static_cast<T>(part_num); center_y /= static_cast<T>(part_num), center_z /= static_cast<T>(part_num);
	return std::array<T, 3> { center_x, center_y, center_z };
}

template <typename T>
void chi_squared(const TestParticles<T> &part, const WoodsSaxon<T> &ws, const Skyrme<T> &skm, int part_per_nucleon) {
	int total = part.protons + part.neutrons;
	T chi_squared_p = 0.0, chi_squared_n = 0.0;
	#pragma omp parallel for reduction(+:chi_squared_p, chi_squared_n)
	for(int i = 0; i < total; i++) {
		int type = (i < part.protons) ? is_proton : is_neutron;
		std::array<T, 3> r_vec;
		copy_particle_pos_to_vector(r_vec, part, i);
		T r = magnitude(r_vec);
		
		T density_p = part.density_p[i];
		T density_n = part.density_n[i];
		
		T v_ws = woods_saxon_potential(ws, r, type);
		T v_skyrme = skyrme_potential(skm, density_p, density_n, type);
		T diff = v_ws - v_skyrme;
		if(type == is_proton)
			chi_squared_p += diff * diff;
		else
			chi_squared_n += diff * diff;
	}
	chi_squared_n /= part_per_nucleon; chi_squared_p /= part_per_nucleon;
	std::printf("CHI SQUARED P %0.2lf CHI SQUARED N %0.2lf\n", chi_squared_p, chi_squared_n);
}

template <typename T>
void relax_woods_saxon(WoodsSaxon<T> &ws, const WoodsSaxon<T> &ws_old, T coef) {
	ws.V0_p = coef * ws.V0_p + (1.0 - coef) * ws_old.V0_p;
	ws.V0_n = coef * ws.V0_n + (1.0 - coef) * ws_old.V0_n;
	ws.R12_p = coef * ws.R12_p + (1.0 - coef) * ws_old.R12_p;
	ws.R12_n = coef * ws.R12_n + (1.0 - coef) * ws_old.R12_n;
	ws.a_p = coef * ws.a_p + (1.0 - coef) * ws_old.a_p;
	ws.a_n = coef * ws.a_n + (1.0 - coef) * ws_old.a_n;
}

template <typename T>
void add_scalar_field_single(ScalarField<T> &field_a, const ScalarField<T> &field_b, const ScalarField<T> &field_c, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < world_size; i++)
		field_a.v[i] = field_b.v[i] + field_c.v[i];
}

template <typename T>
void sub_scalar_field_double(ScalarField<T> &field_a, const ScalarField<T> &field_b, const ScalarField<T> &field_c, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		field_a.v[i] = field_b.v[i] - field_c.v[i];
}

template <typename T>
void copy_scalar_field_double(ScalarField<T> &field_a, const ScalarField<T> &field_b, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		field_a.v[i] = field_b.v[i];
}


template <typename T>
static inline void world_pos_to_vector(std::array<T, 3> &v, const World<T> &world, int idx) {
	int x = world.n[0], y = world.n[1], z = world.n[2];
	int i = idx / (y * z), j = (idx / z) % y, k = idx % z;
	v[0] = world.d_max[0] * (T(2.0) * i / x - T(1.0));
	v[1] = world.d_max[1] * (T(2.0) * j / y - T(1.0));
	v[2] = world.d_max[2] * (T(2.0) * k / z - T(1.0));
}

template <typename T>
static inline void copy_particle_pos_to_vector(std::array<T, 3> &v, const TestParticles<T> &part, int i) {
	v[0] = part.x[i];
	v[1] = part.y[i];
	v[2] = part.z[i];
}

template <typename T>
static inline void copy_particle_vel_to_vector(std::array<T, 3> &v, const TestParticles<T> &part, int i) {
	v[0] = part.kx[i];
	v[1] = part.ky[i];
	v[2] = part.kz[i];
}

template <typename T>
static inline void copy_vector_to_particle_pos(TestParticles<T> &part, const std::array<T, 3> &v, int i) {
	part.x[i] = v[0];
	part.y[i] = v[1];
	part.z[i] = v[2];
}

template <typename T>
static inline void copy_vector_to_particle_vel(TestParticles<T> &part, const std::array<T, 3> &v, int i) {
	part.kx[i] = v[0];
	part.ky[i] = v[1];
	part.kz[i] = v[2];
}

#endif
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

#include <omp.h>
#include <array>
#include <cstdio>
#include <cstdlib>

#include "tools.hpp"
#include "sim_structs.hpp"
#include "math_functions.hpp"
#include "physics_formulas.hpp"
#include "particle_in_cell.hpp"

template <typename T>
T compute_energy(TestParticles<T> &part, const WoodsSaxon<T> &ws, T sigma_k, int z, int i) {
	std::array<T, 3> r_vec, k_vec;
	particle_pos_to_vector(r_vec, part, i);
	particle_vel_to_vector(k_vec, part, i);
	
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
			vector_to_particle_pos(part, r_new, i); i++;
		}
	}
	i = 0;
	while(i < total) {
		std::array<T, 3> k_new;
		random_vec(k_new, k_max<T>);
		if(dot(k_new, k_new) < k_max<T> * k_max<T>) {
			vector_to_particle_vel(part, k_new, i); i++;
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
		
		vector_to_particle_pos(part, r_new, i);
		vector_to_particle_vel(part, k_new, i);
		T energy = compute_energy(part, ws, sigma_k, z, i);
		if(energy < epsilon) {
			r_new = r_new * T(-1.0);
			k_new = k_new * T(-1.0);
			
			vector_to_particle_pos(part, r_new, i + 1);
			vector_to_particle_vel(part, k_new, i + 1);
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
	#pragma omp parallel for reduction(+:part_num, r_sqr)
	for(int i = start; i < end; i++) {
		std::array<T, 3> r_vec;
		particle_pos_to_vector(r_vec, part, i);
		
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
	
	std::array<T, 3> center = { T(0.0) };
	T *center_ptr = center.data();
	#pragma omp parallel for reduction(+:part_num, center_ptr[0 : 3])
	for(int i = start; i < end; i++) {
		std::array<T, 3> r_vec;
		particle_pos_to_vector(r_vec, part, i);
		
		if(r_vec[0] < -d_max_x || r_vec[0] > +d_max_x
		|| r_vec[1] < -d_max_y || r_vec[1] > +d_max_y
		|| r_vec[2] < -d_max_z || r_vec[2] > +d_max_z)
			continue;
		
		center_ptr[0] += r_vec[0];
		center_ptr[1] += r_vec[1];
		center_ptr[2] += r_vec[2];
		part_num++;
	}
	center = center * (T(1.0) / static_cast<T>(part_num));
	return center;
}

template <typename T>
std::array<T, 2> chi_squared(const ScalarField<T> &density, const WoodsSaxon<T> &ws, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	T chi_squared_p = 0.0, chi_squared_n = 0.0;
	#pragma omp parallel for reduction(+:chi_squared_p, chi_squared_n)
	for(int i = 0; i < 2 * world_size; i++) {
		int type = (i < world_size) ? is_proton : is_neutron;
		std::array<T, 3> r_vec;
		world_pos_to_vector(r_vec, world, i % world_size);
		
		T density_p = scatter_scalar_field_cic(density, r_vec, world, is_proton);
		T density_n = scatter_scalar_field_cic(density, r_vec, world, is_neutron);
		
		T r = magnitude(r_vec);
		T v_ws = woods_saxon_potential(ws, r, type);
		T v_skyrme = skyrme_potential(skm, density_p, density_n, type);
		T diff = v_ws - v_skyrme;
		if(type == is_proton)
			chi_squared_p += diff * diff;
		else
			chi_squared_n += diff * diff;
	}
	chi_squared_p /= param.part_per_nucleon; chi_squared_n /= param.part_per_nucleon;
	return std::array<T, 2> { chi_squared_p, chi_squared_n };
}

template <typename T>
void relax_woods_saxon(WoodsSaxon<T> &ws, const WoodsSaxon<T> &ws_old, T coef) {
	ws.V0_p = coef * ws.V0_p + (T(1.0) - coef) * ws_old.V0_p;
	ws.V0_n = coef * ws.V0_n + (T(1.0) - coef) * ws_old.V0_n;
	ws.R12_p = coef * ws.R12_p + (T(1.0) - coef) * ws_old.R12_p;
	ws.R12_n = coef * ws.R12_n + (T(1.0) - coef) * ws_old.R12_n;
	ws.a_p = coef * ws.a_p + (T(1.0) - coef) * ws_old.a_p;
	ws.a_n = coef * ws.a_n + (T(1.0) - coef) * ws_old.a_n;
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
void copy_vector_field_double(VectorField<T> &field_a, const VectorField<T> &field_b, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		field_a.x[i] = field_b.x[i];
		field_a.y[i] = field_b.y[i];
		field_a.z[i] = field_b.z[i];
	}
}

template double compute_energy<double>(TestParticles<double> &part, const WoodsSaxon<double> &ws, double sigma_k, int z, int i);
template void compute_particle_energies<double>(TestParticles<double> &part, const WoodsSaxon<double> &ws, const Parameters<double> &param);
template void generate_random_particles<double>(TestParticles<double> &part, double r_max);
template void generate_checking_particles<double>(TestParticles<double> &part, const WoodsSaxon<double> &ws, const Parameters<double> &param, const Fermi<double> &fermi_levels);
template void compute_coulomb_boundaries<double>(ScalarField<double> &coulomb, const TestParticles<double> &part, const World<double> &world, int z);
template double mean_squared_radius<double>(const TestParticles<double> &part, const World<double> &world, int type);
template std::array<double, 3> center_of_mass<double>(const TestParticles<double> &part, const World<double> &world, int type);
template std::array<double, 2> chi_squared<double>(const ScalarField<double> &density, const WoodsSaxon<double> &ws, const Skyrme<double> &skm, const Parameters<double> &param, const World<double> &world);
template void relax_woods_saxon<double>(WoodsSaxon<double> &ws, const WoodsSaxon<double> &ws_old, double coef);
template void add_scalar_field_single<double>(ScalarField<double> &field_a, const ScalarField<double> &field_b, const ScalarField<double> &field_c, const World<double> &world);
template void sub_scalar_field_double(ScalarField<double> &field_a, const ScalarField<double> &field_b, const ScalarField<double> &field_c, const World<double> &world);
template void copy_scalar_field_double(ScalarField<double> &field_a, const ScalarField<double> &field_b, const World<double> &world);
template void copy_vector_field_double(VectorField<double> &field_a, const VectorField<double> &field_b, const World<double> &world);

template float compute_energy<float>(TestParticles<float> &part, const WoodsSaxon<float> &ws, float sigma_k, int z, int i);
template void compute_particle_energies<float>(TestParticles<float> &part, const WoodsSaxon<float> &ws, const Parameters<float> &param);
template void generate_random_particles<float>(TestParticles<float> &part, float r_max);
template void generate_checking_particles<float>(TestParticles<float> &part, const WoodsSaxon<float> &ws, const Parameters<float> &param, const Fermi<float> &fermi_levels);
template void compute_coulomb_boundaries<float>(ScalarField<float> &coulomb, const TestParticles<float> &part, const World<float> &world, int z);
template float mean_squared_radius<float>(const TestParticles<float> &part, const World<float> &world, int type);
template std::array<float, 3> center_of_mass<float>(const TestParticles<float> &part, const World<float> &world, int type);
template std::array<float, 2> chi_squared<float>(const ScalarField<float> &density, const WoodsSaxon<float> &ws, const Skyrme<float> &skm, const Parameters<float> &param, const World<float> &world);
template void relax_woods_saxon<float>(WoodsSaxon<float> &ws, const WoodsSaxon<float> &ws_old, float coef);
template void add_scalar_field_single<float>(ScalarField<float> &field_a, const ScalarField<float> &field_b, const ScalarField<float> &field_c, const World<float> &world);
template void sub_scalar_field_double(ScalarField<float> &field_a, const ScalarField<float> &field_b, const ScalarField<float> &field_c, const World<float> &world);
template void copy_scalar_field_double(ScalarField<float> &field_a, const ScalarField<float> &field_b, const World<float> &world);
template void copy_vector_field_double(VectorField<float> &field_a, const VectorField<float> &field_b, const World<float> &world);
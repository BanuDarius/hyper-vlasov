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

#ifndef PHYSICS_H
#define PHYSICS_H

#include <array>
#include <cmath>
#include <cstdio>

#include "fit_algorithm.hpp"

template <typename T>
void initialize_particles(TestParticles<T> &part, WoodsSaxon<T> &ws, const Skyrme<T> &skm, Fermi<T> &fermi_levels, const Parameters<T> &param) {
	T total_delta_epsilon, relax_coef = T(0.6);
	int max_part = param.max_test_part, z = param.z, n = param.n, part_per_nucleon = param.part_per_nucleon;
	int total_p = z * part_per_nucleon, total_n = n * part_per_nucleon, it = 0;
	
	TestParticles<T> temp_part(max_part, max_part);
	do {
		generate_random_particles(temp_part, param.r_max);
		compute_particle_energies(temp_part, ws, param);
		int less_p = 0, equal_p = 0, more_p = 0;
		int less_n = 0, equal_n = 0, more_n = 0;
		for(int i = 0; i < max_part; i++) {
			if(temp_part.energy[i] < fermi_levels.epsilon_p)
				equal_p += 2;
			if(temp_part.energy[i] < fermi_levels.epsilon_p + T(0.5))
				more_p += 2;
			if(temp_part.energy[i] < fermi_levels.epsilon_p - T(0.5))
				less_p += 2;
			
			if(temp_part.energy[i + max_part] < fermi_levels.epsilon_n)
				equal_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels.epsilon_n + T(0.5))
				more_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels.epsilon_n - T(0.5))
				less_n += 2;
		}
		T delta_part_n = total_n - equal_n;
		T delta_part_p = total_p - equal_p;
		T delta_epsilon_n = T(0.5) * delta_part_n / (more_n - less_n);
		T delta_epsilon_p = T(0.5) * delta_part_p / (more_p - less_p);
		
		if(std::abs(delta_epsilon_p) > T(0.5)) delta_epsilon_p *= relax_coef;
		if(std::abs(delta_epsilon_n) > T(0.5)) delta_epsilon_n *= relax_coef;
		
		fermi_levels.epsilon_p += delta_epsilon_p;
		fermi_levels.epsilon_n += delta_epsilon_n;
		
		generate_checking_particles(part, ws, param, fermi_levels);
		compute_particle_densities(part, param);
		
		WoodsSaxon<T> ws_old = ws;
		
		minim_woods_saxon(part, ws, skm);
		relax_woods_saxon(ws, ws_old, relax_coef);
		
		total_delta_epsilon = std::abs(delta_epsilon_n) + std::abs(delta_epsilon_p);
		std::printf("------------------------------------\n");
		std::printf("EQUAL P %i EQUAL N %i\n", equal_p, equal_n);
		std::printf("V0 %0.2lf R12 %0.2lf a %0.2lf\n", ws.V0_p, ws.R12_p, ws.a_p);
		std::printf("V0 %0.2lf R12 %0.2lf a %0.2lf\n", ws.V0_n, ws.R12_n, ws.a_n);
		std::printf("FERMI P %0.2lf FERMI N %0.2lf\n", fermi_levels.epsilon_p, fermi_levels.epsilon_n);
		std::printf("DELTA EPSILON %0.2lf\nITERATION %i\n", total_delta_epsilon, it + 1);
		
		it++;
	} while(total_delta_epsilon > delta_epsilon_tolerance<T> && it < max_init_iterations);
	if(it == max_init_iterations)
		fprintf(stderr, "INITIALIZATION DID NOT CONVERGE!\n");
	compute_particle_energies(part, ws, param);
}

template <typename T>
void compute_volumetric_coulomb_potentials_sor(ScalarField<T> &coulomb, const ScalarField<T> &density, const World<T> &world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	T dx = T(2.0) * world.d_max[0] / nx, dy = T(2.0) * world.d_max[1] / ny, dz = T(2.0) * world.d_max[2] / nz;
	T inv_dx2 = T(1.0) / (T(2.0) / (dx * dx) + T(2.0) / (dy * dy) + T(2.0) / (dz * dz)), omega = T(1.50), max_diff;
	for(int it = 0; it < max_sor_iterations; it++) {
		max_diff = T(0.0);
		#pragma omp parallel for collapse(3) reduction(max:max_diff)
		for(int i = 1; i < nx - 1; i++) {
			for(int j = 1; j < ny - 1; j++) {
				for(int k = 1; k < nz - 1; k++) {
					if((i + j + k) % 2 == 0) {
						int idx = grid_idx(i, j, k, nx, ny, nz);
						T rho = density.v[idx];
						
						T phi_x = (coulomb.v[grid_idx(i + 1, j, k, nx, ny, nz)] + coulomb.v[grid_idx(i - 1, j, k, nx, ny, nz)]) / (dx * dx);
						T phi_y = (coulomb.v[grid_idx(i, j + 1, k, nx, ny, nz)] + coulomb.v[grid_idx(i, j - 1, k, nx, ny, nz)]) / (dy * dy);
						T phi_z = (coulomb.v[grid_idx(i, j, k + 1, nx, ny, nz)] + coulomb.v[grid_idx(i, j, k - 1, nx, ny, nz)]) / (dz * dz);
						
						T phi_star = (phi_x + phi_y + phi_z + T(4.0) * pi<T> * T(1.44) * rho) * inv_dx2;
						T phi_old = coulomb.v[idx];
						
						coulomb.v[idx] = (T(1.0) - omega) * phi_old + omega * phi_star;
						T diff = std::abs(coulomb.v[idx] - phi_old);
						if(diff > max_diff)
							max_diff = diff;
					}
				}
			}
		}
		#pragma omp parallel for collapse(3) reduction(max:max_diff)
		for(int i = 1; i < nx - 1; i++) {
			for(int j = 1; j < ny - 1; j++) {
				for(int k = 1; k < nz - 1; k++) {
					if((i + j + k) % 2 != 0) {
						int idx = grid_idx(i, j, k, nx, ny, nz);
						T rho = density.v[idx];
						
						T phi_x = (coulomb.v[grid_idx(i + 1, j, k, nx, ny, nz)] + coulomb.v[grid_idx(i - 1, j, k, nx, ny, nz)]) / (dx * dx);
						T phi_y = (coulomb.v[grid_idx(i, j + 1, k, nx, ny, nz)] + coulomb.v[grid_idx(i, j - 1, k, nx, ny, nz)]) / (dy * dy);
						T phi_z = (coulomb.v[grid_idx(i, j, k + 1, nx, ny, nz)] + coulomb.v[grid_idx(i, j, k - 1, nx, ny, nz)]) / (dz * dz);
						
						T phi_star = (phi_x + phi_y + phi_z + T(4.0) * pi<T> * T(1.44) * rho) * inv_dx2;
						T phi_old = coulomb.v[idx];
						
						coulomb.v[idx] = (T(1.0) - omega) * phi_old + omega * phi_star;
						T diff = std::abs(coulomb.v[idx] - phi_old);
						if(diff > max_diff)
							max_diff = diff;
					}
				}
			}
		}
		if(max_diff < sor_tolerance<T>)
			break;
	}
	if(max_diff > sor_tolerance<T>)
		std::fprintf(stderr, "SOR COULOMB DID NOT CONVERGE!\n");
}

template <typename T>
void compute_volumetric_forces_fdm(VectorField<T> &forces, const ScalarField<T> &potentials, const World<T> &world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	T dx = T(2.0) * world.d_max[0] / nx, dy = T(2.0) * world.d_max[1] / ny, dz = T(2.0) * world.d_max[2] / nz;
	for(int x = 0; x < 2; x++) {
		int offset = (x == 0) ? 0 : world_size;
		#pragma omp parallel for collapse(3)
		for(int i = 0; i < nx; i++) {
			for(int j = 0; j < ny; j++) {
				for(int k = 0; k < nz; k++) {
					int idx = grid_idx(i, j, k, nx, ny, nz);
					std::array<T, 3> gradient = { T(0.0) };
					
					if(i == 0)
						gradient[0] = (potentials.v[grid_idx(1, j, k, nx, ny, nz) + offset] - potentials.v[idx + offset]) / dx;
					else if(i == nx - 1)
						gradient[0] = (potentials.v[idx + offset] - potentials.v[grid_idx(nx - 2, j, k, nx, ny, nz) + offset]) / dx;
					else
						gradient[0] = (potentials.v[grid_idx(i + 1, j, k, nx, ny, nz) + offset] - potentials.v[grid_idx(i - 1, j, k, nx, ny, nz) + offset]) / (T(2.0) * dx);
					
					if(j == 0)
						gradient[1] = (potentials.v[grid_idx(i, 1, k, nx, ny, nz) + offset] - potentials.v[idx + offset]) / dy;
					else if(j == ny - 1)
						gradient[1] = (potentials.v[idx + offset] - potentials.v[grid_idx(i, ny - 2, k, nx, ny, nz) + offset]) / dy;
					else
						gradient[1] = (potentials.v[grid_idx(i, j + 1, k, nx, ny, nz) + offset] - potentials.v[grid_idx(i, j - 1, k, nx, ny, nz) + offset]) / (T(2.0) * dy);
					
					if(k == 0)
						gradient[2] = (potentials.v[grid_idx(i, j, 1, nx, ny, nz) + offset] - potentials.v[idx + offset]) / dz;
					else if(k == nz - 1)
						gradient[2] = (potentials.v[idx + offset] - potentials.v[grid_idx(i, j, nz - 2, nx, ny, nz) + offset]) / dz;
					else
						gradient[2] = (potentials.v[grid_idx(i, j, k + 1, nx, ny, nz) + offset] - potentials.v[grid_idx(i, j, k - 1, nx, ny, nz) + offset]) / (T(2.0) * dz);
					
					forces.x[idx + offset] = -gradient[0];
					forces.y[idx + offset] = -gradient[1];
					forces.z[idx + offset] = -gradient[2];
				}
			}
		}
	}
}

template <typename T>
void compute_volumetric_densities(ScalarField<T> &density, ScalarField<T> &density_temp, const Parameters<T> &param, const World<T> &world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	T sigma_r = param.sigma_r, exp_term = T(1.0) / (T(2.0) * sigma_r * sigma_r);
	T cutoff_squared = T(16.0) * sigma_r * sigma_r;
	
	copy_scalar_field_double(density_temp, density, world);
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		int offset = (i < world_size) ? 0 : world_size;
		std::array<T, 3> r_i, r_j, diff;
		T dist_squared, fact, rho_f = T(0.0);
		world_pos_to_vector(r_i, world, i % world_size);
		
		for(int j = 0; j < world_size; j++) {
			T rho = density_temp.v[j + offset];
			world_pos_to_vector(r_j, world, j);
			
			diff = r_i - r_j;
			dist_squared = dot(diff, diff);
			if(dist_squared > cutoff_squared)
				continue;
			
			fact = std::exp(-dist_squared * exp_term);
			rho_f += rho * fact;
		}
		density.v[i] = rho_f;
	}
	T term = (T(1.0) / param.part_per_nucleon) * (T(1.0) / std::pow(T(2.0) * pi<T> * sigma_r * sigma_r, T(1.5)));
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		density.v[i] *= term;
}

template <typename T>
void compute_current_velocity(VectorField<T> &current_velocity, VectorField<T> &current_velocity_temp, const ScalarField<T> &density, const Parameters<T> &param, const World<T> &world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	T sigma_r = param.sigma_r, exp_term = T(1.0) / (T(2.0) * sigma_r * sigma_r);
	T cutoff_squared = T(16.0) * sigma_r * sigma_r;
	
	copy_vector_field_double(current_velocity_temp, current_velocity, world);
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		int offset = (i < world_size) ? 0 : world_size;
		std::array<T, 3> r_i, r_j, diff, velocity_f = { T(0.0) };
		T dist_squared, fact;
		world_pos_to_vector(r_i, world, i % world_size);
		
		for(int j = 0; j < world_size; j++) {
			world_pos_to_vector(r_j, world, j);
			
			diff = r_i - r_j;
			dist_squared = dot(diff, diff);
			if(dist_squared > cutoff_squared)
				continue;
			
			fact = std::exp(-dist_squared * exp_term);
			velocity_f[0] += fact * current_velocity_temp.x[j + offset];
			velocity_f[1] += fact * current_velocity_temp.y[j + offset];
			velocity_f[2] += fact * current_velocity_temp.z[j + offset];
		}
		current_velocity.x[i] = velocity_f[0];
		current_velocity.y[i] = velocity_f[1];
		current_velocity.z[i] = velocity_f[2];
	}
	T term = (T(1.0) / param.part_per_nucleon) * (T(1.0) / std::pow(T(2.0) * pi<T> * sigma_r * sigma_r, T(1.5)));
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		T rho = density.v[i];
		if(rho > density_tolerance<T>) {
			T final_term = h_bar_c<T> * term / (mc2<T> * rho);
			current_velocity.x[i] *= final_term;
			current_velocity.y[i] *= final_term;
			current_velocity.z[i] *= final_term;
		}
		else {
			current_velocity.x[i] = T(0.0);
			current_velocity.y[i] = T(0.0);
			current_velocity.z[i] = T(0.0);
		}
	}
}

template <typename T>
void center_momentum(TestParticles<T> &part, const World<T> &world) {
	int total = part.protons + part.neutrons, part_num = 0;
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	std::array<T, 3> k_sum = { T(0.0) };
	T *k_sum_ptr = k_sum.data();
	#pragma omp parallel for reduction(+:part_num, k_sum_ptr[0 : 3])
	for(int i = 0; i < total; i++) {
		std::array<T, 3> r_vec;
		particle_pos_to_vector(r_vec, part, i);
		
		if(r_vec[0] < -d_max_x || r_vec[0] > +d_max_x
		|| r_vec[1] < -d_max_y || r_vec[1] > +d_max_y
		|| r_vec[2] < -d_max_z || r_vec[2] > +d_max_z)
			continue;
		
		k_sum_ptr[0] += part.kx[i];
		k_sum_ptr[1] += part.ky[i];
		k_sum_ptr[2] += part.kz[i];
		part_num++;
	}
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		part.kx[i] -= k_sum[0] / static_cast<T>(part_num);
		part.ky[i] -= k_sum[1] / static_cast<T>(part_num);
		part.kz[i] -= k_sum[2] / static_cast<T>(part_num);
	}
}

template <typename T>
void nuclear_excitation(TestParticles<T> &part, const Parameters<T> &param) {
	int protons = part.protons, neutrons = part.neutrons;
	T z = param.z, n = param.n, eta = param.eta_exc;
	#pragma omp parallel for
	for(int i = 0; i < protons; i++)
		part.kz[i] += eta * n / (z + n);
	#pragma omp parallel for
	for(int i = protons; i < protons + neutrons; i++)
		part.kz[i] -= eta * z / (z + n);
}

template <typename T>
void update_momenta_half(TestParticles<T> &part, T dt) {
	int total = part.protons + part.neutrons;
	T fact = dt / (T(2.0) * h_bar_c<T>);
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		part.kx[i] += fact * part.fx[i];
		part.ky[i] += fact * part.fy[i];
		part.kz[i] += fact * part.fz[i];
	}
}

template <typename T>
void update_positions_full(TestParticles<T> &part, T dt) {
	int total = part.protons + part.neutrons;
	T fact = (dt * h_bar_c<T>) / mc2<T>;
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		part.x[i] += fact * part.kx[i];
		part.y[i] += fact * part.ky[i];
		part.z[i] += fact * part.kz[i];
	}
}

template <typename T>
void compute_volumetric_skyrme_potentials(ScalarField<T> &potentials, const ScalarField<T> &density, const Skyrme<T> &skm, const World<T> &world) {
	int x = world.n[0], y = world.n[1], z = world.n[2], world_size = x * y * z;
	#pragma omp parallel for
	for(int i = 0; i < world_size; i++)
		potentials.v[i] = skyrme_potential(skm, density.v[i], density.v[i + world_size], is_proton);
	#pragma omp parallel for
	for(int i = world_size; i < 2 * world_size; i++)
		potentials.v[i] = skyrme_potential(skm, density.v[i - world_size], density.v[i], is_neutron);
}

#endif
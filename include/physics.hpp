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

#include <cmath>
#include <cstdio>

#include "physics_formulas.hpp"
#include "math_functions.hpp"
#include "sim_structs.hpp"
#include "fit_algorithm.hpp"

template <typename T>
void initialize_particles(TestParticles<T> *part, const Parameters<T> &param, WoodsSaxon<T> *ws, const Skyrme<T> &skm, Fermi<T> *fermi_levels) {
	T total_delta_epsilon, relax_coef = T(0.6);
	int max_part = param.max_test_part, z = param.z, n = param.n, part_per_nucleon = param.part_per_nucleon;
	int total_p = z * part_per_nucleon, total_n = n * part_per_nucleon, it = 0;
	
	TestParticles<T> temp_part;
	create_particles(part, total_p, total_n);
	create_particles(&temp_part, max_part, max_part);
	do {
		generate_random_particles(&temp_part, param.r_max);
		compute_particle_energies(&temp_part, ws, param);
		int less_p = 0, equal_p = 0, more_p = 0;
		int less_n = 0, equal_n = 0, more_n = 0;
		for(int i = 0; i < max_part; i++) {
			if(temp_part.energy[i] < fermi_levels->epsilon_p)
				equal_p += 2;
			if(temp_part.energy[i] < fermi_levels->epsilon_p + T(0.5))
				more_p += 2;
			if(temp_part.energy[i] < fermi_levels->epsilon_p - T(0.5))
				less_p += 2;
			
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n)
				equal_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n + T(0.5))
				more_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n - T(0.5))
				less_n += 2;
		}
		T delta_part_n = total_n - equal_n;
		T delta_part_p = total_p - equal_p;
		T delta_epsilon_n = T(0.5) * delta_part_n / (more_n - less_n);
		T delta_epsilon_p = T(0.5) * delta_part_p / (more_p - less_p);
		
		if(std::abs(delta_epsilon_p) > T(0.5)) delta_epsilon_p *= relax_coef;
		if(std::abs(delta_epsilon_n) > T(0.5)) delta_epsilon_n *= relax_coef;
		
		fermi_levels->epsilon_p += delta_epsilon_p;
		fermi_levels->epsilon_n += delta_epsilon_n;
		
		generate_checking_particles(part, ws, param, fermi_levels);
		compute_particle_densities(part, param);
		
		WoodsSaxon<T> ws_old[2];
		ws_old[0] = ws[0]; ws_old[1] = ws[1];
		
		minim_woods_saxon(part, ws, skm);
		relax_woods_saxon(ws, ws_old, relax_coef);
		
		total_delta_epsilon = std::abs(delta_epsilon_n) + std::abs(delta_epsilon_p);
		std::printf("----------------\n");
		std::printf("EQUAL P %i EQUAL N %i\n", equal_p, equal_n);
		std::printf("V0 %0.2lf R12 %0.2lf a %0.2lf\n", ws[0].V0, ws[0].R12, ws[0].a);
		std::printf("V0 %0.2lf R12 %0.2lf a %0.2lf\n", ws[1].V0, ws[1].R12, ws[1].a);
		std::printf("FERMI P %0.2lf FERMI N %0.2lf\n", fermi_levels->epsilon_p, fermi_levels->epsilon_n);
		std::printf("DELTA EPSILON %0.2lf\nITERATION %i\n", total_delta_epsilon, it);
		
		it++;
	} while(total_delta_epsilon > delta_epsilon_tolerance<T> && it < MAX_INIT_ITERATIONS);	
	compute_particle_energies(part, ws, param);
	free_particles(&temp_part);
}

template <typename T>
void compute_volumetric_coulomb_potentials_sor(ScalarField<T> *coulomb, const ScalarField<T> &density, const World<T> &world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	T dx = T(2.0) * world.d_max[0] / nx, dy = T(2.0) * world.d_max[1] / ny, dz = T(2.0) * world.d_max[2] / nz;
	T inv_dx2 = T(1.0) / (T(2.0) / (dx * dx) + T(2.0) / (dy * dy) + T(2.0) / (dz * dz)), omega = T(1.50), max_diff;
	for(int it = 0; it < MAX_SOR_ITERATIONS; it++) {
		max_diff = T(0.0);
		#pragma omp parallel for collapse(3) reduction(max:max_diff)
		for(int i = 1; i < nx - 1; i++) {
			for(int j = 1; j < ny - 1; j++) {
				for(int k = 1; k < nz - 1; k++) {
					if((i + j + k) % 2 == 0) {
						int idx = IDX(i, j, k, nx, ny, nz);
						T rho = density.v[idx];
						
						T phi_x = (coulomb->v[IDX(i + 1, j, k, nx, ny, nz)] + coulomb->v[IDX(i - 1, j, k, nx, ny, nz)]) / (dx * dx);
						T phi_y = (coulomb->v[IDX(i, j + 1, k, nx, ny, nz)] + coulomb->v[IDX(i, j - 1, k, nx, ny, nz)]) / (dy * dy);
						T phi_z = (coulomb->v[IDX(i, j, k + 1, nx, ny, nz)] + coulomb->v[IDX(i, j, k - 1, nx, ny, nz)]) / (dz * dz);
						
						T phi_star = (phi_x + phi_y + phi_z + T(4.0) * pi<T> * T(1.44) * rho) * inv_dx2;
						T phi_old = coulomb->v[idx];
						
						coulomb->v[idx] = (T(1.0) - omega) * phi_old + omega * phi_star;
						T diff = std::abs(coulomb->v[idx] - phi_old);
						if(diff > max_diff) max_diff = diff;
					}
				}
			}
		}
		#pragma omp parallel for collapse(3) reduction(max:max_diff)
		for(int i = 1; i < nx - 1; i++) {
			for(int j = 1; j < ny - 1; j++) {
				for(int k = 1; k < nz - 1; k++) {
					if((i + j + k) % 2 != 0) {
						int idx = IDX(i, j, k, nx, ny, nz);
						T rho = density.v[idx];
						
						T phi_x = (coulomb->v[IDX(i + 1, j, k, nx, ny, nz)] + coulomb->v[IDX(i - 1, j, k, nx, ny, nz)]) / (dx * dx);
						T phi_y = (coulomb->v[IDX(i, j + 1, k, nx, ny, nz)] + coulomb->v[IDX(i, j - 1, k, nx, ny, nz)]) / (dy * dy);
						T phi_z = (coulomb->v[IDX(i, j, k + 1, nx, ny, nz)] + coulomb->v[IDX(i, j, k - 1, nx, ny, nz)]) / (dz * dz);
						
						T phi_star = (phi_x + phi_y + phi_z + T(4.0) * pi<T> * T(1.44) * rho) * inv_dx2;
						T phi_old = coulomb->v[idx];
						
						coulomb->v[idx] = (T(1.0) - omega) * phi_old + omega * phi_star;
						T diff = std::abs(coulomb->v[idx] - phi_old);
						if(diff > max_diff) max_diff = diff;
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
void compute_volumetric_forces_fdm(VectorField<T> *forces, ScalarField<T> potentials, World<T> world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	T dx = T(2.0) * world.d_max[0] / nx, dy = T(2.0) * world.d_max[1] / ny, dz = T(2.0) * world.d_max[2] / nz;
	for(int x = 0; x < 2; x++) {
		int offset = (x == 0) ? 0 : nx * ny * nz;
		#pragma omp parallel for collapse(3)
		for(int i = 0; i < nx; i++) {
			for(int j = 0; j < ny; j++) {
				for(int k = 0; k < nz; k++) {
					T gradient_x, gradient_y, gradient_z;
					
					if(i == 0)
						gradient_x = (potentials.v[IDX(1, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k, nx, ny, nz) + offset]) / dx;
					else if(i == nx - 1)
						gradient_x = (potentials.v[IDX(i, j, k, nx, ny, nz) + offset] - potentials.v[IDX(nx - 2, j, k, nx, ny, nz) + offset]) / dx;
					else
						gradient_x = (potentials.v[IDX(i + 1, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i - 1, j, k, nx, ny, nz) + offset]) / (T(2.0) * dx);
					
					if(j == 0)
						gradient_y = (potentials.v[IDX(i, 1, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k, nx, ny, nz) + offset]) / dy;
					else if(j == ny - 1)
						gradient_y = (potentials.v[IDX(i, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i, ny - 2, k, nx, ny, nz) + offset]) / dy;
					else
						gradient_y = (potentials.v[IDX(i, j + 1, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j - 1, k, nx, ny, nz) + offset]) / (T(2.0) * dy);
					
					if(k == 0)
						gradient_z = (potentials.v[IDX(i, j, 1, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k, nx, ny, nz) + offset]) / dz;
					else if(k == nz - 1)
						gradient_z = (potentials.v[IDX(i, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j, nz - 2, nx, ny, nz) + offset]) / dz;
					else
						gradient_z = (potentials.v[IDX(i, j, k + 1, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k - 1, nx, ny, nz) + offset]) / (T(2.0) * dz);
					
					forces->x[IDX(i, j, k, nx, ny, nz) + offset] = -gradient_x;
					forces->y[IDX(i, j, k, nx, ny, nz) + offset] = -gradient_y;
					forces->z[IDX(i, j, k, nx, ny, nz) + offset] = -gradient_z;
				}
			}
		}
	}
}

template <typename T>
void update_momenta_half(TestParticles<T> *part, T dt) {
	int total = part->protons + part->neutrons;
	T fact = dt / (T(2.0) * h_bar_c<T>);
	#pragma omp parallel for simd
	for(int i = 0; i < total; i++) {
		part->kx[i] += fact * part->fx[i];
		part->ky[i] += fact * part->fy[i];
		part->kz[i] += fact * part->fz[i];
	}
}

template <typename T>
void update_positions_full(TestParticles<T> *part, T dt) {
	int total = part->protons + part->neutrons;
	T fact = (dt * h_bar_c<T>) / mc2<T>;
	#pragma omp parallel for simd
	for(int i = 0; i < total; i++) {
		part->x[i] += fact * part->kx[i];
		part->y[i] += fact * part->ky[i];
		part->z[i] += fact * part->kz[i];
	}
}
template <typename T>
void compute_volumetric_skyrme_potentials(ScalarField<T> *potentials, const ScalarField<T> &density, const Skyrme<T> &skm, const World<T> &world) {
	int x = world.n[0], y = world.n[1], z = world.n[2], world_size = x * y * z;
	#pragma omp parallel for
	for(int i = 0; i < world_size; i++)
		potentials->v[i] = skyrme_potential(skm, density.v[i], density.v[i + world_size], PROTONS);
	#pragma omp parallel for
	for(int i = world_size; i < 2 * world_size; i++)
		potentials->v[i] = skyrme_potential(skm, density.v[i - world_size], density.v[i], NEUTRONS);
}

#endif
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

#include <math.h>
#include <stdio.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"
#include "fit_algorithm.h"

void initialize_particles(TestParticles *part, Parameters param, WoodsSaxon *ws, Skyrme skm, Fermi *fermi_levels) {
	double r_max = param.r_max, total_delta_epsilon, relax_coef = 0.6;
	int max_part = param.max_test_part, z = param.z, n = param.n, part_per_nucleon = param.part_per_nucleon, it = 0;
	int total_p = z * part_per_nucleon, total_n = n * part_per_nucleon;
	
	TestParticles temp_part;
	create_particles(part, total_p, total_n);
	create_particles(&temp_part, max_part, max_part);
	do {
		generate_random_particles(&temp_part, r_max);
		compute_particle_energies(&temp_part, ws, param);
		int check_less_p = 0, check_equal_p = 0, check_more_p = 0;
		int check_less_n = 0, check_equal_n = 0, check_more_n = 0;
		for(int i = 0; i < max_part; i++) {
			if(temp_part.energy[i] < fermi_levels->epsilon_p)
				check_equal_p += 2;
			if(temp_part.energy[i] < fermi_levels->epsilon_p + 0.5)
				check_more_p += 2;
			if(temp_part.energy[i] < fermi_levels->epsilon_p - 0.5)
				check_less_p += 2;
			
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n)
				check_equal_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n + 0.5)
				check_more_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n - 0.5)
				check_less_n += 2;
		}
		double delta_part_n = total_n - check_equal_n;
		double delta_part_p = total_p - check_equal_p;
		double delta_epsilon_n = 0.5 * delta_part_n / (check_more_n - check_less_n);
		double delta_epsilon_p = 0.5 * delta_part_p / (check_more_p - check_less_p);
		
		if(fabs(delta_epsilon_p) > 0.5) delta_epsilon_p *= relax_coef;
		if(fabs(delta_epsilon_n) > 0.5) delta_epsilon_n *= relax_coef;
		
		fermi_levels->epsilon_p += delta_epsilon_p;
		fermi_levels->epsilon_n += delta_epsilon_n;
		
		generate_checking_particles(part, ws, param, fermi_levels);
		compute_particle_densities(part, param);
		
		WoodsSaxon ws_old[2];
		ws_old[0] = ws[0]; ws_old[1] = ws[1];
		
		minim_woods_saxon(part, ws, skm);
		relax_woods_saxon(ws, ws_old, relax_coef);
		
		total_delta_epsilon = fabs(delta_epsilon_n) + fabs(delta_epsilon_p);
		printf("%i %i %i\n", check_less_p, check_equal_p, check_more_p);
		printf("%i %i %i\n", check_less_n, check_equal_n, check_more_n);
		printf("WS PARAM %0.2lf %0.2lf %0.2lf\n", ws[0].V0, ws[0].R12, ws[0].a);
		printf("WS PARAM %0.2lf %0.2lf %0.2lf\n", ws[1].V0, ws[1].R12, ws[1].a);
		printf("FERMI P %0.2lf FERMI N %0.2lf\n", fermi_levels->epsilon_p, fermi_levels->epsilon_n);
		printf("DELTA EPSILON %0.2lf\nITERATION %i\n", total_delta_epsilon, it);
		
		it++;
	} while(total_delta_epsilon > DELTA_EPSILON_TOLERANCE && it < MAX_INIT_ITERATIONS);	
	compute_particle_energies(part, ws, param);
	free_particles(&temp_part);
}

void compute_volumetric_coulomb_potentials_sor(ScalarField *coulomb, ScalarField volume, World world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	double dx = 2.0 * world.d_max[0] / nx, dy = 2.0 * world.d_max[1] / ny, dz = 2.0 * world.d_max[2] / nz;
	double inv_dx2 = 1.0 / (2.0 / (dx * dx) + 2.0 / (dy * dy) + 2.0 / (dz * dz)), omega = 1.80, max_diff;
	for(int it = 0; it < MAX_SOR_ITERATIONS; it++) {
		max_diff = 0.0;
		#pragma omp parallel for collapse(3) reduction(max:max_diff)
		for(int i = 1; i < nx - 1; i++) {
			for(int j = 1; j < ny - 1; j++) {
				for(int k = 1; k < nz - 1; k++) {
					if((i + j + k) % 2 == 0) {
						int idx = IDX(i, j, k, nx, ny, nz);
						double density = volume.v[idx];
						
						double phi_x = (coulomb->v[IDX(i + 1, j, k, nx, ny, nz)] + coulomb->v[IDX(i - 1, j, k, nx, ny, nz)]) / (dx * dx);
						double phi_y = (coulomb->v[IDX(i, j + 1, k, nx, ny, nz)] + coulomb->v[IDX(i, j - 1, k, nx, ny, nz)]) / (dy * dy);
						double phi_z = (coulomb->v[IDX(i, j, k + 1, nx, ny, nz)] + coulomb->v[IDX(i, j, k - 1, nx, ny, nz)]) / (dz * dz);
						
						double phi_star = (phi_x + phi_y + phi_z + 4.0 * M_PI * 1.44 * density) * inv_dx2;
						double phi_old = coulomb->v[idx];
						
						coulomb->v[idx] = (1.0 - omega) * phi_old + omega * phi_star;
						double diff = fabs(coulomb->v[idx] - phi_old);
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
						double density = volume.v[idx];
						
						double phi_x = (coulomb->v[IDX(i + 1, j, k, nx, ny, nz)] + coulomb->v[IDX(i - 1, j, k, nx, ny, nz)]) / (dx * dx);
						double phi_y = (coulomb->v[IDX(i, j + 1, k, nx, ny, nz)] + coulomb->v[IDX(i, j - 1, k, nx, ny, nz)]) / (dy * dy);
						double phi_z = (coulomb->v[IDX(i, j, k + 1, nx, ny, nz)] + coulomb->v[IDX(i, j, k - 1, nx, ny, nz)]) / (dz * dz);
						
						double phi_star = (phi_x + phi_y + phi_z + 4.0 * M_PI * 1.44 * density) * inv_dx2;
						double phi_old = coulomb->v[idx];
						
						coulomb->v[idx] = (1.0 - omega) * phi_old + omega * phi_star;
						double diff = fabs(coulomb->v[idx] - phi_old);
						if(diff > max_diff) max_diff = diff;
					}
				}
			}
		}
		if(max_diff < SOR_TOLERANCE)
			break;
	}
	if(max_diff > SOR_TOLERANCE)
		fprintf(stderr, "SOR COULOMB DID NOT CONVERGE!\n");
}

void compute_volumetric_forces_fdm(VectorField *forces, ScalarField potentials, World world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	double dx = 2.0 * world.d_max[0] / nx, dy = 2.0 * world.d_max[1] / ny, dz = 2.0 * world.d_max[2] / nz;
	for(int x = 0; x < 2; x++) {
		int offset = (x == 0) ? 0 : nx * ny * nz;
		#pragma omp parallel for collapse(3)
		for(int i = 0; i < nx; i++) {
			for(int j = 0; j < ny; j++) {
				for(int k = 0; k < nz; k++) {
					double gradient_x, gradient_y, gradient_z;
					
					if(i == 0)
						gradient_x = (potentials.v[IDX(1, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k, nx, ny, nz) + offset]) / dx;
					else if(i == nx - 1)
						gradient_x = (potentials.v[IDX(i, j, k, nx, ny, nz) + offset] - potentials.v[IDX(nx - 2, j, k, nx, ny, nz) + offset]) / dx;
					else
						gradient_x = (potentials.v[IDX(i + 1, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i - 1, j, k, nx, ny, nz) + offset]) / (2.0 * dx);
					
					if(j == 0)
						gradient_y = (potentials.v[IDX(i, 1, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k, nx, ny, nz) + offset]) / dy;
					else if(j == ny - 1)
						gradient_y = (potentials.v[IDX(i, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i, ny - 2, k, nx, ny, nz) + offset]) / dy;
					else
						gradient_y = (potentials.v[IDX(i, j + 1, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j - 1, k, nx, ny, nz) + offset]) / (2.0 * dy);
					
					if(k == 0)
						gradient_z = (potentials.v[IDX(i, j, 1, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k, nx, ny, nz) + offset]) / dz;
					else if(k == nz - 1)
						gradient_z = (potentials.v[IDX(i, j, k, nx, ny, nz) + offset] - potentials.v[IDX(i, j, nz - 2, nx, ny, nz) + offset]) / dz;
					else
						gradient_z = (potentials.v[IDX(i, j, k + 1, nx, ny, nz) + offset] - potentials.v[IDX(i, j, k - 1, nx, ny, nz) + offset]) / (2.0 * dz);
					
					forces->x[IDX(i, j, k, nx, ny, nz) + offset] = -gradient_x;
					forces->y[IDX(i, j, k, nx, ny, nz) + offset] = -gradient_y;
					forces->z[IDX(i, j, k, nx, ny, nz) + offset] = -gradient_z;
				}
			}
		}
	}
}

void update_momenta_half(TestParticles *part, double dt) {
	int total = part->protons + part->neutrons;
	double fact = dt / (2.0 * H_BAR_C);
	#pragma omp parallel for simd
	for(int i = 0; i < total; i++) {
		part->kx[i] += fact * part->fx[i];
		part->ky[i] += fact * part->fy[i];
		part->kz[i] += fact * part->fz[i];
	}
}

void update_positions_full(TestParticles *part, double dt) {
	int total = part->protons + part->neutrons;
	double fact = dt * (H_BAR_C / MC2);
	#pragma omp parallel for simd
	for(int i = 0; i < total; i++) {
		part->x[i] += fact * part->kx[i];
		part->y[i] += fact * part->ky[i];
		part->z[i] += fact * part->kz[i];
	}
}

void compute_volumetric_skyrme_potentials(ScalarField *potentials, ScalarField volume, Skyrme skm, World world) {
	int x = world.n[0], y = world.n[1], z = world.n[2], world_size = x * y * z;
	#pragma omp parallel for
	for(int i = 0; i < world_size; i++)
		potentials->v[i] = skyrme_potential(skm, volume.v[i], volume.v[i + world_size], PROTONS);
	#pragma omp parallel for
	for(int i = world_size; i < 2 * world_size; i++)
		potentials->v[i] = skyrme_potential(skm, volume.v[i - world_size], volume.v[i], NEUTRONS);
}

double nuclear_radius(unsigned short a) {
	double radius = 1.5 * pow((double)a, 1.0 / 3.0);
	return radius;
}

int max_particles(double r_max, double k_max, int part_per_nucleon) {
	double t = r_max * k_max, ct = 2.0 * M_PI;
	double phase_space_volume = (16.0 / 9.0) * M_PI * M_PI * (t * t * t);
	int max = part_per_nucleon * (int)floor(phase_space_volume / (ct * ct * ct) + 0.5);
	return max;
}
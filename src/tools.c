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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"

void compute_volumetric_density_cic(ScalarField *volume, TestParticles *part, Parameters param, World world) {
	double d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int x = world.n[0], y = world.n[1], z = world.n[2], world_size = x * y * z, total = part->protons + part->neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		volume->v[i] = 0.0;
	
	double *volume_ptr = volume->v;
	#pragma omp parallel for reduction(+: volume_ptr[0 : 2 * world_size])
	for(int i = 0; i < total; i++) {
		double r_vec[3];
		copy_particle_pos_to_vector(r_vec, *part, i);
		
		double cx = (x / 2.0) * (r_vec[0] / d_max_x + 1.0);
		double cy = (y / 2.0) * (r_vec[1] / d_max_y + 1.0);
		double cz = (z / 2.0) * (r_vec[2] / d_max_z + 1.0);
		
		int x0 = (int)cx; int y0 = (int)cy; int z0 = (int)cz;
		
		if(x0 < 0 || y0 < 0 || z0 < 0 || x0 >= x || y0 >= y || z0 >= z)
			continue;
			
		double d_x = cx - x0; double d_y = cy - y0; double d_z = cz - z0;
		double t_x = 1.0 - d_x; double t_y = 1.0 - d_y; double t_z = 1.0 - d_z;
		int x1 = x0 + 1; int y1 = y0 + 1; int z1 = z0 + 1;
		
		int offset = (i < part->protons) ? 0 : world_size;
		int idx000 = x0 * (y * z) + y0 * z + z0 + offset;
		volume_ptr[idx000] += t_x * t_y * t_z;
		if (x1 < x) {
			int idx100 = x1 * (y * z) + y0 * z + z0 + offset;
			volume_ptr[idx100] += d_x * t_y * t_z;
		}
		if (y1 < y) {
			int idx010 = x0 * (y * z) + y1 * z + z0 + offset;
			volume_ptr[idx010] += t_x * d_y * t_z;
		}
		if (x1 < x && y1 < y) {
			int idx110 = x1 * (y * z) + y1 * z + z0 + offset;
			volume_ptr[idx110] += d_x * d_y * t_z;
		}
		if (z1 < z) {
			int idx001 = x0 * (y * z) + y0 * z + z1 + offset;
			volume_ptr[idx001] += t_x * t_y * d_z;
		}
		if (x1 < x && z1 < z) {
			int idx101 = x1 * (y * z) + y0 * z + z1 + offset;
			volume_ptr[idx101] += d_x * t_y * d_z;
		}
		if (y1 < y && z1 < z) {
			int idx011 = x0 * (y * z) + y1 * z + z1 + offset;
			volume_ptr[idx011] += t_x * d_y * d_z;
		}
		if (x1 < x && y1 < y && z1 < z) {
			int idx111 = x1 * (y * z) + y1 * z + z1 + offset;
			volume_ptr[idx111] += d_x * d_y * d_z;
		}
	}
	ScalarField temp_volume;
	create_scalar_field_double(&temp_volume, world);
	copy_scalar_field(&temp_volume, *volume, world);
	
	double sigma_r = param.sigma_r, exp_term = 1.0 / (2.0 * sigma_r * sigma_r);
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		double r_i[3], r_j[3], diff[3];
		double fact, dist_squared, density = 0.0;
		world_pos_to_vector(r_i, world, i % world_size);
		
		for(int j = 0; j < world_size; j++) {
			int offset = (i < world_size) ? 0 : world_size;
			int idx = j + offset;
			double count = temp_volume.v[idx];
			world_pos_to_vector(r_j, world, j);
			
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			fact = exp(-dist_squared * exp_term);
			density += count * fact;
		}
		volume->v[i] = density;
	}
	free_scalar_field(&temp_volume);
	
	double term = (1.0 / param.part_per_nucleon) * (1.0 / (pow(2.0 * M_PI * sigma_r * sigma_r, 1.5)));
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		volume->v[i] *= term;
	}
}

double compute_energy(TestParticles *part, WoodsSaxon *ws, double sigma_k, int z, int i) {
	double r_vec[3], k_vec[3];
	copy_particle_pos_to_vector(r_vec, *part, i);
	copy_particle_vel_to_vector(k_vec, *part, i);
	
	double r = magnitude(r_vec);
	double k = magnitude(k_vec);
	
	double energy = 0.0;
	WoodsSaxon ws_c = (i < part->protons) ? ws[0] : ws[1];
	energy += woods_saxon_potential(ws_c, r);
	energy += (k * k) * kinetic_energy();
	
	if(i < part->protons)
		energy += coulomb_potential(ws_c, (double)z, r);
	
	energy += fluctuation_energy(sigma_k);
	return energy;
}

void compute_particle_energies(TestParticles *part, WoodsSaxon *ws, Parameters param) {
	double sigma_k = param.sigma_k, z = (double)param.z;
	#pragma omp parallel for
	for(int i = 0; i < part->protons + part->neutrons; i++)
		part->energy[i] = compute_energy(part, ws, sigma_k, z, i);
}

void compute_particle_densities(TestParticles *part, Parameters param) {
	int part_per_nucleon = param.part_per_nucleon, p = part->protons, n = part->neutrons, total = p + n;
	
	double sigma_r = param.sigma_r, exp_term = 1.0 / (2.0 * sigma_r * sigma_r);
	double term = (1.0 / part_per_nucleon) * (1.0 / (pow(2.0 * M_PI * sigma_r * sigma_r, 1.5)));
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		double r_i[3], r_j[3], diff[3];
		double fact, dist_squared, density_p = 0.0, density_n = 0.0;
		copy_particle_pos_to_vector(r_i, *part, i);
		
		for(int j = 0; j < total; j++) {
			copy_particle_pos_to_vector(r_j, *part, j);
			
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			fact = exp(-dist_squared * exp_term);
			if(j < p)
				density_p += fact;
			else
				density_n += fact;
		}
		part->density_p[i] = density_p;
		part->density_n[i] = density_n;
	}
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		part->density_p[i] *= term;
		part->density_n[i] *= term;
	}
}

void compute_volumetric_density(ScalarField *volume, ParticleCount part_count, World world_visual, World world_data, Parameters param, int type) {
	int world_size_visual = world_visual.n[0] * world_visual.n[1] * world_visual.n[2], start, end;
	int world_size_data = world_data.n[0] * world_data.n[1] * world_data.n[2];
	
	if(type == PROTONS) { start = 0; end = world_size_data; }
	else if(type == NEUTRONS) { start = world_size_data; end = 2 * world_size_data; }
	else { start = 0; end = 2 * world_size_data; }
	
	double sigma_r = param.sigma_r, exp_term = 1.0 / (2.0 * sigma_r * sigma_r);
	double term = (1.0 / param.part_per_nucleon) * (1.0 / (pow(2.0 * M_PI * sigma_r * sigma_r, 1.5)));
	#pragma omp parallel for
	for(int i = 0; i < world_size_visual; i++) {
		double r_i[3], r_j[3], diff[3];
		double fact, dist_squared, density = 0.0;
		world_pos_to_vector(r_i, world_visual, i);
		
		for(int j = start; j < end; j++) {
			int count = part_count.count[j];
			world_pos_to_vector(r_j, world_data, j % world_size_data);
			
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			fact = exp(-dist_squared * exp_term);
			density += count * fact;
		}
		volume->v[i] += density;
	}
	#pragma omp parallel for
	for(int i = 0; i < world_size_visual; i++)
		volume->v[i] *= term;
}

void scatter_particles(ParticleCount *part_count, TestParticles *part, World world) {
	double d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int x = world.n[0], y = world.n[1], z = world.n[2], world_size = x * y * z, total = part->protons + part->neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		double r_vec[3];
		copy_particle_pos_to_vector(r_vec, *part, i);
		
		int x_idx = (int)(x / 2.0 * (r_vec[0] / d_max_x + 1.0));
		int y_idx = (int)(y / 2.0 * (r_vec[1] / d_max_y + 1.0));
		int z_idx = (int)(z / 2.0 * (r_vec[2] / d_max_z + 1.0));
		
		if(x_idx < 0 || y_idx < 0 || z_idx < 0 || x_idx >= x || y_idx >= y || z_idx >= z)
			continue;
		
		int idx = x_idx * (y * z) + y_idx * z + z_idx;
		if(i >= part->protons)
			idx += world_size;
		
		#pragma omp atomic update
		part_count->count[idx] += 1;
	}
}

void generate_random_particles(TestParticles *part, double r_max) {
	int total = part->protons + part->neutrons, i = 0;
	while(i < total) {
		double r_new[3];
		random_vec(r_new, r_max);
		if(dot(r_new, r_new) < r_max * r_max) {
			copy_vector_to_particle_pos(part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < total) {
		double k_new[3];
		random_vec(k_new, K_MAX);
		if(dot(k_new, k_new) < K_MAX * K_MAX) {
			copy_vector_to_particle_vel(part, k_new, i);
			i++;
		}
	}
}

void generate_checking_particles(TestParticles *part, WoodsSaxon *ws, Parameters param, Fermi *fermi_levels) {
	double r_max = param.r_max, sigma_k = param.sigma_k, z = param.z, epsilon;
	int total = part->protons + part->neutrons, i = 0;
	while(i < total) {
		if(i < part->protons)
			epsilon = fermi_levels->epsilon_p;
		else
			epsilon = fermi_levels->epsilon_n;
		
		double r_new[3], k_new[3];
		random_vec(r_new, r_max);
		random_vec(k_new, K_MAX);
		copy_vector_to_particle_pos(part, r_new, i);
		copy_vector_to_particle_vel(part, k_new, i);
		double energy = compute_energy(part, ws, sigma_k, z, i);
		
		if(energy < epsilon) {
			mult_vec(r_new, r_new, -1.0);
			mult_vec(k_new, k_new, -1.0);
			copy_vector_to_particle_pos(part, r_new, i + 1);
			copy_vector_to_particle_vel(part, k_new, i + 1);
			i+=2;
		}
	}
}

void merge_volumetric_potentials(ScalarField *potential_a, ScalarField potential_b, World world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < world_size; i++) {
		potential_a->v[i] += potential_b.v[i];
	}
}

void copy_scalar_field(ScalarField *volume_a, ScalarField volume_b, World world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		volume_a->v[i] = volume_b.v[i];
}

void chi_squared(TestParticles part, WoodsSaxon *ws, Skyrme skm, int part_per_nucleon) {
	int total = part.protons + part.neutrons;
	
	double chi_squared_p = 0.0, chi_squared_n = 0.0;
	#pragma omp parallel for reduction(+:chi_squared_p, chi_squared_n)
	for(int i = 0; i < total; i++) {
		int type;
		WoodsSaxon ws_c;
		if(i >= part.protons) { type = NEUTRONS; ws_c = ws[1]; }
		else { type = PROTONS; ws_c = ws[0]; }
		
		double r_vec[3];
		copy_particle_pos_to_vector(r_vec, part, i);
		double density_p = part.density_p[i];
		double density_n = part.density_n[i];
		
		double r = magnitude(r_vec);
		double v_ws = woods_saxon_potential(ws_c, r);
		double v_skyrme = skyrme_potential(skm, density_p, density_n, type);
		double diff = v_ws - v_skyrme;
		if(type == PROTONS)
			chi_squared_p += diff * diff;
		else
			chi_squared_n += diff * diff;
	}
	chi_squared_n /= part_per_nucleon; chi_squared_p /= part_per_nucleon;
	printf("CHI SQUARED P %0.2lf CHI SQUARED N %0.2lf\n", chi_squared_p, chi_squared_n);
}

double mean_squared_radius(TestParticles part, int type) {
	int start, end;
	double part_num;
	if(type == PROTONS) { start = 0; end = part.protons; part_num = (double)part.protons;}
	else if(type == NEUTRONS) { start = part.protons; end = part.protons + part.neutrons; part_num = (double)part.neutrons; }
	
	double r_sqr = 0.0;
	#pragma omp parallel for reduction(+:r_sqr)
	for(int i = start; i < end; i++) {
		double r_vec[3], r2;
		
		copy_particle_pos_to_vector(r_vec, part, i);
		r2 = dot(r_vec, r_vec);
		r_sqr += r2;
	}
	return r_sqr / part_num;
}

double kinetic_energy() {
	double hc2 = H_BAR_C * H_BAR_C;
	double e_kin = hc2 / (2.0 * MC2);
	return e_kin;
}

double fluctuation_energy(double sigma_k) {
	double e_fluc = 3.0 * kinetic_energy() * sigma_k * sigma_k;
	return e_fluc;
}

double calc_sigma(double fwhm) {
	double t = 2.0 * sqrt(2.0 * log(2.0));
	double sigma = fwhm / t;
	return sigma;
}

void relax_woods_saxon(WoodsSaxon *ws, WoodsSaxon *ws_old, double coef) {
	for(int i = 0; i < 2; i++) {
		ws[i].V0 = coef * ws[i].V0 + (1.0 - coef) * ws_old[i].V0;
		ws[i].R12 = coef * ws[i].R12 + (1.0 - coef) * ws_old[i].R12;
		ws[i].a = coef * ws[i].a + (1.0 - coef) * ws_old[i].a;
	}
}
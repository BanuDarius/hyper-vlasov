#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"

double compute_energy(struct test_particles *part, struct woods_saxon *ws, double sigma_k, int z, int i) {
	double energy = 0.0, r_vec[3], k_vec[3], r, k;
	struct woods_saxon ws_c;
	if(i < part->protons)
		ws_c = ws[0];
	else
		ws_c = ws[1];
	copy_particle_pos_to_vector(r_vec, *part, i);
	copy_particle_vel_to_vector(k_vec, *part, i);
	
	r = magnitude(r_vec);
	k = magnitude(k_vec);
	
	energy += woods_saxon_potential(ws_c, r);
	energy += (k * k) * kinetic_energy();
	
	if(i < part->protons)
		energy += coulomb_potential(ws_c, (double)z, r);
	
	energy += fluctuation_energy(sigma_k);
	return energy;
}

void compute_particle_energies(struct test_particles *part, struct woods_saxon *ws, struct parameters param) {
	double sigma_k = param.sigma_k, z = param.z;
	for(int i = 0; i < part->protons + part->neutrons; i++)
		part->energy[i] = compute_energy(part, ws, sigma_k, z, i);
}

void compute_particle_densities(struct test_particles *part, struct parameters param) {
	int part_per_nucleon = param.test_part_per_nucleon, p = part->protons, n = part->neutrons, total = p + n;
	
	double sigma_r = param.sigma_r;
	double term = (1.0 / part_per_nucleon) * (1.0 / (pow(2.0 * M_PI * sigma_r, 1.5)));
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		double r_i[3], r_j[3], diff[3];
		double fact, dist_squared, density_p = 0.0, density_n = 0.0;
		copy_particle_pos_to_vector(r_i, *part, i);
		
		for(int j = 0; j < total; j++) {
			copy_particle_pos_to_vector(r_j, *part, j);
			
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			fact = exp(-dist_squared / (2.0 * sigma_r));
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

void compute_volumetric_density(struct volumetric_density *volume, struct particle_count part_count, struct world world_visual, struct world world_data, struct parameters param, int type) {
	int world_size_visual = world_visual.n[0] * world_visual.n[1] * world_visual.n[2], start, end;
	int world_size_data = world_data.n[0] * world_data.n[1] * world_data.n[2], part_per_nucleon = param.test_part_per_nucleon;
	
	double sigma_r = param.sigma_r;
	double term = (1.0 / part_per_nucleon) * (1.0 / (pow(2.0 * M_PI * sigma_r, 1.5)));
	if(type == PROTONS) { start = 0; end = world_size_data; }
	else if(type == NEUTRONS) { start = world_size_data; end = 2 * world_size_data; }
	else { start = 0; end = 2 * world_size_data; }
	
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
			fact = exp(-dist_squared / (2.0 * sigma_r));
			density += count * fact;
		}
		volume->density[i] += density;
	}
	#pragma omp parallel for
	for(int i = 0; i < world_size_visual; i++)
		volume->density[i] *= term;
}

void generate_random_particles(struct test_particles *part, double r_max) {
	double r_new[3], k_new[3];
	int total = part->protons + part->neutrons, i = 0;
	while(i < total) {
		random_vec(r_new, r_max);
		if(dot(r_new, r_new) < r_max * r_max) {
			copy_vector_to_particle_pos(part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < total) {
		random_vec(k_new, K_MAX);
		if(dot(k_new, k_new) < K_MAX * K_MAX) {
			copy_vector_to_particle_vel(part, k_new, i);
			i++;
		}
	}
}

void scatter_particles(struct particle_count *part_count, struct test_particles *part, struct world world) {
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

void generate_checking_particles(struct test_particles *part, struct woods_saxon *ws, struct parameters param, struct fermi *fermi_levels) {
	double r_new[3], k_new[3], energy;
	double r_max = param.r_max, sigma_k = param.sigma_k, z = param.z, epsilon;
	int total = part->protons + part->neutrons, i = 0;
	while(i < total) {
		if(i < part->protons)
			epsilon = fermi_levels->epsilon_p;
		else
			epsilon = fermi_levels->epsilon_n;
		random_vec(r_new, r_max);
		random_vec(k_new, K_MAX);
		copy_vector_to_particle_pos(part, r_new, i);
		copy_vector_to_particle_vel(part, k_new, i);
		energy = compute_energy(part, ws, sigma_k, z, i);
		
		if(energy < epsilon) {
			mult_vec(r_new, r_new, -1.0);
			mult_vec(k_new, k_new, -1.0);
			copy_vector_to_particle_pos(part, r_new, i + 1);
			copy_vector_to_particle_vel(part, k_new, i + 1);
			i+=2;
		}
	}
}

void chi_squared(struct test_particles *part, struct woods_saxon *ws, struct skyrme skm, int part_per_nucleon) {
	struct woods_saxon ws_c;
	int total = part->protons + part->neutrons, type;
	double chi_squared_p = 0.0, chi_squared_n = 0.0, r_vec[3], density_p, density_n, r, v_ws, v_skyrme, diff;
	for(int i = 0; i < total; i++) {
		if(i < part->protons) {
			type = PROTONS; ws_c = ws[0];
		}
		else {
			type = NEUTRONS; ws_c = ws[1];
		}
		copy_particle_pos_to_vector(r_vec, *part, i);
		density_p = part->density_p[i];
		density_n = part->density_n[i];
		
		r = magnitude(r_vec);
		v_ws = woods_saxon_potential(ws_c, r);
		v_skyrme = skyrme_potential(skm, density_p, density_n, type);
		diff = v_ws - v_skyrme;
		if(type == PROTONS)
			chi_squared_p += diff * diff;
		else
			chi_squared_n += diff * diff;
	}
	chi_squared_n /= part_per_nucleon; chi_squared_p /= part_per_nucleon;
	printf("CHI SQUARED P %0.2lf CHI SQUARED N %0.2lf\n", chi_squared_p, chi_squared_n);
}

void relax_woods_saxon(struct woods_saxon *ws, struct woods_saxon *ws_old, double coef) {
	for(int i = 0; i < 2; i++) {
		ws[i].V0 = coef * ws[i].V0 + (1.0 - coef) * ws_old[i].V0;
		ws[i].R12 = coef * ws[i].R12 + (1.0 - coef) * ws_old[i].R12;
		ws[i].a = coef * ws[i].a + (1.0 - coef) * ws_old[i].a;
	}
}

double kinetic_energy() {
	double hc2 = H_BAR_C * H_BAR_C;
	double e_kin = hc2 / (2 * MC2);
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tools.h"
#include "math_tools.h"
#include "sim_structs.h"

double rand_val(double min, double max) {
	double s = (double)rand() / (double)RAND_MAX;
	return min + s * (max - min);
}

void random_vec(double *v, double max) {
	for(int i = 0; i < 3; i++)
		v[i] = rand_val(-max, max);
}

double nuclear_radius(unsigned short a) {
	double radius = 1.5 * pow((double)a, 1.0 / 3.0);
	return radius;
}

int max_particles(double r_max, double k_max, int test_part_per_nucleon) {
	double t = r_max * k_max, ct = 2.0 * M_PI;
	double phase_space_volume = (16.0 / 9.0) * M_PI * M_PI * (t * t * t);
	int max = test_part_per_nucleon * (int)floor(phase_space_volume / (ct * ct * ct) + 0.5);
	return max;
}

double woods_saxon_potential(struct woods_saxon ws, double r) {
	double v = ws.V0 / (1.0 + exp((r - ws.R12) / ws.a));
	return v;
}

double skyrme_potential(struct skyrme skm, double rho) {
	double t = rho / RHO_0;
	double v = skm.A * t + skm.B * pow(t, skm.gamma);
	return v;
}

double coulomb_potential(struct woods_saxon ws, double z, double r) {
	double R12 = ws.R12, v;
	if(r <= ws.R12)
		v = 1.44 * (z - 1.0) / R12 * (1.5 - 0.5 * (r / R12) * (r / R12));
	else
		v = 1.44 * (z - 1.0) / r;
	return v;
}

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
		energy += coulomb_potential(ws_c, z, r);
	
	energy += fluctuation_energy(sigma_k);
	return energy;
}

void compute_particle_energies(struct test_particles *part, struct woods_saxon *ws, struct parameters param) {
	double sigma_k = param.sigma_k, z = param.z;
	for(int i = 0; i < part->protons + part->neutrons; i++)
		part->energy[i] = compute_energy(part, ws, sigma_k, z, i);
}

void compute_particle_densities(struct test_particles *part, double sigma_r, double part_per_nucleon) {
	int p = part->protons, n = part->neutrons, total = p + n;
	for(int i = 0; i < total; i++)
		part->density[i] = 1.0;
	
	double sigma_sqr_4 = 4.0 * sigma_r * sigma_r;
	double term = (1.0 / part_per_nucleon) * (1.0 / (pow(M_PI * sigma_sqr_4, 1.5)));
	double r_i[3], r_j[3], diff[3], dist_squared, fact;
	for(int i = 0; i < total; i++) {
		if(i % 1000 == 0)
			printf("%i\n", i);
		
		copy_particle_pos_to_vector(r_i, *part, i);
		
		for(int j = i + 1; j < total; j++) {
			copy_particle_pos_to_vector(r_j, *part, j);
			
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			fact = exp(-(dist_squared) / sigma_sqr_4);
			
			part->density[i] += fact;
			part->density[j] += fact;
		}
	}
	for(int i = 0; i < total; i++)
		part->density[i] *= term;
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

void scatter_particles(struct test_particles *part, struct world world, int num) {
	
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

void chi_squared(struct test_particles *part, struct woods_saxon ws, struct skyrme skm) {
	int total = part->protons + part->neutrons;
	double chi_squared = 0.0;
	for(int i = 0; i < total; i++) {
		double r_vec[3], density;
		copy_particle_pos_to_vector(r_vec, *part, i);
		density = part->density[i];
		
		double r = magnitude(r_vec);
		double v_ws = woods_saxon_potential(ws, r);
		double v_skyrme = skyrme_potential(skm, density);
		double diff = v_ws - v_skyrme;
		chi_squared += diff * diff;
	}
	printf("CHI SQUARED %lf\n", chi_squared);
}

void copy_accepted_particles(struct test_particles *part_accepted, struct test_particles *part_all, double epsilon, int num) {
	int idx = 0;
	double r[3], r_sym[3], k[3], k_sym[3];
	for(int i = 0; i < num; i++) {
		if(part_all->energy[i] < epsilon) {
			copy_particle_pos_to_vector(r, *part_all, i);
			copy_particle_vel_to_vector(k, *part_all, i);
			mult_vec(r_sym, r, -1.0);
			mult_vec(k_sym, k, -1.0);
			
			copy_vector_to_particle_pos(part_accepted, r, idx);
			copy_vector_to_particle_vel(part_accepted, k, idx);
			idx++;
			copy_vector_to_particle_pos(part_accepted, r_sym, idx);
			copy_vector_to_particle_vel(part_accepted, k_sym, idx);
			idx++;
		}
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

void copy_particle(struct test_particles *part_a, struct test_particles *part_b, int idx, int i) {
	part_a->x[idx] = part_b->x[i];
	part_a->y[idx] = part_b->y[i];
	part_a->z[idx] = part_b->z[i];
	part_a->kx[idx] = part_b->kx[i];
	part_a->ky[idx] = part_b->ky[i];
	part_a->kz[idx] = part_b->kz[i];
	part_a->energy[idx] = part_b->energy[i];
	part_a->density[idx] = part_b->density[i];
}

void copy_particle_pos_to_vector(double *v, struct test_particles part, int i) {
	v[0] = part.x[i];
	v[1] = part.y[i];
	v[2] = part.z[i];
}

void copy_particle_vel_to_vector(double *v, struct test_particles part, int i) {
	v[0] = part.kx[i];
	v[1] = part.ky[i];
	v[2] = part.kz[i];
}

void copy_vector_to_particle_pos(struct test_particles *part, double *v, int i) {
	part->x[i] = v[0];
	part->y[i] = v[1];
	part->z[i] = v[2];
}

void copy_vector_to_particle_vel(struct test_particles *part, double *v, int i) {
	part->kx[i] = v[0];
	part->ky[i] = v[1];
	part->kz[i] = v[2];
}
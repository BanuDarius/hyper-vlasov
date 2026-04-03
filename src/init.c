#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "init.h"
#include "tools.h"
#include "math_tools.h"
#include "sim_structs.h"

#define K_MAX 1.5

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon, double sigma_k, double sigma_r) {
	param->z = z;
	param->n = n;
	param->sigma_k = sigma_k;
	param->sigma_r = sigma_r;
	
	param->test_part_per_nucleon = test_part_per_nucleon;
	param->r_max = nuclear_radius(z + n);
	
	param->max_test_part = max_particles(param->r_max, K_MAX, param->test_part_per_nucleon);
}

void set_woods_saxon(struct woods_saxon *ws, struct parameters param, double V0, double a) {
	ws->V0 = V0;
	ws->a = a;
	ws->R12 = 0.8 * param.r_max;
}

void set_skyrme(struct skyrme *skm, double A, double B, double gamma) {
	skm->A = A;
	skm->B = B;
	skm->gamma = gamma;
}

void set_fermi_levels(struct fermi *fermi, double epsilon_p, double epsilon_n) {
	fermi->epsilon_p = epsilon_p;
	fermi->epsilon_n = epsilon_n;
}

void create_particles(struct test_particles *part, int num) {
	part->x = malloc(num * sizeof(double));
	part->y = malloc(num * sizeof(double));
	part->z = malloc(num * sizeof(double));
	part->kx = malloc(num * sizeof(double));
	part->ky = malloc(num * sizeof(double));
	part->kz = malloc(num * sizeof(double));
	part->energy = malloc(num * sizeof(double));
	part->density = malloc(num * sizeof(double));
}

double compute_energy(struct test_particles *part, struct woods_saxon ws, double sigma_k, int z, int type, int i) {
	double energy = 0.0, r_vec[3], k_vec[3], r, k;
	copy_particle_pos(*part, r_vec, i);
	copy_particle_vel(*part, k_vec, i);
	
	r = magnitude(r_vec);
	k = magnitude(k_vec);
	
	energy += woods_saxon_potential(ws, r);
	energy += (k * k) * kinetic_energy();
	
	if(type == PROTONS)
		energy += coulomb_potential(ws, z, r);
	
	energy += fluctuation_energy(sigma_k);
	return energy;
}

void compute_particle_energies(struct test_particles *part, struct woods_saxon ws, double sigma_k, int z, int type, int num) {
	for(int i = 0; i < num; i++)
		part->energy[i] = compute_energy(part, ws, sigma_k, z, type, i);
}

void compute_particle_densities(struct test_particles *part, double sigma_r, int num) {
	double one = 1.0;
	memcpy(part->density, &one, num * sizeof(double));
}

void generate_random_particles(struct test_particles *part, double r_max, int num) {
	double r_new[3], k_new[3];
	int i = 0;
	while(i < num) {
		random_vec(r_new, r_max);
		if(dot(r_new, r_new) < r_max * r_max) {
			copy_vector_to_particle_pos(*part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < num) {
		random_vec(k_new, K_MAX);
		if(dot(k_new, k_new) < K_MAX * K_MAX) {
			copy_vector_to_particle_vel(*part, k_new, i);
			i++;
		}
	}
}

void generate_checking_particles(struct test_particles *part, struct woods_saxon ws, struct parameters param, double epsilon, int type, int num) {
	double r_new[3], k_new[3], energy;
	double r_max = param.r_max, sigma_k = param.sigma_k, z = param.z;
	int i = 0;
	while(i < num) {
		random_vec(r_new, r_max);
		random_vec(k_new, K_MAX);
		copy_vector_to_particle_pos(*part, r_new, i);
		copy_vector_to_particle_vel(*part, k_new, i);
		energy = compute_energy(part, ws, sigma_k, z, type, i);
		if(energy < epsilon) {
			mult_vec(r_new, r_new, -1.0);
			mult_vec(k_new, k_new, -1.0);
			copy_vector_to_particle_pos(*part, r_new, i + 1);
			copy_vector_to_particle_vel(*part, k_new, i + 1);
			i+=2;
		}
	}
}

void initialize_particles(struct test_particles *part_p, struct test_particles *part_n, struct parameters param, struct woods_saxon ws, struct skyrme skm, struct fermi *fermi_levels) {
	double r_max = param.r_max, sigma_k = param.sigma_k, total_delta_epsilon;
	int max_part = param.max_test_part, z = param.z, n = param.n, part_per_nucleon = param.test_part_per_nucleon, it = 0;
	
	create_particles(part_p,  param.max_test_part);
	create_particles(part_n,  param.max_test_part);
	generate_random_particles(part_p, r_max, max_part);
	generate_random_particles(part_n, r_max, max_part);
	compute_particle_energies(part_p, ws, sigma_k, z, PROTONS, max_part);
	compute_particle_energies(part_n, ws, sigma_k, z, NEUTRONS, max_part);
	
	do {
		int check_less_p = 0, check_equal_p = 0, check_more_p = 0;
		int check_less_n = 0, check_equal_n = 0, check_more_n = 0;
		for(int i = 0; i < max_part; i++) {
			if(part_p->energy[i] < fermi_levels->epsilon_p)
				check_equal_p += 2;
			if(part_p->energy[i] < fermi_levels->epsilon_p + 0.5)
				check_more_p += 2;
			if(part_p->energy[i] < fermi_levels->epsilon_p - 0.5)
				check_less_p += 2;
			
			if(part_n->energy[i] < fermi_levels->epsilon_n)
				check_equal_n += 2;
			if(part_n->energy[i] < fermi_levels->epsilon_n + 0.5)
				check_more_n += 2;
			if(part_n->energy[i] < fermi_levels->epsilon_n - 0.5)
				check_less_n += 2;
		}
		printf("%i %i %i\n", check_less_p, check_equal_p, check_more_p);
		printf("%i %i %i\n", check_less_n, check_equal_n, check_more_n);
		
		double delta_part_n = n * part_per_nucleon - check_equal_n;
		double delta_part_p = z * part_per_nucleon - check_equal_p;
		double delta_epsilon_n = 0.5 * delta_part_n / (check_more_n - check_less_n);
		double delta_epsilon_p = 0.5 * delta_part_p / (check_more_p - check_less_p);
		
		if(fabs(delta_epsilon_p) > 0.5) delta_epsilon_p *= 0.6;
		if(fabs(delta_epsilon_n) > 0.5) delta_epsilon_n *= 0.6;
		
		fermi_levels->epsilon_p += delta_epsilon_p;
		fermi_levels->epsilon_n += delta_epsilon_n;
		
		total_delta_epsilon = fabs(delta_epsilon_n) + fabs(delta_epsilon_p);
		
		printf("%lf\n", total_delta_epsilon);
		it++;
	} while(total_delta_epsilon > DELTA_EPSILON_TOLERANCE && it < MAX_ITERATIONS);
	
	int total_p = z * part_per_nucleon, total_n = n * part_per_nucleon;
	free_particles(part_p);
	free_particles(part_n);
	create_particles(part_p, total_p);
	create_particles(part_n, total_n);
	
	generate_checking_particles(part_p, ws, param, fermi_levels->epsilon_p, PROTONS, total_p);
	generate_checking_particles(part_n, ws, param, fermi_levels->epsilon_n, NEUTRONS, total_n);
	compute_particle_energies(part_p, ws, sigma_k, z, PROTONS, total_p);
	compute_particle_energies(part_n, ws, sigma_k, z, NEUTRONS, total_n);
}

void output_centroids(FILE *out, struct test_particles part, int num) {
	for(int i = 0; i < num; i++) {
		fwrite(&part.x[i], sizeof(double), 1, out);
		fwrite(&part.y[i], sizeof(double), 1, out);
		fwrite(&part.z[i], sizeof(double), 1, out);
		fwrite(&part.kx[i], sizeof(double), 1, out);
		fwrite(&part.ky[i], sizeof(double), 1, out);
		fwrite(&part.kz[i], sizeof(double), 1, out);
	}
}

void free_particles(struct test_particles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->kx); free(part->ky); free(part->kz);
	free(part->energy); free(part->density);
}
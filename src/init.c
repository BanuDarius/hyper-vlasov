#include <omp.h>
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

void set_woods_saxon(struct woods_saxon *ws, double V0, double R12, double a) {
	ws->a = a;
	ws->V0 = V0;
	ws->R12 = R12;
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
	copy_particle_pos_to_vector(r_vec, *part, i);
	copy_particle_vel_to_vector(k_vec, *part, i);
	
	r = magnitude(r_vec);
	k = magnitude(k_vec);
	
	energy += woods_saxon_potential(ws, r);
	energy += (k * k) * kinetic_energy();
	
	if(type == PROTONS)
		energy += coulomb_potential(ws, z, r);
	
	energy += fluctuation_energy(sigma_k);
	return energy;
}

void compute_particle_energies(struct test_particles *part, struct woods_saxon ws, struct parameters param, int type, int num) {
	double sigma_k = param.sigma_k, z = param.z;
	for(int i = 0; i < num; i++)
		part->energy[i] = compute_energy(part, ws, sigma_k, z, type, i);
}

void compute_particle_densities(struct test_particles *part_p, struct test_particles *part_n, double sigma_r, double part_per_nucleon, int total_p, int total_n) {
	int total = total_p + total_n;
	for(int i = 0; i < total_p; i++)
		part_p->density[i] = 1.0;
	for(int i = 0; i < total_n; i++)
		part_n->density[i] = 1.0;
	
	double sigma_sqr_4 = 4.0 * sigma_r * sigma_r;
	double term = (1.0 / part_per_nucleon) * (1.0 / (pow(M_PI * sigma_sqr_4, 1.5)));
	double r_i[3], r_j[3], diff[3], dist_squared, fact;
	for(int i = 0; i < total; i++) {
		if(i % 1000 == 0)
			printf("%i\n", i);
		if(i < total_p)
			copy_particle_pos_to_vector(r_i, *part_p, i);
		else
			copy_particle_pos_to_vector(r_i, *part_n, i - total_p);
		
		for(int j = i + 1; j < total_p + total_n; j++) {
			if(j < total_p)
				copy_particle_pos_to_vector(r_j, *part_p, j);
			else
				copy_particle_pos_to_vector(r_j, *part_n, j - total_p);
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			fact = exp(-(dist_squared) / sigma_sqr_4);
			
			if(i < total_p)
				part_p->density[i] += fact;
			else
				part_n->density[i - total_p] += fact;
			if(j < total_p)
				part_p->density[j] += fact;
			else
				part_n->density[j - total_p] += fact;
		}
	}
	for(int i = 0; i < total_p; i++)
		part_p->density[i] *= term;
	for(int i = 0; i < total_n; i++)
		part_n->density[i] *= term;
}

void generate_random_particles(struct test_particles *part, double r_max, int num) {
	double r_new[3], k_new[3];
	int i = 0;
	while(i < num) {
		random_vec(r_new, r_max);
		if(dot(r_new, r_new) < r_max * r_max) {
			copy_vector_to_particle_pos(part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < num) {
		random_vec(k_new, K_MAX);
		if(dot(k_new, k_new) < K_MAX * K_MAX) {
			copy_vector_to_particle_vel(part, k_new, i);
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
		copy_vector_to_particle_pos(part, r_new, i);
		copy_vector_to_particle_vel(part, k_new, i);
		energy = compute_energy(part, ws, sigma_k, z, type, i);
		if(energy < epsilon) {
			mult_vec(r_new, r_new, -1.0);
			mult_vec(k_new, k_new, -1.0);
			copy_vector_to_particle_pos(part, r_new, i + 1);
			copy_vector_to_particle_vel(part, k_new, i + 1);
			i+=2;
		}
	}
}

void initialize_particles(struct test_particles *part_p, struct test_particles *part_n, struct parameters param, struct woods_saxon *ws, struct skyrme skm, struct fermi *fermi_levels) {
	double r_max = param.r_max, sigma_k = param.sigma_k, sigma_r = param.sigma_r, total_delta_epsilon;
	int max_part = param.max_test_part, z = param.z, n = param.n, part_per_nucleon = param.test_part_per_nucleon, it = 0;
	double delta_part_p, delta_part_n, delta_epsilon_p, delta_epsilon_n;
	
	create_particles(part_p, param.max_test_part);
	create_particles(part_n, param.max_test_part);
	generate_random_particles(part_p, r_max, max_part);
	generate_random_particles(part_n, r_max, max_part);
	
	struct woods_saxon ws_p = *ws, ws_n = *ws;
	do {
		compute_particle_energies(part_p, ws_p, param, PROTONS, max_part);
		compute_particle_energies(part_n, ws_n, param, NEUTRONS, max_part);
		
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
		
		/*struct test_particles part_p_accepted, part_n_accepted;
		create_particles(&part_p_accepted, check_equal_p);
		create_particles(&part_n_accepted, check_equal_n);
		
		copy_accepted_particles(&part_p_accepted, part_p, fermi_levels->epsilon_p, max_part);
		copy_accepted_particles(&part_n_accepted, part_n, fermi_levels->epsilon_n, max_part);
		
		compute_particle_densities(&part_p_accepted, &part_n_accepted, sigma_r, (double)part_per_nucleon, check_equal_p, check_equal_n);

		fit_woods_saxon_param(&ws_p, &part_p_accepted, skm, check_equal_p);
		fit_woods_saxon_param(&ws_n, &part_n_accepted, skm, check_equal_n);;
		
		free_particles(&part_p_accepted);
		free_particles(&part_n_accepted);*/
		
		delta_part_n = n * part_per_nucleon - check_equal_n;
		delta_part_p = z * part_per_nucleon - check_equal_p;
		delta_epsilon_n = 0.5 * delta_part_n / (check_more_n - check_less_n);
		delta_epsilon_p = 0.5 * delta_part_p / (check_more_p - check_less_p);
		
		if(fabs(delta_epsilon_p) > 0.5) delta_epsilon_p *= 0.6;
		if(fabs(delta_epsilon_n) > 0.5) delta_epsilon_n *= 0.6;
		
		fermi_levels->epsilon_p += delta_epsilon_p;
		fermi_levels->epsilon_n += delta_epsilon_n;
		total_delta_epsilon = fabs(delta_epsilon_n) + fabs(delta_epsilon_p);
		
		printf("ITERATION = %i TOTAL DELTA EPSILON = %lf\n", it,  total_delta_epsilon);
		printf("%i %i %i\n", check_less_p, check_equal_p, check_more_p);
		printf("%i %i %i\n", check_less_n, check_equal_n, check_more_n);
		
		it++;
	} while(total_delta_epsilon > DELTA_EPSILON_TOLERANCE && it < MAX_ITERATIONS);
	
	int total_p = z * part_per_nucleon, total_n = n * part_per_nucleon;
	free_particles(part_p);
	free_particles(part_n);
	create_particles(part_p, total_p);
	create_particles(part_n, total_n);
	
	generate_checking_particles(part_p, ws_p, param, fermi_levels->epsilon_p, PROTONS, total_p);
	generate_checking_particles(part_n, ws_n, param, fermi_levels->epsilon_n, NEUTRONS, total_n);
	
	compute_particle_energies(part_p, ws_p, param, PROTONS, total_p);
	compute_particle_energies(part_n, ws_n, param, NEUTRONS, total_n);
	
	printf("WS PARAM %lf %lf %lf\n", ws_p.V0, ws_p.R12, ws_p.a);
	printf("WS PARAM %lf %lf %lf\n", ws_n.V0, ws_n.R12, ws_n.a);
	
	//fit_woods_saxon_param(&ws_p, part_p, skm, total_p);
	//fit_woods_saxon_param(&ws_n, part_n, skm, total_n);
	
	printf("WS PARAM %lf %lf %lf\n", ws_p.V0, ws_p.R12, ws_p.a);
	printf("WS PARAM %lf %lf %lf\n", ws_n.V0, ws_n.R12, ws_n.a);
	
	compute_particle_densities(part_p, part_n, sigma_r, (double)part_per_nucleon, total_p, total_n);
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

void fit_woods_saxon_param(struct woods_saxon *ws, struct test_particles *part, struct skyrme skm, int total_part) {
	double p[3] = {ws->V0, ws->R12, ws->a};
	double h[3] = {0.01, 0.01, 0.01};
	struct woods_saxon ws_new;
	int max_it = 15;
	double tol = 0.01;
	
	for(int it = 0; it < max_it; it++) {
		double A[3][3] = {{0.0}}, b[3] = {0.0};
		double chi_squared = 0.0;
		
		for(int i = 0; i < total_part; i++) {
			double r_vec[3], density, r;
			
			copy_particle_pos_to_vector(r_vec, *part, i);
			density = part->density[i];
			r = magnitude(r_vec);
			
			set_woods_saxon(&ws_new, p[0], p[1], p[2]);
			double v_skyrme = skyrme_potential(skm, density);
			double v_woods_saxon = woods_saxon_potential(ws_new, r);
			
			double t = v_woods_saxon - v_skyrme;
			chi_squared += t * t;
			
			double J[3];
			for(int x = 0; x < 3; x++) {
				struct woods_saxon ws_plus = ws_new;
				struct woods_saxon ws_minus = ws_new;
				if(x == 0) { ws_plus.V0 += h[0]; ws_minus.V0 -= h[0]; }
				if(x == 1) { ws_plus.R12 += h[1]; ws_minus.R12 -= h[1]; }
				if(x == 2) { ws_plus.a += h[2]; ws_minus.a -= h[2]; }
				
				double v_plus = woods_saxon_potential(ws_plus, r);
				double v_minus = woods_saxon_potential(ws_minus, r);
				J[x] = (v_plus - v_minus) / (2.0 * h[x]);
			}
			
			for(int row = 0; row < 3; row++) {
				b[row] -= J[row] * t;
				for(int col = 0; col < 3; col++)
					A[row][col] += J[row] * J[col];
			}
		}
		double delta[3] = {0.0};
		solve_3x3(A, b, delta);
		add_vec(p, p, delta);
		
		if(magnitude(delta) < tol)
			break;
	}
	
	ws->V0 = p[0]; ws->R12 = p[1]; ws->a = p[2];
}

void output_centroids(FILE *out, struct test_particles part, int num) {
	for(int i = 0; i < num; i++) {
		fwrite(&part.x[i], sizeof(double), 1, out);
		fwrite(&part.y[i], sizeof(double), 1, out);
		fwrite(&part.z[i], sizeof(double), 1, out);
		fwrite(&part.kx[i], sizeof(double), 1, out);
		fwrite(&part.ky[i], sizeof(double), 1, out);
		fwrite(&part.kz[i], sizeof(double), 1, out);
		fwrite(&part.energy[i], sizeof(double), 1, out);
		fwrite(&part.density[i], sizeof(double), 1, out);
	}
}

void free_particles(struct test_particles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->kx); free(part->ky); free(part->kz);
	free(part->energy); free(part->density);
}
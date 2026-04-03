#include <stdio.h>
#include <stdlib.h>

#include "tools.h"
#include "math_tools.h"
#include "sim_structs.h"

#define K_MAX 1.5

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon) {
	param->z = z;
	param->n = n;
	param->test_part_per_nucleon = test_part_per_nucleon;
	param->total_test_part = (z + n) * test_part_per_nucleon;
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
}

void generate_random_particles(struct test_particles *part, double r_max, int num) {
	int i = 0;
	while(i < num) {
		double r_new[3];
		for(int j = 0; j < 3; j++)
			r_new[j] = rand_val(-r_max, r_max);
		if(dot(r_new, r_new) < r_max * r_max) {
			copy_vector_to_particle_pos(part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < num) {
		double p_new[3];
		for(int j = 0; j < 3; j++)
			p_new[j] = rand_val(-K_MAX, K_MAX);
		if(dot(p_new, p_new) < K_MAX * K_MAX) {
			copy_vector_to_particle_vel(part, p_new, i);
			i++;
		}
	}
}

void initialize_particles(struct test_particles *part_p, struct test_particles *part_n, struct parameters param, struct woods_saxon ws, struct skyrme skm) {
	generate_random_particles(part_p, param.r_max, param.max_test_part);
	generate_random_particles(part_n, param.r_max, param.max_test_part);
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
}
#include <stdlib.h>
#include <stdio.h>

#include "tools.h"
#include "math_tools.h"
#include "sim_structs.h"

#define MAX_P 1.5

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon) {
	param->z = z;
	param->n = n;
	param->test_part_per_nucleon = test_part_per_nucleon;
	param->total_test_part = (z + n) * test_part_per_nucleon;
	param->r = nuclear_radius(z + n);
}

void create_particles(struct test_particles *part, struct parameters param) {
	int num = param.total_test_part;
	part->x = malloc(num * sizeof(double));
	part->y = malloc(num * sizeof(double));
	part->z = malloc(num * sizeof(double));
	part->px = malloc(num * sizeof(double));
	part->py = malloc(num * sizeof(double));
	part->pz = malloc(num * sizeof(double));
}


void generate_random_particles(struct test_particles *part, struct parameters param) {
	double r = param.r;
	double num = param.total_test_part;
	int i = 0;
	while(i < num) {
		double r_new[3];
		for(int j = 0; j < 3; j++)
			r_new[j] = rand_val(-r, r);
		if(dot(r_new, r_new) < r * r) {
			copy_vector_to_particle_pos(part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < num) {
		double p_new[3];
		for(int j = 0; j < 3; j++)
			p_new[j] = rand_val(-MAX_P, MAX_P);
		if(dot(p_new, p_new) < MAX_P * MAX_P) {
			copy_vector_to_particle_vel(part, p_new, i);
			i++;
		}
	}
}

void output_centroids(FILE *out, struct test_particles *part, struct parameters param) {
	int num = param.total_test_part;
	for(int i = 0; i < num; i++) {
		fwrite(&part->x[i], sizeof(double), 1, out);
		fwrite(&part->y[i], sizeof(double), 1, out);
		fwrite(&part->z[i], sizeof(double), 1, out);
		fwrite(&part->px[i], sizeof(double), 1, out);
		fwrite(&part->py[i], sizeof(double), 1, out);
		fwrite(&part->pz[i], sizeof(double), 1, out);
	}
}

void output_centroid_positions(FILE *out, struct test_particles *part, struct parameters param) {
	int num = param.total_test_part;
	for(int i = 0; i < num; i++) {
		fwrite(&part->x[i], sizeof(double), 1, out);
		fwrite(&part->y[i], sizeof(double), 1, out);
		fwrite(&part->z[i], sizeof(double), 1, out);
	} //This is bad for the structure of arrays data type but it works for now
}

void output_centroid_velocities(FILE *out, struct test_particles *part, struct parameters param) {
	int num = param.total_test_part;
	for(int i = 0; i < num; i++) {
		fwrite(&part->px[i], sizeof(double), 1, out);
		fwrite(&part->py[i], sizeof(double), 1, out);
		fwrite(&part->pz[i], sizeof(double), 1, out);
	} //This is also bad
}

void free_particles(struct test_particles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->px); free(part->py); free(part->pz);
	free(part);
}
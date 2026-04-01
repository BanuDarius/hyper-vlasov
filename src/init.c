#include <stdlib.h>
#include <stdio.h>

#include "tools.h"
#include "math_tools.h"
#include "sim_structs.h"

void create_particles(struct test_particles *part, struct parameters param) {
	printf("here");
	int num = param.num_test_part;
	part->x = malloc(num * sizeof(double));
	part->y = malloc(num * sizeof(double));
	part->z = malloc(num * sizeof(double));
	part->px = malloc(num * sizeof(double));
	part->py = malloc(num * sizeof(double));
	part->pz = malloc(num * sizeof(double));
}


void generate_random_particles(struct test_particles *part, struct parameters param) {
	double r = param.r;
	double num = param.num_test_part;
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
}

void free_particles(struct test_particles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->px); free(part->py); free(part->pz);
	free(part);
}
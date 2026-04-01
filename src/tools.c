#include <math.h>
#include <stdlib.h>

#include "sim_structs.h"

double rand_val(double min, double max) {
	double s = (double)rand() / (double)RAND_MAX;
	return min + s * (max - min);
}

double nuclear_radius(unsigned short a) {
	double radius = 1.5 * pow((double)a, 1.0 / 3.0);
	return radius;
}

void copy_vector_to_particle_pos(struct test_particles *part, double *v, int i) {
	part->x[i] = v[0];
	part->y[i] = v[1];
	part->z[i] = v[2];
}

void copy_vector_to_particle_vel(struct test_particles *part, double *v, int i) {
	part->px[i] = v[0];
	part->py[i] = v[1];
	part->pz[i] = v[2];
}
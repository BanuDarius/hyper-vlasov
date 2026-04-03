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

int max_particles(double r_max, double k_max, int total_test_part) {
	double t = r_max * k_max, ct = 2.0 * M_PI;
	double phase_space_volume = 64.0 * t * t * t;
	int max = total_test_part * (int)floor(phase_space_volume / (ct * ct * ct) + 0.5);
	return max;
}

double woods_saxon_potential(struct woods_saxon ws, double r) {
	double v = ws.V0 / (1.0 + exp((r - ws.R12) / ws.a));
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
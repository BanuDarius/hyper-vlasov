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

int max_particles(double r_max, double k_max, int test_part_per_nucleon) {
	double t = r_max * k_max, ct = 2.0 * M_PI;
	double phase_space_volume = 64.0 * t * t * t;
	int max = test_part_per_nucleon * (int)floor(phase_space_volume / (ct * ct * ct) + 0.5);
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

double kinetic_energy() {
	double hc2 = H_BAR_C * H_BAR_C;
	double e_kin = hc2 / (2 * MC2);
	return e_kin;
}

double fluctuation_energy(double sigma_k) {
	double e_fluc = 3.0 * kinetic_energy() * sigma_k * sigma_k;
	return e_fluc;
}

double calc_sigma_k(double k_fwhm) {
	double t = 2.0 * sqrt(2.0 * log(2.0));
	double sigma_k = k_fwhm / t;
	return sigma_k;
}

void copy_particle_pos(struct test_particles part, double *v, int i) {
	v[0] = part.x[i];
	v[1] = part.y[i];
	v[2] = part.z[i];
}

void copy_particle_vel(struct test_particles part, double *v, int i) {
	v[0] = part.kx[i];
	v[1] = part.ky[i];
	v[2] = part.kz[i];
}

void copy_vector_to_particle_pos(struct test_particles part, double *v, int i) {
	part.x[i] = v[0];
	part.y[i] = v[1];
	part.z[i] = v[2];
}

void copy_vector_to_particle_vel(struct test_particles part, double *v, int i) {
	part.kx[i] = v[0];
	part.ky[i] = v[1];
	part.kz[i] = v[2];
}
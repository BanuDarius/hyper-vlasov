#include <math.h>
#include <stdlib.h>

#include "tools.h"
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

void solve_3x3(double A[3][3], double b[3], double x[3]) {
	double detA = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
	
	double invA[3][3];
	invA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]) / detA;
	invA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]) / detA;
	invA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]) / detA;
	invA[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]) / detA;
	invA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]) / detA;
	invA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]) / detA;
	invA[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]) / detA;
	invA[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]) / detA;
	invA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]) / detA;
	
	for (int i = 0; i < 3; i++) {
		x[i] = 0.0;
		for (int j = 0; j < 3; j++) {
			x[i] += invA[i][j] * b[j];
		}
	}
}
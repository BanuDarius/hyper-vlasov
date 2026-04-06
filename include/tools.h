#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>

#include "sim_structs.h"

double compute_energy(struct test_particles *part, struct woods_saxon *ws, double sigma_k, int z, int i);
void compute_particle_energies(struct test_particles *part, struct woods_saxon *ws, struct parameters param);
void compute_particle_densities(struct test_particles *part, struct parameters param);
void compute_volumetric_density(struct volumetric_density *volume, struct particle_count part_count, struct world world_visual, struct world world_data, struct parameters param, int type);
void apply_constant_to_density(struct volumetric_density *volume, struct world world, struct parameters param);
void generate_random_particles(struct test_particles *part, double r_max);
void scatter_particles(struct particle_count *part_count, struct test_particles *part, struct world world);
void generate_checking_particles(struct test_particles *part, struct woods_saxon *ws, struct parameters param, struct fermi *fermi_levels);
void chi_squared(struct test_particles *part, struct woods_saxon *ws, struct skyrme skm, int part_per_nucleon);
void relax_woods_saxon(struct woods_saxon *ws, struct woods_saxon *ws_old, double coef);
double kinetic_energy();
double fluctuation_energy(double sigma_k);
double calc_sigma(double fwhm);

static inline double rand_val(double min, double max) {
	double s = (double)rand() / (double)RAND_MAX;
	return min + s * (max - min);
}

static inline void random_vec(double *v, double max) {
	for(int i = 0; i < 3; i++)
		v[i] = rand_val(-max, max);
}

static inline void world_pos_to_vector(double *v, struct world world, int idx) {
	int x = world.n[0], y = world.n[1], z = world.n[2];
	int i = idx / (y * z), j = (idx / z) % y, k = idx % z;
	v[0] = world.d_max[0] * (2.0 * (double)i / x - 1.0);
	v[1] = world.d_max[1] * (2.0 * (double)j / y - 1.0);
	v[2] = world.d_max[2] * (2.0 * (double)k / z - 1.0);
}

static inline void copy_particle_pos_to_vector(double *v, struct test_particles part, int i) {
	v[0] = part.x[i];
	v[1] = part.y[i];
	v[2] = part.z[i];
}

static inline void copy_particle_vel_to_vector(double *v, struct test_particles part, int i) {
	v[0] = part.kx[i];
	v[1] = part.ky[i];
	v[2] = part.kz[i];
}

static inline void copy_vector_to_particle_pos(struct test_particles *part, double *v, int i) {
	part->x[i] = v[0];
	part->y[i] = v[1];
	part->z[i] = v[2];
}

static inline void copy_vector_to_particle_vel(struct test_particles *part, double *v, int i) {
	part->kx[i] = v[0];
	part->ky[i] = v[1];
	part->kz[i] = v[2];
}

#endif
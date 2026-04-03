#ifndef TOOLS_H
#define TOOLS_H

#include "sim_structs.h"

double rand_val(double min, double max);
void random_vec(double *v, double max);
double nuclear_radius(unsigned short a);
int max_particles(double r_max, double k_max, int total_test_part);
double woods_saxon_potential(struct woods_saxon ws, double r);
double skyrme_potential(struct skyrme skm, double rho);
double coulomb_potential(struct woods_saxon ws, double z, double r);
double kinetic_energy();
double fluctuation_energy(double sigma_k);
double calc_sigma(double fwhm);
void copy_particle(struct test_particles *part_a, struct test_particles *part_b, int idx, int i);
void copy_particle_pos_to_vector(double *v, struct test_particles part, int i);
void copy_particle_vel_to_vector(double *v, struct test_particles part, int i);
void copy_vector_to_particle_pos(struct test_particles *part, double *v, int i);
void copy_vector_to_particle_vel(struct test_particles *part, double *v, int i);
void solve_3x3(double A[3][3], double b[3], double x[3]);

#endif
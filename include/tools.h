#ifndef TOOLS_H
#define TOOLS_H

#include "sim_structs.h"

double rand_val(double min, double max);
double nuclear_radius(unsigned short a);
int max_particles(double r_max, double k_max, int total_test_part);
double woods_saxon_potential(struct woods_saxon ws, double r);
double coulomb_potential(struct woods_saxon ws, double z, double r);
double kinetic_energy();
double fluctuation_energy(double sigma_k);
double calc_sigma_k(double k_fwhm);
void copy_particle_pos(struct test_particles part, double *v, int i);
void copy_particle_vel(struct test_particles part, double *v, int i);
void copy_vector_to_particle_pos(struct test_particles part, double *v, int i);
void copy_vector_to_particle_vel(struct test_particles part, double *v, int i);

#endif
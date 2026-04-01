#ifndef TOOLS_H
#define TOOLS_H

#include "sim_structs.h"

double rand_val(double min, double max);
double nuclear_radius(unsigned short a);
void copy_vector_to_particle_pos(struct test_particles *part, double *v, int i);
void copy_vector_to_particle_vel(struct test_particles *part, double *v, int i);

#endif
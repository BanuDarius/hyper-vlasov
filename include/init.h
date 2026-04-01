#ifndef INIT_H
#define INIT_H

#include "sim_structs.h"

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon);
void create_particles(struct test_particles *part, struct parameters param);
void generate_random_particles(struct test_particles *part, struct parameters param);
void output_centroids(FILE *out, struct test_particles *part, struct parameters param);
void output_centroid_positions(FILE *out, struct test_particles *part, struct parameters param);
void output_centroid_velocities(FILE *out, struct test_particles *part, struct parameters param);
void free_particles(struct test_particles *part);

#endif
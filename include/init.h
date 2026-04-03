#ifndef INIT_H
#define INIT_H

#include "sim_structs.h"

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon);
void set_woods_saxon(struct woods_saxon *ws, struct parameters param, double V0, double a);
void set_skyrme(struct skyrme *skm, double A, double B, double gamma);
void create_particles(struct test_particles *part, int num);
void initialize_particles(struct test_particles *part, struct parameters param, struct woods_saxon ws, struct skyrme skm);
void generate_random_particles(struct test_particles *part, double r_max, int num);
void output_centroids(FILE *out, struct test_particles *part, struct parameters param);
void free_particles(struct test_particles *part);

#endif
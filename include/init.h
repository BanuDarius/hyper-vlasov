#ifndef INIT_H
#define INIT_H

#include "sim_structs.h"

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon, double sigma_k);
void set_woods_saxon(struct woods_saxon *ws, struct parameters param, double V0, double a);
void set_skyrme(struct skyrme *skm, double A, double B, double gamma);
void set_fermi_levels(struct fermi *fermi, double epsilon_p, double epsilon_n);
void create_particles(struct test_particles *part, int num);
void generate_random_particles(struct test_particles *part, double r_max, int num);
void compute_particle_energies(struct test_particles *part, struct woods_saxon ws, double sigma_k, int num, int z, int type);
void initialize_particles(struct test_particles *part_p, struct test_particles *part_n, struct parameters param, struct woods_saxon ws, struct skyrme skm, struct fermi *fermi_levels);
void output_centroids(FILE *out, struct test_particles part, int num);
void free_particles(struct test_particles *part);

#endif
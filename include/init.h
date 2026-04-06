#ifndef INIT_H
#define INIT_H

#include <stdint.h>

#include "sim_structs.h"

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon, double sigma_k, double sigma_r);
void set_woods_saxon(struct woods_saxon *ws, double V0, double R12, double a);
void set_skyrme(struct skyrme *skm, double A, double B, double C, double gamma);
void set_fermi_levels(struct fermi *fermi, double epsilon_p, double epsilon_n);
void set_world(struct world *world, double d_max, int n);
void create_particle_count(struct particle_count *part_count, struct world world);
void create_volumetric_density(struct volumetric_density *volume, struct world world);
void create_particles(struct test_particles *part, int protons, int neutrons);
void output_centroids(FILE *out, struct test_particles part, int type);
void output_particle_count(FILE *out, struct particle_count particle_count, struct world world);
void output_volumetric_density(FILE *out, struct volumetric_density volume, struct world world);
void output_vtk_header_count(FILE *out, struct world world);
void output_vtk_header_volumetric(FILE *out, struct world world);
void free_particles(struct test_particles *part);
void free_particle_count(struct particle_count *part_count);
void free_volumetric_density(struct volumetric_density *volume);
static inline uint64_t swap_endian(double v);

#endif
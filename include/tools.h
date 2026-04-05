#ifndef TOOLS_H
#define TOOLS_H

#include "sim_structs.h"

double rand_val(double min, double max);
void random_vec(double *v, double max);
double compute_energy(struct test_particles *part, struct woods_saxon *ws, double sigma_k, int z, int i);
void compute_particle_energies(struct test_particles *part, struct woods_saxon *ws, struct parameters param);
void compute_particle_densities(struct test_particles *part, struct parameters param);
void compute_volumetric_density(struct volumetric_density *volume_dens, struct particle_count part_count, struct world world_visual, struct world world_data, struct parameters param);
void apply_constant_to_density(struct volumetric_density *volume_dens, struct world world, struct parameters param);
void generate_random_particles(struct test_particles *part, double r_max);
void scatter_particles(struct particle_count *part_count, struct test_particles *part, struct world world, int type);
void generate_checking_particles(struct test_particles *part, struct woods_saxon *ws, struct parameters param, struct fermi *fermi_levels);
void chi_squared(struct test_particles *part, struct woods_saxon *ws, struct skyrme skm, int part_per_nucleon);
void relax_woods_saxon(struct woods_saxon *ws, struct woods_saxon *ws_old, double coef);
double kinetic_energy();
double fluctuation_energy(double sigma_k);
double calc_sigma(double fwhm);
void world_pos_to_vector(double *v, struct world world, int idx);
void copy_particle_pos_to_vector(double *v, struct test_particles part, int i);
void copy_particle_vel_to_vector(double *v, struct test_particles part, int i);
void copy_vector_to_particle_pos(struct test_particles *part, double *v, int i);
void copy_vector_to_particle_vel(struct test_particles *part, double *v, int i);

#endif
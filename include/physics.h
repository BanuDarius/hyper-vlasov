#ifndef PHYSICS_H
#define PHYSICS_H

void initialize_particles(struct test_particles *part, struct parameters param, struct woods_saxon *ws, struct skyrme skm, struct fermi *fermi_levels);
double nuclear_radius(unsigned short a);
int max_particles(double r_max, double k_max, int total_test_part);
double woods_saxon_potential(struct woods_saxon ws, double r);
double skyrme_potential(struct skyrme skm, double rho_p, double rho_n, int type);
double coulomb_potential(struct woods_saxon ws, double z, double r);

#endif
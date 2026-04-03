#ifndef SIM_STRUCTS_H
#define SIM_STRUCTS_H

#define PROTONS 0
#define NEUTRONS 1
#define MAX_ITERATIONS 30
#define DELTA_EPSILON_TOLERANCE 0.01

#define MC2 935.0
#define H_BAR_C 197.33
#define RHO_0 0.16

struct parameters {
	int test_part_per_nucleon, max_test_part, z, n;
	double sigma_k, sigma_r, r_max;
};

struct test_particles {
	double *x, *y, *z, *kx, *ky, *kz;
	double *energy, *density;
};

struct woods_saxon {
	double V0, R12, a;
};

struct skyrme {
	double A, B, gamma;
};

struct fermi {
	double epsilon_p, epsilon_n;
};

#endif
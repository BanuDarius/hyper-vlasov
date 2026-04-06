#ifndef SIM_STRUCTS_H
#define SIM_STRUCTS_H

#define PROTONS 0
#define NEUTRONS 1
#define PROTONS_AND_NEUTRONS 2
#define MAX_ITERATIONS 32
#define DELTA_EPSILON_TOLERANCE 0.1

#define INPUT_FILE_COUNT 14

#define MC2 935.0
#define H_BAR_C 197.33
#define RHO_0 0.16
#define K_MAX 1.5

struct parameters {
	double sigma_k, sigma_r, r_max;
	int test_part_per_nucleon, max_test_part, z, n;
};

struct test_particles {
	double *energy, *density_p, *density_n;
	double *x, *y, *z, *kx, *ky, *kz;
	int protons, neutrons;
};

struct particle_count {
	int *count;
};

struct volumetric_density {
	double *density;
};

struct woods_saxon {
	double V0, R12, a;
};

struct skyrme {
	double A, B, C, gamma;
};

struct fermi {
	double epsilon_p, epsilon_n;
};

struct world {
	int n[3];
	double d_max[3];
};

#endif
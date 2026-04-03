#ifndef SIM_STRUCTS_H
#define SIM_STRUCTS_H

#define PROTONS 0
#define NEUTRONS 1

#define H_BAR_C 197.33
#define MC2 935.0

struct parameters {
	int total_test_part, test_part_per_nucleon, max_test_part, z, n;
	double r_max, sigma_k;
};

struct test_particles {
	double *x, *y, *z, *kx, *ky, *kz;
	double *energy;
};

struct energies {
	double *e_gauss;
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
#ifndef SIM_STRUCTS_H
#define SIM_STRUCTS_H

struct parameters {
	int total_test_part, test_part_per_nucleon, max_test_part, z, n;
	double r_max;
};

struct test_particles {
	double *x, *y, *z, *kx, *ky, *kz;
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
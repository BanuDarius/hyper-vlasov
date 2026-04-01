#ifndef SIM_STRUCTS_H
#define SIM_STRUCTS_H

struct parameters {
	int total_test_part, test_part_per_nucleon, z, n;
	double r;
};

struct test_particles {
	double *x, *y, *z, *px, *py, *pz;
};

#endif
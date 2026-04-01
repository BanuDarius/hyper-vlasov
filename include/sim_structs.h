#ifndef SIM_STRUCTS_H
#define SIM_STRUCTS_H

struct parameters {
	double num_z, num_p, num_test_part;
	double r;
};

struct test_particles {
	double *x, *y, *z, *px, *py, *pz;
};

#endif
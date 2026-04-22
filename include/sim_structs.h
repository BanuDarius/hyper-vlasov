/* MIT License

Copyright (c) 2026 Banu Darius-Matei

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */

#ifndef SIM_STRUCTS_H
#define SIM_STRUCTS_H

#define PROTONS 0
#define NEUTRONS 1
#define PROTONS_AND_NEUTRONS 2
#define MAX_ITERATIONS 32
#define DELTA_EPSILON_TOLERANCE 0.1

#define INPUT_FILE_COUNT 16

#define MC2 935.0
#define H_BAR_C 197.33
#define RHO_0 0.16
#define K_MAX 1.5

struct parameters {
	double sigma_k, sigma_r, r_max, t_f;
	int part_per_nucleon, max_test_part, steps, z, n;
};

struct test_particles {
	int protons, neutrons;
	double *x, *y, *z, *kx, *ky, *kz;
	double *energy, *density_p, *density_n;
};

struct particle_count {
	int *count;
};

struct scalar_field {
	double *v;
};

struct vector_field {
	double *x, *y, *z;
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
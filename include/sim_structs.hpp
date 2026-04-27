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

#define MAX_INIT_ITERATIONS 32
#define MAX_SOR_ITERATIONS 128
#define DELTA_EPSILON_TOLERANCE 0.1

#define INPUT_FILE_COUNT 17

#define IDX(i, j, k, nx, ny, nz) (((i) * (ny) * (nz)) + ((j) * (nz)) + (k))

template <typename T>
constexpr T mc2 = T(935.0);

template <typename T>
constexpr T rho_0 = T(0.16);

template <typename T>
constexpr T h_bar_c = T(197.33);

template <typename T>
constexpr T k_max = T(1.5);

template <typename T>
constexpr T sor_tolerance = T(1e-4);

template <typename T>
struct Parameters {
	T sigma_k, sigma_r, r_max, t_f;
	int part_per_nucleon, max_test_part, substeps, steps, z, n;
};

template <typename T>
struct TestParticles {
	int protons, neutrons;
	T *x, *y, *z, *kx, *ky, *kz;
	T *energy, *density_p, *density_n, *fx, *fy, *fz;
};

template <typename T>
struct ScalarField {
	T *v;
};

template <typename T>
struct VectorField {
	T *x, *y, *z;
};

template <typename T>
struct World {
	int n[3];
	T d_max[3];
};

template <typename T>
struct WoodsSaxon {
	T V0, R12, a;
};

template <typename T>
struct Skyrme {
	T A, B, C, gamma;
};

template <typename T>
struct Fermi {
	T epsilon_p, epsilon_n;
};

/*struct ParticleCount {
	int *count;
};*/

#endif
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

#include <array>
#include <vector>
#include <numbers>

#define PROTONS 0
#define NEUTRONS 1
#define PROTONS_AND_NEUTRONS 2

#define RESET_STEPS 10
#define MAX_INIT_ITERATIONS 32
#define MAX_SOR_ITERATIONS 512

#define STRING_SIZE 128
#define INPUT_FILE_COUNT 21

#define IDX(i, j, k, nx, ny, nz) (((i) * (ny) * (nz)) + ((j) * (nz)) + (k))

template <typename T> constexpr T mc2 = T(935.0);
template <typename T> constexpr T k_max = T(1.5);
template <typename T> constexpr T rho_0 = T(0.16);
template <typename T> constexpr T h_bar_c = T(197.33);
template <typename T> constexpr T sor_tolerance = T(1e-4);
template <typename T> constexpr T pi = std::numbers::pi_v<T>;
template <typename T> constexpr T delta_epsilon_tolerance = T(0.1);

template <typename T>
struct Parameters {
	bool use_gpu;
	T sigma_k, sigma_r, r_max, t_f, t_exc, eta_exc, d_max_scale;
	int part_per_nucleon, max_test_part, substeps, steps, z, n;
};

template <typename T>
struct TestParticles {
	int protons, neutrons;
	std::vector<T> energy, density_p, density_n;
	std::vector<T> x, y, z, kx, ky, kz, fx, fy, fz;
	TestParticles(int p, int n) {
		protons = p;
		neutrons = n;
		int total = protons + neutrons;
		x.resize(total); y.resize(total); z.resize(total);
		kx.resize(total); ky.resize(total); kz.resize(total);
		fx.resize(total); fy.resize(total); fz.resize(total);
		energy.resize(total); density_p.resize(total); density_n.resize(total);
	}
};

template <typename T>
struct ScalarField {
	std::vector<T> v;
	ScalarField(int field_size) {
		v.resize(field_size);
	}
};

template <typename T>
struct VectorField {
	std::vector<T> x, y, z;
	VectorField(int field_size) {
		x.resize(field_size); y.resize(field_size); z.resize(field_size);
	}
};

template <typename T>
struct World {
	std::array<int, 3> n;
	std::array<T, 3> d_max;
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
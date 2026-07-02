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
#include <memory>
#include <cassert>

constexpr int is_proton = 0;
constexpr int is_neutron = 1;
constexpr int is_proton_or_neutron = 2;

constexpr int reset_steps = 10;
constexpr int max_init_iterations = 32;
constexpr int max_sor_iterations = 512;

constexpr int string_size = 128;
constexpr int input_file_count = 23;

constexpr int grid_idx(int i, int j, int k, int nx, int ny, int nz) noexcept {
	(void)nx;
	return (i * ny * nz) + (j * nz) + k;
}

template <typename T> constexpr T mc2 = T(935.0);
template <typename T> constexpr T k_max = T(1.5);
template <typename T> constexpr T rho_0 = T(0.16);
template <typename T> constexpr T h_bar_c = T(197.33);
template <typename T> constexpr T sor_tolerance = T(1e-4);
template <typename T> constexpr T density_tolerance = T(0.01);
template <typename T> constexpr T pi = T(3.14159265358979323846);
template <typename T> constexpr T delta_epsilon_tolerance = T(0.1);

template <typename T>
struct Parameters {
	bool use_gpu;
	T sigma_k, sigma_r, r_max, t_f, t_exc, eta_exc, d_max_scale, sample_position;
	int part_per_nucleon, max_test_part, density_samples, substeps, steps, z, n;
};

template <typename T>
struct TestParticles {
	int protons, neutrons;
	std::unique_ptr<T[]> x, y, z, kx, ky, kz, fx, fy, fz, energy;
	TestParticles(int protons_new, int neutrons_new) : protons(protons_new), neutrons(neutrons_new) {
		int total = protons + neutrons;
		x = std::unique_ptr<T[]>(new T[total]); y = std::unique_ptr<T[]>(new T[total]); z = std::unique_ptr<T[]>(new T[total]);
		kx = std::unique_ptr<T[]>(new T[total]); ky = std::unique_ptr<T[]>(new T[total]); kz = std::unique_ptr<T[]>(new T[total]);
		fx = std::unique_ptr<T[]>(new T[total]); fy = std::unique_ptr<T[]>(new T[total]); fz = std::unique_ptr<T[]>(new T[total]);
		energy = std::unique_ptr<T[]>(new T[total]);
		#pragma omp parallel for simd schedule(static)
		for(int i = 0; i < total; i++) {
			x[i] = T(0.0); y[i] = T(0.0); z[i] = T(0.0);
			kx[i] = T(0.0); ky[i] = T(0.0); kz[i] = T(0.0);
			fx[i] = T(0.0); fy[i] = T(0.0); fz[i] = T(0.0);
			energy[i] = T(0.0);
		}
	}
};

template <typename T>
struct ScalarField {
	std::size_t size;
	std::unique_ptr<T[]> v;
	ScalarField(int field_size) : size(field_size) {
		v = std::unique_ptr<T[]>(new T[size]);
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++)
			v[i] = T(0.0);
	}
	ScalarField(const ScalarField &other) : size(other.size) {
		v = std::unique_ptr<T[]>(new T[size]);
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++)
			v[i] = other.v[i];
	}
	ScalarField &operator=(const ScalarField &other) {
		if(this == &other) return *this;
		assert(size == other.size && "FIELD SIZES DO NOT MATCH!");
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++)
			v[i] = other.v[i];
		return *this;
	}
	ScalarField &operator+=(const ScalarField &other) {
		assert(size == other.size && "FIELD SIZES DO NOT MATCH!");
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++)
			v[i] += other.v[i];
		return *this;
	}
	ScalarField &operator-=(const ScalarField &other) {
		assert(size == other.size && "FIELD SIZES DO NOT MATCH!");
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++)
			v[i] -= other.v[i];
		return *this;
	}
	ScalarField(ScalarField &&other) noexcept = default;
	ScalarField &operator=(ScalarField &&other) noexcept = default;
	~ScalarField() = default;
};

template <typename T>
struct VectorField {
	std::size_t size;
	std::unique_ptr<T[]> x, y, z;
	VectorField(int field_size) : size(field_size) {
		x = std::unique_ptr<T[]>(new T[size]);
		y = std::unique_ptr<T[]>(new T[size]);
		z = std::unique_ptr<T[]>(new T[size]);
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++) {
			x[i] = T(0.0); y[i] = T(0.0); z[i] = T(0.0); 
		}
	}
	VectorField &operator=(const VectorField &other) {
		if(this == &other) return *this;
		assert(size == other.size && "FIELD SIZES DO NOT MATCH!");
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++) {
			x[i] = other.x[i]; y[i] = other.y[i]; z[i] = other.z[i]; 
		}
		return *this;
	}
	VectorField &operator+=(const VectorField &other) {
		assert(size == other.size && "FIELD SIZES DO NOT MATCH!");
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++) {
			x[i] += other.x[i]; y[i] += other.y[i]; z[i] += other.z[i]; 
		}
		return *this;
	}
	VectorField &operator-=(const VectorField &other) {
		assert(size == other.size && "FIELD SIZES DO NOT MATCH!");
		#pragma omp parallel for simd schedule(static)
		for(std::size_t i = 0; i < size; i++) {
			x[i] -= other.x[i]; y[i] -= other.y[i]; z[i] -= other.z[i]; 
		}
		return *this;
	}
	VectorField(VectorField &&other) noexcept = default;
	VectorField &operator=(VectorField &&other) noexcept = default;
	~VectorField() = default;
};

template <typename T>
struct World {
	std::array<int, 3> n;
	std::array<T, 3> d_max;
	World() = default;
	World(int n_new, T d_max_new) {
		for(int i = 0; i < 3; i++) {
			n[i] = n_new;
			d_max[i] = d_max_new;
		}
	}
};

template <typename T>
struct WoodsSaxon {
	T V0_p, R12_p, a_p;
	T V0_n, R12_n, a_n;
	WoodsSaxon() = default;
	WoodsSaxon(T V0, T R12, T a) {
		V0_p = V0; V0_n = V0;
		R12_p = R12, R12_n = R12;
		a_p = a; a_n = a;
	}
};

template <typename T>
struct Skyrme {
	T A, B, C, gamma;
	Skyrme() = default;
	Skyrme(T A_new, T B_new, T C_new, T gamma_new) {
		A = A_new; B = B_new; C = C_new; gamma = gamma_new;
	}
};

template <typename T>
struct Fermi {
	T epsilon_p, epsilon_n;
	Fermi() = default;
	Fermi(T epsilon_p_new, T epsilon_n_new) {
		epsilon_p = epsilon_p_new; epsilon_n = epsilon_n_new;
	}
};

#endif
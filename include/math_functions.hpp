// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <array>
#include <cmath>
#include <random>
#include <concepts>

template <std::floating_point T>
inline T rand_val(T min, T max) noexcept {
	static std::mt19937 mt(128);
	T s = static_cast<T>(mt()) / static_cast<T>(mt.max());
	return min + s * (max - min);
}

template <std::floating_point T>
inline void random_vec(std::array<T, 3> &v, T max) noexcept {
	for(int i = 0; i < 3; i++)
		v[i] = rand_val(-max, max);
}

constexpr int grid_idx(int i, int j, int k, int nx, int ny, int nz) noexcept {
	(void)nx;
	return (i * ny * nz) + (j * nz) + k;
}

template <std::floating_point T>
inline T dot(const std::array<T, 3> a, const std::array<T, 3> b) noexcept {
	T x = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return x;
}

template <std::floating_point T>
inline T magnitude(const std::array<T, 3> a) noexcept {
	T x = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	return x;
}

template <std::floating_point T>
inline std::array<T, 3> operator+(std::array<T, 3> lhs, const std::array<T, 3> rhs) noexcept {
	std::array<T, 3> x = { lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2] };
	return x;
}

template <std::floating_point T>
inline std::array<T, 3> operator-(std::array<T, 3> lhs, const std::array<T, 3> rhs) noexcept {
	std::array<T, 3> x = { lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2] };
	return x;
}

template <std::floating_point T>
inline std::array<T, 3> &operator+=(std::array<T, 3> &lhs, const std::array<T, 3> &rhs) noexcept {
	lhs[0] += rhs[0]; lhs[1] += rhs[1]; lhs[2] += rhs[2];
	return lhs;
}

template <std::floating_point T>
inline std::array<T, 3> operator*(std::array<T, 3> lhs, T b) noexcept {
	std::array<T, 3> x = { lhs[0] *= b, lhs[1] *= b, lhs[2] *= b };
	return x;
}

#endif
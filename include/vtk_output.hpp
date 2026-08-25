// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include <fstream>
#include <cstring>
#include <concepts>

#include "sim_structs.hpp"

inline uint32_t swap_endian(float v) noexcept {
	uint32_t data;
	std::memcpy(&data, &v, 4);
	return __builtin_bswap32(data);
}

template <std::floating_point T> void output_vtk_header_start(std::ofstream &output_file, World<T> world);
void output_vtk_header_scalar_next(std::ofstream &output_file, const char *name, int type);
void output_vtk_header_vector_next(std::ofstream &output_file, const char *name, int type);
template <std::floating_point T> void output_density_samples_positions(std::ofstream &output_file, const Parameters<T> &param, const World<T> &world);
template <std::floating_point T> void output_density_samples(std::ofstream &output_file, const float *samples, const Parameters<T> &param);
template <std::floating_point T> void output_scalar_field(std::ofstream &output_file, const ScalarField<T> &field, const World<T> &world, const char *name);
template <std::floating_point T> void output_vector_field(std::ofstream &output_file, const VectorField<T> &field, const World<T> &world, const char *name);

#endif
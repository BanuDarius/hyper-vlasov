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
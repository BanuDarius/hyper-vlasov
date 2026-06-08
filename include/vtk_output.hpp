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

#include <cstdio>
#include <cstring>

#include "sim_structs.hpp"

inline uint32_t swap_endian(float v) noexcept {
	uint32_t data;
	std::memcpy(&data, &v, 4);
	return __builtin_bswap32(data);
}

template <typename T> void output_vtk_header_start(std::FILE *out, World<T> world);
void output_vtk_header_scalar_next(std::FILE *out, const char *name, int type);
void output_vtk_header_vector_next(std::FILE *out, const char *name, int type);
template <typename T> void output_density_samples_positions(std::FILE *out, const Parameters<T> &param, const World<T> &world);
template <typename T> void output_density_samples(std::FILE *out, const float *samples, const Parameters<T> &param);
template <typename T> void output_scalar_field(std::FILE *out, const ScalarField<T> &field, const World<T> &world, const char *name);
template <typename T> void output_vector_field(std::FILE *out, const VectorField<T> &field, const World<T> &world, const char *name);

#endif
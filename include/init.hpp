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

#ifndef INIT_H
#define INIT_H

#include <cstdio>
#include <cstdint>
#include <cstring>

#include "sim_structs.hpp"

static inline uint32_t swap_endian(float v) {
	uint32_t data;
	std::memcpy(&data, &v, 4);
	return __builtin_bswap32(data);
}

template <typename T> void set_parameters(Parameters<T> &param, int z, int n, int part_per_nucleon, int steps, int substeps, int density_samples, bool use_gpu, T sample_position, T sigma_k, T sigma_r, T t_f, T t_exc, T eta_exc);
template <typename T> void output_vtk_header_start(FILE *out, World<T> world);
void output_vtk_header_scalar_next(FILE *out, const char *name, int type);
void output_vtk_header_vector_next(FILE *out, const char *name, int type);
template <typename T> void output_density_samples_positions(FILE *out, const Parameters<T> &param, const World<T> &world);
template <typename T> void output_density_samples(FILE *out, const float *samples, const Parameters<T> &param);
template <typename T> void output_scalar_field(FILE *out, const ScalarField<T> &field, const World<T> &world, const char *name);
template <typename T> void output_vector_field(FILE *out, const VectorField<T> &field, const World<T> &world, const char *name);
template <typename T> void read_input_file(FILE *in, Skyrme<T> &skm, World<T> &world, Fermi<T> &fermi_levels, Parameters<T> &param, WoodsSaxon<T> &ws);

#endif
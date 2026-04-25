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

#include <stdint.h>
#include <string.h>

#include "sim_structs.h"

void set_parameters(Parameters *param, int z, int n, int part_per_nucleon, double sigma_k, double sigma_r, double t_f, int steps);
void set_woods_saxon(WoodsSaxon *ws, double V0, double R12, double a);
void set_skyrme(Skyrme *skm, double A, double B, double C, double gamma);
void set_fermi_levels(Fermi *fermi, double epsilon_p, double epsilon_n);
void set_world(World *world, double d_max, int n);
void create_scalar_field_single(ScalarField *field, World world);
void create_scalar_field_double(ScalarField *field, World world);
void create_vector_field(VectorField *field, World world);
void create_particles(TestParticles *part, int protons, int neutrons);
void output_centroids(FILE *out, TestParticles part, int type);
void output_scalar_field(FILE *out, ScalarField volume, World world, char *name);
void output_vector_field(FILE *out, VectorField field, World world, char *name);
void output_vtk_header_start(FILE *out, World world);
void output_vtk_header_scalar_next(FILE *out, char *name, int type);
void output_vtk_header_vector_next(FILE *out, char *name, int type);
void free_particles(TestParticles *part);
void free_particle_count(ParticleCount *part_count);
void free_vector_field(VectorField *field);
void free_scalar_field(ScalarField *field);
void read_input_file(FILE *in, Skyrme *skm, World *world, Fermi *fermi_levels, Parameters *param, WoodsSaxon *ws);
void set_output_filename(char *output_filename, char *output_directory, int i);
//void create_particle_count(ParticleCount *part_count, World world);
//void output_particle_count(FILE *out, ParticleCount particle_count, World world);

static inline uint64_t swap_endian(double v) {
	uint64_t data;
	memcpy(&data, &v, sizeof(double));
	return __builtin_bswap64(data);
}

#endif
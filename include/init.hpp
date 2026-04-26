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

#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include "math_functions.hpp"
#include "sim_structs.hpp"
#include "physics_formulas.hpp"

template <typename T>
static inline auto swap_endian(T v) {
	if constexpr(sizeof(T) == 8) {
		uint64_t data;
		std::memcpy(&data, &v, 8);
		return __builtin_bswap64(data);
	} else {
		uint32_t data;
		std::memcpy(&data, &v, 4);
		return __builtin_bswap32(data);
	}
}

template <typename T>
void set_parameters(Parameters<T> *param, int z, int n, int part_per_nucleon, T sigma_k, T sigma_r, T t_f, int steps) {
	param->z = z;
	param->n = n;
	param->t_f = t_f;
	param->steps = steps;
	param->sigma_k = sigma_k;
	param->sigma_r = sigma_r;
	
	param->part_per_nucleon = part_per_nucleon;
	param->r_max = nuclear_radius<T>(z + n);
	
	param->max_test_part = max_particles(param->r_max, K_MAX, param->part_per_nucleon);
}

template <typename T>
void set_woods_saxon(WoodsSaxon<T> *ws, T V0, T R12, T a) {
	ws->a = a;
	ws->V0 = V0;
	ws->R12 = R12;
}

template <typename T>
void set_skyrme(Skyrme<T> *skm, T A, T B, T C, T gamma) {
	skm->A = A;
	skm->B = B;
	skm->C = C;
	skm->gamma = gamma;
}

template <typename T>
void set_fermi_levels(Fermi<T> *fermi, T epsilon_p, T epsilon_n) {
	fermi->epsilon_p = epsilon_p;
	fermi->epsilon_n = epsilon_n;
}

template <typename T>
void set_world(World<T> *world, T d_max, int n) {
	for(int i = 0; i < 3; i++) {
		world->n[i] = n;
		world->d_max[i] = d_max;
	}
}

template <typename T>
void create_scalar_field_single(ScalarField<T> *field, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	field->v = (T*)malloc(world_size * sizeof(T));
	
	if(field->v == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for simd
	for(int i = 0; i < world_size; i++)
		field->v[i] = T(0.0);
}

template <typename T>
void create_scalar_field_double(ScalarField<T> *field, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	field->v = (T*)malloc(2 * world_size * sizeof(T));
	
	if(field->v == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for simd
	for(int i = 0; i < 2 * world_size; i++)
		field->v[i] = T(0.0);
}

template <typename T>
void create_vector_field(VectorField<T> *field, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	field->x = (T*)malloc(2 * world_size * sizeof(T));
	field->y = (T*)malloc(2 * world_size * sizeof(T));
	field->z = (T*)malloc(2 * world_size * sizeof(T));
	
	if(field->x == NULL || field->y == NULL || field->z == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for simd
	for(int i = 0; i < 2 * world_size; i++) {
		field->x[i] = T(0.0);
		field->y[i] = T(0.0);
		field->z[i] = T(0.0);
	}
}

template <typename T>
void create_particles(TestParticles<T> *part, int protons, int neutrons) {
	int total = protons + neutrons;
	part->protons = protons;
	part->neutrons = neutrons;
	part->x = (T*)malloc(total * sizeof(T));
	part->y = (T*)malloc(total * sizeof(T));
	part->z = (T*)malloc(total * sizeof(T));
	part->kx = (T*)malloc(total * sizeof(T));
	part->ky = (T*)malloc(total * sizeof(T));
	part->kz = (T*)malloc(total * sizeof(T));
	part->fx = (T*)malloc(total * sizeof(T));
	part->fy = (T*)malloc(total * sizeof(T));
	part->fz = (T*)malloc(total * sizeof(T));
	part->energy = (T*)malloc(total * sizeof(T));
	part->density_p = (T*)malloc(total * sizeof(T));
	part->density_n = (T*)malloc(total * sizeof(T));
	
	if(part->x == NULL || part->y == NULL || part->z == NULL
	|| part->kx == NULL || part->ky == NULL || part->kz == NULL
	|| part->fx == NULL || part->fy == NULL || part->fz == NULL
	|| part->energy == NULL || part->density_p == NULL || part->density_n == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
}

template <typename T>
void output_vtk_header_start(FILE *out, World<T> world) {
	std::fprintf(out, "# vtk DataFile Version 3.0\n");
	std::fprintf(out, "Volumetric data\n");
	std::fprintf(out, "BINARY\n");
	std::fprintf(out, "DATASET STRUCTURED_POINTS\n");
	std::fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	std::fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	std::fprintf(out, "SPACING %lf %lf %lf\n", 2.0 * world.d_max[0] / world.n[0], 2.0 * world.d_max[1] / world.n[1], 2.0 * world.d_max[2] / world.n[2]);
	std::fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
}

void output_vtk_header_scalar_next(FILE *out, const char *name, int type) {
	char tag;
	if(type == PROTONS) tag = 'p';
	else if(type == NEUTRONS) tag = 'n';
	else tag = 't';
	std::fprintf(out, "SCALARS %s_%c double 1\n", name, tag);
	std::fprintf(out, "LOOKUP_TABLE default\n");
}

void output_vtk_header_vector_next(FILE *out, const char *name, int type) {
	char tag;
	if(type == PROTONS) tag = 'p';
	else if(type == NEUTRONS) tag = 'n';
	std::fprintf(out, "VECTORS %s_%c double\n", name, tag);
}

template <typename T>
void output_scalar_field(FILE *out, const ScalarField<T> &field, const World<T> &world, const char *name) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	uint64_t *vtk_density_p = (uint64_t*)malloc(world_size * sizeof(uint64_t));
	uint64_t *vtk_density_n = (uint64_t*)malloc(world_size * sizeof(uint64_t));
	uint64_t *vtk_density_t = (uint64_t*)malloc(world_size * sizeof(uint64_t));
	if(vtk_density_p == NULL || vtk_density_n == NULL || vtk_density_t == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for collapse(3)
	for(int k = 0; k < nz; k++) {
		for(int j = 0; j < ny; j++) {
			for(int i = 0; i < nx; i++) {
				int idx = IDX(i, j, k, nx, ny, nz);
				int write_idx = (k * ny * nx) + (j * nx) + i;
				
				vtk_density_p[write_idx] = swap_endian(field.v[idx]);
				vtk_density_n[write_idx] = swap_endian(field.v[idx + world_size]);
				vtk_density_t[write_idx] = swap_endian(field.v[idx] + field.v[idx + world_size]);
			}
		}
	}
	output_vtk_header_scalar_next(out, name, PROTONS);
	fwrite(vtk_density_p, sizeof(uint64_t), world_size, out);
	
	output_vtk_header_scalar_next(out, name, NEUTRONS);
	fwrite(vtk_density_n, sizeof(uint64_t), world_size, out);
	
	output_vtk_header_scalar_next(out, name, PROTONS_AND_NEUTRONS);
	fwrite(vtk_density_t, sizeof(uint64_t), world_size, out);
	
	free(vtk_density_p); free(vtk_density_n); free(vtk_density_t);
}

template <typename T>
void output_vector_field(FILE *out, const VectorField<T> &field, const World<T> &world, const char *name) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	uint64_t *vtk_force_p = (uint64_t*)malloc(3 * world_size * sizeof(uint64_t));
	uint64_t *vtk_force_n = (uint64_t*)malloc(3 * world_size * sizeof(uint64_t));
	
	if(vtk_force_p == NULL || vtk_force_n == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for collapse(3)
	for(int k = 0; k < nz; k++) {
		for(int j = 0; j < ny; j++) {
			for(int i = 0; i < nx; i++) {
				int idx = IDX(i, j, k, nx, ny, nz);
				int write_idx = (k * ny * nx) + (j * nx) + i;
				
				vtk_force_p[3 * write_idx] = swap_endian(field.x[idx]);
				vtk_force_p[3 * write_idx + 1] = swap_endian(field.y[idx]);
				vtk_force_p[3 * write_idx + 2] = swap_endian(field.z[idx]);
				
				vtk_force_n[3 * write_idx] = swap_endian(field.x[idx + world_size]);
				vtk_force_n[3 * write_idx + 1] = swap_endian(field.y[idx + world_size]);
				vtk_force_n[3 * write_idx + 2] = swap_endian(field.z[idx + world_size]);
			}
		}
	}
	output_vtk_header_vector_next(out, name, PROTONS);
	fwrite(vtk_force_p, sizeof(uint64_t), 3 * world_size, out);
	
	output_vtk_header_vector_next(out, name, NEUTRONS);
	fwrite(vtk_force_n, sizeof(uint64_t), 3 * world_size, out);
	
	free(vtk_force_p); free(vtk_force_n);
}

template <typename T>
void free_particles(TestParticles<T> *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->kx); free(part->ky); free(part->kz);
	free(part->fx); free(part->fy); free(part->fz);
	free(part->density_p); free(part->density_n);
	free(part->energy);
}

template <typename T>
void free_vector_field(VectorField<T> *field) {
	free(field->x); free(field->y); free(field->z);
}

template <typename T>
void free_scalar_field(ScalarField<T> *field) {
	free(field->v);
}

template <typename T>
void read_input_file(FILE *in, Skyrme<T> *skm, World<T> *world, Fermi<T> *fermi_levels, Parameters<T> *param, WoodsSaxon<T> *ws) {
	double V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, t_f;
	int i = 0, num_test_part, steps, nx, z, n;
	char current[32];
	
	while(std::fscanf(in, "%s", current) != EOF) {
		if(!std::strcmp(current, "V0"))
			i += std::fscanf(in, "%lf", &V0);
		else if(!std::strcmp(current, "a"))
			i += std::fscanf(in, "%lf", &a);
		else if(!std::strcmp(current, "A"))
			i += std::fscanf(in, "%lf", &A);
		else if(!std::strcmp(current, "B"))
			i += std::fscanf(in, "%lf", &B);
		else if(!std::strcmp(current, "C"))
			i += std::fscanf(in, "%lf", &C);
		else if(!std::strcmp(current, "gamma"))
			i += std::fscanf(in, "%lf", &gamma);
		else if(!std::strcmp(current, "epsilon_p"))
			i += std::fscanf(in, "%lf", &epsilon_p);
		else if(!std::strcmp(current, "epsilon_n"))
			i += std::fscanf(in, "%lf", &epsilon_n);
		else if(!std::strcmp(current, "k_fwhm"))
			i += std::fscanf(in, "%lf", &k_fwhm);
		else if(!std::strcmp(current, "r_fwhm"))
			i += std::fscanf(in, "%lf", &r_fwhm);
		else if(!std::strcmp(current, "t_f"))
			i += std::fscanf(in, "%lf", &t_f);
		else if(!std::strcmp(current, "nx"))
			i += std::fscanf(in, "%i", &nx);
		else if(!std::strcmp(current, "num_test_part"))
			i += std::fscanf(in, "%i", &num_test_part);
		else if(!std::strcmp(current, "steps"))
			i += std::fscanf(in, "%i", &steps);
		else if(!std::strcmp(current, "n"))
			i += std:: fscanf(in, "%i", &n);
		else if(!std::strcmp(current, "z"))
			i += std::fscanf(in, "%i", &z);
	}
	if(i != INPUT_FILE_COUNT) {
		std::fprintf(stderr, "Error: Invalid input file.\n");
		exit(1);
	}
	
	T sigma_k = calc_sigma(T(k_fwhm)), sigma_r = calc_sigma(T(r_fwhm));
	T d_max = T(1.35) * nuclear_radius<T>(z + n);
	
	set_skyrme(skm, T(A), T(B), T(C), gamma);
	set_world(world, T(d_max), nx);
	set_fermi_levels(fermi_levels, T(epsilon_p), T(epsilon_n));
	set_parameters(param, z, n, num_test_part, T(sigma_k), T(sigma_r), T(t_f), steps);
	set_woods_saxon(&ws[0], T(V0), T(0.8) * T(param->r_max), a);
	set_woods_saxon(&ws[1], T(V0), T(0.8) * T(param->r_max), a);
}

void set_output_filename(char *output_filename, char *output_directory, int i) {
	std::sprintf(output_filename, "%sout-%04d.vtk", output_directory, i);
}

template <typename T>
void output_vtk_header_count(FILE *out, const World<T> &world) {
	std::fprintf(out, "# vtk DataFile Version 3.0\n");
	std::fprintf(out, "Volumetric count\n");
	std::fprintf(out, "BINARY\n");
	std::fprintf(out, "DATASET STRUCTURED_POINTS\n");
	std::fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	std::fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	std::fprintf(out, "SPACING %lf %lf %lf\n", 2.0 * world.d_max[0] / world.n[0], 2.0 * world.d_max[1] / world.n[1], 2.0 * world.d_max[2] / world.n[2]);
	std::fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
	std::fprintf(out, "SCALARS count int 1\n");
	std::fprintf(out, "LOOKUP_TABLE default\n");
}

/*void output_particle_count(FILE *out, ParticleCount<T> particle_count, World<T> world) {
	output_vtk_header_count(out, world);
	int total = world.n[0] * world.n[1] * world.n[2];
	uint32_t *vtk_count = malloc(total * sizeof(uint32_t));
	if(vtk_count == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for
	for(int i = 0; i < total; i++)
		vtk_count[i] = __builtin_bswap32(particle_count.count[i]);
	
	fwrite(vtk_count, sizeof(uint32_t), total, out);
	free(vtk_count);
}

void create_particle_count(ParticleCount<T> *part_count, World<T> world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	part_count->count = malloc(2 * world_size * sizeof(int));
	if(part_count->count == NULL) {
		std::fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for simd
	for(int i = 0; i < 2 * world_size; i++)
		part_count->count[i] = 0;
}

void output_centroids(FILE *out, TestParticles<T> part, int type) {
	int start, end;
	if(type == PROTONS) { start = 0; end = part.protons; }
	else if(type == NEUTRONS) { start = 0; end = part.protons + part.neutrons; }
	else { start = 0; end = part.protons + part.neutrons; }
	
	for(int i = start; i < end; i++) {
		fwrite(&part.x[i], sizeof(double), 1, out);
		fwrite(&part.y[i], sizeof(double), 1, out);
		fwrite(&part.z[i], sizeof(double), 1, out);
		fwrite(&part.kx[i], sizeof(double), 1, out);
		fwrite(&part.ky[i], sizeof(double), 1, out);
		fwrite(&part.kz[i], sizeof(double), 1, out);
		fwrite(&part.energy[i], sizeof(double), 1, out);
		fwrite(&part.density_p[i], sizeof(double), 1, out);
		fwrite(&part.density_n[i], sizeof(double), 1, out);
	}
}
template <typename T>
void free_particle_count(ParticleCount<T> *part_count) {
	free(part_count->count);
}*/

//void create_particle_count(ParticleCount<T> *part_count, World<T> world);
//void output_particle_count(FILE *out, ParticleCount<T> particle_count, World<T> world);
//void output_centroids(FILE *out, TestParticles<T> part, int type);

#endif
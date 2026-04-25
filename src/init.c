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

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"

void set_parameters(Parameters *param, int z, int n, int part_per_nucleon, double sigma_k, double sigma_r, double t_f, int steps) {
	param->z = z;
	param->n = n;
	param->t_f = t_f;
	param->steps = steps;
	param->sigma_k = sigma_k;
	param->sigma_r = sigma_r;
	
	param->part_per_nucleon = part_per_nucleon;
	param->r_max = nuclear_radius(z + n);
	
	param->max_test_part = max_particles(param->r_max, K_MAX, param->part_per_nucleon);
}

void set_woods_saxon(WoodsSaxon *ws, double V0, double R12, double a) {
	ws->a = a;
	ws->V0 = V0;
	ws->R12 = R12;
}

void set_skyrme(Skyrme *skm, double A, double B, double C, double gamma) {
	skm->A = A;
	skm->B = B;
	skm->C = C;
	skm->gamma = gamma;
}

void set_fermi_levels(Fermi *fermi, double epsilon_p, double epsilon_n) {
	fermi->epsilon_p = epsilon_p;
	fermi->epsilon_n = epsilon_n;
}

void set_world(World *world, double d_max, int n) {
	for(int i = 0; i < 3; i++) {
		world->n[i] = n;
		world->d_max[i] = d_max;
	}
}

void create_particle_count(ParticleCount *part_count, World world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	part_count->count = malloc(2 * world_size * sizeof(int));
	if(part_count->count == NULL) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		part_count->count[i] = 0;
}

void create_scalar_field(ScalarField *field, World world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	field->v = malloc(2 * world_size * sizeof(double));
	if(field->v == NULL) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		field->v[i] = 0.0;
}

void create_vector_field(VectorField *field, World world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	field->x = malloc(2 * world_size * sizeof(double));
	field->y = malloc(2 * world_size * sizeof(double));
	field->z = malloc(2 * world_size * sizeof(double));
	if(field->x == NULL || field->y == NULL || field->z == NULL) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		field->x[i] = 0.0;
		field->y[i] = 0.0;
		field->z[i] = 0.0;
	}
}

void create_particles(TestParticles *part, int protons, int neutrons) {
	int total = protons + neutrons;
	part->protons = protons;
	part->neutrons = neutrons;
	part->x = malloc(total * sizeof(double));
	part->y = malloc(total * sizeof(double));
	part->z = malloc(total * sizeof(double));
	part->kx = malloc(total * sizeof(double));
	part->ky = malloc(total * sizeof(double));
	part->kz = malloc(total * sizeof(double));
	part->fx = malloc(total * sizeof(double));
	part->fy = malloc(total * sizeof(double));
	part->fz = malloc(total * sizeof(double));
	part->energy = malloc(total * sizeof(double));
	part->density_p = malloc(total * sizeof(double));
	part->density_n = malloc(total * sizeof(double));
	
	if(part->x == NULL || part->y == NULL || part->z == NULL
	|| part->kx == NULL || part->ky == NULL || part->kz == NULL
	|| part->fx == NULL || part->fy == NULL || part->fz == NULL
	|| part->energy == NULL || part->density_p == NULL || part->density_n == NULL) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
}

void output_centroids(FILE *out, TestParticles part, int type) {
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

void output_scalar_field(FILE *out, ScalarField field, World world, char *name) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	uint64_t *vtk_density_p = malloc(world_size * sizeof(uint64_t));
	uint64_t *vtk_density_n = malloc(world_size * sizeof(uint64_t));
	uint64_t *vtk_density_t = malloc(world_size * sizeof(uint64_t));
	if(vtk_density_p == NULL || vtk_density_n == NULL || vtk_density_t == NULL) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
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

void output_vector_field(FILE *out, VectorField field, World world, char *name) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	uint64_t *vtk_force_p = malloc(3 * world_size * sizeof(uint64_t));
	uint64_t *vtk_force_n = malloc(3 * world_size * sizeof(uint64_t));
	
	if(vtk_force_p == NULL || vtk_force_n == NULL) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
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

void output_vtk_header_start(FILE *out, World world) {
	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "Volumetric data\n");
	fprintf(out, "BINARY\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	fprintf(out, "SPACING %lf %lf %lf\n", 2.0 * world.d_max[0] / world.n[0], 2.0 * world.d_max[1] / world.n[1], 2.0 * world.d_max[2] / world.n[2]);
	fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
}

void output_vtk_header_scalar_next(FILE *out, char *name, int type) {
	char tag;
	if(type == PROTONS) tag = 'p';
	else if(type == NEUTRONS) tag = 'n';
	else tag = 't';
	fprintf(out, "SCALARS %s_%c double 1\n", name, tag);
	fprintf(out, "LOOKUP_TABLE default\n");
}

void output_vtk_header_vector_next(FILE *out, char *name, int type) {
	char tag;
	if(type == PROTONS) tag = 'p';
	else if(type == NEUTRONS) tag = 'n';
	fprintf(out, "VECTORS %s_%c double\n", name, tag);
}

void free_particles(TestParticles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->kx); free(part->ky); free(part->kz);
	free(part->fx); free(part->fy); free(part->fz);
	free(part->density_p); free(part->density_n);
	free(part->energy);
}

void free_particle_count(ParticleCount *part_count) {
	free(part_count->count);
}

void free_vector_field(VectorField *field) {
	free(field->x); free(field->y); free(field->z);
}

void free_scalar_field(ScalarField *field) {
	free(field->v);
}

void read_input_file(FILE *in, Skyrme *skm, World *world, World *world_visual, Fermi *fermi_levels, Parameters *param, WoodsSaxon *ws) {
	double V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, t_f;
	char current[32];
	int i = 0, num_test_part, steps, nx, z, n;
	
	while(fscanf(in, "%s", current) != EOF) {
		if(!strcmp(current, "V0"))
			i += fscanf(in, "%lf", &V0);
		else if(!strcmp(current, "a"))
			i += fscanf(in, "%lf", &a);
		else if(!strcmp(current, "A"))
			i += fscanf(in, "%lf", &A);
		else if(!strcmp(current, "B"))
			i += fscanf(in, "%lf", &B);
		else if(!strcmp(current, "C"))
			i += fscanf(in, "%lf", &C);
		else if(!strcmp(current, "gamma"))
			i += fscanf(in, "%lf", &gamma);
		else if(!strcmp(current, "epsilon_p"))
			i += fscanf(in, "%lf", &epsilon_p);
		else if(!strcmp(current, "epsilon_n"))
			i += fscanf(in, "%lf", &epsilon_n);
		else if(!strcmp(current, "k_fwhm"))
			i += fscanf(in, "%lf", &k_fwhm);
		else if(!strcmp(current, "r_fwhm"))
			i += fscanf(in, "%lf", &r_fwhm);
		else if(!strcmp(current, "t_f"))
			i += fscanf(in, "%lf", &t_f);
		else if(!strcmp(current, "nx"))
			i += fscanf(in, "%i", &nx);
		else if(!strcmp(current, "num_test_part"))
			i += fscanf(in, "%i", &num_test_part);
		else if(!strcmp(current, "steps"))
			i += fscanf(in, "%i", &steps);
		else if(!strcmp(current, "n"))
			i += fscanf(in, "%i", &n);
		else if(!strcmp(current, "z"))
			i += fscanf(in, "%i", &z);
	}
	if(i != INPUT_FILE_COUNT) {
		fprintf(stderr, "Error: Invalid input file.\n");
		exit(1);
	}
	
	double sigma_k = calc_sigma(k_fwhm), sigma_r = calc_sigma(r_fwhm);
	double d_max = 1.2 * nuclear_radius(z + n);
	
	set_skyrme(skm, A, B, C, gamma);
	set_world(world, d_max, nx);
	set_world(world_visual, d_max, 4 * nx);
	set_fermi_levels(fermi_levels, epsilon_p, epsilon_n);
	set_parameters(param, z, n, num_test_part, sigma_k, sigma_r, t_f, steps);
	set_woods_saxon(&ws[0], V0, 0.8 * param->r_max, a);
	set_woods_saxon(&ws[1], V0, 0.8 * param->r_max, a);

}

void set_output_filename(char *output_filename, char *output_directory, int i) {
	sprintf(output_filename, "%sout-%04d.vtk", output_directory, i);
}

void output_vtk_header_count(FILE *out, World world) {
	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "Volumetric count\n");
	fprintf(out, "BINARY\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	fprintf(out, "SPACING %lf %lf %lf\n", 2.0 * world.d_max[0] / world.n[0], 2.0 * world.d_max[1] / world.n[1], 2.0 * world.d_max[2] / world.n[2]);
	fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
	fprintf(out, "SCALARS count int 1\n");
	fprintf(out, "LOOKUP_TABLE default\n");
}

void output_particle_count(FILE *out, ParticleCount particle_count, World world) {
	output_vtk_header_count(out, world);
	int total = world.n[0] * world.n[1] * world.n[2];
	uint32_t *vtk_count = malloc(total * sizeof(uint32_t));
	if(vtk_count == NULL) {
		fprintf(stderr, "ERROR ALLOCATING MEMORY!\n");
		exit(1);
	}
	#pragma omp parallel for
	for(int i = 0; i < total; i++)
		vtk_count[i] = __builtin_bswap32(particle_count.count[i]);
	
	fwrite(vtk_count, sizeof(uint32_t), total, out);
	free(vtk_count);
}
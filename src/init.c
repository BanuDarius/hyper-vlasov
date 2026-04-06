#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"

void set_parameters(struct parameters *param, int z, int n, int test_part_per_nucleon, double sigma_k, double sigma_r) {
	param->z = z;
	param->n = n;
	param->sigma_k = sigma_k;
	param->sigma_r = sigma_r;
	
	param->test_part_per_nucleon = test_part_per_nucleon;
	param->r_max = nuclear_radius(z + n);
	
	param->max_test_part = max_particles(param->r_max, K_MAX, param->test_part_per_nucleon);
}

void set_woods_saxon(struct woods_saxon *ws, double V0, double R12, double a) {
	ws->a = a;
	ws->V0 = V0;
	ws->R12 = R12;
}

void set_skyrme(struct skyrme *skm, double A, double B, double C, double gamma) {
	skm->A = A;
	skm->B = B;
	skm->C = C;
	skm->gamma = gamma;
}

void set_fermi_levels(struct fermi *fermi, double epsilon_p, double epsilon_n) {
	fermi->epsilon_p = epsilon_p;
	fermi->epsilon_n = epsilon_n;
}

void set_world(struct world *world, double d_max, int n) {
	for(int i = 0; i < 3; i++) {
		world->n[i] = n;
		world->d_max[i] = d_max;
	}
}

void create_particle_count(struct particle_count *part_count, struct world world) {
	int world_size = world.n[0] * world.n[1] * world.n[2] ;
	part_count->count = malloc(2 * world_size * sizeof(int));
	for(int i = 0; i < 2 * world_size; i++)
		part_count->count[i] = 0;
}

void create_volumetric_density(struct volumetric_density *volume, struct world world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	volume->density = malloc(world_size * sizeof(double));
	for(int i = 0; i < world_size; i++)
		volume->density[i] = 0.0;
}

void create_particles(struct test_particles *part, int protons, int neutrons) {
	int total = protons + neutrons;
	part->protons = protons;
	part->neutrons = neutrons;
	part->x = malloc(total * sizeof(double));
	part->y = malloc(total * sizeof(double));
	part->z = malloc(total * sizeof(double));
	part->kx = malloc(total * sizeof(double));
	part->ky = malloc(total * sizeof(double));
	part->kz = malloc(total * sizeof(double));
	part->energy = malloc(total * sizeof(double));
	part->density_p = malloc(total * sizeof(double));
	part->density_n = malloc(total * sizeof(double));
}

void output_centroids(FILE *out, struct test_particles part, int type) {
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

void output_particle_count(FILE *out, struct particle_count particle_count, struct world world) {
	output_vtk_header_count(out, world);
	int total = world.n[0] * world.n[1] * world.n[2];
	uint32_t *vtk_count = malloc(total * sizeof(uint32_t));
	for(int i = 0; i < total; i++)
		vtk_count[i] = __builtin_bswap32(particle_count.count[i]);
	fwrite(vtk_count, sizeof(uint32_t), total, out);
	free(vtk_count);
}

void output_volumetric_density(FILE *out, struct volumetric_density volume, struct world world) {
	output_vtk_header_volumetric(out, world);
	int total = world.n[0] * world.n[1] * world.n[2];
	uint64_t *vtk_density = malloc(total * sizeof(uint64_t));
	for(int i = 0; i < total; i++)
		vtk_density[i] = swap_endian(volume.density[i]);
	fwrite(vtk_density, sizeof(uint64_t), total, out);
	free(vtk_density);
}

void output_vtk_header_count(FILE *out, struct world world) {
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

void output_vtk_header_volumetric(FILE *out, struct world world) {
	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "Volumetric density\n");
	fprintf(out, "BINARY\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	fprintf(out, "SPACING %lf %lf %lf\n", 2.0 * world.d_max[0] / world.n[0], 2.0 * world.d_max[1] / world.n[1], 2.0 * world.d_max[2] / world.n[2]);
	fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
	fprintf(out, "SCALARS density double 1\n");
	fprintf(out, "LOOKUP_TABLE default\n");
}

void free_particles(struct test_particles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->kx); free(part->ky); free(part->kz);
	free(part->density_p); free(part->density_n);
	free(part->energy);
}

void free_particle_count(struct particle_count *part_count) {
	free(part_count->count);
}

void free_volumetric_density(struct volumetric_density *volume) {
	free(volume->density);
}

static inline uint64_t swap_endian(double v) {
	uint64_t data;
	memcpy(&data, &v, sizeof(double));
	return __builtin_bswap64(data);
}
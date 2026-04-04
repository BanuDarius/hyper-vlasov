#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "init.h"
#include "tools.h"
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
	int world_size = n * n * n;
	for(int i = 0; i < 3; i++) {
		world->n[i] = n;
		world->d_max[i] = d_max;
	}
}

void create_volumetric_density(struct volumetric_density *dens, struct world world) {
	int world_size = world.n[0] *world.n[1] * world.n[2] ;
	dens->density = malloc(world_size * sizeof(int));
	for(int i = 0; i < world_size; i++)
		dens->density[i] = 0;
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
	if(type == PROTONS) {
		start = 0; end = part.protons;
	}
	else {
		start = part.protons; end = part.protons + part.neutrons;
	}
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

void output_volumetric_density(FILE *out, struct volumetric_density dens, struct world world) {
	output_vtk_header(out, world);
	int total = world.n[0] * world.n[1] * world.n[2];
	uint32_t *vtk_density = malloc(total * sizeof(uint32_t));
	for(int i = 0; i < total; i++)
		vtk_density[i] = __builtin_bswap32(dens.density[i]);
	fwrite(vtk_density, sizeof(uint32_t), total, out);
	free(vtk_density);
}

void output_vtk_header(FILE *out, struct world world) {
	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "Volumetric density\n");
	fprintf(out, "BINARY\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	fprintf(out, "SPACING %lf %lf %lf\n", 2.0 * world.d_max[0] / world.n[0], 2.0 * world.d_max[1] / world.n[1], 2.0 * world.d_max[2] / world.n[2]);
	fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
	fprintf(out, "SCALARS density int 1\n");
	fprintf(out, "LOOKUP_TABLE default\n");
}

void free_particles(struct test_particles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->kx); free(part->ky); free(part->kz);
	free(part->density_p); free(part->density_n);
	free(part->energy);
}

void free_volumetric_density(struct volumetric_density *dens) {
	free(dens->density);
}
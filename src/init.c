#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void set_skyrme(struct skyrme *skm, double A, double B, double gamma) {
	skm->A = A;
	skm->B = B;
	skm->gamma = gamma;
}

void set_fermi_levels(struct fermi *fermi, double epsilon_p, double epsilon_n) {
	fermi->epsilon_p = epsilon_p;
	fermi->epsilon_n = epsilon_n;
}

void set_world(struct world *world, double d_max, int n) {
	for(int i = 0; i < 3; i++) {
		world->d_max[i] = d_max;
		world->n[i] = n;
	}
}

void create_particles(struct test_particles *part, int protons, int neutrons) {
	int num = protons + neutrons;
	part->protons = protons;
	part->neutrons = neutrons;
	part->x = malloc(num * sizeof(double));
	part->y = malloc(num * sizeof(double));
	part->z = malloc(num * sizeof(double));
	part->kx = malloc(num * sizeof(double));
	part->ky = malloc(num * sizeof(double));
	part->kz = malloc(num * sizeof(double));
	part->energy = malloc(num * sizeof(double));
	part->density = malloc(num * sizeof(double));
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
		fwrite(&part.density[i], sizeof(double), 1, out);
	}
}

void free_particles(struct test_particles *part) {
	free(part->x); free(part->y); free(part->z);
	free(part->kx); free(part->ky); free(part->kz);
	free(part->energy); free(part->density);
}
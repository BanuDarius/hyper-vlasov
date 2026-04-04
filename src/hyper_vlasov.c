#include <stdio.h>
#include <stdlib.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "sim_structs.h"

int main(int argc, char **argv) {
	srand(128);
	if(argc != 2) {
		printf("BAD!\n");
		return 1;
	}
	
	FILE *out = fopen(argv[1], "wb");
	int z = 10, n = 12, num_test_part = 3000, nx = 16;
	double V0 = -50.0, a = 0.66;
	double A = -356.0, B = 303.0, C = 32.0, gamma = 7.0 / 6.0;
	double epsilon_p = -8.0, epsilon_n = -12.0;
	double k_fwhm = 0.346, r_fwhm = 1.444;
	double sigma_k = calc_sigma(k_fwhm), sigma_r = calc_sigma(r_fwhm);
	double d_max = nuclear_radius(z + n);
	
	struct skyrme skm;
	struct world world;
	struct parameters param;
	struct fermi fermi_levels;
	struct woods_saxon ws[2];
	struct test_particles part;
	struct volumetric_density dens;
	
	set_skyrme(&skm, A, B, C, gamma);
	set_world(&world, d_max, nx);
	set_fermi_levels(&fermi_levels, epsilon_p, epsilon_n);
	set_parameters(&param, z, n, num_test_part, sigma_k, sigma_r);
	set_woods_saxon(&ws[0], V0, 0.8 * param.r_max, a);
	set_woods_saxon(&ws[1], V0, 0.8 * param.r_max, a);
	
	printf("MAX TEST PART %i\n", param.max_test_part);
	
	initialize_particles(&part, param, ws, skm, &fermi_levels);
	
	//output_centroids(out, part, NEUTRONS);
	create_volumetric_density(&dens, world);
	scatter_particles(&dens, &part, world, NEUTRONS);
	output_volumetric_density(out, dens, world);
	
	printf("FERMI P %lf FERMI N %lf\n", fermi_levels.epsilon_p, fermi_levels.epsilon_n);
	printf("Done\n");
	
	free_particles(&part);
	free_volumetric_density(&dens);
	fclose(out);
	return 0;
}
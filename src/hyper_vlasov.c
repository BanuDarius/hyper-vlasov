#include <stdio.h>
#include <stdlib.h>

#include "init.h"
#include "tools.h"
#include "sim_structs.h"

int main(int argc, char **argv) {
	srand(128);
	if(argc != 2) {
		printf("BAD!\n");
		return 1;
	}
	
	FILE *out = fopen(argv[1], "wb");
	int z = 50, n = 82, num_test_part = 200;
	double V0 = -50.0, a = 0.66;
	double A = -356.0, B = 303.0, gamma = 7.0 / 6.0;
	double epsilon_p = -8.0, epsilon_n = -12.0;
	double k_fwhm = 0.346, r_fwhm = 1.444;
	double sigma_k = calc_sigma(k_fwhm), sigma_r = calc_sigma(r_fwhm);
	
	struct skyrme skm;
	struct woods_saxon ws;
	struct parameters param;
	struct fermi fermi_levels;
	struct test_particles part_p, part_n;
	
	set_fermi_levels(&fermi_levels, epsilon_p, epsilon_n);
	set_parameters(&param, z, n, num_test_part, sigma_k, sigma_r);
	set_woods_saxon(&ws, param, V0, a);
	set_skyrme(&skm, A, B, gamma);
	
	printf("%i\n", param.max_test_part);
	printf("%lf %lf\n", param.sigma_k, param.sigma_r);
	
	initialize_particles(&part_p, &part_n, param, ws, skm, &fermi_levels);
	
	output_centroids(out, part_n, n * num_test_part);
	//output_centroids(out, part_n, param.max_test_part);
	
	printf("%lf %lf\n", fermi_levels.epsilon_p, fermi_levels.epsilon_n);
	
	free_particles(&part_p);
	free_particles(&part_n);
	printf("Done\n");
	fclose(out);
	return 0;
}
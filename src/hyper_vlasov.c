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
	int z = 10, n = 20, num_test_part = 100;
	double V0 = -50.0, a = 0.66;
	double A = -356.0, B = 303.0, gamma = 7.0 / 6.0;
	double epsilon_p = -8.0, epsilon_n = -12.0;
	
	struct skyrme skm;
	struct woods_saxon ws;
	struct parameters param;
	struct fermi fermi_init;
	struct test_particles part_init_p, part_init_n;
	struct energies en_p, en_n;
	
	set_fermi_levels(&fermi_init, epsilon_p, epsilon_n);
	set_parameters(&param, z, n, num_test_part);
	set_woods_saxon(&ws, param, V0, a);
	set_skyrme(&skm, A, B, gamma);
	
	printf("%i\n", param.max_test_part);
	
	create_particles(&part_init_p,  param.max_test_part);
	create_particles(&part_init_n,  param.max_test_part);
	initialize_particles(&part_init_p, &part_init_n, param, ws, skm);
	
	output_centroids(out, part_init_p, param.max_test_part);
	output_centroids(out, part_init_n, param.max_test_part);
	
	free_particles(&part_init_p);
	free_particles(&part_init_n);
	printf("Done\n");
	fclose(out);
	return 0;
}
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
	int z = 10, n = 10, num_test_part = 1000;
	double V0 = -50.0, a = 0.66;
	double A = -356.0, B = 303.0, gamma = 7.0 / 6.0;
	
	struct test_particles *part_init;
	struct parameters param;
	struct woods_saxon ws;
	struct skyrme skm;
	
	part_init = malloc(sizeof(struct test_particles));
	
	set_parameters(&param, z, n, num_test_part);
	set_woods_saxon(&ws, param, V0, a);
	set_skyrme(&skm, A, B, gamma);
	
	create_particles(part_init, param.max_test_part);
	initialize_particles(part_init, param, ws, skm);
	
	output_centroids(out, part_init, param);
	
	free_particles(part_init);
	printf("Done\n");
	fclose(out);
	return 0;
}
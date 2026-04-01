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
	
	struct test_particles *part;
	struct parameters param;
	
	part = malloc(sizeof(struct test_particles));
	
	set_parameters(&param, z, n, num_test_part);
	create_particles(part, param);
	generate_random_particles(part, param);
	
	//output_centroid_positions(out, part, param);
	//output_centroid_velocities(out, part, param);
	output_centroids(out, part, param);
	
	free_particles(part);
	printf("Done\n");
	fclose(out);
	return 0;
}
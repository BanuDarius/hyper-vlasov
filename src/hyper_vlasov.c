#include <stdio.h>
#include <stdlib.h>

#include "init.h"
#include "tools.h"
#include "sim_structs.h"

int main(int argc, char **argv) {
	printf("here");
	srand(128);
	struct test_particles *part;
	struct parameters param;
	
	param.num_test_part = 10 * 10;
	
	part = malloc(sizeof(struct test_particles));
	
	create_particles(part, param);
	
	generate_random_particles(part, param);
	
	for(int i = 0; i < 100; i++) {
		double r = nuclear_radius(i);
		printf("%lf\n", r);
	}
	
	free_particles(part);
	printf("Done\n");
	return 0;
}
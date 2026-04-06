#include <omp.h>
#include <math.h>
#include <stdio.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "sim_structs.h"

int main(int argc, char **argv) {
	srand(128);
	if(argc != 3) {
		printf("BAD!\n");
		return 1;
	}
	double start_time = omp_get_wtime();
	FILE *in = fopen(argv[1], "r"), *out = fopen(argv[2], "wb");
	
	struct skyrme skm;
	struct parameters param;
	struct woods_saxon ws[2];
	struct fermi fermi_levels;
	struct test_particles part;
	struct particle_count part_count;
	struct world world, world_visual;
	struct volumetric_density volume;
	
	read_input_file(in, &skm, &world, &world_visual, &fermi_levels, &param, ws);
	
	create_particle_count(&part_count, world);
	create_volumetric_density(&volume, world_visual);
	printf("MAX TEST PART %i\n", param.max_test_part);
	
	printf("Simulation started.\n");
	initialize_particles(&part, param, ws, skm, &fermi_levels);
	
	chi_squared(part, ws, skm, param.test_part_per_nucleon);
	double msr_p = mean_squared_radius(part, PROTONS);
	double msr_n = mean_squared_radius(part, NEUTRONS);

	scatter_particles(&part_count, &part, world);
	compute_volumetric_density(&volume, part_count, world_visual, world, param, PROTONS_AND_NEUTRONS);
	
	//output_centroids(out, part, PROTONS);
	//output_particle_count(out, part_count, world);
	output_volumetric_density(out, volume, world_visual);
	
	printf("FERMI P %0.2lf FERMI N %0.2lf\n", fermi_levels.epsilon_p, fermi_levels.epsilon_n);
	printf("RADIUS N %0.2lf RADIUS P %0.2lf\n", sqrt(msr_n), sqrt(msr_p));
	printf("Simulation ended.\n");
	printf("Time taken: %0.3lfs\n", omp_get_wtime() - start_time);
	
	fclose(out); fclose(in);
	free_particles(&part);
	free_particle_count(&part_count);
	free_volumetric_density(&volume);
	return 0;
}
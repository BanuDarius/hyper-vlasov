#include <omp.h>
#include <math.h>
#include <stdio.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "sim_structs.h"

void simulate(FILE *out, struct test_particles *part, struct woods_saxon *ws, struct skyrme skm, struct parameters param, struct world world, struct world world_visual) {
	chi_squared(*part, ws, skm, param.test_part_per_nucleon);
	double msr_p = mean_squared_radius(*part, PROTONS);
	double msr_n = mean_squared_radius(*part, NEUTRONS);
	
	//struct particle_count part_count;
	struct volumetric_density volume;
	
	//create_particle_count(&part_count, world);
	create_volumetric_density(&volume, world_visual);
	
	//scatter_particles(&part_count, part, world);
	//compute_volumetric_density(&volume, part_count, world_visual, world, param, PROTONS_AND_NEUTRONS);
	distribute_particles_cic(&volume, part, world, PROTONS);
	
	double x = 0.0;
	for(int i = 0; i < world.n[0] * world.n[1] * world.n[2]; i++) {
		x += volume.density[i];
	}
	printf("TOTAL SUM %0.2lf\n", x);
	
	//output_centroids(out, part, PROTONS);
	//output_particle_count(out, part_count, world);
	output_volumetric_density(out, volume, world);
	printf("RADIUS N %0.2lf RADIUS P %0.2lf\n", sqrt(msr_n), sqrt(msr_p));
	
	//free_particle_count(&part_count);
	free_volumetric_density(&volume);
}

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
	struct world world, world_visual;
	
	read_input_file(in, &skm, &world, &world_visual, &fermi_levels, &param, ws);
	printf("MAX TEST PART %i\n", param.max_test_part);
	
	printf("Simulation started.\n");
	initialize_particles(&part, param, ws, skm, &fermi_levels);
	simulate(out, &part, ws, skm, param, world, world_visual);
	printf("Simulation ended.\n");
	
	printf("FERMI P %0.2lf FERMI N %0.2lf\n", fermi_levels.epsilon_p, fermi_levels.epsilon_n);
	printf("Time taken: %0.3lfs\n", omp_get_wtime() - start_time);
	
	fclose(out); fclose(in);
	free_particles(&part);
	return 0;
}
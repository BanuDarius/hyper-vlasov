/* MIT License

Copyright (c) 2026 Banu Darius-Matei

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */

#include <math.h>
#include <stdio.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "sim_structs.h"

void simulate(FILE *out, TestParticles *part, WoodsSaxon *ws, Skyrme skm, Parameters param, World world, World world_visual) {
	chi_squared(*part, ws, skm, param.part_per_nucleon);
	double msr_p = mean_squared_radius(*part, PROTONS);
	double msr_n = mean_squared_radius(*part, NEUTRONS);
	
	ScalarField volume;
	create_volumetric_density(&volume, world);
	compute_volumetric_density_cic(&volume, part, param, world);
	
	output_volumetric_density(out, volume, world);
	printf("RADIUS N %0.2lf RADIUS P %0.2lf\n", sqrt(msr_n), sqrt(msr_p));
	
	free_scalar_field(&volume);
	//ParticleCount part_count;
	//create_particle_count(&part_count, world);
	//scatter_particles(&part_count, part, world);
	//free_particle_count(&part_count);
	//compute_volumetric_density(&volume, part_count, world_visual, world, param, PROTONS_AND_NEUTRONS);
	//output_centroids(out, part, PROTONS);
	//output_particle_count(out, part_count, world);
}

int main(int argc, char **argv) {
	srand(128);
	if(argc != 3) {
		printf("BAD!\n");
		return 1;
	}
	double start_time = omp_get_wtime();
	FILE *in = fopen(argv[1], "r"), *out = fopen(argv[2], "wb");
	
	Skyrme skm;
	Parameters param;
	WoodsSaxon ws[2];
	Fermi fermi_levels;
	TestParticles part;
	World world, world_visual;
	
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
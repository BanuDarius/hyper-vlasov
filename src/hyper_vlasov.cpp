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

#include <omp.h>
#include <cmath>
#include <cstdio>

#include "init.hpp"
#include "tools.hpp"
#include "physics.hpp"
#include "sim_structs.hpp"

template <typename T>
void simulate(char *output_directory, TestParticles<T> *part, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world) {
	VectorField<T> forces;
	create_vector_field(&forces, world);
	
	ScalarField<T> volume, potentials, coulomb;
	create_scalar_field_single(&coulomb, world);
	create_scalar_field_double(&volume, world);
	create_scalar_field_double(&potentials, world);
	set_initial_coulomb_boundaries(&coulomb, world, param.z);
	
	compute_volumetric_density_cic(&volume, part, param, world);
	
	compute_volumetric_skyrme_potentials(&potentials, volume, skm, world);
	compute_volumetric_coulomb_potentials_sor(&coulomb, volume, world);
	merge_volumetric_potentials(&potentials, coulomb, world);
	compute_volumetric_forces_fdm(&forces, potentials, world);
	
	distribute_forces_to_particles_cic(part, forces, world);
	char stats_filename[32];
	std::sprintf(stats_filename, "%sstats.txt", output_directory);
	FILE *stats = fopen(stats_filename, "w");
	
	for(int step = 1; step <= param.steps; step++) {
		char output_filename[32];
		set_output_filename(output_filename, output_directory, step);
		FILE *out = fopen(output_filename, "wb");
		if(out != NULL) {
			output_vtk_header_start(out, world);
			output_scalar_field(out, volume, world, "density");
			output_scalar_field(out, potentials, world, "potentials");
			output_vector_field(out, forces, world, "forces");
			fclose(out);
		}
		else {
			std::fprintf(stderr, "CANNOT OPEN FILE!\n");
			exit(1);
		}
		T dt = param.t_f / param.steps;
		update_momenta_half(part, dt);
		update_positions_full(part, dt);
		
		compute_volumetric_density_cic(&volume, part, param, world);
		
		compute_volumetric_skyrme_potentials(&potentials, volume, skm, world);
		compute_volumetric_coulomb_potentials_sor(&coulomb, volume, world);
		merge_volumetric_potentials(&potentials, coulomb, world);
		compute_volumetric_forces_fdm(&forces, potentials, world);
		
		distribute_forces_to_particles_cic(part, forces, world);
		
		update_momenta_half(part, dt);
		
		std::printf("TIME STEP %i/%i\n", step, param.steps);
		T msr_p = mean_squared_radius(*part, PROTONS);
		T msr_n = mean_squared_radius(*part, NEUTRONS);
		std::fprintf(stats, "%0.4lf %0.4lf %0.4lf\n", std::sqrt(msr_n), std::sqrt(msr_p), step * dt);
	}
	free_vector_field(&forces);
	free_scalar_field(&volume);
	free_scalar_field(&coulomb);
	free_scalar_field(&potentials);
	fclose(stats);
	//ParticleCount<T> part_count;
	//create_particle_count(&part_count, world);
	//scatter_particles(&part_count, part, world);
	//compute_volumetric_density(&volume, part_count, world_visual, world, param, PROTONS_AND_NEUTRONS);
	//output_centroids(out, part, PROTONS);
	//output_particle_count(out, part_count, world);
	//free_particle_count(&part_count);
}

int main(int argc, char **argv) {
	srand(128);
	if(argc != 3) {
		std::fprintf(stderr, "NEED 3 ARGUMENTS!\n");
		return 1;
	}
	
	Skyrme<double> skm;
	World<double> world;
	Parameters<double> param;
	WoodsSaxon<double> ws[2];
	Fermi<double> fermi_levels;
	TestParticles<double> part;
	
	FILE *in = fopen(argv[1], "r");
	read_input_file(in, &skm, &world, &fermi_levels, &param, ws);
	std::printf("MAX TEST PART %i\n", param.max_test_part);
	
	std::printf("Simulation started.\n");
	double start_time = omp_get_wtime();
	initialize_particles(&part, param, ws, skm, &fermi_levels);
	chi_squared(part, ws, skm, param.part_per_nucleon);
	
	simulate(argv[2], &part, skm, param, world);
	std::printf("Simulation ended.\n");
	std::printf("Time taken: %0.3lfs\n", omp_get_wtime() - start_time);
	
	free_particles(&part);
	fclose(in);
	return 0;
}
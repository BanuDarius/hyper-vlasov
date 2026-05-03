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
void cpu_simulate(const char *output_directory, TestParticles<T> *part, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world) {
	bool excited_nucleus = false;
	T dt = param.t_f / param.steps;
	
	char stats_filename[STRING_SIZE];
	set_stats_filename(stats_filename, output_directory);
	FILE *stats = fopen(stats_filename, "w");
	if(stats == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN STATS FILE!\n"); exit(1);
	}
	
	VectorField<T> forces;
	create_vector_field_double(&forces, world);
	
	ScalarField<T> potentials, coulomb, temp_density, density;
	create_scalar_field_single(&coulomb, world);
	create_scalar_field_double(&density, world);
	create_scalar_field_double(&potentials, world);
	create_scalar_field_double(&temp_density, world);
	
	distribute_volumetric_particles_cic(&density, part, world);
	compute_volumetric_densities(&density, &temp_density, param, world);
	
	set_initial_coulomb_boundaries(&coulomb, world, param.z);
	compute_volumetric_skyrme_potentials(&potentials, density, skm, world);
	compute_volumetric_coulomb_potentials_sor(&coulomb, density, world);
	merge_volumetric_potentials(&potentials, coulomb, world);
	compute_volumetric_forces_fdm(&forces, potentials, world);
	
	distribute_forces_to_particles_cic(part, forces, world);
	for(int step = 0; step < param.steps; step++) {
		if(step % param.substeps == 0) {
			T x_p = center_of_mass(*part, PROTONS);
			T x_n = center_of_mass(*part, NEUTRONS);
			T msr_p = mean_squared_radius(*part, world, PROTONS);
			T msr_n = mean_squared_radius(*part, world, NEUTRONS);
			std::fprintf(stats, "%e %e %e %e %e\n",
			step * dt, std::sqrt(msr_p), std::sqrt(msr_n), x_p, x_n);
			
			char output_filename[STRING_SIZE];
			set_output_filename(output_filename, output_directory, step / param.substeps);
			FILE *out = fopen(output_filename, "wb");
			if(out == nullptr) {
				std::fprintf(stderr, "CANNOT OPEN OUTPUT FILE!\n");
				exit(1);
			}
			output_vtk_header_start(out, world);
			output_vector_field(out, forces, world, "forces");
			output_scalar_field(out, density, world, "density");
			output_scalar_field(out, potentials, world, "potentials");
			std::printf("Processed step: %i/%i.\n", step, param.steps);
			fclose(out);
		}
		if(step * dt >= param.t_exc && !excited_nucleus) {
			excited_nucleus = true;
			nuclear_excitation(part, param);
		}
		update_momenta_half(part, dt);
		update_positions_full(part, dt);
		
		distribute_volumetric_particles_cic(&density, part, world);
		compute_volumetric_densities(&density, &temp_density, param, world);
		
		compute_volumetric_skyrme_potentials(&potentials, density, skm, world);
		compute_volumetric_coulomb_potentials_sor(&coulomb, density, world);
		merge_volumetric_potentials(&potentials, coulomb, world);
		compute_volumetric_forces_fdm(&forces, potentials, world);
		
		distribute_forces_to_particles_cic(part, forces, world);
		
		update_momenta_half(part, dt);
	}
	fclose(stats);
	free_vector_field(&forces);
	free_scalar_field(&density);
	free_scalar_field(&coulomb);
	free_scalar_field(&potentials);
	free_scalar_field(&temp_density);
	/*ParticleCount<T> part_count;
	create_particle_count(&part_count, world);
	scatter_particles(&part_count, part, world);
	compute_volumetric_density(&density, part_count, world_visual, world, param, PROTONS_AND_NEUTRONS);
	output_centroids(out, part, PROTONS);
	output_particle_count(out, part_count, world);
	free_particle_count(&part_count);*/
}

/*template <typename T>
void gpu_simulate(const char *output_directory, TestParticles<T> *part, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world) {
	return;
}*/ //Will add this GPU function after all the CPU functions have been implemented

template <typename T>
void run_simulation(const char *input_filename, const char *output_filename) {
	Skyrme<T> skm;
	World<T> world;
	Parameters<T> param;
	WoodsSaxon<T> ws[2];
	Fermi<T> fermi_levels;
	TestParticles<T> part;
	
	FILE *in = fopen(input_filename, "r");
	if(in == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN INPUT FILE!\n"); exit(1);
	}
	read_input_file(in, &skm, &world, &fermi_levels, &param, ws);
	std::printf("MAX TEST PART %i\n", param.max_test_part);
	
	initialize_particles(&part, ws, skm, &fermi_levels, param);
	chi_squared(part, ws, skm, param.part_per_nucleon);
	if(param.use_gpu == true)
		//gpu_simulate(output_filename, &part, skm, param, world);
		return;
	else
		cpu_simulate(output_filename, &part, skm, param, world);
	
	free_particles(&part);
	fclose(in);
}

int main(int argc, char **argv) {
	srand(128);
	if(argc < 3) {
		std::fprintf(stderr, "BAD ARGUMENTS!\n");
		return 1;
	}
	std::printf("Simulation started.\n");
	double start_time = omp_get_wtime();
	
	if(!std::strcmp(argv[1], "--float"))
		run_simulation<float>(argv[2], argv[3]);
	else
		run_simulation<double>(argv[1], argv[2]);
	
	std::printf("Simulation ended.\n");
	std::printf("Time taken: %0.3lfs\n", omp_get_wtime() - start_time);
	return 0;
}
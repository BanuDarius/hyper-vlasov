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
void cpu_simulate(const char *output_directory, TestParticles<T> &part, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world) {
	bool excited_nucleus = false;
	T dt = param.t_f / param.steps;
	std::vector<float> samples(2 * param.density_samples);
	
	char stats_filename[string_size], density_samples_filename[string_size], density_samples_diff_filename[string_size];
	std::sprintf(stats_filename, "%s/%s", output_directory, "stats.txt");
	std::sprintf(density_samples_filename, "%s/%s", output_directory, "density_samples.bin");
	std::sprintf(density_samples_diff_filename, "%s/%s", output_directory, "density_samples_diff.bin");
	
	FILE *out_stats = fopen(stats_filename, "w");
	FILE *out_samples = fopen(density_samples_filename, "wb");
	FILE *out_samples_diff = fopen(density_samples_diff_filename, "wb");
	if(out_stats == nullptr || out_samples == nullptr || out_samples_diff == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN STATS FILES!\n"); exit(1);
	}
	output_density_samples_positions(out_samples, param, world);
	output_density_samples_positions(out_samples_diff, param, world);
	
	int world_size = world.n[0] * world.n[1] * world.n[2];
	VectorField<T> forces(2 * world_size);
	ScalarField<T> coulomb(world_size), potentials(2 * world_size), density_temp(2 * world_size), density_before(2 * world_size), density(2 * world_size);
	
	distribute_volumetric_particles_cic(density, part, world);
	compute_volumetric_densities(density, density_temp, param, world);
	
	compute_coulomb_boundaries(coulomb, part, world, param.z);
	compute_volumetric_skyrme_potentials(potentials, density, skm, world);
	compute_volumetric_coulomb_potentials_sor(coulomb, density, world);
	add_scalar_field_single(potentials, potentials, coulomb, world);
	compute_volumetric_forces_fdm(forces, potentials, world);
	
	distribute_forces_to_particles_cic(part, forces, world);
	for(int step = 0; step < param.steps; step++) {
		if(step % param.substeps == 0) {
			std::printf("Processed step: %i/%i.\n", step, param.steps);
			std::array<T, 3> x_p = center_of_mass(part, world, is_proton);
			std::array<T, 3> x_n = center_of_mass(part, world, is_neutron);
			T msr_p = mean_squared_radius(part, world, is_proton);
			T msr_n = mean_squared_radius(part, world, is_neutron);
			std::fprintf(out_stats, "%e %e %e %e %e\n",
			step * dt, std::sqrt(msr_p), std::sqrt(msr_n), x_p[2], x_n[2]);
			
			compute_density_samples_cic(samples, density, param, world);
			output_density_samples(out_samples, samples, param);
			
			char volume_filename[string_size];
			std::sprintf(volume_filename, "%s/out-%04d.vtk", output_directory, step / param.substeps);
			FILE *out_volume = fopen(volume_filename, "wb");
			if(out_volume == nullptr) {
				std::fprintf(stderr, "CANNOT OPEN VOLUME FILE!\n"); exit(1);
			}
			output_vtk_header_start(out_volume, world);
			output_vector_field(out_volume, forces, world, "forces");
			output_scalar_field(out_volume, density, world, "density");
			output_scalar_field(out_volume, potentials, world, "potentials");
			if(excited_nucleus) {
				sub_scalar_field_double(density_temp, density, density_before, world);
				output_scalar_field(out_volume, density_temp, world, "density_difference");
				
				compute_density_samples_cic(samples, density_temp, param, world);
				output_density_samples(out_samples_diff, samples, param);
			}
			fclose(out_volume);
		}
		if(step * dt >= param.t_exc && !excited_nucleus) {
			excited_nucleus = true;
			nuclear_excitation(part, param);
			copy_scalar_field_double(density_before, density, world);
		}
		if(step % reset_steps == 0)
			center_momentum(part, world);
		
		update_momenta_half(part, dt);
		update_positions_full(part, dt);
		
		distribute_volumetric_particles_cic(density, part, world);
		compute_volumetric_densities(density, density_temp, param, world);
		
		compute_coulomb_boundaries(coulomb, part, world, param.z);
		compute_volumetric_skyrme_potentials(potentials, density, skm, world);
		compute_volumetric_coulomb_potentials_sor(coulomb, density, world);
		add_scalar_field_single(potentials, potentials, coulomb, world);
		compute_volumetric_forces_fdm(forces, potentials, world);
		
		distribute_forces_to_particles_cic(part, forces, world);
		
		update_momenta_half(part, dt);
	}
	fclose(out_stats);
	fclose(out_samples);
	/*ParticleCount<T> part_count;
	create_particle_count(&part_count, world);
	scatter_particles(&part_count, part, world);
	compute_volumetric_density(&density, part_count, world_visual, world, param, is_proton_AND_is_neutron);
	output_centroids(out, part, is_proton);
	output_particle_count(out, part_count, world);
	free_particle_count(&part_count);*/
}

template <typename T>
void run_simulation(const char *input_filename, const char *output_filename) {
	FILE *in = fopen(input_filename, "r");
	if(in == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN INPUT FILE!\n"); exit(1);
	}
	Skyrme<T> skm;
	World<T> world;
	Parameters<T> param;
	WoodsSaxon<T> ws;
	Fermi<T> fermi_levels;
	read_input_file(in, skm, world, fermi_levels, param, ws);
	fclose(in);
	
	TestParticles<T> part(param.z * param.part_per_nucleon, param.n * param.part_per_nucleon);
	std::printf("MAX TEST PART %i\n", param.max_test_part);
	
	initialize_particles(part, ws, skm, fermi_levels, param);
	chi_squared(part, ws, skm, param.part_per_nucleon);
	if(param.use_gpu == true)
		return;
	else
		cpu_simulate(output_filename, part, skm, param, world);
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
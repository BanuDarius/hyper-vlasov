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

#include <cmath>
#include <cstdlib>

#include "simulation.hpp"
#include "init.hpp"
#include "tools.hpp"
#include "physics.hpp"
#include "particle_in_cell.hpp"

template <typename T>
void cpu_simulate(const char *output_directory, TestParticles<T> &part, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world) {
	bool excited_nucleus = false;
	T dt = param.t_f / param.steps;
	std::unique_ptr<float[]> samples(new float[2 * param.density_samples]);
	#pragma omp parallel for simd schedule(static)
	for(int i = 0; i < 2 * param.density_samples; i++)
		samples[i] = 0.0f;
	
	char stats_filename[string_size], density_samples_filename[string_size], density_samples_diff_filename[string_size];
	std::sprintf(stats_filename, "%s/%s", output_directory, "stats.txt");
	std::sprintf(density_samples_filename, "%s/%s", output_directory, "density_samples.bin");
	std::sprintf(density_samples_diff_filename, "%s/%s", output_directory, "density_samples_diff.bin");
	
	FILE *out_stats = fopen(stats_filename, "w");
	FILE *out_samples = fopen(density_samples_filename, "wb");
	FILE *out_samples_diff = fopen(density_samples_diff_filename, "wb");
	if(out_stats == nullptr || out_samples == nullptr || out_samples_diff == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN STATS FILES!\n"); std::exit(1);
	}
	output_density_samples_positions(out_samples, param, world);
	output_density_samples_positions(out_samples_diff, param, world);
	
	int world_size = world.n[0] * world.n[1] * world.n[2];
	VectorField<T> forces(2 * world_size), current_velocity(2 * world_size), temp_vector_field(2 * world_size);
	ScalarField<T> coulomb(2 * world_size), potentials(2 * world_size), density_before(2 * world_size), density(2 * world_size), temp_scalar_field(2 * world_size);
	
	distribute_volumetric_particles_cic(density, part, world);
	compute_volumetric_densities(density, temp_scalar_field, param, world);
	
	compute_coulomb_boundaries(coulomb, part, world, param.z);
	compute_volumetric_skyrme_potentials(potentials, density, skm, world);
	compute_volumetric_coulomb_potentials_sor(coulomb, density, world);
	potentials += coulomb;
	compute_volumetric_forces_fdm(forces, potentials, world);
	
	distribute_forces_to_particles_cic(part, forces, world);
	for(int step = 0; step < param.steps; step++) {
		update_momenta_half(part, dt);
		update_positions_full(part, dt);
		
		distribute_volumetric_particles_cic(density, part, world);
		compute_volumetric_densities(density, temp_scalar_field, param, world);
		
		compute_coulomb_boundaries(coulomb, part, world, param.z);
		compute_volumetric_skyrme_potentials(potentials, density, skm, world);
		compute_volumetric_coulomb_potentials_sor(coulomb, density, world);
		potentials += coulomb;
		compute_volumetric_forces_fdm(forces, potentials, world);
		
		distribute_forces_to_particles_cic(part, forces, world);
		
		update_momenta_half(part, dt);
		if(step % param.substeps == 0) {
			std::printf("Processed step: %i/%i.\n", step, param.steps);
			std::array<T, 3> x_p = center_of_mass(part, world, is_proton);
			std::array<T, 3> x_n = center_of_mass(part, world, is_neutron);
			T msr_p = mean_squared_radius(part, world, is_proton);
			T msr_n = mean_squared_radius(part, world, is_neutron);
			std::fprintf(out_stats, "%e %e %e %e %e\n",
			step * dt, std::sqrt(msr_p), std::sqrt(msr_n), x_p[2], x_n[2]);
			
			distribute_volumetric_momenta_cic(current_velocity, part, world);
			compute_current_velocity(current_velocity, temp_vector_field, density, param, world);
			
			compute_density_samples_cic(samples.get(), density, param, world);
			output_density_samples(out_samples, samples.get(), param);
			
			char volume_filename[string_size];
			std::sprintf(volume_filename, "%s/out-%04d.vtk", output_directory, step / param.substeps);
			FILE *out_volume = fopen(volume_filename, "wb");
			if(out_volume == nullptr) {
				std::fprintf(stderr, "CANNOT OPEN VOLUME FILE!\n"); std::exit(1);
			}
			output_vtk_header_start(out_volume, world);
			output_vector_field(out_volume, forces, world, "forces");
			output_scalar_field(out_volume, density, world, "density");
			output_scalar_field(out_volume, potentials, world, "potentials");
			output_vector_field(out_volume, current_velocity, world, "current");
			if(excited_nucleus) {
				temp_scalar_field = density;
				temp_scalar_field -= density_before;
				output_scalar_field(out_volume, temp_scalar_field, world, "density_difference");
				
				compute_density_samples_cic(samples.get(), temp_scalar_field, param, world);
				output_density_samples(out_samples_diff, samples.get(), param);
			}
			fclose(out_volume);
		}
		if(step * dt >= param.t_exc && !excited_nucleus) {
			excited_nucleus = true;
			nuclear_excitation(part, param);
			density_before = density;
		}
		if(step % reset_steps == 0)
			center_momentum(part, world);
	}
	fclose(out_stats);
	fclose(out_samples);
	fclose(out_samples_diff);
}

template void cpu_simulate<double>(const char *output_directory, TestParticles<double> &part, const Skyrme<double> &skm, const Parameters<double> &param, const World<double> &world);

template void cpu_simulate<float>(const char *output_directory, TestParticles<float> &part, const Skyrme<float> &skm, const Parameters<float> &param, const World<float> &world);
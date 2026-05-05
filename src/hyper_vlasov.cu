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

#include "init.cuh"
#include "physics.cuh"

template <typename T>
void gpu_simulate(const char *output_directory, TestParticles<T> *h_part, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world) {
	TestParticles<T> d_part;
	gpu_create_particles(&h_part, param.part_per_nucleon * param.z, param.part_per_nucleon * param.n);
	copy_particles_cpu_to_gpu(&h_part, d_part);
	bool excited_nucleus = false;
	T dt = param.t_f / param.steps;
	
	char stats_filename[STRING_SIZE];
	set_stats_filename(stats_filename, output_directory);
	FILE *stats = fopen(stats_filename, "w");
	if(stats == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN STATS FILE!\n"); exit(1);
	}
	VectorField<T> d_forces;
	gpu_create_vector_field_double(&d_forces, world);
	
	ScalarField<T> d_potentials, d_coulomb, d_temp_density, d_density;
	gpu_create_scalar_field_single(&d_coulomb, world);
	gpu_create_scalar_field_double(&d_density, world);
	gpu_create_scalar_field_double(&d_potentials, world);
	gpu_create_scalar_field_double(&d_temp_density, world);
	
	VectorField<T> h_forces;
	create_vector_field_double(&h_forces, world);
	
	ScalarField<T> h_potentials, h_density;
	create_scalar_field_double(&h_density, world);
	create_scalar_field_double(&h_potentials, world);
	
	distribute_volumetric_particles_cic(&density, part, world);
	compute_volumetric_densities(&density, &temp_density, param, world);
	
	compute_coulomb_boundaries(&coulomb, *part, world, param.z);
	compute_volumetric_skyrme_potentials(&potentials, density, skm, world);
	compute_volumetric_coulomb_potentials_sor(&coulomb, density, world);
	merge_volumetric_potentials(&potentials, coulomb, world);
	compute_volumetric_forces_fdm(&forces, potentials, world);
	
	distribute_forces_to_particles_cic(part, forces, world);
	for(int step = 0; step < param.steps; step++) {
		if(step % param.substeps == 0) {
			T x_p = center_of_mass(*part, world, PROTONS);
			T x_n = center_of_mass(*part, world, NEUTRONS);
			T msr_p = mean_squared_radius(*part, world, PROTONS);
			T msr_n = mean_squared_radius(*part, world, NEUTRONS);
			std::fprintf(stats, "%e %e %e %e %e\n",
			step * dt, std::sqrt(msr_p), std::sqrt(msr_n), x_p, x_n);
			
			char output_filename[STRING_SIZE];
			set_output_filename(output_filename, output_directory, step / param.substeps);
			FILE *out = fopen(output_filename, "wb");
			if(out == nullptr) {
				std::fprintf(stderr, "CANNOT OPEN OUTPUT FILE!\n"); exit(1);
			}
			copy_scalar_field_gpu_to_cpu(&h_density, &d_density);
			copy_scalar_field_gpu_to_cpu(&h_potentials, &d_potentials);
			copy_vector_field_gpu_to_cpu(&h_forces, &d_forces);
			
			output_vtk_header_start(out, world);
			output_vector_field(out, h_forces, world, "forces");
			output_scalar_field(out, h_density, world, "density");
			output_scalar_field(out, h_potentials, world, "potentials");
			std::printf("Processed step: %i/%i.\n", step, param.steps);
			fclose(out);
		}
		if(step * dt >= param.t_exc && !excited_nucleus) {
			excited_nucleus = true;
			nuclear_excitation(part, param);
		}
		if(step % RESET_STEPS == 0)
			gpu_center_momentum(part);
		
		gpu_update_momenta_half(part, dt);
		gpu_update_positions_full(part, dt);
		
		distribute_volumetric_particles_cic(&density, part, world);
		compute_volumetric_densities(&density, &temp_density, param, world);
		
		compute_coulomb_boundaries(&coulomb, *part, world, param.z);
		compute_volumetric_skyrme_potentials(&potentials, density, skm, world);
		compute_volumetric_coulomb_potentials_sor(&coulomb, density, world);
		merge_volumetric_potentials(&potentials, coulomb, world);
		compute_volumetric_forces_fdm(&forces, potentials, world);
		
		distribute_forces_to_particles_cic(part, forces, world);
		
		gpu_update_momenta_half(part, dt);
	}
	fclose(stats);
	gpu_free_particles(&h_part);
	gpu_free_vector_field(&forces);
	gpu_free_scalar_field(&density);
	gpu_free_scalar_field(&coulomb);
	gpu_free_scalar_field(&potentials);
	gpu_free_scalar_field(&temp_density);
}
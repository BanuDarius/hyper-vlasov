// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <concepts>

#include "init.hpp"
#include "physics.hpp"
#include "simulation.hpp"

template <std::floating_point T>
void run_simulation(const char *input_filename, const char *output_filename) {
	Skyrme<T> skm;
	World<T> world;
	WoodsSaxon<T> ws;
	Parameters<T> param;
	Fermi<T> fermi_levels;
	std::FILE *in = fopen(input_filename, "r");
	if(in == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN INPUT FILE!\n"); std::exit(1);
	}
	read_input_file(in, skm, world, fermi_levels, param, ws);
	std::fclose(in);
	
	TestParticles<T> part(param.z * param.part_per_nucleon, param.n * param.part_per_nucleon);
	std::printf("MAX TEST PART %i\n", param.max_test_part);
	
	initialize_particles(part, ws, skm, fermi_levels, param, world);
	if(param.use_gpu)
		return;
	else
		cpu_simulate(output_filename, part, skm, param, world);
}

int main(int argc, char **argv) {
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
	std::printf("Time taken: %0.3lfs.\n", omp_get_wtime() - start_time);
	return 0;
}
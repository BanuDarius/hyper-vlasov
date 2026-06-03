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
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "init.hpp"
#include "physics.hpp"
#include "simulation.hpp"

template <typename T>
void run_simulation(const char *input_filename, const char *output_filename) {
	Skyrme<T> skm;
	World<T> world;
	WoodsSaxon<T> ws;
	Parameters<T> param;
	Fermi<T> fermi_levels;
	FILE *in = fopen(input_filename, "r");
	if(in == nullptr) {
		std::fprintf(stderr, "CANNOT OPEN INPUT FILE!\n"); std::exit(1);
	}
	read_input_file(in, skm, world, fermi_levels, param, ws);
	fclose(in);
	
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
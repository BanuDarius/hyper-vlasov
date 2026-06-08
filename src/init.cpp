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

#include <cstring>
#include <cstdlib>

#include "init.hpp"
#include "physics_formulas.hpp"

template <typename T>
void set_parameters(Parameters<T> &param, int z, int n, int part_per_nucleon, int steps, int substeps, int density_samples, bool use_gpu, T sample_position, T sigma_k, T sigma_r, T t_f, T t_exc, T eta_exc) {
	param.z = z;
	param.n = n;
	param.t_f = t_f;
	param.t_exc = t_exc;
	param.steps = steps;
	param.use_gpu = use_gpu;
	param.sigma_k = sigma_k;
	param.sigma_r = sigma_r;
	param.eta_exc = eta_exc;
	param.substeps = substeps;
	param.r_max = nuclear_radius<T>(z + n);
	param.density_samples = density_samples;
	param.sample_position = sample_position;
	param.part_per_nucleon = part_per_nucleon;
	param.max_test_part = max_particles(param.r_max, k_max<T>, param.part_per_nucleon);
}

template <typename T>
void read_input_file(std::FILE *in, Skyrme<T> &skm, World<T> &world, Fermi<T> &fermi_levels, Parameters<T> &param, WoodsSaxon<T> &ws) {
	double V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, t_f, t_exc, eta_exc, d_max_scale, sample_position;
	int i = 0, num_test_part, use_gpu, density_samples, substeps, steps, nx, z, n;
	char current[string_size];
	
	while(std::fscanf(in, "%s", current) != EOF) {
		if(!std::strcmp(current, "V0"))
			i += std::fscanf(in, "%lf", &V0);
		else if(!std::strcmp(current, "a"))
			i += std::fscanf(in, "%lf", &a);
		else if(!std::strcmp(current, "A"))
			i += std::fscanf(in, "%lf", &A);
		else if(!std::strcmp(current, "B"))
			i += std::fscanf(in, "%lf", &B);
		else if(!std::strcmp(current, "C"))
			i += std::fscanf(in, "%lf", &C);
		else if(!std::strcmp(current, "gamma"))
			i += std::fscanf(in, "%lf", &gamma);
		else if(!std::strcmp(current, "epsilon_p"))
			i += std::fscanf(in, "%lf", &epsilon_p);
		else if(!std::strcmp(current, "epsilon_n"))
			i += std::fscanf(in, "%lf", &epsilon_n);
		else if(!std::strcmp(current, "k_fwhm"))
			i += std::fscanf(in, "%lf", &k_fwhm);
		else if(!std::strcmp(current, "r_fwhm"))
			i += std::fscanf(in, "%lf", &r_fwhm);
		else if(!std::strcmp(current, "t_f"))
			i += std::fscanf(in, "%lf", &t_f);
		else if(!std::strcmp(current, "t_exc"))
			i += std::fscanf(in, "%lf", &t_exc);
		else if(!std::strcmp(current, "eta_exc"))
			i += std::fscanf(in, "%lf", &eta_exc);
		else if(!std::strcmp(current, "d_max_scale"))
			i += std::fscanf(in, "%lf", &d_max_scale);
		else if(!std::strcmp(current, "sample_position"))
			i += std::fscanf(in, "%lf", &sample_position);
		else if(!std::strcmp(current, "nx"))
			i += std::fscanf(in, "%i", &nx);
		else if(!std::strcmp(current, "num_test_part"))
			i += std::fscanf(in, "%i", &num_test_part);
		else if(!std::strcmp(current, "density_samples"))
			i += std::fscanf(in, "%i", &density_samples);
		else if(!std::strcmp(current, "steps"))
			i += std::fscanf(in, "%i", &steps);
		else if(!std::strcmp(current, "n"))
			i += std::fscanf(in, "%i", &n);
		else if(!std::strcmp(current, "z"))
			i += std::fscanf(in, "%i", &z);
		else if(!std::strcmp(current, "substeps"))
			i += std::fscanf(in, "%i", &substeps);
		else if(!std::strcmp(current, "use_gpu"))
			i += std::fscanf(in, "%i", &use_gpu);
	}
	if(i != input_file_count) {
		std::fprintf(stderr, "INVALID INPUT FILE!\n"); std::exit(1);
	}
	T sigma_k = calc_sigma(T(k_fwhm)), sigma_r = calc_sigma(T(r_fwhm));
	T d_max = T(d_max_scale) * nuclear_radius<T>(z + n);
	
	set_parameters(param, z, n, num_test_part, steps, substeps, density_samples, static_cast<bool>(use_gpu), T(sample_position), T(sigma_k), T(sigma_r), T(t_f), T(t_exc), T(eta_exc));
	world = World(nx, T(d_max));
	skm = Skyrme(T(A), T(B), T(C), T(gamma));
	fermi_levels = Fermi(T(epsilon_p), T(epsilon_n));
	ws = WoodsSaxon(T(V0), T(0.8) * param.r_max, T(a));
}

template void set_parameters<double>(Parameters<double> &param, int z, int n, int part_per_nucleon, int steps, int substeps, int density_samples, bool use_gpu, double sample_position, double sigma_k, double sigma_r, double t_f, double t_exc, double eta_exc);
template void read_input_file<double>(std::FILE *in, Skyrme<double> &skm, World<double> &world, Fermi<double> &fermi_levels, Parameters<double> &param, WoodsSaxon<double> &ws);

template void set_parameters<float>(Parameters<float> &param, int z, int n, int part_per_nucleon, int steps, int substeps, int density_samples, bool use_gpu, float sample_position, float sigma_k, float sigma_r, float t_f, float t_exc, float eta_exc);
template void read_input_file<float>(std::FILE *in, Skyrme<float> &skm, World<float> &world, Fermi<float> &fermi_levels, Parameters<float> &param, WoodsSaxon<float> &ws);
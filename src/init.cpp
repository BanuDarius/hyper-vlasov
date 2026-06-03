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

#include <vector>
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
void output_vtk_header_start(FILE *out, World<T> world) {
	std::fprintf(out, "# vtk DataFile Version 3.0\n");
	std::fprintf(out, "Volumetric data\n");
	std::fprintf(out, "BINARY\n");
	std::fprintf(out, "DATASET STRUCTURED_POINTS\n");
	std::fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	std::fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	std::fprintf(out, "SPACING %lf %lf %lf\n", 2.0 * world.d_max[0] / world.n[0], 2.0 * world.d_max[1] / world.n[1], 2.0 * world.d_max[2] / world.n[2]);
	std::fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
}

void output_vtk_header_scalar_next(FILE *out, const char *name, int type) {
	char tag;
	if(type == is_proton) tag = 'p';
	else if(type == is_neutron) tag = 'n';
	else tag = 't';
	std::fprintf(out, "SCALARS %s_%c float 1\n", name, tag);
	std::fprintf(out, "LOOKUP_TABLE default\n");
}

void output_vtk_header_vector_next(FILE *out, const char *name, int type) {
	char tag;
	if(type == is_proton) tag = 'p';
	else if(type == is_neutron) tag = 'n';
	std::fprintf(out, "VECTORS %s_%c float\n", name, tag);
}

template <typename T>
void output_density_samples_positions(FILE *out, const Parameters<T> &param, const World<T> &world) {
	T d_max_z = world.d_max[2];
	int samples = param.density_samples;
	std::vector<float> positions(samples);
	#pragma omp parallel for
	for(int i = 0; i < samples; i++)
		positions[i] = d_max_z * 2.0 * i / samples - d_max_z;
	
	fwrite(positions.data(), sizeof(float), samples, out);
}

template <typename T>
void output_density_samples(FILE *out, const float *samples_ptr, const Parameters<T> &param) {
	fwrite(samples_ptr, sizeof(float), 2 * param.density_samples, out);
}

template <typename T>
void output_scalar_field(FILE *out, const ScalarField<T> &field, const World<T> &world, const char *name) {
	std::size_t nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	std::unique_ptr<uint32_t[]> vtk_scalar_p(new uint32_t[world_size]);
	std::unique_ptr<uint32_t[]> vtk_scalar_n(new uint32_t[world_size]);
	std::unique_ptr<uint32_t[]> vtk_scalar_t(new uint32_t[world_size]);
	#pragma omp parallel for collapse(3)
	for(std::size_t k = 0; k < nz; k++) {
		for(std::size_t j = 0; j < ny; j++) {
			for(std::size_t i = 0; i < nx; i++) {
				int idx = grid_idx(i, j, k, nx, ny, nz);
				int write_idx = (k * ny * nx) + (j * nx) + i;
				
				vtk_scalar_p[write_idx] = swap_endian(static_cast<float>(field.v[idx]));
				vtk_scalar_n[write_idx] = swap_endian(static_cast<float>(field.v[idx + world_size]));
				vtk_scalar_t[write_idx] = swap_endian(static_cast<float>(field.v[idx] + field.v[idx + world_size]));
			}
		}
	}
	output_vtk_header_scalar_next(out, name, is_proton);
	fwrite(vtk_scalar_p.get(), sizeof(uint32_t), world_size, out);
	
	output_vtk_header_scalar_next(out, name, is_neutron);
	fwrite(vtk_scalar_n.get(), sizeof(uint32_t), world_size, out);
	
	output_vtk_header_scalar_next(out, name, is_proton_or_neutron);
	fwrite(vtk_scalar_t.get(), sizeof(uint32_t), world_size, out);
}

template <typename T>
void output_vector_field(FILE *out, const VectorField<T> &field, const World<T> &world, const char *name) {
	size_t nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	std::unique_ptr<uint32_t[]> vtk_vector_p(new uint32_t[3 * world_size]);
	std::unique_ptr<uint32_t[]> vtk_vector_n(new uint32_t[3 * world_size]);
	#pragma omp parallel for collapse(3)
	for(size_t k = 0; k < nz; k++) {
		for(size_t j = 0; j < ny; j++) {
			for(size_t i = 0; i < nx; i++) {
				int idx = grid_idx(i, j, k, nx, ny, nz);
				int write_idx = (k * ny * nx) + (j * nx) + i;
				
				vtk_vector_p[3 * write_idx] = swap_endian(static_cast<float>(field.x[idx]));
				vtk_vector_p[3 * write_idx + 1] = swap_endian(static_cast<float>(field.y[idx]));
				vtk_vector_p[3 * write_idx + 2] = swap_endian(static_cast<float>(field.z[idx]));
				
				vtk_vector_n[3 * write_idx] = swap_endian(static_cast<float>(field.x[idx + world_size]));
				vtk_vector_n[3 * write_idx + 1] = swap_endian(static_cast<float>(field.y[idx + world_size]));
				vtk_vector_n[3 * write_idx + 2] = swap_endian(static_cast<float>(field.z[idx + world_size]));
			}
		}
	}
	output_vtk_header_vector_next(out, name, is_proton);
	fwrite(vtk_vector_p.get(), sizeof(uint32_t), 3 * world_size, out);
	
	output_vtk_header_vector_next(out, name, is_neutron);
	fwrite(vtk_vector_n.get(), sizeof(uint32_t), 3 * world_size, out);
}

template <typename T>
void read_input_file(FILE *in, Skyrme<T> &skm, World<T> &world, Fermi<T> &fermi_levels, Parameters<T> &param, WoodsSaxon<T> &ws) {
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
		std::fprintf(stderr, "Error: Invalid input file.\n"); std::exit(1);
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
template void output_vtk_header_start<double>(FILE *out, World<double> world);
template void output_density_samples_positions<double>(FILE *out, const Parameters<double> &param, const World<double> &world);
template void output_density_samples<double>(FILE *out, const float *samples, const Parameters<double> &param);
template void output_scalar_field<double>(FILE *out, const ScalarField<double> &field, const World<double> &world, const char *name);
template void output_vector_field<double>(FILE *out, const VectorField<double> &field, const World<double> &world, const char *name);
template void read_input_file<double>(FILE *in, Skyrme<double> &skm, World<double> &world, Fermi<double> &fermi_levels, Parameters<double> &param, WoodsSaxon<double> &ws);

template void set_parameters<float>(Parameters<float> &param, int z, int n, int part_per_nucleon, int steps, int substeps, int density_samples, bool use_gpu, float sample_position, float sigma_k, float sigma_r, float t_f, float t_exc, float eta_exc);
template void output_vtk_header_start<float>(FILE *out, World<float> world);
template void output_density_samples_positions<float>(FILE *out, const Parameters<float> &param, const World<float> &world);
template void output_density_samples<float>(FILE *out, const float *samples, const Parameters<float> &param);
template void output_scalar_field<float>(FILE *out, const ScalarField<float> &field, const World<float> &world, const char *name);
template void output_vector_field<float>(FILE *out, const VectorField<float> &field, const World<float> &world, const char *name);
template void read_input_file<float>(FILE *in, Skyrme<float> &skm, World<float> &world, Fermi<float> &fermi_levels, Parameters<float> &param, WoodsSaxon<float> &ws);
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

#include <cstdio>
#include <vector>

#include "vtk_output.hpp"

template <typename T>
void output_vtk_header_start(std::FILE *out, World<T> world) {
	std::fprintf(out, "# vtk DataFile Version 3.0\n");
	std::fprintf(out, "Volumetric data\n");
	std::fprintf(out, "BINARY\n");
	std::fprintf(out, "DATASET STRUCTURED_POINTS\n");
	std::fprintf(out, "DIMENSIONS %d %d %d\n", world.n[0], world.n[1], world.n[2]);
	std::fprintf(out, "ORIGIN %lf %lf %lf\n", -world.d_max[0], -world.d_max[1], -world.d_max[2]);
	std::fprintf(out, "SPACING %lf %lf %lf\n", T(2.0) * world.d_max[0] / world.n[0], T(2.0) * world.d_max[1] / world.n[1], T(2.0) * world.d_max[2] / world.n[2]);
	std::fprintf(out, "POINT_DATA %d\n", world.n[0] * world.n[1] * world.n[2]);
}

void output_vtk_header_scalar_next(std::FILE *out, const char *name, int type) {
	char tag;
	if(type == is_proton) tag = 'p';
	else if(type == is_neutron) tag = 'n';
	else tag = 't';
	std::fprintf(out, "SCALARS %s_%c float 1\n", name, tag);
	std::fprintf(out, "LOOKUP_TABLE default\n");
}

void output_vtk_header_vector_next(std::FILE *out, const char *name, int type) {
	char tag;
	if(type == is_proton) tag = 'p';
	else if(type == is_neutron) tag = 'n';
	std::fprintf(out, "VECTORS %s_%c float\n", name, tag);
}

template <typename T>
void output_density_samples_positions(std::FILE *out, const Parameters<T> &param, const World<T> &world) {
	T d_max_z = world.d_max[2];
	int samples = param.density_samples;
	std::vector<float> positions(samples);
	#pragma omp parallel for
	for(int i = 0; i < samples; i++)
		positions[i] = d_max_z * 2.0 * i / samples - d_max_z;
	
	std::fwrite(positions.data(), sizeof(float), samples, out);
}

template <typename T>
void output_density_samples(std::FILE *out, const float *samples_ptr, const Parameters<T> &param) {
	std::fwrite(samples_ptr, sizeof(float), 2 * param.density_samples, out);
}

template <typename T>
void output_scalar_field(std::FILE *out, const ScalarField<T> &field, const World<T> &world, const char *name) {
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
	std::fwrite(vtk_scalar_p.get(), sizeof(uint32_t), world_size, out);
	
	output_vtk_header_scalar_next(out, name, is_neutron);
	std::fwrite(vtk_scalar_n.get(), sizeof(uint32_t), world_size, out);
	
	output_vtk_header_scalar_next(out, name, is_proton_or_neutron);
	std::fwrite(vtk_scalar_t.get(), sizeof(uint32_t), world_size, out);
}

template <typename T>
void output_vector_field(std::FILE *out, const VectorField<T> &field, const World<T> &world, const char *name) {
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
	std::fwrite(vtk_vector_p.get(), sizeof(uint32_t), 3 * world_size, out);
	
	output_vtk_header_vector_next(out, name, is_neutron);
	std::fwrite(vtk_vector_n.get(), sizeof(uint32_t), 3 * world_size, out);
}

template void output_vtk_header_start<double>(std::FILE *out, World<double> world);
template void output_density_samples_positions<double>(std::FILE *out, const Parameters<double> &param, const World<double> &world);
template void output_density_samples<double>(std::FILE *out, const float *samples_ptr, const Parameters<double> &param);
template void output_scalar_field<double>(std::FILE *out, const ScalarField<double> &field, const World<double> &world, const char *name);
template void output_vector_field<double>(std::FILE *out, const VectorField<double> &field, const World<double> &world, const char *name);

template void output_vtk_header_start<float>(std::FILE *out, World<float> world);
template void output_density_samples_positions<float>(std::FILE *out, const Parameters<float> &param, const World<float> &world);
template void output_density_samples<float>(std::FILE *out, const float *samples_ptr, const Parameters<float> &param);
template void output_scalar_field<float>(std::FILE *out, const ScalarField<float> &field, const World<float> &world, const char *name);
template void output_vector_field<float>(std::FILE *out, const VectorField<float> &field, const World<float> &world, const char *name);
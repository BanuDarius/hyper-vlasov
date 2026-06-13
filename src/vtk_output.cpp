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
void output_vtk_header_start(std::ofstream &output_file, World<T> world) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	T r_max_x = world.d_max[0], r_max_y = world.d_max[1], r_max_z = world.d_max[2];
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "Volumetric data\n";
	output_file << "BINARY\n";
	output_file << "DATASET STRUCTURED_POINTS\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " "  << nz << "\n";
	output_file << "ORIGIN " << -r_max_x << " " << -r_max_y << " " <<  -r_max_z << "\n";
	output_file << "SPACING " << T(2.0) * r_max_x / nx << " " << T(2.0) * r_max_y / ny << " " << T(2.0) * r_max_z / nz << "\n";
	output_file << "POINT_DATA " << nx * ny * nz << "\n";
}

void output_vtk_header_scalar_next(std::ofstream &output_file, const char *name, int type) {
	char tag;
	if(type == is_proton) tag = 'p';
	else if(type == is_neutron) tag = 'n';
	else tag = 't';
	output_file << "SCALARS " << name << "_" << tag << " float 1\n";
	output_file << "LOOKUP_TABLE default\n";
}

void output_vtk_header_vector_next(std::ofstream &output_file, const char *name, int type) {
	char tag;
	if(type == is_proton) tag = 'p';
	else if(type == is_neutron) tag = 'n';
	output_file << "VECTORS " << name << "_" << tag << " float\n";
}

template <typename T>
void output_density_samples_positions(std::ofstream &output_file, const Parameters<T> &param, const World<T> &world) {
	T d_max_z = world.d_max[2];
	int samples = param.density_samples;
	std::vector<float> positions(samples);
	#pragma omp parallel for
	for(int i = 0; i < samples; i++)
		positions[i] = d_max_z * 2.0 * i / samples - d_max_z;
	
	output_file.write(reinterpret_cast<const char*>(positions.data()), samples * sizeof(float));
}

template <typename T>
void output_density_samples(std::ofstream &output_file, const float *samples_ptr, const Parameters<T> &param) {
	output_file.write(reinterpret_cast<const char*>(samples_ptr), 2 * param.density_samples * sizeof(float));
}

template <typename T>
void output_scalar_field(std::ofstream &output_file, const ScalarField<T> &field, const World<T> &world, const char *name) {
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
	output_vtk_header_scalar_next(output_file, name, is_proton);
	output_file.write(reinterpret_cast<const char*>(vtk_scalar_p.get()), world_size * sizeof(uint32_t));
	
	output_vtk_header_scalar_next(output_file, name, is_neutron);
	output_file.write(reinterpret_cast<const char*>(vtk_scalar_n.get()), world_size * sizeof(uint32_t));
	
	output_vtk_header_scalar_next(output_file, name, is_proton_or_neutron);
	output_file.write(reinterpret_cast<const char*>(vtk_scalar_t.get()), world_size * sizeof(uint32_t));
}

template <typename T>
void output_vector_field(std::ofstream &output_file, const VectorField<T> &field, const World<T> &world, const char *name) {
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
	output_vtk_header_vector_next(output_file, name, is_proton);
	output_file.write(reinterpret_cast<const char*>(vtk_vector_p.get()), 3 * world_size * sizeof(uint32_t));
	
	output_vtk_header_vector_next(output_file, name, is_neutron);
	output_file.write(reinterpret_cast<const char*>(vtk_vector_n.get()), 3 * world_size * sizeof(uint32_t));
}

template void output_vtk_header_start<double>(std::ofstream &output_file, World<double> world);
template void output_density_samples_positions<double>(std::ofstream &output_file, const Parameters<double> &param, const World<double> &world);
template void output_density_samples<double>(std::ofstream &output_file, const float *samples_ptr, const Parameters<double> &param);
template void output_scalar_field<double>(std::ofstream &output_file, const ScalarField<double> &field, const World<double> &world, const char *name);
template void output_vector_field<double>(std::ofstream &output_file, const VectorField<double> &field, const World<double> &world, const char *name);

template void output_vtk_header_start<float>(std::ofstream &output_file, World<float> world);
template void output_density_samples_positions<float>(std::ofstream &output_file, const Parameters<float> &param, const World<float> &world);
template void output_density_samples<float>(std::ofstream &output_file, const float *samples_ptr, const Parameters<float> &param);
template void output_scalar_field<float>(std::ofstream &output_file, const ScalarField<float> &field, const World<float> &world, const char *name);
template void output_vector_field<float>(std::ofstream &output_file, const VectorField<float> &field, const World<float> &world, const char *name);
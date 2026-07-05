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

#ifndef PARTICLE_IN_CELL_H
#define PARTICLE_IN_CELL_H

#include <concepts>

#include "sim_structs.hpp"

template <std::floating_point T>
inline T scatter_scalar_field_cic(const ScalarField<T> &field, std::array<T, 3> r_vec, const World<T> &world, int type) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	
	T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
	T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
	T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
	
	if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz)
		return T(0.0);
	
	int x0 = static_cast<int>(cx), y0 = static_cast<int>(cy), z0 = static_cast<int>(cz);
	int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
	T d_x = cx - x0, d_y = cy - y0, d_z = cz - z0;
	T t_x = T(1.0) - d_x, t_y = T(1.0) - d_y, t_z = T(1.0) - d_z;
	
	int offset = (type == is_proton) ? 0 : world_size;
	int idx = grid_idx(x0, y0, z0, nx, ny, nz) + offset;
	T v = t_x * t_y * t_z * field.v[idx];
	
	if(x1 < nx) {
		idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
		v += d_x * t_y * t_z * field.v[idx];
	}
	if(y1 < ny) {
		idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
		v += t_x * d_y * t_z * field.v[idx];
	}
	if(x1 < nx && y1 < ny) {
		idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
		v += d_x * d_y * t_z * field.v[idx];
	}
	if(z1 < nz) {
		idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
		v += t_x * t_y * d_z * field.v[idx];
	}
	if(x1 < nx && z1 < nz) {
		idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
		v += d_x * t_y * d_z * field.v[idx];
	}
	if(y1 < ny && z1 < nz) {
		idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
		v += t_x * d_y * d_z * field.v[idx];
	}
	if(x1 < nx && y1 < ny && z1 < nz) {
		idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
		v += d_x * d_y * d_z * field.v[idx];
	}
	return v;
}

template <std::floating_point T> void distribute_volumetric_particles_cic(ScalarField<T> &density, const TestParticles<T> &part, const World<T> &world);
template <std::floating_point T> void distribute_forces_to_particles_cic(TestParticles<T> &part, const VectorField<T> &forces, const World<T> &world);
template <std::floating_point T> void distribute_volumetric_momenta_cic(VectorField<T> &current_density, const TestParticles<T> &part, const World<T> &world);
template <std::floating_point T> void compute_density_samples_cic(float *density_samples, const ScalarField<T> &density, const Parameters<T> &param, const World<T> &world);

#endif
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

#include <omp.h>
#include <array>

template <typename T>
void distribute_volumetric_particles_cic(ScalarField<T> &density, const TestParticles<T> &part, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part.protons + part.neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		density.v[i] = T(0.0);
	
	T *density_ptr = density.v.data();
	#pragma omp parallel for reduction(+:density_ptr[0 : 2 * world_size])
	for(int i = 0; i < total; i++) {
		std::array<T, 3> r_vec;
		particle_pos_to_vector(r_vec, part, i);
		
		T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
		T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
		T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
		
		if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz)
			continue;
		
		int x0 = static_cast<int>(cx), y0 = static_cast<int>(cy), z0 = static_cast<int>(cz);
		int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
		T d_x = cx - x0, d_y = cy - y0, d_z = cz - z0;
		T t_x = T(1.0) - d_x, t_y = T(1.0) - d_y, t_z = T(1.0) - d_z;
		
		int offset = (i < part.protons) ? 0 : world_size;
		int idx = grid_idx(x0, y0, z0, nx, ny, nz) + offset;
		density_ptr[idx] += t_x * t_y * t_z;
		
		if (x1 < nx) {
			idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
			density_ptr[idx] += d_x * t_y * t_z;
		}
		if (y1 < ny) {
			idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
			density_ptr[idx] += t_x * d_y * t_z;
		}
		if (x1 < nx && y1 < ny) {
			idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
			density_ptr[idx] += d_x * d_y * t_z;
		}
		if (z1 < nz) {
			idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
			density_ptr[idx] += t_x * t_y * d_z;
		}
		if (x1 < nx && z1 < nz) {
			idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
			density_ptr[idx] += d_x * t_y * d_z;
		}
		if (y1 < ny && z1 < nz) {
			idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
			density_ptr[idx] += t_x * d_y * d_z;
		}
		if (x1 < nx && y1 < ny && z1 < nz) {
			idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
			density_ptr[idx] += d_x * d_y * d_z;
		}
	}
}

template <typename T>
void distribute_forces_to_particles_cic(TestParticles<T> &part, const VectorField<T> &forces, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part.protons + part.neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		std::array<T, 3> r_vec;
		particle_pos_to_vector(r_vec, part, i);
		
		T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
		T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
		T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
		
		if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz) {
			part.fx[i] = T(0.0); part.fy[i] = T(0.0); part.fz[i] = T(0.0);
			continue;
		}
		
		int x0 = static_cast<int>(cx), y0 = static_cast<int>(cy), z0 = static_cast<int>(cz);
		int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
		T d_x = cx - x0, d_y = cy - y0, d_z = cz - z0;
		T t_x = T(1.0) - d_x, t_y = T(1.0) - d_y, t_z = T(1.0) - d_z;
		
		int offset = (i < part.protons) ? 0 : world_size;
		int idx = grid_idx(x0, y0, z0, nx, ny, nz) + offset;
		
		T w = t_x * t_y * t_z;
		T fx = w * forces.x[idx], fy = w * forces.y[idx], fz = w * forces.z[idx];
		
		if(x1 < nx) {
			idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
			w = d_x * t_y * t_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(y1 < ny) {
			idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
			w = t_x * d_y * t_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && y1 < ny) {
			idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
			w = d_x * d_y * t_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(z1 < nz) {
			idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
			w = t_x * t_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && z1 < nz) {
			idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
			w = d_x * t_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(y1 < ny && z1 < nz) {
			idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
			w = t_x * d_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && y1 < ny && z1 < nz) {
			idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
			w = d_x * d_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		part.fx[i] = fx;
		part.fy[i] = fy;
		part.fz[i] = fz;
	}
}

template <typename T>
void distribute_volumetric_momenta_cic(ScalarField<T> &current_density, const TestParticles<T> &part, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part.protons + part.neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		current_density.x[i] = T(0.0); current_density.y[i] = T(0.0); current_density.z[i] = T(0.0);
	}
	
	T *current_density_ptr_x = current_density.x.data(), *current_density_ptr_y = current_density.y.data(), *current_density_ptr_z = current_density.z.data();
	#pragma omp parallel for reduction(+:current_density_ptr_x[0 : 2 * world_size], current_density_ptr_y[0 : 2 * world_size], current_density_ptr_z[0 : 2 * world_size])
	for(int i = 0; i < total; i++) {
		std::array<T, 3> r_vec, k_vec;
		particle_pos_to_vector(r_vec, part, i);
		particle_vel_to_vector(k_vec, part, i);
		
		T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
		T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
		T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
		
		if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz)
			continue;
		
		int x0 = static_cast<int>(cx), y0 = static_cast<int>(cy), z0 = static_cast<int>(cz);
		int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
		T d_x = cx - x0, d_y = cy - y0, d_z = cz - z0;
		T t_x = T(1.0) - d_x, t_y = T(1.0) - d_y, t_z = T(1.0) - d_z;
		
		int offset = (i < part.protons) ? 0 : world_size;
		int idx = grid_idx(x0, y0, z0, nx, ny, nz) + offset;
		
		T w = t_x * t_y * t_z;
		current_density_ptr_x[idx] += w * part.kx[i];
		current_density_ptr_y[idx] += w * part.ky[i];
		current_density_ptr_z[idx] += w * part.kz[i];
		
		if (x1 < nx) {
			w = d_x * t_y * t_z;
			idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
			current_density_ptr_x[idx] += w * part.kx[i];
			current_density_ptr_y[idx] += w * part.ky[i];
			current_density_ptr_z[idx] += w * part.kz[i];
		}
		if (y1 < ny) {
			w = t_x * d_y * t_z;
			idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
			current_density_ptr_x[idx] += w * part.kx[i];
			current_density_ptr_y[idx] += w * part.ky[i];
			current_density_ptr_z[idx] += w * part.kz[i];
		}
		if (x1 < nx && y1 < ny) {
			w = d_x * d_y * t_z;
			idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
			current_density_ptr_x[idx] += w * part.kx[i];
			current_density_ptr_y[idx] += w * part.ky[i];
			current_density_ptr_z[idx] += w * part.kz[i];
		}
		if (z1 < nz) {
			w = t_x * t_y * d_z;
			idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
			current_density_ptr_x[idx] += w * part.kx[i];
			current_density_ptr_y[idx] += w * part.ky[i];
			current_density_ptr_z[idx] += w * part.kz[i];
		}
		if (x1 < nx && z1 < nz) {
			w = d_x * t_y * d_z;
			idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
			current_density_ptr_x[idx] += w * part.kx[i];
			current_density_ptr_y[idx] += w * part.ky[i];
			current_density_ptr_z[idx] += w * part.kz[i];
		}
		if (y1 < ny && z1 < nz) {
			w = t_x * d_y * d_z;
			idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
			current_density_ptr_x[idx] += w * part.kx[i];
			current_density_ptr_y[idx] += w * part.ky[i];
			current_density_ptr_z[idx] += w * part.kz[i];
		}
		if (x1 < nx && y1 < ny && z1 < nz) {
			w = d_x * d_y * d_z;
			idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
			current_density_ptr_x[idx] += w * part.kx[i];
			current_density_ptr_y[idx] += w * part.ky[i];
			current_density_ptr_z[idx] += w * part.kz[i];
		}
	}
}

template <typename T>
void compute_density_samples_cic(std::vector<float> &density_samples, const ScalarField<T> &density, const Parameters<T> &param, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz;
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * param.density_samples; i++) {
		int i_new = (i < param.density_samples) ? i : i - param.density_samples;
		T z = 2.0 * d_max_z * i_new / param.density_samples - d_max_z;
		T y = d_max_y * param.sample_position;
		std::array<T, 3> r_vec = { T(0.0), y, z };
		
		T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
		T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
		T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
		
		if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz) {
			density_samples[i] = 0.0f;
			continue;
		}
		
		int x0 = static_cast<int>(cx), y0 = static_cast<int>(cy), z0 = static_cast<int>(cz);
		int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
		T d_x = cx - x0, d_y = cy - y0, d_z = cz - z0;
		T t_x = T(1.0) - d_x, t_y = T(1.0) - d_y, t_z = T(1.0) - d_z;
		
		int offset = (i < param.density_samples) ? 0 : world_size;
		int idx = grid_idx(x0, y0, z0, nx, ny, nz) + offset;
		T rho = t_x * t_y * t_z * density.v[idx];
		
		if(x1 < nx) {
			idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
			rho += d_x * t_y * t_z * density.v[idx];
		}
		if(y1 < ny) {
			idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
			rho += t_x * d_y * t_z * density.v[idx];
		}
		if(x1 < nx && y1 < ny) {
			idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
			rho += d_x * d_y * t_z * density.v[idx];
		}
		if(z1 < nz) {
			idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
			rho += t_x * t_y * d_z * density.v[idx];
		}
		if(x1 < nx && z1 < nz) {
			idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
			rho += d_x * t_y * d_z * density.v[idx];
		}
		if(y1 < ny && z1 < nz) {
			idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
			rho += t_x * d_y * d_z * density.v[idx];
		}
		if(x1 < nx && y1 < ny && z1 < nz) {
			idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
			rho += d_x * d_y * d_z * density.v[idx];
		}
		density_samples[i] = static_cast<float>(rho);
	}
}

template <typename T>
void distribute_volumetric_momenta_ngp(ScalarField<T> &current_density, const TestParticles<T> &part, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part.protons + part.neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		current_density.x[i] = T(0.0); current_density.y[i] = T(0.0); current_density.z[i] = T(0.0);
	}
	
	T *current_density_ptr_x = current_density.x.data(), *current_density_ptr_y = current_density.y.data(), *current_density_ptr_z = current_density.z.data();
	#pragma omp parallel for reduction(+:current_density_ptr_x[0 : 2 * world_size], current_density_ptr_y[0 : 2 * world_size], current_density_ptr_z[0 : 2 * world_size])
	for(int i = 0; i < total; i++) {
		std::array<T, 3> r_vec, k_vec;
		particle_pos_to_vector(r_vec, part, i);
		particle_vel_to_vector(k_vec, part, i);
		
		T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
		T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
		T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
		
		if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz)
			continue;
		
		int x = static_cast<int>(cx), y = static_cast<int>(cy), z = static_cast<int>(cz);
		
		int offset = (i < part.protons) ? 0 : world_size;
		int idx = grid_idx(x, y, z, nx, ny, nz) + offset;
		
		current_density_ptr_x[idx] += part.kx[i];
		current_density_ptr_y[idx] += part.ky[i];
		current_density_ptr_z[idx] += part.kz[i];
	}
}

#endif
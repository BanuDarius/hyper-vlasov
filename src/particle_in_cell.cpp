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

#include "particle_in_cell.hpp"
#include "tools.hpp"

template <typename T>
void distribute_volumetric_particles_cic(ScalarField<T> &density, const TestParticles<T> &part, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part.protons + part.neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		density.v[i] = T(0.0);
	
	#pragma omp parallel for
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
		#pragma omp atomic
		density.v[idx] += t_x * t_y * t_z;
		
		if (x1 < nx) {
			idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
			#pragma omp atomic
			density.v[idx] += d_x * t_y * t_z;
		}
		if (y1 < ny) {
			idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
			#pragma omp atomic
			density.v[idx] += t_x * d_y * t_z;
		}
		if (x1 < nx && y1 < ny) {
			idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
			#pragma omp atomic
			density.v[idx] += d_x * d_y * t_z;
		}
		if (z1 < nz) {
			idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			density.v[idx] += t_x * t_y * d_z;
		}
		if (x1 < nx && z1 < nz) {
			idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			density.v[idx] += d_x * t_y * d_z;
		}
		if (y1 < ny && z1 < nz) {
			idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			density.v[idx] += t_x * d_y * d_z;
		}
		if (x1 < nx && y1 < ny && z1 < nz) {
			idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			density.v[idx] += d_x * d_y * d_z;
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
		T w = t_x * t_y * t_z;
		
		int offset = (i < part.protons) ? 0 : world_size;
		int idx = grid_idx(x0, y0, z0, nx, ny, nz) + offset;
		
		T fx = w * forces.x[idx], fy = w * forces.y[idx], fz = w * forces.z[idx];
		
		if(x1 < nx) {
			w = d_x * t_y * t_z;
			idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(y1 < ny) {
			w = t_x * d_y * t_z;
			idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && y1 < ny) {
			w = d_x * d_y * t_z;
			idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(z1 < nz) {
			w = t_x * t_y * d_z;
			idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && z1 < nz) {
			w = d_x * t_y * d_z;
			idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(y1 < ny && z1 < nz) {
			w = t_x * d_y * d_z;
			idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && y1 < ny && z1 < nz) {
			w = d_x * d_y * d_z;
			idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		part.fx[i] = fx;
		part.fy[i] = fy;
		part.fz[i] = fz;
	}
}

template <typename T>
void distribute_volumetric_momenta_cic(VectorField<T> &current_density, const TestParticles<T> &part, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part.protons + part.neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++) {
		current_density.x[i] = T(0.0); current_density.y[i] = T(0.0); current_density.z[i] = T(0.0);
	}
	
	#pragma omp parallel for
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
		T w = t_x * t_y * t_z;
		
		int offset = (i < part.protons) ? 0 : world_size;
		int idx = grid_idx(x0, y0, z0, nx, ny, nz) + offset;
		
		#pragma omp atomic
		current_density.x[idx] += w * k_vec[0];
		#pragma omp atomic
		current_density.y[idx] += w * k_vec[1];
		#pragma omp atomic
		current_density.z[idx] += w * k_vec[2];
		
		if (x1 < nx) {
			w = d_x * t_y * t_z;
			idx = grid_idx(x1, y0, z0, nx, ny, nz) + offset;
			#pragma omp atomic
			current_density.x[idx] += w * k_vec[0];
			#pragma omp atomic
			current_density.y[idx] += w * k_vec[1];
			#pragma omp atomic
			current_density.z[idx] += w * k_vec[2];
		}
		if (y1 < ny) {
			w = t_x * d_y * t_z;
			idx = grid_idx(x0, y1, z0, nx, ny, nz) + offset;
			#pragma omp atomic
			current_density.x[idx] += w * k_vec[0];
			#pragma omp atomic
			current_density.y[idx] += w * k_vec[1];
			#pragma omp atomic
			current_density.z[idx] += w * k_vec[2];
		}
		if (x1 < nx && y1 < ny) {
			w = d_x * d_y * t_z;
			idx = grid_idx(x1, y1, z0, nx, ny, nz) + offset;
			#pragma omp atomic
			current_density.x[idx] += w * k_vec[0];
			#pragma omp atomic
			current_density.y[idx] += w * k_vec[1];
			#pragma omp atomic
			current_density.z[idx] += w * k_vec[2];
		}
		if (z1 < nz) {
			w = t_x * t_y * d_z;
			idx = grid_idx(x0, y0, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			current_density.x[idx] += w * k_vec[0];
			#pragma omp atomic
			current_density.y[idx] += w * k_vec[1];
			#pragma omp atomic
			current_density.z[idx] += w * k_vec[2];
		}
		if (x1 < nx && z1 < nz) {
			w = d_x * t_y * d_z;
			idx = grid_idx(x1, y0, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			current_density.x[idx] += w * k_vec[0];
			#pragma omp atomic
			current_density.y[idx] += w * k_vec[1];
			#pragma omp atomic
			current_density.z[idx] += w * k_vec[2];
		}
		if (y1 < ny && z1 < nz) {
			w = t_x * d_y * d_z;
			idx = grid_idx(x0, y1, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			current_density.x[idx] += w * k_vec[0];
			#pragma omp atomic
			current_density.y[idx] += w * k_vec[1];
			#pragma omp atomic
			current_density.z[idx] += w * k_vec[2];
		}
		if (x1 < nx && y1 < ny && z1 < nz) {
			w = d_x * d_y * d_z;
			idx = grid_idx(x1, y1, z1, nx, ny, nz) + offset;
			#pragma omp atomic
			current_density.x[idx] += w * k_vec[0];
			#pragma omp atomic
			current_density.y[idx] += w * k_vec[1];
			#pragma omp atomic
			current_density.z[idx] += w * k_vec[2];
		}
	}
}

template <typename T>
void compute_density_samples_cic(float *density_samples, const ScalarField<T> &density, const Parameters<T> &param, const World<T> &world) {
	T d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	
	#pragma omp parallel for
	for(int i = 0; i < 2 * param.density_samples; i++) {
		int i_new = (i < param.density_samples) ? i : i - param.density_samples;
		T z = 2.0 * d_max_z * i_new / param.density_samples - d_max_z;
		T y = d_max_y * param.sample_position;
		std::array<T, 3> r_vec = { T(0.0), y, z };
		
		int type = (i < param.density_samples) ? is_proton : is_neutron;
		density_samples[i] = scatter_scalar_field_cic(density, r_vec, world, type);
	}
}

template void distribute_volumetric_particles_cic<double>(ScalarField<double> &density, const TestParticles<double> &part, const World<double> &world);
template void distribute_forces_to_particles_cic<double>(TestParticles<double> &part, const VectorField<double> &forces, const World<double> &world);
template void distribute_volumetric_momenta_cic<double>(VectorField<double> &current_density, const TestParticles<double> &part, const World<double> &world);
template void compute_density_samples_cic<double>(float *density_samples, const ScalarField<double> &density, const Parameters<double> &param, const World<double> &world);

template void distribute_volumetric_particles_cic<float>(ScalarField<float> &density, const TestParticles<float> &part, const World<float> &world);
template void distribute_forces_to_particles_cic<float>(TestParticles<float> &part, const VectorField<float> &forces, const World<float> &world);
template void distribute_volumetric_momenta_cic<float>(VectorField<float> &current_density, const TestParticles<float> &part, const World<float> &world);
template void compute_density_samples_cic<float>(float *density_samples, const ScalarField<float> &density, const Parameters<float> &param, const World<float> &world);
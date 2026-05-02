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

#ifndef TOOLS_H
#define TOOLS_H

#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "math_functions.hpp"
#include "sim_structs.hpp"
#include "physics_formulas.hpp"

template <typename T>
void distribute_volumetric_particles_cic(ScalarField<T> *density, TestParticles<T> *part, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part->protons + part->neutrons;
	
	#pragma omp parallel for simd
	for(int i = 0; i < 2 * world_size; i++)
		density->v[i] = 0.0;
	
	T *density_ptr = density->v;
	#pragma omp parallel for reduction(+: density_ptr[0 : 2 * world_size])
	for(int i = 0; i < total; i++) {
		T r_vec[3];
		copy_particle_pos_to_vector(r_vec, *part, i);
		
		T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
		T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
		T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
		
		if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz)
			continue;
		
		int x0 = (int)cx, y0 = (int)cy, z0 = (int)cz;
		int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
		T d_x = cx - x0, d_y = cy - y0, d_z = cz - z0;
		T t_x = T(1.0) - d_x, t_y = T(1.0) - d_y, t_z = T(1.0) - d_z;
		
		int offset = (i < part->protons) ? 0 : world_size;
		
		int idx000 = x0 * (ny * nz) + y0 * nz + z0 + offset;
		density_ptr[idx000] += t_x * t_y * t_z;
		
		if (x1 < nx) {
			int idx100 = x1 * (ny * nz) + y0 * nz + z0 + offset;
			density_ptr[idx100] += d_x * t_y * t_z;
		}
		if (y1 < ny) {
			int idx010 = x0 * (ny * nz) + y1 * nz + z0 + offset;
			density_ptr[idx010] += t_x * d_y * t_z;
		}
		if (x1 < nx && y1 < ny) {
			int idx110 = x1 * (ny * nz) + y1 * nz + z0 + offset;
			density_ptr[idx110] += d_x * d_y * t_z;
		}
		if (z1 < nz) {
			int idx001 = x0 * (ny * nz) + y0 * nz + z1 + offset;
			density_ptr[idx001] += t_x * t_y * d_z;
		}
		if (x1 < nx && z1 < nz) {
			int idx101 = x1 * (ny * nz) + y0 * nz + z1 + offset;
			density_ptr[idx101] += d_x * t_y * d_z;
		}
		if (y1 < ny && z1 < nz) {
			int idx011 = x0 * (ny * nz) + y1 * nz + z1 + offset;
			density_ptr[idx011] += t_x * d_y * d_z;
		}
		if (x1 < nx && y1 < ny && z1 < nz) {
			int idx111 = x1 * (ny * nz) + y1 * nz + z1 + offset;
			density_ptr[idx111] += d_x * d_y * d_z;
		}
	}
}

template <typename T>
void distribute_forces_to_particles_cic(TestParticles<T> *part, const VectorField<T> &forces, const World<T> &world) {
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int nx = world.n[0], ny = world.n[1], nz = world.n[2], world_size = nx * ny * nz, total = part->protons + part->neutrons;
	
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		T r_vec[3];
		copy_particle_pos_to_vector(r_vec, *part, i);
		
		T cx = (nx / T(2.0)) * (r_vec[0] / d_max_x + T(1.0));
		T cy = (ny / T(2.0)) * (r_vec[1] / d_max_y + T(1.0));
		T cz = (nz / T(2.0)) * (r_vec[2] / d_max_z + T(1.0));
		
		if(cx < T(0.0) || cy < T(0.0) || cz < T(0.0) || cx >= nx || cy >= ny || cz >= nz) {
			part->fx[i] = T(0.0); part->fy[i] = T(0.0); part->fz[i] = T(0.0);
			continue;
		}
		
		int x0 = (int)cx, y0 = (int)cy, z0 = (int)cz;
		int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
		T d_x = cx - x0, d_y = cy - y0, d_z = cz - z0;
		T t_x = T(1.0) - d_x, t_y = T(1.0) - d_y, t_z = T(1.0) - d_z;
		
		int offset = (i < part->protons) ? 0 : world_size;
		T fx = T(0.0), fy = T(0.0), fz = T(0.0);
		
		int idx = IDX(x0, y0, z0, nx, ny, nz) + offset;
		T w = t_x * t_y * t_z;
		fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		
		if(x1 < nx) {
			idx = IDX(x1, y0, z0, nx, ny, nz) + offset;
			w = d_x * t_y * t_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(y1 < ny) {
			idx = IDX(x0, y1, z0, nx, ny, nz) + offset;
			w = t_x * d_y * t_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && y1 < ny) {
			idx = IDX(x1, y1, z0, nx, ny, nz) + offset;
			w = d_x * d_y * t_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(z1 < nz) {
			idx = IDX(x0, y0, z1, nx, ny, nz) + offset;
			w = t_x * t_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && z1 < nz) {
			idx = IDX(x1, y0, z1, nx, ny, nz) + offset;
			w = d_x * t_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(y1 < ny && z1 < nz) {
			idx = IDX(x0, y1, z1, nx, ny, nz) + offset;
			w = t_x * d_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		if(x1 < nx && y1 < ny && z1 < nz) {
			idx = IDX(x1, y1, z1, nx, ny, nz) + offset;
			w = d_x * d_y * d_z;
			fx += w * forces.x[idx]; fy += w * forces.y[idx]; fz += w * forces.z[idx];
		}
		part->fx[i] = fx;
		part->fy[i] = fy;
		part->fz[i] = fz;
	}
}

template <typename T>
void compute_particle_densities(TestParticles<T> *part, const Parameters<T> &param) {
	int part_per_nucleon = param.part_per_nucleon, protons = part->protons, neutrons = part->neutrons, total = protons + neutrons;
	T sigma_r = param.sigma_r, exp_term = T(1.0) / (T(2.0) * sigma_r * sigma_r);
	T cutoff_squared = T(16.0) * sigma_r * sigma_r;
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		T r_i[3], r_j[3], diff[3];
		T fact, dist_squared, density_p = T(0.0), density_n = T(0.0);
		copy_particle_pos_to_vector(r_i, *part, i);
		
		for(int j = 0; j < total; j++) {
			copy_particle_pos_to_vector(r_j, *part, j);
			
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			if(dist_squared > cutoff_squared)
				continue;
			fact = exp(-dist_squared * exp_term);
			if(j < protons)
				density_p += fact;
			else
				density_n += fact;
		}
		part->density_p[i] = density_p;
		part->density_n[i] = density_n;
	}
	T term = (T(1.0) / part_per_nucleon) * (T(1.0) / std::pow(T(2.0) * pi<T> * sigma_r * sigma_r, T(1.5)));
	#pragma omp parallel for simd
	for(int i = 0; i < total; i++) {
		part->density_p[i] *= term;
		part->density_n[i] *= term;
	}
}

template <typename T>
T compute_energy(TestParticles<T> *part, const WoodsSaxon<T> *ws, T sigma_k, int z, int i) {
	T r_vec[3], k_vec[3];
	copy_particle_pos_to_vector(r_vec, *part, i);
	copy_particle_vel_to_vector(k_vec, *part, i);
	
	T r = magnitude(r_vec);
	T k = magnitude(k_vec);
	
	T energy = T(0.0);
	WoodsSaxon<T> ws_c = (i < part->protons) ? ws[0] : ws[1];
	energy += woods_saxon_potential(ws_c, r);
	energy += (k * k) * kinetic_energy<T>();
	
	if(i < part->protons)
		energy += coulomb_potential<T>(ws_c, z, r);
	
	energy += fluctuation_energy(sigma_k);
	return energy;
}

template <typename T>
void compute_particle_energies(TestParticles<T> *part, const WoodsSaxon<T> *ws, const Parameters<T> &param) {
	T sigma_k = param.sigma_k, z = param.z;
	#pragma omp parallel for
	for(int i = 0; i < part->protons + part->neutrons; i++)
		part->energy[i] = compute_energy(part, ws, sigma_k, z, i);
}

template <typename T>
void generate_random_particles(TestParticles<T> *part, T r_max) {
	int total = part->protons + part->neutrons, i = 0;
	while(i < total) {
		T r_new[3];
		random_vec(r_new, r_max);
		if(dot(r_new, r_new) < r_max * r_max) {
			copy_vector_to_particle_pos(part, r_new, i);
			i++;
		}
	}
	i = 0;
	while(i < total) {
		T k_new[3];
		random_vec(k_new, k_max<T>);
		if(dot(k_new, k_new) < k_max<T> * k_max<T>) {
			copy_vector_to_particle_vel(part, k_new, i);
			i++;
		}
	}
}

template <typename T>
void generate_checking_particles(TestParticles<T> *part, const WoodsSaxon<T> *ws, const Parameters<T> &param, const Fermi<T> *fermi_levels) {
	T r_max = param.r_max, sigma_k = param.sigma_k, z = param.z, epsilon;
	int total = part->protons + part->neutrons, i = 0;
	while(i < total) {
		if(i < part->protons)
			epsilon = fermi_levels->epsilon_p;
		else
			epsilon = fermi_levels->epsilon_n;
		
		T r_new[3], k_new[3];
		random_vec(r_new, r_max);
		random_vec(k_new, k_max<T>);
		copy_vector_to_particle_pos(part, r_new, i);
		copy_vector_to_particle_vel(part, k_new, i);
		T energy = compute_energy(part, ws, sigma_k, z, i);
		if(energy < epsilon) {
			mult_vec(r_new, r_new, T(-1.0));
			mult_vec(k_new, k_new, T(-1.0));
			copy_vector_to_particle_pos(part, r_new, i + 1);
			copy_vector_to_particle_vel(part, k_new, i + 1);
			i+=2;
		}
	}
}

template <typename T>
void set_initial_coulomb_boundaries(ScalarField<T> *coulomb, const World<T> &world, int z) {
	int nx = world.n[0], ny = world.n[1], nz = world.n[2];
	#pragma omp parallel for collapse(3)
	for(int i = 0; i < nx; i++) {
		for(int j = 0; j < ny; j++) {
			for(int k = 0; k < nz; k++) {
				if(i == 0 || j == 0 || k == 0 || i == nx - 1 || j == ny - 1 || k == nz - 1) {
					int idx = IDX(i, j, k, nx, ny, nz);
					T r_vec[3];
					world_pos_to_vector(r_vec, world, idx);
					T r = magnitude(r_vec);
					coulomb->v[idx] = T(1.44) * z / r;
				}
			}
		}
	}
}

template <typename T>
void chi_squared(const TestParticles<T> &part, const WoodsSaxon<T> *ws, Skyrme<T> skm, int part_per_nucleon) {
	int total = part.protons + part.neutrons;
	T chi_squared_p = 0.0, chi_squared_n = 0.0;
	#pragma omp parallel for reduction(+:chi_squared_p, chi_squared_n)
	for(int i = 0; i < total; i++) {
		int type;
		WoodsSaxon<T> ws_c;
		if(i >= part.protons) { type = NEUTRONS; ws_c = ws[1]; }
		else { type = PROTONS; ws_c = ws[0]; }
		
		T r_vec[3];
		copy_particle_pos_to_vector(r_vec, part, i);
		T density_p = part.density_p[i];
		T density_n = part.density_n[i];
		T r = magnitude(r_vec);
		T v_ws = woods_saxon_potential(ws_c, r);
		T v_skyrme = skyrme_potential(skm, density_p, density_n, type);
		T diff = v_ws - v_skyrme;
		if(type == PROTONS)
			chi_squared_p += diff * diff;
		else
			chi_squared_n += diff * diff;
	}
	chi_squared_n /= part_per_nucleon; chi_squared_p /= part_per_nucleon;
	std::printf("CHI SQUARED P %0.2lf CHI SQUARED N %0.2lf\n", chi_squared_p, chi_squared_n);
}

template <typename T>
T mean_squared_radius(const TestParticles<T> &part, const World<T> &world, int type) {
	int start, end, part_num = 0;
	T d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2], r_sqr = T(0.0);
	
	if(type == PROTONS) { start = 0; end = part.protons; }
	else if(type == NEUTRONS) { start = part.protons; end = part.protons + part.neutrons; }
	
	#pragma omp parallel for reduction(+:r_sqr, part_num)
	for(int i = start; i < end; i++) {
		T r_vec[3];
		copy_particle_pos_to_vector(r_vec, part, i);
		
		if(r_vec[0] < -d_max_x || r_vec[0] > +d_max_x
		|| r_vec[1] < -d_max_y || r_vec[1] > +d_max_y
		|| r_vec[2] < -d_max_z || r_vec[2] > +d_max_z)
			continue;
		
		T r2 = dot(r_vec, r_vec);
		r_sqr += r2;
		part_num++;
	}
	return r_sqr / (T)part_num;
}

template <typename T>
void relax_woods_saxon(WoodsSaxon<T> *ws, WoodsSaxon<T> *ws_old, T coef) {
	for(int i = 0; i < 2; i++) {
		ws[i].V0 = coef * ws[i].V0 + (1.0 - coef) * ws_old[i].V0;
		ws[i].R12 = coef * ws[i].R12 + (1.0 - coef) * ws_old[i].R12;
		ws[i].a = coef * ws[i].a + (1.0 - coef) * ws_old[i].a;
	}
}

template <typename T>
void merge_volumetric_potentials(ScalarField<T> *potential_a, const ScalarField<T> &potential_b, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < world_size; i++) {
		potential_a->v[i] += potential_b.v[i];
	}
}

template <typename T>
void copy_scalar_field(ScalarField<T> *density_a, const ScalarField<T> &density_b, const World<T> &world) {
	int world_size = world.n[0] * world.n[1] * world.n[2];
	#pragma omp parallel for
	for(int i = 0; i < 2 * world_size; i++)
		density_a->v[i] = density_b.v[i];
}

template <typename T>
static inline void world_pos_to_vector(T *v, const World<T> &world, int idx) {
	int x = world.n[0], y = world.n[1], z = world.n[2];
	int i = idx / (y * z), j = (idx / z) % y, k = idx % z;
	v[0] = world.d_max[0] * (T(2.0) * i / x - T(1.0));
	v[1] = world.d_max[1] * (T(2.0) * j / y - T(1.0));
	v[2] = world.d_max[2] * (T(2.0) * k / z - T(1.0));
}

template <typename T>
static inline void copy_particle_pos_to_vector(T *v, const TestParticles<T> &part, int i) {
	v[0] = part.x[i];
	v[1] = part.y[i];
	v[2] = part.z[i];
}

template <typename T>
static inline void copy_particle_vel_to_vector(T *v, const TestParticles<T> &part, int i) {
	v[0] = part.kx[i];
	v[1] = part.ky[i];
	v[2] = part.kz[i];
}

template <typename T>
static inline void copy_vector_to_particle_pos(const TestParticles<T> *part, const T *v, int i) {
	part->x[i] = v[0];
	part->y[i] = v[1];
	part->z[i] = v[2];
}

template <typename T>
static inline void copy_vector_to_particle_vel(const TestParticles<T> *part, const T *v, int i) {
	part->kx[i] = v[0];
	part->ky[i] = v[1];
	part->kz[i] = v[2];
}

/*void compute_volumetric_density(ScalarField<T> *density, ParticleCount<T> part_count, World<T> world_visual, World<T> world_data, Parameters<T> param, int type) {
	int world_size_visual = world_visual.n[0] * world_visual.n[1] * world_visual.n[2], start, end;
	int world_size_data = world_data.n[0] * world_data.n[1] * world_data.n[2];
	
	if(type == PROTONS) { start = 0; end = world_size_data; }
	else if(type == NEUTRONS) { start = world_size_data; end = 2 * world_size_data; }
	else { start = 0; end = 2 * world_size_data; }
	
	double sigma_r = param.sigma_r, exp_term = 1.0 / (2.0 * sigma_r * sigma_r);
	double term = (1.0 / param.part_per_nucleon) * (1.0 / (pow(2.0 * pi<T> * sigma_r * sigma_r, 1.5)));
	#pragma omp parallel for
	for(int i = 0; i < world_size_visual; i++) {
		double r_i[3], r_j[3], diff[3];
		double fact, dist_squared, density = 0.0;
		world_pos_to_vector(r_i, world_visual, i);
		
		for(int j = start; j < end; j++) {
			int count = part_count.count[j];
			world_pos_to_vector(r_j, world_data, j % world_size_data);
			
			sub_vec(diff, r_i, r_j);
			dist_squared = dot(diff, diff);
			fact = exp(-dist_squared * exp_term);
			density += count * fact;
		}
		density->v[i] += density;
	}
	#pragma omp parallel for
	for(int i = 0; i < world_size_visual; i++)
		density->v[i] *= term;
}

void scatter_particles(ParticleCount<T> *part_count, TestParticles<T> *part, World<T> world) {
	double d_max_x = world.d_max[0], d_max_y = world.d_max[1], d_max_z = world.d_max[2];
	int x = world.n[0], y = world.n[1], z = world.n[2], world_size = x * y * z, total = part->protons + part->neutrons;
	#pragma omp parallel for
	for(int i = 0; i < total; i++) {
		double r_vec[3];
		copy_particle_pos_to_vector(r_vec, *part, i);
		
		int x_idx = (int)(x / 2.0 * (r_vec[0] / d_max_x + 1.0));
		int y_idx = (int)(y / 2.0 * (r_vec[1] / d_max_y + 1.0));
		int z_idx = (int)(z / 2.0 * (r_vec[2] / d_max_z + 1.0));
		
		if(x_idx < 0 || y_idx < 0 || z_idx < 0 || x_idx >= x || y_idx >= y || z_idx >= z)
			continue;
		
		int idx = x_idx * (y * z) + y_idx * z + z_idx;
		if(i >= part->protons)
			idx += world_size;
		
		#pragma omp atomic update
		part_count->count[idx] += 1;
	}
}*/

#endif
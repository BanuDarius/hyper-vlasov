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

#include <concepts>

#include "sim_structs.hpp"

template <std::floating_point T>
inline void world_pos_to_vector(std::array<T, 3> &v, const World<T> &world, int idx) {
	int x = world.n[0], y = world.n[1], z = world.n[2];
	int i = idx / (y * z), j = (idx / z) % y, k = idx % z;
	v[0] = world.d_max[0] * (T(2.0) * i / x - T(1.0));
	v[1] = world.d_max[1] * (T(2.0) * j / y - T(1.0));
	v[2] = world.d_max[2] * (T(2.0) * k / z - T(1.0));
}

template <std::floating_point T>
inline void particle_pos_to_vector(std::array<T, 3> &v, const TestParticles<T> &part, int i) {
	v[0] = part.x[i];
	v[1] = part.y[i];
	v[2] = part.z[i];
}

template <std::floating_point T>
inline void particle_vel_to_vector(std::array<T, 3> &v, const TestParticles<T> &part, int i) {
	v[0] = part.kx[i];
	v[1] = part.ky[i];
	v[2] = part.kz[i];
}

template <std::floating_point T>
inline void vector_to_particle_pos(TestParticles<T> &part, const std::array<T, 3> &v, int i) {
	part.x[i] = v[0];
	part.y[i] = v[1];
	part.z[i] = v[2];
}

template <std::floating_point T>
inline void vector_to_particle_vel(TestParticles<T> &part, const std::array<T, 3> &v, int i) {
	part.kx[i] = v[0];
	part.ky[i] = v[1];
	part.kz[i] = v[2];
}

template <std::floating_point T> T compute_energy(TestParticles<T> &part, const WoodsSaxon<T> &ws, T sigma_k, int z, int i);
template <std::floating_point T> void compute_particle_energies(TestParticles<T> &part, const WoodsSaxon<T> &ws, const Parameters<T> &param);
template <std::floating_point T> void generate_random_particles(TestParticles<T> &part, T r_max);
template <std::floating_point T> void generate_checking_particles(TestParticles<T> &part, const WoodsSaxon<T> &ws, const Parameters<T> &param, const Fermi<T> &fermi_levels);
template <std::floating_point T> void compute_coulomb_boundaries(ScalarField<T> &coulomb, const TestParticles<T> &part, const World<T> &world, int z);
template <std::floating_point T> T mean_squared_radius(const TestParticles<T> &part, const World<T> &world, int type);
template <std::floating_point T> std::array<T, 3> center_of_mass(const TestParticles<T> &part, const World<T> &world, int type);
template <std::floating_point T> std::array<T, 2> chi_squared(const ScalarField<T> &density, const WoodsSaxon<T> &ws, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world);
template <std::floating_point T> void relax_woods_saxon(WoodsSaxon<T> &ws, const WoodsSaxon<T> &ws_old, T coef);

#endif
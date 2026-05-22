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

#ifndef PHYSICS_H
#define PHYSICS_H

#include "sim_structs.hpp"

template <typename T> void initialize_particles(TestParticles<T> &part, WoodsSaxon<T> &ws, const Skyrme<T> &skm, Fermi<T> &fermi_levels, const Parameters<T> &param);
template <typename T> void compute_volumetric_coulomb_potentials_sor(ScalarField<T> &coulomb, const ScalarField<T> &density, const World<T> &world);
template <typename T> void compute_volumetric_forces_fdm(VectorField<T> &forces, const ScalarField<T> &potentials, const World<T> &world);
template <typename T> void compute_volumetric_densities(ScalarField<T> &density, ScalarField<T> &density_temp, const Parameters<T> &param, const World<T> &world);
template <typename T> void compute_current_velocity(VectorField<T> &current_velocity, VectorField<T> &current_velocity_temp, const ScalarField<T> &density, const Parameters<T> &param, const World<T> &world);
template <typename T> void center_momentum(TestParticles<T> &part, const World<T> &world);
template <typename T> void nuclear_excitation(TestParticles<T> &part, const Parameters<T> &param);
template <typename T> void update_momenta_half(TestParticles<T> &part, T dt);
template <typename T> void update_positions_full(TestParticles<T> &part, T dt);
template <typename T> void compute_volumetric_skyrme_potentials(ScalarField<T> &potentials, const ScalarField<T> &density, const Skyrme<T> &skm, const World<T> &world);

#endif
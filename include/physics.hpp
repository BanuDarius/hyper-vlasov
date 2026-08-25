// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#ifndef PHYSICS_H
#define PHYSICS_H

#include <concepts>

#include "sim_structs.hpp"

template <std::floating_point T> void initialize_particles(TestParticles<T> &part, WoodsSaxon<T> &ws, const Skyrme<T> &skm, Fermi<T> &fermi_levels, const Parameters<T> &param, const World<T> &world);
template <std::floating_point T> void compute_volumetric_coulomb_potentials_sor(ScalarField<T> &coulomb, const ScalarField<T> &density, const World<T> &world);
template <std::floating_point T> void compute_volumetric_forces_fdm(VectorField<T> &forces, const ScalarField<T> &potentials, const World<T> &world);
template <std::floating_point T> void compute_volumetric_densities(ScalarField<T> &density, ScalarField<T> &density_temp, const Parameters<T> &param, const World<T> &world);
template <std::floating_point T> void compute_current_velocity(VectorField<T> &current_velocity, VectorField<T> &current_velocity_temp, const ScalarField<T> &density, const Parameters<T> &param, const World<T> &world);
template <std::floating_point T> void center_momentum(TestParticles<T> &part, const World<T> &world);
template <std::floating_point T> void nuclear_excitation(TestParticles<T> &part, const Parameters<T> &param);
template <std::floating_point T> void update_momenta_half(TestParticles<T> &part, T dt);
template <std::floating_point T> void update_positions_full(TestParticles<T> &part, T dt);
template <std::floating_point T> void compute_volumetric_skyrme_potentials(ScalarField<T> &potentials, const ScalarField<T> &density, const Skyrme<T> &skm, const World<T> &world);

#endif
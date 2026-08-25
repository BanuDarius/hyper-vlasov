// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#ifndef INIT_H
#define INIT_H

#include <cstdio>
#include <concepts>

#include "sim_structs.hpp"

template <std::floating_point T> void set_parameters(Parameters<T> &param, int z, int n, int part_per_nucleon, int steps, int substeps, int density_samples, bool use_gpu, T sample_position, T sigma_k, T sigma_r, T t_f, T t_exc, T eta_exc);
template <std::floating_point T> void read_input_file(std::FILE *in, Skyrme<T> &skm, World<T> &world, Fermi<T> &fermi_levels, Parameters<T> &param, WoodsSaxon<T> &ws);

#endif
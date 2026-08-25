// Copyright (c) 2026 Banu Darius-Matei
// SPDX-License-Identifier: MIT

#ifndef SIMULATION_H
#define SIMULATION_H

#include <concepts>

#include "sim_structs.hpp"

template <std::floating_point T> void cpu_simulate(const char *output_directory, TestParticles<T> &part, const Skyrme<T> &skm, const Parameters<T> &param, const World<T> &world);

#endif
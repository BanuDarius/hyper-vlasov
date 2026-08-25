# Copyright (c) 2026 Banu Darius-Matei
# SPDX-License-Identifier: MIT

import scripts.sim_init as sim_init
import scripts.programs as programs
import scripts.plotting as plotting

use_gpu = False
use_floats = False

num_test_part = 3000
z = 50
n = 82
nx = 16
t_f = 1200.0
t_exc = 200.0
eta_exc = 0.20
steps = 1200
substeps = 4
d_max_scale = 1.3
density_samples = 256
sample_position = 0.0

V0 = -50.0
a = 0.66
A = -356.8
B = 303.9
C = 32.0
k_fwhm = 0.346
r_fwhm = 1.444
epsilon_p = -8.0
epsilon_n = -12.0
gamma = 7.0 / 6.0

if __name__ == "__main__":
    sim_parameters = sim_init.SimParameters(num_test_part, z, n, nx, t_f, t_exc, eta_exc, steps, substeps, d_max_scale, density_samples, sample_position,
    V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, use_floats, use_gpu)
    
    programs.run_simulation(sim_parameters)
    
    programs.compute_energy_spectrum(sim_parameters)
    
    plotting.plot_radius()
    plotting.plot_center_of_mass()
    plotting.plot_dipole(sim_parameters)
    plotting.plot_energy_spectrum(sim_parameters)
    plotting.plot_density_samples(sim_parameters)
    plotting.plot_density_samples_differences(sim_parameters)
    #plotting.plot_density_samples_differences_lines(sim_parameters)
    
    print("Hyper-Vlasov finished!\a")
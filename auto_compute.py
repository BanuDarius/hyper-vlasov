'''MIT License

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
SOFTWARE.'''

# ---------------------------------------------------------- #

import scripts.sim_init as sim_init
import scripts.programs as programs
import scripts.plotting as plotting

# ---------------------------------------------------------- #

use_gpu = False
use_floats = False

num_test_part = 5000
z = 10
n = 16
nx = 16
t_f = 800.0
t_exc = 100.0
eta_exc = 0.05
steps = 1600
substeps = 4

d_max_scale = 1.3
V0 = -50.0
a = 0.66
A = -356.8
B = 303.9
C = 32.0
gamma = 7.0 / 6.0
epsilon_p = -8.0
epsilon_n = -12.0
k_fwhm = 0.346
r_fwhm = 1.444

# ---------------------------------------------------------- #

if __name__ == "__main__":
    sim_parameters = sim_init.SimParameters(num_test_part, z, n, nx, t_f, t_exc, eta_exc,
    steps, substeps, d_max_scale, V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, use_floats, use_gpu)
    
    programs.run_simulation(sim_parameters)
    
    programs.compute_energy_spectrum(sim_parameters)
    
    plotting.plot_radius()
    plotting.plot_center_of_mass()
    plotting.plot_energy_spectrum()
    
    print("Hyper-Vlasov finished!\a")

# ---------------------------------------------------------- #
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

# ---------------------------------------------------------- #

num_test_part = 1000
z = 16
n = 20
nx = 16
t_f = 600
steps = 1200
substeps = 4
use_floats = False

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
    sim_parameters = sim_init.SimParameters(num_test_part, z, n, nx, t_f, steps, substeps, V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, use_floats)
    
    programs.run_simulation(sim_parameters)
    
    print("Simulation finished.\a")

# ---------------------------------------------------------- #
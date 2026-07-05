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

from pathlib import Path

MAIN_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = MAIN_DIR.parent
INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_DIR = PROJECT_ROOT / "output"

# ---------------------------------------------------------- #

class SimParameters():
    def __init__(self, num_test_part, z, n, nx, t_f, t_exc, eta_exc, steps, substeps, d_max_scale, density_samples, sample_position, V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, use_floats, use_gpu):
        self.a = a
        self.A = A
        self.B = B
        self.C = C
        self.z = z
        self.n = n
        self.nx = nx
        self.V0 = V0
        self.t_f = t_f
        self.t_exc = t_exc
        self.steps = steps
        self.gamma = gamma
        self.r_fwhm = r_fwhm
        self.k_fwhm = k_fwhm
        self.use_gpu = use_gpu
        self.eta_exc = eta_exc
        self.substeps = substeps
        self.epsilon_p = epsilon_p
        self.epsilon_n = epsilon_n
        self.use_floats = use_floats
        self.d_max_scale = d_max_scale
        self.num_test_part = num_test_part
        self.output_directory = OUTPUT_DIR
        self.density_samples = density_samples
        self.sample_position = sample_position
        self.input_file = INPUT_DIR / "input.txt"
        
# ---------------------------------------------------------- #
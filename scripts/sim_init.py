from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_DIR = PROJECT_ROOT / "output"

# ---------------------------------------------------------- #

class SimParameters():
    def __init__(self, num_test_part, z, n, nx, t_f, steps, substeps, V0, a, A, B, C, gamma, epsilon_p, epsilon_n, k_fwhm, r_fwhm, use_floats):
        self.a = a
        self.A = A
        self.B = B
        self.C = C
        self.z = z
        self.n = n
        self.nx = nx
        self.V0 = V0
        self.t_f = t_f
        self.steps = steps
        self.gamma = gamma
        self.r_fwhm = r_fwhm
        self.k_fwhm = k_fwhm
        self.substeps = substeps
        self.epsilon_p = epsilon_p
        self.epsilon_n = epsilon_n
        self.use_floats = use_floats
        self.num_test_part = num_test_part
        self.output_directory = OUTPUT_DIR
        self.input_file = INPUT_DIR / "input.txt"
        
# ---------------------------------------------------------- #
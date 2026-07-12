
from types import SimpleNamespace
import numpy as np
from read_data import processing_data
from diffusion_mechanism_emcha import run_emcha_symmetry_analysis

config = SimpleNamespace(**processing_data())

_, dec = run_emcha_symmetry_analysis(
    config,
    voltage_values=np.linspace(0.1, 1.0, 10),
    n_repeats=3,
    random_seed=1,
    save_csv_prefix='results/emcha_fast',
    make_plot=True,
)

print(dec[dec['chirality'] == 0][['voltage_abs', 'emcha_total', 'emcha_spin']])


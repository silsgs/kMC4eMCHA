import os
import time
from types import SimpleNamespace

import imageio
import matplotlib.pyplot as plt
import numpy as np

from diffusion_mechanism import apply_diffusion_mechanism
from plotting import plotting
from read_data import processing_data


def build_initial_state(config):
    alpha_state_matrix = np.zeros((config.alpha_spins, config.n_steps), dtype=int)
    beta_state_matrix = np.zeros((config.beta_spins, config.n_steps), dtype=int)

    alpha_state_matrix[:, 0] = config.alpha_init_position
    beta_state_matrix[:, 0] = config.beta_init_position

    return alpha_state_matrix, beta_state_matrix


def run_single_simulation(config, ciss_effect=None, dV=None):
    alpha_state_matrix, beta_state_matrix = build_initial_state(config)
    return apply_diffusion_mechanism(
        config,
        alpha_state_matrix,
        beta_state_matrix,
        config.ciss_effect if ciss_effect is None else ciss_effect,
        config.dV if dV is None else dV,
    )


def run_qciss_frame_sweep(config):
    ciss_values = np.linspace(0, 1, 40)
    frames = []

    for ciss_effect in ciss_values:
        _, df_alpha, df_beta, *_ = run_single_simulation(config, ciss_effect=ciss_effect)
        plotting(config, df_alpha, df_beta, 0, 0, ciss_effect)
        frames.append(os.path.join("results_" + config.simulation_name, "temp_dir", f"drain_electrons_qciss_{ciss_effect:.2f}.png"))

    with imageio.get_writer(os.path.join("results_" + config.simulation_name, "drained_electrons_qciss.gif"), mode="I", duration=0.5) as writer:
        for frame in frames:
            writer.append_data(imageio.imread(frame))

    for frame in frames:
        os.remove(frame)
    os.rmdir(os.path.join("results_" + config.simulation_name, "temp_dir"))


def run_voltage_sweep(config):
    voltage_magnitude = np.linspace(-0.10, 0.10, 50)   ## hardcoded voltage range
    diff_alpha = np.zeros(len(voltage_magnitude))
    diff_beta = np.zeros(len(voltage_magnitude))
    diff_total = np.zeros(len(voltage_magnitude))

    for i, voltage in enumerate(voltage_magnitude):
        df_summary, *_ = run_single_simulation(config, dV=voltage / config.positions)
        diff_alpha[i] = df_summary.loc[0, "I difference between drain & source"]
        diff_beta[i] = df_summary.loc[1, "I difference between drain & source"]
        diff_total[i] = (diff_alpha[i] + diff_beta[i]) / 2

    plotting(config, diff_alpha, diff_beta, diff_total, voltage_magnitude, config.ciss_effect)
    
    output_path = os.path.join("results_" + config.simulation_name, "voltage_sweep_results.txt")
    with open(output_path, "w", encoding="utf-8") as out_f:
        out_f.write("#Voltage (V) \t I difference alpha (A) \t I difference beta (A) \t Total difference (A)\n")
        for v, da, db, dt in zip(voltage_magnitude, diff_alpha, diff_beta, diff_total):
            out_f.write(f"{v:.4f}, {da:.6e}, {db:.6e}, {dt:.6e}\n")


def run_qciss_current_sweep(config):
    ciss_values = np.linspace(0, 1, 40)
    alpha_current = []
    beta_current = []

    for ciss_effect in ciss_values:
        df_summary, *_ = run_single_simulation(config, ciss_effect=ciss_effect)
        alpha_current.append(df_summary.loc[0, "I difference between drain & source"])
        beta_current.append(df_summary.loc[1, "I difference between drain & source"])

    plt.plot(ciss_values, alpha_current, label="alpha spins")
    plt.plot(ciss_values, beta_current, label="beta spins")
    plt.xlabel("CISS effect")
    plt.ylabel("I difference between drain & source")
    plt.legend()
    plt.show()


def main():
    start_time = time.time()
    config = SimpleNamespace(**processing_data())

    if config.setting_seed == 1:
        np.random.seed(0)

    if config.save_results == 1:
        os.makedirs("results_" + config.simulation_name, exist_ok=True)

    print("Starting simulation...")
    print(f"Temperature = {config.Temperature} K")
    print(f"Simulation type = {config.simulation_type}")

    if config.simulation_type == 0:
        _, df_alpha, df_beta, *_ = run_single_simulation(config)
        plotting(config, df_alpha, df_beta, 0, 0, config.ciss_effect)
    elif config.simulation_type == 1:
        run_qciss_frame_sweep(config)
    elif config.simulation_type == 5:
        run_voltage_sweep(config)
    elif config.simulation_type == 10:
        run_qciss_current_sweep(config)
    elif config.simulation_type in (6, 7):
        _, _, _, alpha_state_matrix, beta_state_matrix, *_ = run_single_simulation(config)
        plotting(config, alpha_state_matrix, beta_state_matrix, 0, 0, config.ciss_effect)
    else:
        raise ValueError(
            "main.py currently supports simulation_type 0, 1, 5, 6, 7, and 10. "
            f"Received {config.simulation_type}."
        )

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} s")


if __name__ == "__main__":
    main()

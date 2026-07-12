import os
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


kB = 8.617333262e-5  # eV/K
amperec = 6.241509074e18  # electrons/s -> A


def build_initial_state(config):
    alpha_state_matrix = np.zeros((config.alpha_spins, config.n_steps), dtype=int)
    beta_state_matrix = np.zeros((config.beta_spins, config.n_steps), dtype=int)
    alpha_state_matrix[:, 0] = config.alpha_init_position
    beta_state_matrix[:, 0] = config.beta_init_position
    return alpha_state_matrix, beta_state_matrix


def diffusion_mechanism(
    spin_type,
    config,
    state_matrix,
    ciss_effect,
    dV,
    chirality=None,
    magnetic_field=10.0,
    temperature=None,
):
    """
    eMCHA-oriented hopping model.

    The drive contains an ordinary voltage term plus a second-order,
    chirality-dependent correction:

        drive = dV + ciss_effect * chirality * spin_type * B * dV**2

    The second term is even in voltage, odd in chirality, and spin dependent.
    """
    if chirality is None:
        chirality = getattr(config, "helix_twisting", 1)
    if temperature is None:
        temperature = getattr(config, "Temperature", 300)

    drained_spins = np.zeros(config.n_steps, dtype=int)
    sourced_spins = np.zeros(config.n_steps, dtype=int)

    emcha_drive = dV + ciss_effect * chirality * spin_type * magnetic_field * dV**2
    right_mov_probability = (
        1 / (1 + np.exp(-(emcha_drive / (kB * temperature))))
    ) * config.diff_coefficient
    left_mov_probability = config.diff_coefficient - right_mov_probability

    max_position = config.positions
    diff_coeff = config.diff_coefficient

    for i in range(1, config.n_steps):
        random_numbers = np.random.rand(len(state_matrix))
        prev_state = state_matrix[:, i - 1]

        right_mask = random_numbers <= right_mov_probability
        right_not_max_mask = right_mask & (prev_state != max_position)
        right_max_mask = right_mask & (prev_state == max_position)

        left_mask = (random_numbers > right_mov_probability) & (random_numbers <= diff_coeff)
        left_not_zero_mask = left_mask & (prev_state != 0)
        left_zero_mask = left_mask & (prev_state == 0)

        state_matrix[right_not_max_mask, i] = prev_state[right_not_max_mask] + 1
        state_matrix[right_max_mask, i] = 0
        drained_spins[i] += np.sum(right_max_mask)

        state_matrix[left_not_zero_mask, i] = prev_state[left_not_zero_mask] - 1
        state_matrix[left_zero_mask, i] = max_position
        sourced_spins[i] += np.sum(left_zero_mask)

        stay_mask = ~(
            right_not_max_mask
            | right_max_mask
            | left_not_zero_mask
            | left_zero_mask
        )
        state_matrix[stay_mask, i] = prev_state[stay_mask]

    total_drained_spins = int(np.sum(drained_spins))
    total_sourced_spins = int(np.sum(sourced_spins))
    I_difference = (total_drained_spins - total_sourced_spins) / config.n_steps

    time_step_length = 1e-14
    I_difference_A = (I_difference / time_step_length) / amperec

    if spin_type == -1:
        df = pd.DataFrame({
            "Step": np.arange(config.n_steps),
            "Drained alpha spins": drained_spins,
            "Sourced alpha spins": sourced_spins,
            "Right Probability (alpha)": right_mov_probability,
            "Left Probability (alpha)": left_mov_probability,
        })
    elif spin_type == 1:
        df = pd.DataFrame({
            "Step": np.arange(config.n_steps),
            "Drained beta spins": drained_spins,
            "Sourced beta spins": sourced_spins,
            "Right Probability (beta)": right_mov_probability,
            "Left Probability (beta)": left_mov_probability,
        })
    else:
        raise ValueError("spin_type must be -1 for alpha or 1 for beta")

    return (
        state_matrix,
        total_drained_spins,
        total_sourced_spins,
        df,
        float(np.mean(right_mov_probability)),
        float(np.mean(left_mov_probability)),
        I_difference_A,
        emcha_drive,
    )


def apply_diffusion_mechanism(
    config,
    alpha_state_matrix,
    beta_state_matrix,
    ciss_effect,
    dV,
    chirality=None,
    magnetic_field=10.0,
    temperature=None,
):
    alpha_state_matrix, total_drained_alpha_spins, total_sourced_alpha_spins, df_alpha, \
        r_prob_mean_alpha, l_prob_mean_alpha, I_diff_alpha, drive_alpha = (
            diffusion_mechanism(
                -1,
                config,
                alpha_state_matrix,
                ciss_effect,
                dV,
                chirality=chirality,
                magnetic_field=magnetic_field,
                temperature=temperature,
            )
        )

    beta_state_matrix, total_drained_beta_spins, total_sourced_beta_spins, df_beta, \
        r_prob_mean_beta, l_prob_mean_beta, I_diff_beta, drive_beta = (
            diffusion_mechanism(
                1,
                config,
                beta_state_matrix,
                ciss_effect,
                dV,
                chirality=chirality,
                magnetic_field=magnetic_field,
                temperature=temperature,
            )
        )

    df_summary = pd.DataFrame({
        "Spin Type": ["Spin alpha", "Spin beta"],
        "Total Drained electrons": [total_drained_alpha_spins, total_drained_beta_spins],
        "Total Sourced electrons": [total_sourced_alpha_spins, total_sourced_beta_spins],
        "Average Right Probability": [r_prob_mean_alpha, r_prob_mean_beta],
        "Average Left Probability": [l_prob_mean_alpha, l_prob_mean_beta],
        "I difference between drain & source": [I_diff_alpha, I_diff_beta],
        "eMCHA drive": [drive_alpha, drive_beta],
    })

    return (
        df_summary,
        df_alpha,
        df_beta,
        alpha_state_matrix,
        beta_state_matrix,
        drive_alpha,
        drive_beta,
    )


def run_single_simulation(
    config,
    voltage,
    chirality,
    ciss_effect=None,
    magnetic_field=10.0,
    temperature=None,
):
    alpha_state_matrix, beta_state_matrix = build_initial_state(config)
    df_summary, *_ = apply_diffusion_mechanism(
        config,
        alpha_state_matrix,
        beta_state_matrix,
        config.ciss_effect if ciss_effect is None else ciss_effect,
        voltage / config.positions,
        chirality=chirality,
        magnetic_field=magnetic_field,
        temperature=temperature,
    )
    alpha_current = df_summary.loc[0, "I difference between drain & source"]
    beta_current = df_summary.loc[1, "I difference between drain & source"]
    total_current = (alpha_current + beta_current) / 2
    spin_current = (alpha_current - beta_current) / 2
    return alpha_current, beta_current, total_current, spin_current


def _simulate_current_only(
    spin_type,
    config,
    ciss_effect,
    voltage,
    chirality,
    magnetic_field,
    temperature,
):
    dV = voltage / config.positions
    drive = dV + ciss_effect * chirality * spin_type * magnetic_field * dV**2
    right_probability = (
        1 / (1 + np.exp(-(drive / (kB * temperature))))
    ) * config.diff_coefficient
    diff_coeff = config.diff_coefficient

    n_spins = config.alpha_spins if spin_type == -1 else config.beta_spins
    init_position = config.alpha_init_position if spin_type == -1 else config.beta_init_position
    positions = np.full(n_spins, init_position, dtype=np.int16)

    drained_spins = 0
    sourced_spins = 0
    max_position = config.positions

    for _ in range(1, config.n_steps):
        random_numbers = np.random.rand(n_spins)

        right_mask = random_numbers <= right_probability
        right_max_mask = right_mask & (positions == max_position)
        right_not_max_mask = right_mask & ~right_max_mask

        left_mask = (random_numbers > right_probability) & (random_numbers <= diff_coeff)
        left_zero_mask = left_mask & (positions == 0)
        left_not_zero_mask = left_mask & ~left_zero_mask

        positions[right_not_max_mask] += 1
        positions[right_max_mask] = 0
        drained_spins += int(np.sum(right_max_mask))

        positions[left_not_zero_mask] -= 1
        positions[left_zero_mask] = max_position
        sourced_spins += int(np.sum(left_zero_mask))

    current_steps = (drained_spins - sourced_spins) / config.n_steps
    return (current_steps / 1e-14) / amperec


def run_single_simulation_fast(
    config,
    voltage,
    chirality,
    ciss_effect=None,
    magnetic_field=10.0,
    temperature=None,
):
    if temperature is None:
        temperature = getattr(config, "Temperature", 300)
    if ciss_effect is None:
        ciss_effect = config.ciss_effect

    alpha_current = _simulate_current_only(
        -1,
        config,
        ciss_effect,
        voltage,
        chirality,
        magnetic_field,
        temperature,
    )
    beta_current = _simulate_current_only(
        1,
        config,
        ciss_effect,
        voltage,
        chirality,
        magnetic_field,
        temperature,
    )
    total_current = (alpha_current + beta_current) / 2
    spin_current = (alpha_current - beta_current) / 2
    return alpha_current, beta_current, total_current, spin_current


def run_emcha_symmetry_analysis(
    config,
    voltage_values=None,
    ciss_effect=None,
    magnetic_field=10.0,
    temperature=None,
    n_repeats=10,
    random_seed=None,
    save_csv_prefix=None,
    make_plot=False,
    fast=True,
):
    """
    Measure the second-order eMCHA candidate without changing main.py.

    The analysis runs paired simulations at +V and -V for both chiralities.
    It returns:
      - raw_df: one row per Monte Carlo repeat and symmetry condition
      - decomposed_df: odd/even voltage parts and the chirality-odd eMCHA signal

    The eMCHA observable is:

        eMCHA(V) = [I_even(chirality=+1) - I_even(chirality=-1)] / 2

    where I_even = [I(+V) + I(-V)] / 2.
    """
    if isinstance(config, dict):
        config = SimpleNamespace(**config)

    if voltage_values is None:
        voltage_values = np.linspace(0.05, abs(config.voltage_magnitude), 20)
    voltage_values = np.asarray(voltage_values, dtype=float)
    voltage_values = voltage_values[voltage_values > 0]

    rng = np.random.default_rng(random_seed)
    raw_rows = []
    simulation_runner = run_single_simulation_fast if fast else run_single_simulation

    for repeat in range(n_repeats):
        repeat_seed = None if random_seed is None else int(rng.integers(0, 2**32 - 1))

        for voltage_abs in voltage_values:
            for chirality in (-1, 1):
                for voltage_sign in (-1, 1):
                    if repeat_seed is not None:
                        seed_offset = (voltage_sign + 1) + 2 * (chirality + 1)
                        np.random.seed((repeat_seed + seed_offset) % (2**32 - 1))

                    voltage = voltage_sign * voltage_abs
                    alpha_current, beta_current, total_current, spin_current = simulation_runner(
                        config,
                        voltage,
                        chirality,
                        ciss_effect=ciss_effect,
                        magnetic_field=magnetic_field,
                        temperature=temperature,
                    )
                    raw_rows.append({
                        "repeat": repeat,
                        "voltage_abs": voltage_abs,
                        "voltage": voltage,
                        "voltage_sign": voltage_sign,
                        "chirality": chirality,
                        "alpha_current": alpha_current,
                        "beta_current": beta_current,
                        "total_current": total_current,
                        "spin_current": spin_current,
                    })

    raw_df = pd.DataFrame(raw_rows)
    mean_df = (
        raw_df
        .groupby(["voltage_abs", "voltage_sign", "chirality"], as_index=False)
        [["alpha_current", "beta_current", "total_current", "spin_current"]]
        .mean()
    )

    decomposed_rows = []
    for voltage_abs in voltage_values:
        by_chirality = {}
        for chirality in (-1, 1):
            plus = mean_df[
                (mean_df["voltage_abs"] == voltage_abs)
                & (mean_df["chirality"] == chirality)
                & (mean_df["voltage_sign"] == 1)
            ].iloc[0]
            minus = mean_df[
                (mean_df["voltage_abs"] == voltage_abs)
                & (mean_df["chirality"] == chirality)
                & (mean_df["voltage_sign"] == -1)
            ].iloc[0]

            total_odd = (plus["total_current"] - minus["total_current"]) / 2
            total_even = (plus["total_current"] + minus["total_current"]) / 2
            spin_odd = (plus["spin_current"] - minus["spin_current"]) / 2
            spin_even = (plus["spin_current"] + minus["spin_current"]) / 2

            by_chirality[chirality] = {
                "total_odd": total_odd,
                "total_even": total_even,
                "spin_odd": spin_odd,
                "spin_even": spin_even,
            }

            decomposed_rows.append({
                "voltage_abs": voltage_abs,
                "chirality": chirality,
                "total_odd": total_odd,
                "total_even": total_even,
                "spin_odd": spin_odd,
                "spin_even": spin_even,
                "emcha_total": np.nan,
                "emcha_spin": np.nan,
            })

        emcha_total = (by_chirality[1]["total_even"] - by_chirality[-1]["total_even"]) / 2
        emcha_spin = (by_chirality[1]["spin_even"] - by_chirality[-1]["spin_even"]) / 2

        decomposed_rows.append({
            "voltage_abs": voltage_abs,
            "chirality": 0,
            "total_odd": np.nan,
            "total_even": np.nan,
            "spin_odd": np.nan,
            "spin_even": np.nan,
            "emcha_total": emcha_total,
            "emcha_spin": emcha_spin,
        })

    decomposed_df = pd.DataFrame(decomposed_rows)

    if save_csv_prefix:
        output_dir = os.path.dirname(save_csv_prefix)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        raw_df.to_csv(f"{save_csv_prefix}_raw.csv", index=False)
        decomposed_df.to_csv(f"{save_csv_prefix}_decomposed.csv", index=False)

    if make_plot:
        emcha_df = decomposed_df[decomposed_df["chirality"] == 0]
        plt.figure()
        plt.plot(emcha_df["voltage_abs"] ** 2, emcha_df["emcha_total"], "o-")
        plt.xlabel("Voltage^2")
        plt.ylabel("eMCHA total current candidate (A)")
        plt.title("eMCHA scaling check")
        if save_csv_prefix:
            plt.savefig(f"{save_csv_prefix}_scaling.png", dpi=150, bbox_inches="tight")
        plt.show()

    return raw_df, decomposed_df


apply_diffussion_mechanism = apply_diffusion_mechanism

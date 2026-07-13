import numpy as np
import pandas as pd

####  cte definitions
kB = 8.617333262*10**(-5) # eV/K
amperec = 6.241509074 * 10 **(18) ## from e/s --> ampere

def diffusion_mechanism(spin_type, config, state_matrix, emcha_effect, dV):
    """
    Simulate one spin population diffusing on a 1D lattice.

    Particles move right, move left, or stay in place at each step. A particle
    crossing the right boundary is counted as drained and re-enters at 0; a
    particle crossing the left boundary is counted as sourced and re-enters at
    the maximum position.
    """
    drained_spins = np.zeros(config.n_steps, dtype=int)
    sourced_spins = np.zeros(config.n_steps, dtype=int)

    ## Testing
    B = 10  # units
    T = 300  # Kelvin
    helix_twisting  = 1 
    emcha_contribution = dV + emcha_effect * helix_twisting  * spin_type * B * dV**2
    D = config.diff_coefficient

    right_mov_probability = (
        1/(1+np.exp(-(emcha_contribution/(kB*T))))) * D

    #emcha_contribution = emcha_effect * spin_type * np.tanh(10 * dV)
    #right_mov_probability = ( 1 / (1 + np.exp( - ( (dV) / (kB * 300) ) * (1 + emcha_contribution))) ) * config.diff_coefficient
    #right_mov_probability = ( 1 / (1 + np.exp( -( 10* dV ) * (1 + emcha_contribution))) ) * config.diff_coefficient
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

        stay_mask = ~( right_not_max_mask | right_max_mask | left_not_zero_mask | left_zero_mask )
        state_matrix[stay_mask, i] = prev_state[stay_mask]

    total_drained_spins = int(np.sum(drained_spins))
    total_sourced_spins = int(np.sum(sourced_spins))
    
    I_drain = total_drained_spins / config.n_steps
    I_source = total_sourced_spins / config.n_steps
    I_difference = I_drain - I_source

    time_step_length = 10**(-14) # seconds

    I_drain_s = I_drain / time_step_length
    I_source_s = I_source / time_step_length
    I_diff_s = I_difference / time_step_length
    
    I_drain_A = I_drain_s / amperec
    I_source_A = I_source_s / amperec
    I_diff_A = I_diff_s / amperec

    I_difference = I_diff_A   # Current Intensity difference in Ampere

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
        I_difference,
        emcha_contribution,
    )


def apply_diffusion_mechanism(config, alpha_state_matrix, beta_state_matrix, emcha_effect, dV):
    alpha_state_matrix, total_drained_alpha_spins, total_sourced_alpha_spins, df_alpha, \
        r_prob_mean_alpha, l_prob_mean_alpha, I_diff_alpha, emcha_contribution_alpha = (
            diffusion_mechanism(-1, config, alpha_state_matrix, emcha_effect, dV)
        )

    beta_state_matrix, total_drained_beta_spins, total_sourced_beta_spins, df_beta, \
        r_prob_mean_beta, l_prob_mean_beta, I_diff_beta, emcha_contribution_beta = (
            diffusion_mechanism(1, config, beta_state_matrix, emcha_effect, dV)
        )

    df_summary = pd.DataFrame({
        "Spin Type": ["Spin alpha", "Spin beta"],
        "Total Drained electrons": [total_drained_alpha_spins, total_drained_beta_spins],
        "Total Sourced electrons": [total_sourced_alpha_spins, total_sourced_beta_spins],
        "Average Right Probability": [r_prob_mean_alpha, r_prob_mean_beta],
        "Average Left Probability": [l_prob_mean_alpha, l_prob_mean_beta],
        "I difference between drain & source": [I_diff_alpha, I_diff_beta],
    })

    return (
        df_summary,
        df_alpha,
        df_beta,
        alpha_state_matrix,
        beta_state_matrix,
        emcha_contribution_alpha,
        emcha_contribution_beta,
    )


# Backwards-compatible name used by the original scripts.
apply_diffussion_mechanism = apply_diffusion_mechanism

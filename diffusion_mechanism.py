import numpy as np
import pandas as pd

####  cte definitions
kB = 8.617333262*10**(-5) # eV/K
amperec = 6.241509074 * 10 **(18) ## from e/s --> ampere


def diffusion_mechanism(spin_type, config, state_matrix, ciss_effect, dV):
    """
    Simulate one spin population diffusing on a 1D lattice.

    Particles move right, move left, or stay in place at each step. A particle
    crossing the right boundary is counted as drained and re-enters at 0; a
    particle crossing the left boundary is counted as sourced and re-enters at
    the maximum position.
    """
    ciss_effect = config.ciss_effect
    helix_twisting  = config.helix_twisting  # signo : signo del campo magnetico que genera el e- pasando por la helice

    drained_spins = np.zeros(config.n_steps, dtype=int)
    sourced_spins = np.zeros(config.n_steps, dtype=int)
     
    ciss_contribution = np.exp( ciss_effect * spin_type * helix_twisting * np.abs(np.tanh( dV )) ) 
    
    """
    #### v 1.2 --> not valid yet
    barrier = 0.1 # eV 
    bias = dV
    T = 300 # harcoded for testing
    #print(T)

    ciss_contribution = np.exp( ciss_effect * spin_type * helix_twisting * np.abs(np.tanh( dV )) ) 
    
    E_right = barrier - bias
    E_left = barrier + bias
    E_stay = 0
    fund = min(E_right, E_stay, E_left)

    norm_E_right = E_right - (fund)
    norm_E_left = E_left - (fund)
    norm_E_stay = E_stay - (fund)

    bltz_E_right = ( np.exp(-(norm_E_right / (kB * T) ))) * ciss_contribution
    bltz_E_left = ( np.exp(-(norm_E_left / (kB * T) ))) / ciss_contribution
    bltz_E_stay = np.exp( (norm_E_stay / ( kB * T) ) )

    #bltz_E_right = np.exp( (norm_E_right / ( kB * T) ) ) * ( 1+ np.exp(ciss_contribution) ) 
    #bltz_E_left = np.exp( (norm_E_left / ( kB * T) ) ) / ( 1+ np.exp(ciss_contribution) ) 
    #bltz_E_stay = np.exp( (norm_E_stay / ( kB * T) ) )

    p_bltz_E_right = bltz_E_right / (bltz_E_right + bltz_E_left + bltz_E_stay)
    p_bltz_E_left = bltz_E_left / (bltz_E_right + bltz_E_left + bltz_E_stay)

    right_mov_probability = p_bltz_E_right 
    left_mov_probability = p_bltz_E_left

    print(f"Right  probability: {right_mov_probability}")
    print(f"Left  probability: {left_mov_probability}")
    print(f"CISS contribution: {ciss_contribution}")

    diff_coeff = right_mov_probability + left_mov_probability
    #### end v 1.2
    """

    #### v 1.1 --> stable
    D = 0.5
    #T = 300 # harcoded for testing
    T = config.Temperature
    ciss_contribution = ciss_effect * spin_type * helix_twisting * np.tanh( dV )
    
    right_mov_probability = (
        1/ ( 1 + np.exp( -( ( dV / (kB * T) ) * ( np.exp(ciss_contribution) ) ) ) ) ) * D 
    
    left_mov_probability = D - right_mov_probability

    diff_coeff = config.diff_coefficient
    #### end v 1.1

    max_position = config.positions

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

    time_step_length = 1.0 * 10**(-17) # seconds

    I_diff_s = I_difference / time_step_length
    
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
        ciss_contribution,
    )


def apply_diffusion_mechanism(config, alpha_state_matrix, beta_state_matrix, ciss_effect, dV):
    alpha_state_matrix, total_drained_alpha_spins, total_sourced_alpha_spins, df_alpha, \
        r_prob_mean_alpha, l_prob_mean_alpha, I_diff_alpha, ciss_contribution_alpha = (
            diffusion_mechanism(-1, config, alpha_state_matrix, ciss_effect, dV)
        )

    beta_state_matrix, total_drained_beta_spins, total_sourced_beta_spins, df_beta, \
        r_prob_mean_beta, l_prob_mean_beta, I_diff_beta, ciss_contribution_beta = (
            diffusion_mechanism(1, config, beta_state_matrix, ciss_effect, dV)
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
        ciss_contribution_alpha,
        ciss_contribution_beta,
    )

#def magnetic_field():
    

# Backwards-compatible name used by the original scripts.
apply_diffussion_mechanism = apply_diffusion_mechanism

#%%
"""
Data used:
    - First simulations

# mu_B/k_B
mu_b = 0.67171388  # (K/T)

# magnetic field in T
B_arr = np.arange(0, 1, 0.1)

# spin 1/2
spin=0.5

# landau g factor
g_factor = 2.00232

# temperature in K
temp = 8 

## bias
V = 0.0765

# current from simulations
alpha_current_p = 0.16695
beta_current_p = 0.21643

alpha_current_n = -0.21626
beta_current_n = -0.16705


## bias
V = 0.15

# current from simulations
alpha_current_p = 0.26959 
beta_current_p = 0.41686

alpha_current_n = -0.41667
beta_current_n = -0.26949


    - Second simulations qEMCHA

qEMCHA = [0.90, 0.95, 0.99]

## bias
V = 0.25

## qEMCHA = 0.90

# current from simulations
# 0.2500, 3.415491e-01, 5.865533e-01, 4.640512e-01
alpha_current_p = 3.415491e-01
beta_current_p = 5.865533e-01

## -0.2500, -5.863181e-01, -3.413764e-01, -4.638473e-01
alpha_current_n = -5.863181e-01
beta_current_n = -3.413764e-01



"""

#%%

### Post analysis script

import numpy as np
import matplotlib.pyplot as plt

#%%
#####
# params
#####
# Bohr magneton/boltzmann's constant
# mu_B/k_B
mu_b = 0.67171388  # (K/T)
# magnetic field in T
B_arr = np.arange(0, 1, 0.1)
# spin 1/2
spin=0.5
# landau g factor
g_factor = 2.00232
# temperature in K
temp = 8 

#Zeeman energy
Ezee_arr = g_factor * mu_b * B_arr * spin

# delta Ezee
delta_Ezee_arr = Ezee_arr*2

# boltzmann populations
gs_pop_arr = 1/(1+np.exp(-delta_Ezee_arr/temp))
exc_pop_arr = ((np.exp(-delta_Ezee_arr/temp) / (1+np.exp(-delta_Ezee_arr/temp))))

# spin ground state
gs = 0 # ground state alpha = 0, beta = 1

# popultions
alpha_pop_arr, beta_pop_arr = (gs_pop_arr, exc_pop_arr) if gs == 0 else (exc_pop_arr, gs_pop_arr)


#####
# Voltage sweep
#####
## bias
V = 0.0765

# current from simulations
alpha_current_p = 0.16695
beta_current_p = 0.21643

alpha_current_n = -0.21626
beta_current_n = -0.16705


# total current
total_current_p_arr = (alpha_pop_arr * alpha_current_p) + (beta_pop_arr * beta_current_p)
total_current_n_arr = (alpha_pop_arr * alpha_current_n) + (beta_pop_arr * beta_current_n)

# R
R_arr_p = V / total_current_p_arr
R_arr_n = -V / total_current_n_arr

#### delta R
delta_R_arr = R_arr_n - R_arr_p

delta_R_prime = delta_R_arr - delta_R_arr[0]
deltaRR = delta_R_prime/0.4
gamma_arr = deltaRR/0.192/2 

print(gamma_arr)


#%%
#####
# Voltage sweep
#####
## bias
V = 0.15

# current from simulations
alpha_current_p = 0.26959
beta_current_p = 0.41686

alpha_current_n = -0.41667
beta_current_n = -0.26949

# total current
total_current_p_arr = (alpha_pop_arr * alpha_current_p) + (beta_pop_arr * beta_current_p)
total_current_n_arr = (alpha_pop_arr * alpha_current_n) + (beta_pop_arr * beta_current_n)

# R
R_arr_p = V / total_current_p_arr
R_arr_n = -V / total_current_n_arr

#### delta R
delta_R_arr = R_arr_n - R_arr_p
delta_R_prime = delta_R_arr - delta_R_arr[0]
deltaRR = delta_R_prime/0.4
gamma_arr = deltaRR/0.192/2 




#%%
# =============================================================================
# qEMCHA vs gamma
# =============================================================================

#####
# params
#####
# Bohr magneton/boltzmann's constant
# mu_B/k_B
mu_b = 0.67171388  # (K/T)
# magnetic field in T
B_arr = np.arange(0, 1, 0.1)
# spin 1/2
spin=0.5
# landau g factor
g_factor = 2.00232
# temperature in K
temp = 8 

#Zeeman energy
Ezee_arr = g_factor * mu_b * B_arr * spin

# delta Ezee
delta_Ezee_arr = Ezee_arr*2

# boltzmann populations
gs_pop_arr = 1/(1+np.exp(-delta_Ezee_arr/temp))
exc_pop_arr = ((np.exp(-delta_Ezee_arr/temp) / (1+np.exp(-delta_Ezee_arr/temp))))

# spin ground state
gs = 0 # ground state alpha = 0, beta = 1

# popultions
alpha_pop_arr, beta_pop_arr = (gs_pop_arr, exc_pop_arr) if gs == 0 else (exc_pop_arr, gs_pop_arr)


#####
# Voltage sweep
#####
## bias
V = 0.25

# current from simulations
# qEMCHA = 0.90 : 0.2500, 3.415491e-01, 5.865533e-01, 4.640512e-01
# qEMCHA = 0.95 : 0.2500, 3.304630e-01, 5.898823e-01, 4.601727e-01
# qEMCHA = 0.99 : 0.2500, 3.213287e-01, 5.924558e-01, 4.568922e-01
alpha_current_p = 3.415491e-01
beta_current_p = 5.865533e-01

# qEMCHA = 0.90 : -0.2500, -5.863181e-01, -3.413764e-01, -4.638473e-01
# qEMCHA = 0.95 : -0.2500, -5.896593e-01, -3.302704e-01, -4.599649e-01
# qEMCHA = 0.99 : -0.2500, -5.922247e-01, -3.211707e-01, -4.566977e-01
alpha_current_n = -5.922247e-01
beta_current_n = -3.211707e-01


# total current
total_current_p_arr = (alpha_pop_arr * alpha_current_p) + (beta_pop_arr * beta_current_p)
total_current_n_arr = (alpha_pop_arr * alpha_current_n) + (beta_pop_arr * beta_current_n)

# R
R_arr_p = V / total_current_p_arr
R_arr_n = -V / total_current_n_arr

#### delta R
delta_R_arr = R_arr_n - R_arr_p

delta_R_prime = delta_R_arr - delta_R_arr[0]
deltaRR = delta_R_prime/0.4
gamma_arr = deltaRR/0.192/2 

print(gamma_arr)


#


#%%
# =============================================================================
# TESTING
# =============================================================================

MU_B = 0.67171388  # mu_B / k_B (K/T)
G_FACTOR = 2.00232
SPIN = 0.5


def spin_populations(B, T, gs=0,
                     g_factor=G_FACTOR,
                     spin=SPIN,
                     mu_b=MU_B):
    """
    Calculate alpha and beta spin populations.

    Parameters
    ----------
    B : float or array-like
        Magnetic field (T).
    T : float
        Temperature (K).
    gs : int
        Ground-state spin:
            0 -> alpha
            1 -> beta

    Returns
    -------
    alpha_pop, beta_pop
    """
    B = np.asarray(B)

    delta_E = 2 * g_factor * mu_b * spin * B

    gs_pop = 1 / (1 + np.exp(-delta_E / T))
    exc_pop = np.exp(-delta_E / T) / (1 + np.exp(-delta_E / T))

    return (gs_pop, exc_pop) if gs == 0 else (exc_pop, gs_pop)

def compute_gamma(
    B,
    T,
    V,
    alpha_current_p,
    beta_current_p,
    alpha_current_n,
    beta_current_n,
    gs=0,
):
    """
    Compute resistance and gamma vs magnetic field.
    """

    alpha_pop, beta_pop = spin_populations(B, T, gs)

    current_p = (
        alpha_pop * alpha_current_p
        + beta_pop * beta_current_p
    )

    current_n = (
        alpha_pop * alpha_current_n
        + beta_pop * beta_current_n
    )

    R_p = V / current_p
    R_n = -V / current_n

    delta_R = R_n - R_p
    delta_R_prime = delta_R - delta_R[0]

    deltaRR = delta_R_prime / 0.4
    gamma = deltaRR / 0.192 / 2

    return {
        "gamma": gamma,
        "delta_R": delta_R,
        "R_p": R_p,
        "R_n": R_n,
        "current_p": current_p,
        "current_n": current_n,
        "alpha_pop": alpha_pop,
        "beta_pop": beta_pop,
    }

#%%
## barrido de campos

B = np.arange(0, 1.1, 0.1)
T = 8
V = 0.0765

results = compute_gamma(
    B=B,
    T=T,
    V=V,
    alpha_current_p=0.16695,
    beta_current_p=0.21643,
    alpha_current_n=-0.21626,
    beta_current_n=-0.16705,
)

gamma_Bs = results["gamma"]
print(gamma_Bs)

#%%
## barrido de temperaturas
temps = [300,250,200,150,100,50,25,10,8,5,4,3,2,1, 0.001]

gammas_temps = []

for i in temps:
    B = np.arange(0, 1.1, 0.1)
    T = i
    V = 0.0765

    results = compute_gamma(
        B=B,
        T=T,
        V=V,
        alpha_current_p=0.16695,
        beta_current_p=0.21643,
        alpha_current_n=-0.21626,
        beta_current_n=-0.16705,
    )

    gamma = results["gamma"]
    gammas_temps.append(gamma[-1])
    
print(gammas_temps)


#%%
# =============================================================================
# FIG
# =============================================================================
## init fig
plt.rcParams.update({
    "font.family": "arial",
    "font.size": 9,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

color1 = "#1f77b4"
color2 = "#b22222"

fig, ax = plt.subplots(1,2, layout='tight', figsize=(6.75,4))

##
ax[0].plot(
    temps,
    gammas_temps,
    color=color1,
    lw=1.4,
    marker="o",
    ms=5.5,
    #mfc="white",
    mec=color1,
    #mew=1.0,
)

ax[0].set_xlabel("Temperature, $T$(K)")
ax[0].set_ylabel("$\gamma*T$")

##
ax[1].plot(
    B,
    gamma_Bs,
    color=color2,
    lw=1.4,
    marker="s",
    ms=5.5,
    #mfc="white",
    mec=color2,
    #mew=1.0,
)

ax[1].set_xlabel("Magnetic field, $B$(T)")
ax[1].set_ylabel("$\gamma$")

for i in range(2):
    ax[i].grid(True)
    
    ax[i].spines['right'].set_visible(False)
    ax[i].spines['left'].set_visible(True)
    ax[i].spines['bottom'].set_visible(True)
    ax[i].spines['top'].set_visible(False)


plt.savefig(r"C:\Users\silvi\OneDrive - Trinity College Dublin\Projects\eMCHA-STOSS\paper\gammavsTB.png", 
            dpi = 300)



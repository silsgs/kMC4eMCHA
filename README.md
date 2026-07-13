# kMC4eMCHA code

## Code overview

This repository contains the code for the kMC4eMCHA project, which is designed to simulate spin-dependent charge transport in a one dimensional molecular system, motivated by electric mangnetochiral anisotropy effects in molecular systems. 
The simulation follows two electron spin populations, spin-up and spin-down labeled as alpha and beta electrons, as they hop between dicrete positions along a molecule under the combined effect of voltage bias, temerature, molecular handedness, diffusion strengthm and a phenomenological eMCHA term. 

The workflow is implemented as follows: 

1. The user specifies the simulation parameters in a configuration file "user_configurations.yaml".
2. Initial spin-position alpha/beta matrices are built. The user can define the initial spin-position distribution in the configuration file.
3. Simulation of stochastic hopping processes, source events, and drain events. 
4. Update and record final spin-position matrices, and final intensities of the spin populations.
5. Output the results, generate plots, write spin-position matrices and final intensities of the spin populations.

## Main files

### `user_configurations.yaml`

This file contains the user-facing simulation parameters. It controls the simulation name, system size, initial spin positions, physical/model parameters, random seed behavior, and output mode.

Important parameters include:

- `simulation_name`: name used to create the output directory `results_<simulation_name>`.
- `simulation_type`: selects the type of simulation or plot to run.
- `n_steps`: number of Monte Carlo time steps.
- `positions`: number of discrete positions along the molecule. Internally, positions run from `0` to `positions`.
- `number_spins`: total number of simulated spins/electrons.
- `spin_ratio`: controls the relative number of alpha and beta spins. A value of `1` gives equal populations.
- `alpha_init_position`, `beta_init_position`: initial positions for the two spin populations.
- `helix_twisting`: chirality/sign parameter of the molecular helix, usually `-1` or `1`.
- `emcha_effect`: strength of the spin/chirality contribution.
- `diff_coefficient`: total probability scale for attempted hopping.
- `Temperature`: simulation temperature in Kelvin.
- `voltage_magnitude`: voltage used to define the default voltage step `dV`.
- `setting_seed`: if set to `1`, NumPy's random seed is fixed for reproducible simulations.
- `save_results`: controls whether the result directory is explicitly created at startup.

### `read_data.py`

`read_data.py` handles configuration parsing and validation.

The function `read_config_file()` reads the configuration file as a lightweight YAML-like key/value file. It removes comments, ignores blank lines, and automatically casts values to `bool`, `int`, or `float` when possible.

The function `validate_config()` checks that all required configuration keys are present and verifies basic constraints, such as positive values for `n_steps`, `positions`, and `number_spins`, and a valid range for `diff_coefficient`.

### `main.py`

`main.py` is the main entry point of the simulator.

It first loads the configuration using `processing_data()`, optionally fixes the random seed, creates the results directory, and dispatches the selected simulation mode.

The helper function `build_initial_state()` creates two state matrices:

- `alpha_state_matrix`
- `beta_state_matrix`

Each matrix has shape:

```python
(number_of_spins, n_steps)
```

The first column stores the initial position of each spin, while later columns are updated during the simulation. 

The helper function `run_single_simulation()` initializes the spin matrices and passes them to the core simulation framework implemented in `diffussion_mechanism.py`. 

Supported `simulation_type` values include:
- 0: run a single simulation and plot cummulative drained/sourced electrons. 
- 1: sweep the EMCHA parameter and generate frames/GIFs of drained electrons. 
- 5: perform a voltage sweep and plot spin-resolved current response. 
- 6: generate a spin-polarization heatmap animation.
- 7:generate an animation of alpha/beta spin distributions along the molecule.
- 10: sweep the EMCHA parameter and plot current differences interactively. 

For voltage sweeps, the current version uses a hardcoded voltage range from -0.10 V to 0.10 V.

### `diffusion_mechanism.py`

This file contains the central transport model.

The simulator treats the molecule as a one-dimensional lattice. At each time step, every spin can:
- move one position to the right,
- move one position to the left,
- remain in place.

Boundary crossing events are interpreted as transport events:

- A spin moving right from the maximum position exits through the drain and re-enters at position 0.
- A spin moving left from position 0 exits through the source and re-enters at the maximum position.

The function `diffusion_mechanism()` simulates one spin population at a time. It uses the following inputs:
- spin_type = -1 for alpha spins,
- spin_type = 1 for beta spins,
- the configuration object,
- the state matrix,
- the EMCHA parameter,
- the voltage step dV.

The hopping probabilities are spin-dependent. eMCHA contribution is coded as follows:

```python
emcha_contribution = emcha_effect * spin_type * helix_twisting * tanh(dV)
```
This contribution modifies the right-moving probability through a voltage-, temperature-, spin-, and chirality-dependent expression. The left-moving probability is then defined as:

```python
left_probability = D - right_probability
```

where D = 0.5 in the current stable implementation. The actual hopping threshold uses `config.diff_coefficient`, allowing the user to control the total hopping probability scale.

For each spin and time step, a random number determines whether the spin moves right, moves left, or stays in place. The function records the number of drained and sourced spins at every step.

At the end of the simulation, the code computes:
- total drained spins,
- total sourced spins,
- average right/left hopping probabilities,
- current difference between drain and source.

The current conversion uses a fixed time step of 1e-17 s and converts electron counts per second into amperes using the constant:
```python
amperec = 6.241509074e18   
```

The wrapper function `apply_diffusion_mechanism()` runs the diffusion model separately for alpha and beta spin populations, then combines their summary statistics into a single DataFrame.

### `plotting.py`

plotting.py contains all visualization routines. It uses Matplotlib, Seaborn, Pillow, and ImageIO.

The main function is:

```{python}
plotting(config, df_alpha, df_beta, I_total, voltage_vector, qEMCHA)
```

Its behavior depends on `config.simulation_type`.

For `simulation_type == 0`, it plots cumulative drained and sourced electrons for alpha and beta spins, and also saves an animated GIF of the drain evolution.
For `simulation_type == 1`, it saves temporary frames showing drained-electron evolution for different EMCHA values. These frames are later combined into a GIF by main.py.
For `simulation_type == 5`, it plots the voltage-dependent current response for alpha spins, beta spins, and the unpolarized average.
For `simulation_type == 6`, it generates a heatmap animation of spin polarization along the molecule, where the local polarization is computed from the difference between alpha and beta occupation counts.
For `simulation_type == 7`, it creates an animation showing the spatial distribution of alpha and beta spins over time.

All plots are saved inside:
```{python}
results_<simulation_name>/
```

### Simulation Outputs

Depending on the selected simulation type, the code can produce:

- PNG plots of drained and sourced spin populations.
- GIF animations of spin transport over time.
- Voltage-current sweep plots.
- Text files containing voltage sweep data.
- Spin-polarization heatmap animations.
- Spatial spin-distribution animations.
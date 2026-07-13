import argparse

import numpy as np


REQUIRED_CONFIG_KEYS = (
    "n_strands",
    "n_turns_step",
    "n_steps",
    "helix_twisting",
    "molecule_length",
    "positions",
    "n_cycles",
    "number_spins",
    "electron_charge",
    "emcha_effect",
    "spin_ratio",
    "type_of_pulse",
    "voltage_magnitude",
    "voltage_frequency",
    "Temperature",
    "relaxation_time",
    "alpha_init_position",
    "beta_init_position",
    "init_polarization_fraction",
    "magnetochiral_anisotropy",
    "starting_mode",
    "setting_seed",
    "simulation_type",
    "save_results",
    "Model_type",
    "diff_coefficient",
)


def read_config_file(filename):
    """
    Read a simple key-value configuration file.

    The project config uses YAML-like ``key: value`` lines with comments. This
    parser keeps that lightweight format while casting basic bool, int, and
    float values.
    """
    config = {}

    with open(filename, "r", encoding="utf-8") as file:
        for line in file:
            line = line.split("#")[0].strip()
            if not line or ":" not in line:
                continue

            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()

            if value.lower() in ("true", "false"):
                value = value.lower() == "true"
            else:
                try:
                    if "." in value or "e" in value.lower():
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass

            config[key] = value

    return config


def validate_config(config_data):
    missing = [key for key in REQUIRED_CONFIG_KEYS if key not in config_data]
    if missing:
        raise KeyError(f"Missing required configuration values: {', '.join(missing)}")

    for key in ("n_steps", "positions", "number_spins"):
        if int(config_data[key]) <= 0:
            raise ValueError(f"{key} must be a positive integer")

    diff_coefficient = float(config_data["diff_coefficient"])
    if diff_coefficient <= 0 or diff_coefficient > 1:
        raise ValueError("diff_coefficient must be greater than 0 and less than or equal to 1")

    if float(config_data["spin_ratio"]) < 0:
        raise ValueError("spin_ratio must be greater than or equal to 0")

    for key in ("alpha_init_position", "beta_init_position"):
        value = int(config_data[key])
        if value < 0 or value > int(config_data["positions"]):
            raise ValueError(f"{key} must be between 0 and positions")


def processing_data():
    """
    Read simulation parameters and derive the values used by the simulation.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config",
        nargs="?",
        default="user_configurations.yaml",
        help="Path to configuration file (default: user_configurations.yaml)",
    )
    args = parser.parse_args()

    config_data = read_config_file(args.config)
    validate_config(config_data)

    config = {}

    config["n_strands"] = int(config_data["n_strands"])
    config["n_turns_step"] = int(config_data["n_turns_step"])
    config["n_steps"] = int(config_data["n_steps"])
    config["helix_twisting"] = int(config_data["helix_twisting"])
    config["molecule_length"] = float(config_data["molecule_length"])
    config["positions"] = int(config_data["positions"])
    config["n_cycles"] = int(config_data["n_cycles"])
    config["number_spins"] = int(config_data["number_spins"])
    config["electron_charge"] = float(config_data["electron_charge"])
    config["emcha_effect"] = float(config_data["emcha_effect"])
    config["spin_ratio"] = float(config_data["spin_ratio"])

    config["type_of_pulse"] = int(config_data["type_of_pulse"])
    config["voltage_magnitude"] = float(config_data["voltage_magnitude"])
    config["voltage_frequency"] = float(config_data["voltage_frequency"])

    config["Temperature"] = float(config_data["Temperature"])
    config["relaxation_time"] = float(config_data["relaxation_time"])
    config["alpha_init_position"] = int(config_data["alpha_init_position"])
    config["beta_init_position"] = int(config_data["beta_init_position"])
    config["init_polarization_fraction"] = float(config_data["init_polarization_fraction"])
    config["magnetochiral_anisotropy"] = float(config_data["magnetochiral_anisotropy"])

    config["starting_mode"] = float(config_data["starting_mode"]) / 100.0
    config["setting_seed"] = int(config_data["setting_seed"])
    config["simulation_type"] = int(config_data["simulation_type"])
    config["save_results"] = int(config_data["save_results"])
    config["Model_type"] = int(config_data["Model_type"])
    config["diff_coefficient"] = float(config_data["diff_coefficient"])

    config["k"] = 11604.525
    config["dx"] = config["molecule_length"] / config["positions"]
    config["dt"] = 2e-10 / (6 * config["diff_coefficient"] * config["positions"])
    config["dV"] = config["voltage_magnitude"] / config["positions"]
    config["Area"] = config["molecule_length"] * config["positions"] * config["number_spins"]
    config["total_strands"] = config["n_strands"] * 2
    config["relaxation_rate"] = (
        1 / config["relaxation_time"] if config["relaxation_time"] != 0 else np.nan
    )
    config["alpha_spins"] = int(config["number_spins"] / (config["spin_ratio"] + 1))
    config["beta_spins"] = config["number_spins"] - config["alpha_spins"]

    return config

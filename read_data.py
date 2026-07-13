import argparse


REQUIRED_CONFIG_KEYS = (
    "simulation_name",
    "n_steps",
    "helix_twisting",
    "positions",
    "number_spins",
    "emcha_effect",
    "spin_ratio",
    "voltage_magnitude",
    "Temperature",
    "alpha_init_position",
    "beta_init_position",
    "setting_seed",
    "simulation_type",
    "save_results",
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
    
    if float(config_data["Temperature"]) <= 0:
        raise ValueError("Temperature must be greater than 0")

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
    config["simulation_name"] = config_data["simulation_name"]
    config["n_steps"] = int(config_data["n_steps"])
    config["helix_twisting"] = int(config_data["helix_twisting"])
    config["positions"] = int(config_data["positions"])
    config["number_spins"] = int(config_data["number_spins"])
    config["emcha_effect"] = float(config_data["emcha_effect"])
    config["spin_ratio"] = float(config_data["spin_ratio"])

    config["voltage_magnitude"] = float(config_data["voltage_magnitude"])

    config["Temperature"] = float(config_data["Temperature"])
    config["alpha_init_position"] = int(config_data["alpha_init_position"])
    config["beta_init_position"] = int(config_data["beta_init_position"])

    config["setting_seed"] = int(config_data["setting_seed"])
    config["simulation_type"] = int(config_data["simulation_type"])
    config["save_results"] = int(config_data["save_results"])
    config["diff_coefficient"] = float(config_data["diff_coefficient"])

    config["dV"] = config["voltage_magnitude"] / config["positions"]
    config["alpha_spins"] = int(config["number_spins"] / (config["spin_ratio"] + 1))
    config["beta_spins"] = config["number_spins"] - config["alpha_spins"]

    return config

import os

import imageio
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from PIL import Image


sns.set_theme(style="whitegrid")

beta_spins_color = "#b7b1f2"
alpha_spins_color = "#e78895"

font = {"family": "Serif", "weight": "bold", "size": 15}
plt.rc("font", **font)
plt.rc("legend", fontsize=10)
plt.rcParams["lines.linewidth"] = 2.5
plt.rcParams["lines.markersize"] = 6
plt.rcParams["figure.figsize"] = [6, 4]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.dpi"] = 150


def _results_dir(config, *parts):
    name = config.simulation_name if hasattr(config, "simulation_name") else "results"
    path = os.path.join("results_" + name, *parts)
    os.makedirs(path, exist_ok=True)
    return path


def _cumulative(df, column):
    return df[column].cumsum().to_numpy()


def plotting(config, df_alpha, df_beta, I_total, voltage_vector, qCISS):
    results_dir = _results_dir(config)

    if config.simulation_type == 0:
        x = np.arange(config.n_steps)

        alpha_drain_history = _cumulative(df_alpha, "Drained alpha spins")
        beta_drain_history = _cumulative(df_beta, "Drained beta spins")
        alpha_source_history = _cumulative(df_alpha, "Sourced alpha spins")
        beta_source_history = _cumulative(df_beta, "Sourced beta spins")

        max_value_d = max(alpha_drain_history.max(), beta_drain_history.max(), 1)
        max_value_s = max(alpha_source_history.max(), beta_source_history.max(), 1)

        plt.figure()
        sns.lineplot(x=x, y=alpha_drain_history, color=alpha_spins_color, alpha=0.7, label=r"$\alpha$ spins")
        sns.lineplot(x=x, y=beta_drain_history, color=beta_spins_color, alpha=0.7, label=r"$\beta$ spins")
        plt.title(
            "%d Spins Drain Evolution at %d K, D %.2f, %.2f V & Qciss = %.1f"
            % (
                config.number_spins,
                config.Temperature,
                config.diff_coefficient,
                config.voltage_magnitude,
                config.ciss_effect,
            )
        )
        plt.xlabel("Steps")
        plt.xlim([0, config.n_steps])
        plt.ylim([0, max_value_d])
        plt.ylabel("Drained Electrons")
        plt.legend(title="Spin Type", loc="upper left")
        plt.savefig(os.path.join(results_dir, f"spin_drain_evolution_{config.voltage_magnitude}.png"))
        plt.show()

        plt.figure()
        sns.lineplot(x=x, y=alpha_source_history, color=alpha_spins_color, label=r"$\alpha$ spins")
        sns.lineplot(x=x, y=beta_source_history, color=beta_spins_color, label=r"$\beta$ spins")
        plt.title(
            "%d Spins Source Evolution at %d K, D %.2f, %.2f V & Qciss = %.1f"
            % (
                config.number_spins,
                config.Temperature,
                config.diff_coefficient,
                config.voltage_magnitude,
                config.ciss_effect,
            )
        )
        plt.xlabel("Steps")
        plt.xlim([0, config.n_steps])
        plt.ylim([0, max_value_s])
        plt.ylabel("Sourced Electrons")
        plt.legend(title="Spin Type", loc="upper left")
        plt.savefig(os.path.join(results_dir, f"spin_source_evolution_{config.voltage_magnitude}.png"))
        plt.show()

        fig, ax = plt.subplots()
        ax.set_title(
            f"{config.number_spins} Spins Drain Evolution at {config.Temperature} K, "
            f"D {config.diff_coefficient:.2f}, Qciss = {config.ciss_effect} & "
            f"{config.voltage_magnitude} V"
        )
        ax.set_xlabel("Steps")
        ax.set_ylabel("Drained Electrons")
        ax.set_xlim([0, config.n_steps])
        ax.set_ylim([0, max_value_d])

        line_alpha, = ax.plot([], [], label=r"$\alpha$ spins", color=alpha_spins_color)
        line_beta, = ax.plot([], [], label=r"$\beta$ spins", color=beta_spins_color)
        ax.legend(title="Spin Type", loc="upper left")

        animation_frames = np.linspace(
            0,
            config.n_steps - 1,
            min(config.n_steps, 60),
            dtype=int,
        )

        def update(frame):
            line_alpha.set_data(x[: frame + 1], alpha_drain_history[: frame + 1])
            line_beta.set_data(x[: frame + 1], beta_drain_history[: frame + 1])
            return line_alpha, line_beta

        ani = animation.FuncAnimation(fig, update, frames=animation_frames, interval=20, blit=True)
        ani.save(os.path.join(results_dir, f"spin_drain_evolution_{config.voltage_magnitude}.gif"), writer="pillow")
        plt.show()

    elif config.simulation_type == 1:
        temp_dir = _results_dir(config, "temp_dir")
        x = np.arange(config.n_steps)
        alpha_drain_history = _cumulative(df_alpha, "Drained alpha spins")
        beta_drain_history = _cumulative(df_beta, "Drained beta spins")
        max_value_d = max(alpha_drain_history.max(), beta_drain_history.max(), 1)

        plt.figure(figsize=(6, 4))
        sns.lineplot(x=x, y=alpha_drain_history, color=alpha_spins_color, alpha=0.7, label=r"$\alpha$ spins")
        sns.lineplot(x=x, y=beta_drain_history, color=beta_spins_color, alpha=0.7, label=r"$\beta$ spins")
        plt.title(f"Drain Electrons Evolution (q_CISS = {qCISS:.2f})")
        plt.xlabel("Time Steps")
        plt.ylabel("Drained Electrons")
        plt.xlim([0, config.n_steps])
        plt.ylim([0, max_value_d])
        plt.legend(title="Spin Type", loc="upper left")
        plt.savefig(os.path.join(temp_dir, f"drain_electrons_qciss_{qCISS:.2f}.png"))
        plt.close()

    elif config.simulation_type == 5:
        x = voltage_vector
        sns.lineplot(x=x, y=df_alpha, color=alpha_spins_color, label=r"$\alpha$ spins")
        sns.lineplot(x=x, y=df_beta, color=beta_spins_color, label=r"$\beta$ spins")
        sns.lineplot(x=x, y=I_total, color="black", label=r"Unpolarized")
        plt.xlabel("Voltage magnitude (V)")
        plt.xlim([x[0], x[-1]])
        plt.ylabel("Intensity (A)")
        plt.savefig(os.path.join(results_dir, "IV_alpha_beta.png"), dpi=150, bbox_inches="tight")
        plt.show()

    elif config.simulation_type == 6:
        frames_folder = _results_dir(config, "heatmap_frames")
        frame_files = []
        for t in range(config.n_steps):
            alpha_counts = np.bincount(df_alpha[:, t], minlength=config.positions + 1)
            beta_counts = np.bincount(df_beta[:, t], minlength=config.positions + 1)
            total_counts = alpha_counts + beta_counts
            polarization = (alpha_counts - beta_counts) / np.where(total_counts > 0, total_counts, 1)
            polarization_2d = np.nan_to_num(polarization, nan=0).reshape((1, -1))

            plt.figure(figsize=(8, 2))
            sns.heatmap(polarization_2d, cmap="coolwarm", cbar=True, vmin=-1, vmax=1, xticklabels=False, yticklabels=False)
            plt.title(f"Spin Polarization at t = {t}")
            frame_path = os.path.join(frames_folder, f"frame_{t:03d}.png")
            plt.savefig(frame_path, dpi=300)
            plt.close()
            frame_files.append(frame_path)

        images = [Image.open(f) for f in frame_files]
        images[0].save(os.path.join(results_dir, "heatmap.gif"), save_all=True, append_images=images[1:], duration=100, loop=0)

    elif config.simulation_type == 7:
        frames = []
        temp_dir = _results_dir(config, "temp_frames")
        for t in range(config.n_steps):
            plt.figure(figsize=(10, 5))
            alpha_counts = np.bincount(df_alpha[:, t], minlength=config.positions + 1)
            beta_counts = np.bincount(df_beta[:, t], minlength=config.positions + 1)

            plt.subplot(2, 1, 1)
            markerline, stemlines, _ = plt.stem(np.arange(len(alpha_counts)), alpha_counts, basefmt=" ")
            plt.setp(markerline, color=alpha_spins_color, marker="o")
            plt.setp(stemlines, color=alpha_spins_color)
            plt.xlabel("Position")
            plt.ylabel("Number of alpha spins")
            plt.title(f"Distribution of alpha spins. Nstep {t}")

            plt.subplot(2, 1, 2)
            markerline, stemlines, _ = plt.stem(np.arange(len(beta_counts)), beta_counts, basefmt=" ")
            plt.setp(markerline, color=beta_spins_color, marker="o")
            plt.setp(stemlines, color=beta_spins_color)
            plt.xlabel("Position")
            plt.ylabel("Number of beta spins")
            plt.title(f"Distribution of beta spins. Nstep {t}")

            plt.tight_layout()
            frame_path = os.path.join(temp_dir, f"frame_{t}.png")
            plt.savefig(frame_path)
            plt.close()
            frames.append(frame_path)

        with imageio.get_writer(os.path.join(results_dir, "spins_evolution.gif"), mode="I", duration=0.1) as writer:
            for frame in frames:
                writer.append_data(imageio.imread(frame))

        for frame in frames:
            os.remove(frame)
        os.rmdir(temp_dir)

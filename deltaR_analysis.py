from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


MU_B_OVER_K_B = 0.67171388  # K/T, value used in the spreadsheet


@dataclass(frozen=True)
class ZeemanInputs:
    temperature: float = 8.0
    spin: float = 0.5
    lande_g: float = 2.00232
    field_start: float = 0.0
    field_stop: float = 1.0
    field_step: float = 0.1


@dataclass(frozen=True)
class CurrentBlock:
    name: str
    voltage: float
    positive_alpha_current: float
    positive_beta_current: float
    negative_alpha_current: float
    negative_beta_current: float
    ground_state: int = 0
    reference_resistance: float | None = None
    gamma_denominator: float | None = None


DEFAULT_BLOCKS = (
    CurrentBlock(
        name="V=0.15",
        voltage=0.15,
        positive_alpha_current=0.2695899,
        positive_beta_current=0.4168643,
        negative_alpha_current=-0.4166726,
        negative_beta_current=-0.2694855,
    ),
    CurrentBlock(
        name="V=0.0765",
        voltage=0.0765,
        positive_alpha_current=0.1669507,
        positive_beta_current=0.21643,
        negative_alpha_current=-0.216258,
        negative_beta_current=-0.1670535,
        reference_resistance=0.4,
        gamma_denominator=0.192 * 2,
    ),
)


def inclusive_range(start: float, stop: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("field_step must be greater than 0")
    n_steps = int(round((stop - start) / step))
    return [round(start + i * step, 12) for i in range(n_steps + 1)]


def zeeman_energy(field: float, inputs: ZeemanInputs) -> float:
    return inputs.lande_g * MU_B_OVER_K_B * field * inputs.spin


def delta_zeeman_energy(field: float, inputs: ZeemanInputs) -> float:
    return 2.0 * zeeman_energy(field, inputs)


def ground_population(delta_zeeman: float, temperature: float) -> float:
    if temperature <= 0:
        raise ValueError("temperature must be greater than 0")
    return 1.0 / (1.0 + math.exp(-delta_zeeman / temperature))


def excited_population(delta_zeeman: float, temperature: float) -> float:
    return 1.0 - ground_population(delta_zeeman, temperature)


def alpha_population(delta_zeeman: float, inputs: ZeemanInputs, ground_state: int) -> float:
    ground_pop = ground_population(delta_zeeman, inputs.temperature)
    excited_pop = 1.0 - ground_pop
    return ground_pop if ground_state == 1 else excited_pop


def resistance(
    voltage: float,
    alpha_pop: float,
    alpha_current: float,
    beta_current: float,
) -> tuple[float, float, float]:
    beta_pop = 1.0 - alpha_pop
    total_current = alpha_pop * alpha_current + beta_pop * beta_current
    if total_current == 0:
        raise ZeroDivisionError("total current is zero, so resistance is undefined")
    return beta_pop, total_current, voltage / total_current


def analyze_block(block: CurrentBlock, inputs: ZeemanInputs) -> list[dict[str, float | str]]:
    fields = inclusive_range(inputs.field_start, inputs.field_stop, inputs.field_step)
    rows = []

    for field in fields:
        zeeman = zeeman_energy(field, inputs)
        delta_zeeman = 2.0 * zeeman
        ground_pop = ground_population(delta_zeeman, inputs.temperature)
        excited_pop = 1.0 - ground_pop
        alpha_pop = ground_pop if block.ground_state == 1 else excited_pop

        beta_pop_pos, total_current_pos, resistance_pos = resistance(
            block.voltage,
            alpha_pop,
            block.positive_alpha_current,
            block.positive_beta_current,
        )
        beta_pop_neg, total_current_neg, resistance_neg = resistance(
            -block.voltage,
            alpha_pop,
            block.negative_alpha_current,
            block.negative_beta_current,
        )

        rows.append(
            {
                "block": block.name,
                "field_T": field,
                "zeeman_E_K": zeeman,
                "delta_zeeman_E_K": delta_zeeman,
                "ground_pop": ground_pop,
                "excited_pop": excited_pop,
                "alpha_pop": alpha_pop,
                "beta_pop": beta_pop_pos,
                "positive_total_current": total_current_pos,
                "positive_R": resistance_pos,
                "negative_total_current": total_current_neg,
                "negative_R": resistance_neg,
                "deltaR": resistance_pos - resistance_neg,
            }
        )

    baseline_delta_r = float(rows[0]["deltaR"])
    for row in rows:
        delta_r_prime = float(row["deltaR"]) - baseline_delta_r
        row["deltaR_prime"] = delta_r_prime

        if block.reference_resistance is not None:
            row["deltaR_over_R"] = delta_r_prime / block.reference_resistance
        else:
            row["deltaR_over_R"] = ""

        if block.reference_resistance is not None and block.gamma_denominator is not None:
            row["gamma"] = (delta_r_prime / block.reference_resistance) / block.gamma_denominator
        else:
            row["gamma"] = ""

        # Keep this explicit in case future blocks use different spin-state assumptions.
        row["ground_state_flag"] = block.ground_state
        row["positive_alpha_current"] = block.positive_alpha_current
        row["positive_beta_current"] = block.positive_beta_current
        row["negative_alpha_current"] = block.negative_alpha_current
        row["negative_beta_current"] = block.negative_beta_current

    return rows


def run_analysis(inputs: ZeemanInputs = ZeemanInputs()) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for block in DEFAULT_BLOCKS:
        rows.extend(analyze_block(block, inputs))
    return rows


def write_csv(rows: list[dict[str, float | str]], output_path: Path) -> None:
    if not rows:
        raise ValueError("no rows to write")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the deltaR / gamma equations from the currentEMCHA workbook."
    )
    parser.add_argument("--temperature", type=float, default=8.0)
    parser.add_argument("--spin", type=float, default=0.5)
    parser.add_argument("--lande-g", type=float, default=2.00232)
    parser.add_argument("--field-start", type=float, default=0.0)
    parser.add_argument("--field-stop", type=float, default=1.0)
    parser.add_argument("--field-step", type=float, default=0.1)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("deltaR_analysis_results.csv"),
        help="CSV output path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = ZeemanInputs(
        temperature=args.temperature,
        spin=args.spin,
        lande_g=args.lande_g,
        field_start=args.field_start,
        field_stop=args.field_stop,
        field_step=args.field_step,
    )
    rows = run_analysis(inputs)
    write_csv(rows, args.output)
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()

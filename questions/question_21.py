"""Question 21: sluice flow, normal depth, and force on hydraulic-jump blocks."""

if __package__ is None or __package__ == "":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from open_channel import (
    RectangularSection,
    normal_depth,
    rectangular_discharge_from_depths,
    rectangular_velocity,
    sluice_gate_force,
)


def solve(
    width: float = 3.0,
    upstream_sluice_depth: float = 2.0,
    downstream_sluice_depth: float = 0.3,
    roughness: float = 0.014,
    bed_slope: float = 1.0 / 1000.0,
) -> dict[str, float]:
    discharge = rectangular_discharge_from_depths(
        width,
        upstream_sluice_depth,
        downstream_sluice_depth,
    )
    section = RectangularSection(width)
    normal = normal_depth(section, discharge, roughness, bed_slope)
    block_force = sluice_gate_force(
        discharge,
        width,
        downstream_sluice_depth,
        normal,
    )

    return {
        "discharge": discharge,
        "normal_depth": normal,
        "velocity_before_jump": rectangular_velocity(
            discharge, width, downstream_sluice_depth
        ),
        "velocity_after_jump": rectangular_velocity(discharge, width, normal),
        "force_on_blocks": block_force,
    }


def _input_float(prompt: str, default: float) -> float:
    value = input(f"{prompt} [{default}]: ").strip()
    return default if value == "" else float(value)


def main() -> None:
    print("Question 21: sluice flow, normal depth, and force on hydraulic-jump blocks")
    print("Press Enter to use the value stated in the question.\n")

    width = _input_float("Channel width (m)", 3.0)
    upstream_sluice_depth = _input_float("Upstream sluice depth (m)", 2.0)
    downstream_sluice_depth = _input_float("Downstream sluice depth (m)", 0.3)
    roughness = _input_float("Manning roughness n", 0.014)
    bed_slope = _input_float("Bed slope S0", 1.0 / 1000.0)

    results = solve(
        width=width,
        upstream_sluice_depth=upstream_sluice_depth,
        downstream_sluice_depth=downstream_sluice_depth,
        roughness=roughness,
        bed_slope=bed_slope,
    )

    print("\n--- Results ---")
    for name, value in results.items():
        unit = (
            "m^3/s"
            if name == "discharge"
            else "N" if "force" in name else "m" if "depth" in name else "m/s"
        )
        print(f"{name}: {value:.6f} {unit}")


if __name__ == "__main__":
    main()

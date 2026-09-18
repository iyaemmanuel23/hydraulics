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


def solve() -> dict[str, float]:
    width = 3.0
    upstream_sluice_depth = 2.0
    downstream_sluice_depth = 0.3
    roughness = 0.014
    bed_slope = 1.0 / 1000.0

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


if __name__ == "__main__":
    for name, value in solve().items():
        unit = (
            "m^3/s"
            if name == "discharge"
            else "N" if "force" in name else "m" if "depth" in name else "m/s"
        )
        print(f"{name}: {value:.6f} {unit}")

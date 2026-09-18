"""Question 15: flow controlled by an undershot sluice gate."""

if __package__ is None or __package__ == "":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from open_channel import (
    rectangular_depth_from_specific_energy,
    rectangular_froude,
    rectangular_specific_energy,
    rectangular_velocity,
    sluice_gate_force,
)


def solve() -> dict[str, float]:
    width = 0.8
    discharge = 0.9
    downstream_depth = 0.25

    total_head = rectangular_specific_energy(discharge, width, downstream_depth)
    upstream_depth = rectangular_depth_from_specific_energy(
        discharge,
        width,
        total_head,
        regime="subcritical",
    )
    velocity_upstream = rectangular_velocity(discharge, width, upstream_depth)
    velocity_downstream = rectangular_velocity(discharge, width, downstream_depth)

    return {
        "total_head": total_head,
        "upstream_depth": upstream_depth,
        "upstream_velocity": velocity_upstream,
        "downstream_velocity": velocity_downstream,
        "upstream_froude": rectangular_froude(discharge, width, upstream_depth),
        "downstream_froude": rectangular_froude(discharge, width, downstream_depth),
        "gate_force": sluice_gate_force(
            discharge, width, upstream_depth, downstream_depth
        ),
    }


if __name__ == "__main__":
    for name, value in solve().items():
        unit = (
            "N"
            if name == "gate_force"
            else (
                "m"
                if "depth" in name or name == "total_head"
                else "m/s" if "velocity" in name else ""
            )
        )
        print(f"{name}: {value:.6f} {unit}".rstrip())

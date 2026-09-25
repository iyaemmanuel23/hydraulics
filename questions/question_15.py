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


def solve(
    width: float = 0.8,
    discharge: float = 0.9,
    downstream_depth: float = 0.25,
) -> dict[str, float]:
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


def _input_float(prompt: str, default: float) -> float:
    value = input(f"{prompt} [{default}]: ").strip()
    return default if value == "" else float(value)


def main() -> None:
    print("Question 15: flow controlled by an undershot sluice gate")
    print("Press Enter to use the value stated in the question.\n")

    width = _input_float("Channel width (m)", 0.8)
    discharge = _input_float("Discharge Q (m^3/s)", 0.9)
    downstream_depth = _input_float("Downstream depth (m)", 0.25)

    results = solve(
        width=width,
        discharge=discharge,
        downstream_depth=downstream_depth,
    )

    print("\n--- Results ---")
    for name, value in results.items():
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


if __name__ == "__main__":
    main()

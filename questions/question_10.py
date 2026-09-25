"""Question 10: broad-crested weir and downstream hydraulic jump."""

if __package__ is None or __package__ == "":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from open_channel import (
    broad_crested_weir_head,
    critical_depth,
    rectangular_conjugate_depth,
    supercritical_reach_length,
    wide_rectangular_normal_depth,
    RectangularSection,
)


def solve(
    discharge_per_width: float = 0.5,
    bed_slope: float = 2.0e-5,
    roughness: float = 0.01,
    weir_height: float = 0.7,
) -> dict[str, float]:
    section = RectangularSection(width=1.0)

    normal = wide_rectangular_normal_depth(
        discharge_per_width, roughness, bed_slope
    )
    weir_critical = critical_depth(section, discharge_per_width)
    head_over_weir = broad_crested_weir_head(discharge_per_width)
    jump_upstream = rectangular_conjugate_depth(discharge_per_width, normal)
    jump_distance = supercritical_reach_length(
        discharge_per_width,
        upstream_depth=weir_critical,
        downstream_depth=jump_upstream,
        roughness=roughness,
        slope=bed_slope,
    )

    return {
        "normal_depth": normal,
        "depth_over_weir": head_over_weir,
        "critical_depth_at_crest": weir_critical,
        "downstream_depth_of_jump": normal,
        "upstream_depth_of_jump": jump_upstream,
        "weir_height": weir_height,
        "jump_distance_downstream_of_weir": jump_distance,
    }


def _input_float(prompt: str, default: float) -> float:
    value = input(f"{prompt} [{default}]: ").strip()
    return default if value == "" else float(value)


def main() -> None:
    print("Question 10: broad-crested weir and downstream hydraulic jump")
    print("Press Enter to use the value stated in the question.\n")

    discharge_per_width = _input_float(
        "Discharge per unit width q (m^2/s)", 0.5
    )
    bed_slope = _input_float("Bed slope S0", 2.0e-5)
    roughness = _input_float("Manning roughness n", 0.01)
    weir_height = _input_float("Weir height (m)", 0.7)

    results = solve(
        discharge_per_width=discharge_per_width,
        bed_slope=bed_slope,
        roughness=roughness,
        weir_height=weir_height,
    )

    print("\n--- Results ---")
    for name, value in results.items():
        print(f"{name}: {value:.6f} m")


if __name__ == "__main__":
    main()

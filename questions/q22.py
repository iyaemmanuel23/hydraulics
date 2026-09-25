"""Solution for Open-Channel Flow Question 22.

A rectangular channel expands abruptly from one width to another and a
hydraulic jump occurs through the expansion.  Run this file directly to enter
problem data in the terminal.
"""

from pathlib import Path
import sys

# Allow ``python questions/question_22.py`` to import the project-level
# open_channel.py when this file is launched from the project root.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from open_channel import (
    GRAVITY,
    abrupt_expansion_jump_depth,
    rectangular_froude,
)


def solve_question_22(
    discharge: float,
    upstream_width: float,
    downstream_width: float,
    upstream_depth: float,
    gravity: float = GRAVITY,
) -> dict[str, float | str]:
    """Solve Question 22 using the momentum equation across the expansion.

    Assumptions:
    - steady, incompressible flow;
    - rectangular channel sections;
    - hydrostatic pressure distributions at the two sections;
    - the expansion is sufficiently short that bed friction and weight effects
      across the transition are negligible;
    - momentum is therefore conserved through the abrupt expansion;
    - the physically relevant downstream solution is the subcritical root.
    """
    if discharge <= 0.0:
        raise ValueError("Discharge must be positive.")
    if upstream_width <= 0.0 or downstream_width <= 0.0:
        raise ValueError("Channel widths must be positive.")
    if upstream_depth <= 0.0:
        raise ValueError("Upstream depth must be positive.")
    if downstream_width <= upstream_width:
        raise ValueError("Downstream width must be greater than upstream width.")

    downstream_depth = abrupt_expansion_jump_depth(
        discharge=discharge,
        upstream_width=upstream_width,
        downstream_width=downstream_width,
        upstream_depth=upstream_depth,
        gravity=gravity,
    )

    upstream_froude = rectangular_froude(
        discharge, upstream_width, upstream_depth, gravity
    )
    downstream_froude = rectangular_froude(
        discharge, downstream_width, downstream_depth, gravity
    )

    return {
        "downstream_depth_m": downstream_depth,
        "upstream_froude": upstream_froude,
        "downstream_froude": downstream_froude,
        "upstream_flow_regime": (
            "supercritical" if upstream_froude > 1.0 else "subcritical"
        ),
        "downstream_flow_regime": (
            "supercritical" if downstream_froude > 1.0 else "subcritical"
        ),
    }


def main() -> None:
    print("Question 22 — Abrupt Expansion / Hydraulic Jump")
    print("Press Enter for the stated value in brackets.\n")

    def value(prompt: str, default: float) -> float:
        raw = input(f"{prompt} [{default}]: ").strip()
        return default if not raw else float(raw)

    discharge = value("Discharge Q (m^3/s)", 10.0)
    upstream_width = value("Upstream width b1 (m)", 4.0)
    downstream_width = value("Downstream width b2 (m)", 8.0)
    upstream_depth = value("Upstream depth y1 (m)", 0.5)

    result = solve_question_22(
        discharge=discharge,
        upstream_width=upstream_width,
        downstream_width=downstream_width,
        upstream_depth=upstream_depth,
    )

    print("\n--- Results ---")
    print(
        f"Upstream Froude number = {result['upstream_froude']:.4f} "
        f"({result['upstream_flow_regime']})"
    )
    print(f"Downstream depth y2 = {result['downstream_depth_m']:.4f} m")
    print(
        f"Downstream Froude number = {result['downstream_froude']:.4f} "
        f"({result['downstream_flow_regime']})"
    )
    print(
        "The expansion changes the flow from supercritical upstream to "
        "subcritical downstream, consistent with a hydraulic jump through "
        "the abrupt expansion."
    )


if __name__ == "__main__":
    main()

"""Solution for Open-Channel Flow Question 34.

Run this file directly to enter the problem data in the terminal.
"""

if __package__ is None or __package__ == "":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))




from open_channel import (
    GRAVITY,
    broad_crested_weir_transition,
    critical_depth,
    rectangular_conjugate_depth,
    rectangular_froude,
    rectangular_specific_energy,
    manning_discharge,
    gvf_one_step_jump_distance,
    RectangularSection,
)


def solve_question_34(
    width: float,
    slope: float,
    roughness: float,
    undisturbed_depth: float,
    weir_height: float,
    gravity: float = GRAVITY,
) -> dict[str, float | str]:
    """Solve all numerical parts of Question 34.

    The broad-crested weir is treated as a frictionless control with a
    discharge coefficient of 1.0, as implied by the question. The hydraulic
    jump is estimated one GVF step downstream of the weir.
    """
    if width <= 0 or slope <= 0 or roughness <= 0 or undisturbed_depth <= 0:
        raise ValueError("Width, slope, roughness, and depth must be positive.")
    if weir_height <= 0:
        raise ValueError("Weir height must be positive.")

    # (a) Undisturbed discharge from Manning's equation.
    area = width * undisturbed_depth
    perimeter = width + 2.0 * undisturbed_depth
    radius = area / perimeter
    discharge = manning_discharge(
        section=RectangularSection(width),
        depth=undisturbed_depth,
        roughness=roughness,
        slope=slope,
    )
    discharge_per_width = discharge / width

    # (b) Critical depth and slope classification.
    yc = critical_depth(
        RectangularSection(width), discharge, gravity
    )
    froude_normal = rectangular_froude(
        discharge, width, undisturbed_depth, gravity
    )
    channel_slope_type = "mild" if undisturbed_depth > yc else "steep"

    # (c) Weir control: upstream subcritical and downstream supercritical.
    weir = broad_crested_weir_transition(
        discharge=discharge,
        width=width,
        weir_height=weir_height,
        gravity=gravity,
    )

    # The jump occurs when the rising supercritical depth reaches the
    # conjugate depth of the normal subcritical flow.
    jump_supercritical_depth = rectangular_conjugate_depth(
        discharge_per_width, undisturbed_depth, gravity
    )
    jump_distance = gvf_one_step_jump_distance(
        discharge_per_width=discharge_per_width,
        initial_depth=weir["downstream_depth"],
        target_depth=jump_supercritical_depth,
        roughness=roughness,
        slope=slope,
        gravity=gravity,
    )

    return {
        "discharge_m3_s": discharge,
        "discharge_per_width_m2_s": discharge_per_width,
        "hydraulic_radius_m": radius,
        "critical_depth_m": yc,
        "normal_froude": froude_normal,
        "slope_type": channel_slope_type,
        "weir_critical_depth_m": weir["critical_depth"],
        "weir_upstream_depth_m": weir["upstream_depth"],
        "weir_downstream_depth_m": weir["downstream_depth"],
        "weir_upstream_froude": weir["upstream_froude"],
        "weir_downstream_froude": weir["downstream_froude"],
        "jump_supercritical_depth_m": jump_supercritical_depth,
        "jump_distance_downstream_m": jump_distance,
    }


def main() -> None:
    print("Question 34 — Open Channel Flow")
    print("Press Enter for the stated value in brackets.\n")

    def value(prompt: str, default: float) -> float:
        raw = input(f"{prompt} [{default}]: ").strip()
        return default if not raw else float(raw)

    width = value("Channel width b (m)", 7.0)
    slope = value("Channel slope S0", 0.005)
    roughness = value("Manning roughness n", 0.035)
    undisturbed_depth = value("Undisturbed depth y_n (m)", 1.6)
    weir_height = value("Broad-crested weir height (m)", 1.8)

    result = solve_question_34(
        width=width,
        slope=slope,
        roughness=roughness,
        undisturbed_depth=undisturbed_depth,
        weir_height=weir_height,
    )

    print("\n--- Results ---")
    print(f"(a) Flow rate Q = {result['discharge_m3_s']:.4f} m^3/s")
    print(f"    Unit discharge q = {result['discharge_per_width_m2_s']:.4f} m^2/s")
    print(f"(b) Critical depth yc = {result['critical_depth_m']:.4f} m")
    print(
        f"    Since yn = {undisturbed_depth:.4f} m > yc, the channel is "
        f"{result['slope_type']} (Fr = {result['normal_froude']:.4f})."
    )
    print(
        f"(c) Critical depth over the weir = "
        f"{result['weir_critical_depth_m']:.4f} m"
    )
    print(
        f"    Depth just upstream of weir = "
        f"{result['weir_upstream_depth_m']:.4f} m "
        f"(Fr = {result['weir_upstream_froude']:.4f})"
    )
    print(
        f"    Depth just downstream of weir = "
        f"{result['weir_downstream_depth_m']:.4f} m "
        f"(Fr = {result['weir_downstream_froude']:.4f})"
    )
    print(
        "    Therefore the weir changes the flow from subcritical upstream "
        "to supercritical downstream."
    )
    print(
        f"(d) Supercritical depth immediately downstream = "
        f"{result['weir_downstream_depth_m']:.4f} m"
    )
    print(
        f"    Conjugate depth for the normal depth {undisturbed_depth:.4f} m = "
        f"{result['jump_supercritical_depth_m']:.4f} m"
    )
    print(
        f"    One-step GVF estimate of hydraulic-jump position = "
        f"{result['jump_distance_downstream_m']:.2f} m downstream of the weir"
    )


if __name__ == "__main__":
    main()

"""Reusable open-channel flow calculations in SI units."""

from dataclasses import dataclass
from math import acos, sin, sqrt
from typing import Callable

GRAVITY = 9.81


class Section:
    """Interface required by the flow solvers."""

    def area(self, depth: float) -> float:
        raise NotImplementedError

    def wetted_perimeter(self, depth: float) -> float:
        raise NotImplementedError

    def top_width(self, depth: float) -> float:
        raise NotImplementedError


@dataclass(frozen=True)
class RectangularSection(Section):
    width: float

    def area(self, depth: float) -> float:
        return self.width * depth

    def wetted_perimeter(self, depth: float) -> float:
        return self.width + 2.0 * depth

    def top_width(self, depth: float) -> float:
        return self.width


@dataclass(frozen=True)
class TrapezoidalSection(Section):
    bottom_width: float
    side_slope: float

    def area(self, depth: float) -> float:
        return depth * (self.bottom_width + self.side_slope * depth)

    def wetted_perimeter(self, depth: float) -> float:
        return self.bottom_width + 2.0 * depth * sqrt(1.0 + self.side_slope**2)

    def top_width(self, depth: float) -> float:
        return self.bottom_width + 2.0 * self.side_slope * depth


@dataclass(frozen=True)
class CircularSection(Section):
    diameter: float

    def _angle(self, depth: float) -> float:
        radius = self.diameter / 2.0
        return 2.0 * acos((radius - depth) / radius)

    def area(self, depth: float) -> float:
        radius = self.diameter / 2.0
        angle = self._angle(depth)
        return radius**2 * (angle - sin(angle)) / 2.0

    def wetted_perimeter(self, depth: float) -> float:
        return self._angle(depth) * self.diameter / 2.0

    def top_width(self, depth: float) -> float:
        radius = self.diameter / 2.0
        return 2.0 * sqrt(max(0.0, radius**2 - (radius - depth) ** 2))


def _bisect(function: Callable[[float], float], lower: float, upper: float) -> float:
    lower_value = function(lower)
    upper_value = function(upper)
    if lower_value * upper_value > 0.0:
        raise ValueError("The root is not bracketed by the supplied bounds.")
    for _ in range(100):
        midpoint = (lower + upper) / 2.0
        midpoint_value = function(midpoint)
        if abs(midpoint_value) < 1e-11 or upper - lower < 1e-10:
            return midpoint
        if lower_value * midpoint_value <= 0.0:
            upper, upper_value = midpoint, midpoint_value
        else:
            lower, lower_value = midpoint, midpoint_value
    return (lower + upper) / 2.0


def hydraulic_radius(section: Section, depth: float) -> float:
    return section.area(depth) / section.wetted_perimeter(depth)


def manning_discharge(
    section: Section, depth: float, roughness: float, slope: float
) -> float:
    radius = hydraulic_radius(section, depth)
    return (1.0 / roughness) * section.area(depth) * radius ** (2.0 / 3.0) * sqrt(slope)


def normal_depth(
    section: Section,
    discharge: float,
    roughness: float,
    slope: float,
    upper_depth: float | None = None,
) -> float:
    """Solve Manning's equation for normal depth."""
    if discharge <= 0.0 or roughness <= 0.0 or slope <= 0.0:
        raise ValueError("Discharge, roughness, and slope must be positive.")
    upper = upper_depth or 1.0
    while manning_discharge(section, upper, roughness, slope) < discharge:
        upper *= 2.0
    return _bisect(
        lambda depth: manning_discharge(section, depth, roughness, slope) - discharge,
        1e-9,
        upper,
    )


def wide_rectangular_normal_depth(
    discharge_per_width: float,
    roughness: float,
    slope: float,
) -> float:
    """Solve Manning's equation for a wide rectangular channel per unit width."""
    if discharge_per_width <= 0.0 or roughness <= 0.0 or slope <= 0.0:
        raise ValueError("Discharge, roughness, and slope must be positive.")
    return (discharge_per_width * roughness / sqrt(slope)) ** (3.0 / 5.0)


def critical_depth(
    section: Section, discharge: float, gravity: float = GRAVITY
) -> float:
    """Solve the critical-flow condition Fr = 1."""
    if discharge <= 0.0:
        raise ValueError("Discharge must be positive.")
    upper = 1.0

    def critical_function(depth: float) -> float:
        area = section.area(depth)
        return gravity * area**3 - discharge**2 * section.top_width(depth)

    while critical_function(upper) < 0.0:
        upper *= 2.0
    return _bisect(critical_function, 1e-9, upper)


def broad_crested_weir_head(
    discharge_per_width: float,
    discharge_coefficient: float = 1.0,
    gravity: float = GRAVITY,
) -> float:
    """Return the upstream energy head above a broad-crested weir crest."""
    if discharge_per_width <= 0.0 or discharge_coefficient <= 0.0:
        raise ValueError("Discharge and discharge coefficient must be positive.")
    critical = (discharge_per_width**2 / gravity) ** (1.0 / 3.0)
    return 1.5 * critical / discharge_coefficient ** (2.0 / 3.0)


def rectangular_conjugate_depth(
    discharge_per_width: float,
    downstream_depth: float,
    gravity: float = GRAVITY,
) -> float:
    """Return the supercritical conjugate depth for a rectangular jump."""
    froude_downstream = discharge_per_width / (
        downstream_depth * sqrt(gravity * downstream_depth)
    )
    return 0.5 * downstream_depth * (sqrt(1.0 + 8.0 * froude_downstream**2) - 1.0)


def rectangular_velocity(discharge: float, width: float, depth: float) -> float:
    """Return mean velocity in a rectangular channel."""
    if discharge <= 0.0 or width <= 0.0 or depth <= 0.0:
        raise ValueError("Discharge, width, and depth must be positive.")
    return discharge / (width * depth)


def rectangular_specific_energy(
    discharge: float,
    width: float,
    depth: float,
    gravity: float = GRAVITY,
) -> float:
    """Return specific energy relative to the channel bed."""
    velocity = rectangular_velocity(discharge, width, depth)
    return depth + velocity**2 / (2.0 * gravity)


def rectangular_discharge_from_depths(
    width: float,
    depth_a: float,
    depth_b: float,
    gravity: float = GRAVITY,
) -> float:
    """Return discharge from equal specific energy at two rectangular depths."""
    if width <= 0.0 or depth_a <= 0.0 or depth_b <= 0.0:
        raise ValueError("Width and depths must be positive.")
    if depth_a == depth_b:
        raise ValueError("The two depths must be different.")
    numerator = 2.0 * gravity * width**2 * abs(depth_a - depth_b)
    denominator = abs(1.0 / depth_b**2 - 1.0 / depth_a**2)
    return sqrt(numerator / denominator)


def rectangular_depth_from_specific_energy(
    discharge: float,
    width: float,
    specific_energy: float,
    regime: str = "subcritical",
    gravity: float = GRAVITY,
) -> float:
    """Solve rectangular-channel depth for a specified specific energy."""
    if specific_energy <= 0.0:
        raise ValueError("Specific energy must be positive.")
    critical = critical_depth(RectangularSection(width), discharge, gravity)
    critical_energy = rectangular_specific_energy(discharge, width, critical, gravity)
    if specific_energy < critical_energy:
        raise ValueError("Specific energy is below the minimum critical energy.")
    if regime not in {"subcritical", "supercritical"}:
        raise ValueError("regime must be 'subcritical' or 'supercritical'.")
    if regime == "subcritical":
        lower, upper = critical, specific_energy
    else:
        lower, upper = 1e-9, critical
    return _bisect(
        lambda depth: rectangular_specific_energy(discharge, width, depth, gravity)
        - specific_energy,
        lower,
        upper,
    )


def rectangular_froude(
    discharge: float,
    width: float,
    depth: float,
    gravity: float = GRAVITY,
) -> float:
    """Return the Froude number for rectangular-channel flow."""
    velocity = rectangular_velocity(discharge, width, depth)
    return velocity / sqrt(gravity * depth)


def sluice_gate_force(
    discharge: float,
    width: float,
    upstream_depth: float,
    downstream_depth: float,
    density: float = 1000.0,
    gravity: float = GRAVITY,
) -> float:
    """Return force exerted by water on a vertical sluice gate.

    Positive force acts in the downstream direction. Hydrostatic pressure and
    momentum flux are included; channel friction and gate thickness are omitted.
    """
    velocity_upstream = rectangular_velocity(discharge, width, upstream_depth)
    velocity_downstream = rectangular_velocity(discharge, width, downstream_depth)
    pressure_difference = (
        gravity * width * (upstream_depth**2 - downstream_depth**2) / 2.0
    )
    momentum_change = discharge * (velocity_downstream - velocity_upstream)
    return density * (pressure_difference - momentum_change)


def supercritical_reach_length(
    discharge_per_width: float,
    upstream_depth: float,
    downstream_depth: float,
    roughness: float,
    slope: float,
    steps: int = 20000,
    gravity: float = GRAVITY,
) -> float:
    """Estimate distance from a control section to a jump by direct integration.

    The flow is assumed rectangular, steady, gradually varied, and supercritical.
    ``upstream_depth`` is the depth at the control section and should exceed
    ``downstream_depth`` for the profile used here.
    """
    if not 0.0 < downstream_depth < upstream_depth:
        raise ValueError("Expected 0 < downstream_depth < upstream_depth.")
    if steps < 2:
        raise ValueError("steps must be at least 2.")
    increment = (upstream_depth - downstream_depth) / steps

    def integrand(depth: float) -> float:
        area = depth
        radius = depth
        friction_slope = (
            roughness * discharge_per_width / (area * radius ** (2.0 / 3.0))
        ) ** 2
        froude_squared = discharge_per_width**2 / (gravity * depth**3)
        return (1.0 - froude_squared) / (slope - friction_slope)

    total = 0.0
    for index in range(steps):
        left = downstream_depth + index * increment
        right = left + increment
        total += 0.5 * (integrand(left) + integrand(right)) * increment
    return abs(total)


def broad_crested_weir_transition(
    discharge: float,
    width: float,
    weir_height: float,
    gravity: float = GRAVITY,
) -> dict[str, float]:
    """Calculate ideal broad-crested-weir control depths.

    The crest is assumed to be high enough to force critical flow and
    head losses over the control are neglected.  The returned upstream
    depth is the subcritical alternate depth and the downstream depth is
    the supercritical alternate depth, both referenced to the original
    channel bed.
    """
    if discharge <= 0.0 or width <= 0.0 or weir_height < 0.0:
        raise ValueError("Discharge and width must be positive; weir height must be non-negative.")

    section = RectangularSection(width)
    unit_discharge = discharge / width
    critical = (unit_discharge**2 / gravity) ** (1.0 / 3.0)
    critical_energy = 1.5 * critical
    required_total_energy = weir_height + critical_energy

    upstream_depth = rectangular_depth_from_specific_energy(
        discharge=discharge,
        width=width,
        specific_energy=required_total_energy,
        regime="subcritical",
        gravity=gravity,
    )
    downstream_depth = rectangular_depth_from_specific_energy(
        discharge=discharge,
        width=width,
        specific_energy=required_total_energy,
        regime="supercritical",
        gravity=gravity,
    )

    return {
        "critical_depth": critical,
        "critical_energy_above_crest": critical_energy,
        "required_total_energy": required_total_energy,
        "upstream_depth": upstream_depth,
        "downstream_depth": downstream_depth,
        "upstream_froude": rectangular_froude(discharge, width, upstream_depth, gravity),
        "downstream_froude": rectangular_froude(discharge, width, downstream_depth, gravity),
    }


def gvf_one_step_jump_distance(
    discharge_per_width: float,
    initial_depth: float,
    target_depth: float,
    roughness: float,
    slope: float,
    gravity: float = GRAVITY,
) -> float:
    """Estimate a supercritical hydraulic-jump position with one GVF step.

    The gradually-varied-flow equation is evaluated at the mean depth
    between the initial supercritical depth and the target (conjugate)
    depth.  This gives a single-step estimate of the horizontal distance.
    """
    if discharge_per_width <= 0.0 or initial_depth <= 0.0 or target_depth <= 0.0:
        raise ValueError("Discharge per width and depths must be positive.")
    if roughness <= 0.0 or slope <= 0.0:
        raise ValueError("Roughness and slope must be positive.")
    if target_depth <= initial_depth:
        raise ValueError("For this downstream supercritical profile, target_depth must exceed initial_depth.")

    mean_depth = 0.5 * (initial_depth + target_depth)
    friction_slope = (
        roughness * discharge_per_width / mean_depth ** (5.0 / 3.0)
    ) ** 2
    froude_squared = discharge_per_width**2 / (gravity * mean_depth**3)

    denominator = slope - friction_slope
    numerator = 1.0 - froude_squared
    if abs(denominator) < 1e-14:
        raise ValueError("GVF denominator is too close to zero at the mean depth.")

    dx_dy = numerator / denominator
    delta_y = target_depth - initial_depth
    return abs(dx_dy * delta_y)


def abrupt_expansion_jump_depth(
    discharge: float,
    upstream_width: float,
    downstream_width: float,
    upstream_depth: float,
    gravity: float = GRAVITY,
) -> float:
    """Solve the downstream depth after an abrupt rectangular expansion.

    The calculation treats the expansion as a short hydraulic transition and
    applies conservation of momentum (specific force) between the two sections.
    Hydrostatic pressure distributions are assumed at both sections, while
    bed-friction and the weight component over the very short transition are
    neglected.  The physically relevant solution is the *subcritical* root
    downstream of the expansion.

    For a rectangular channel the momentum function is

        M = Q^2 / (g A) + A*y/2

    where A = b*y.  The downstream depth is found from M2 = M1.
    """
    if discharge <= 0.0:
        raise ValueError("Discharge must be positive.")
    if upstream_width <= 0.0 or downstream_width <= 0.0:
        raise ValueError("Channel widths must be positive.")
    if upstream_depth <= 0.0:
        raise ValueError("Upstream depth must be positive.")
    if gravity <= 0.0:
        raise ValueError("Gravity must be positive.")
    if downstream_width <= upstream_width:
        raise ValueError("This function requires a downstream expansion in width.")

    upstream_area = upstream_width * upstream_depth
    upstream_momentum = (
        discharge**2 / (gravity * upstream_area)
        + upstream_area * upstream_depth / 2.0
    )

    downstream_section = RectangularSection(downstream_width)
    downstream_critical = critical_depth(
        downstream_section, discharge, gravity
    )

    def momentum_difference(depth: float) -> float:
        area = downstream_width * depth
        momentum = discharge**2 / (gravity * area) + area * depth / 2.0
        return momentum - upstream_momentum

    # The larger root is the subcritical downstream solution.  Start at the
    # critical depth and expand the upper bracket until the root is enclosed.
    lower = downstream_critical
    upper = max(2.0 * lower, upstream_depth)
    while momentum_difference(upper) < 0.0:
        upper *= 2.0

    return _bisect(momentum_difference, lower, upper)

"""Footprint / coverage geometry of a satellite over a spherical Earth.

For a satellite at altitude h above a spherical Earth of radius Re,
the horizon as seen from the satellite subtends an Earth central
angle

    lambda_max = acos(Re / (Re + h))

The ground footprint is a spherical cap of that angular radius; its
radius measured along the surface is Re * lambda_max.  The nadir
angle (half-angle of the Earth as seen from the satellite) is

    eta = asin(Re / (Re + h)).

For a minimum elevation mask el_min at the ground station, the
central angle of the visible cap is

    lambda = acos( Re / (Re + h) * cos(el_min) ) - el_min

(Vallado, "Fundamentals of Astrodynamics and Applications", ch. 5 /
coverage geometry; the same construction appears in most satellite
communications textbooks, e.g. for link-budget coverage analysis.)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import DEG2RAD, EARTH_RADIUS_KM, RAD2DEG


@dataclass(frozen=True)
class Footprint:
    """Coverage geometry for one altitude (and optional elevation mask)."""

    altitude_km: float
    min_elevation_deg: float = 0.0
    central_angle_deg: float = 0.0  # Earth central angle of the cap
    footprint_radius_km: float = 0.0  # along the surface
    nadir_angle_deg: float = 0.0  # half-angle of Earth seen from satellite
    footprint_area_km2: float = 0.0  # spherical-cap area

    @property
    def diameter_km(self) -> float:
        """Footprint diameter along the surface, km."""
        return 2.0 * self.footprint_radius_km


def compute_footprint(altitude_km: float,
                      min_elevation_deg: float = 0.0) -> Footprint:
    """Footprint geometry for a satellite at ``altitude_km``.

    ``min_elevation_deg`` is the ground-station elevation mask; 0
    means the geometric horizon.
    """
    if altitude_km <= 0:
        raise ValueError(f"altitude must be positive, got {altitude_km}")
    if not 0.0 <= min_elevation_deg < 90.0:
        raise ValueError(f"elevation mask must be in [0, 90), got {min_elevation_deg}")
    re = EARTH_RADIUS_KM
    r = re + altitude_km
    el = min_elevation_deg * DEG2RAD
    lam = math.acos((re / r) * math.cos(el)) - el
    nadir = math.asin(re / r)
    radius = re * lam
    area = 2.0 * math.pi * re * re * (1.0 - math.cos(lam))
    return Footprint(
        altitude_km=altitude_km,
        min_elevation_deg=min_elevation_deg,
        central_angle_deg=lam * RAD2DEG,
        footprint_radius_km=radius,
        nadir_angle_deg=nadir * RAD2DEG,
        footprint_area_km2=area,
    )


def max_slant_range_km(altitude_km: float, min_elevation_deg: float = 0.0) -> float:
    """Longest possible slant range (satellite at the mask elevation)."""
    from .coordinates import slant_range_km
    return slant_range_km(altitude_km, min_elevation_deg)


def ground_track_speed_km_s(altitude_km: float) -> float:
    """Ground-track speed of a circular orbit at ``altitude_km``.

    v_ground = v_orbit * Re / (Re + h), with v_orbit = sqrt(mu / r).
    """
    from .constants import MU_EARTH_KM3_S2
    r = EARTH_RADIUS_KM + altitude_km
    v_orbit = math.sqrt(MU_EARTH_KM3_S2 / r)
    return v_orbit * EARTH_RADIUS_KM / r

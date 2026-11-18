"""Physical and numerical constants for satellite-pass-predictor.

Sources
-------
* WGS84 ellipsoid: National Geospatial-Intelligence Agency,
  "Department of Defense World Geodetic System 1984" (NIMA TR8350.2,
  3rd ed., 2000). a = 6378137.0 m, 1/f = 298.257223563.
* Earth gravitational parameter and J2: Vallado, "Fundamentals of
  Astrodynamics and Applications", 4th ed. (2013), Appendix D.
  mu = 398600.4418 km^3/s^2, J2 = 0.0010826269.
* Speed of light: exact SI value, 299792458 m/s.
"""

from __future__ import annotations

import math

# --- Earth (WGS84 / Vallado) ---------------------------------------------
MU_EARTH_KM3_S2: float = 398600.4418  # km^3 / s^2
EARTH_RADIUS_KM: float = 6378.137  # WGS84 semi-major axis, km
EARTH_RADIUS_M: float = 6378137.0  # WGS84 semi-major axis, m
EARTH_FLATTENING: float = 1.0 / 298.257223563
EARTH_ECC2: float = EARTH_FLATTENING * (2.0 - EARTH_FLATTENING)
EARTH_J2: float = 0.0010826269
EARTH_ROTATION_RATE_RAD_S: float = 7.2921159e-5  # rad/s (WGS84)
SECONDS_PER_DAY: float = 86400.0
SPEED_OF_LIGHT_M_S: float = 299792458.0

# --- Time ------------------------------------------------------------------
J2000_JD: float = 2451545.0  # JD of 2000-01-01 12:00 TT

# --- Numerical ---------------------------------------------------------------
DEG2RAD: float = math.pi / 180.0
RAD2DEG: float = 180.0 / math.pi


def wrap_degrees(angle_deg: float) -> float:
    """Wrap an angle in degrees to [0, 360)."""
    return angle_deg % 360.0


def wrap_radians(angle_rad: float) -> float:
    """Wrap an angle in radians to [0, 2*pi)."""
    return angle_rad % (2.0 * math.pi)

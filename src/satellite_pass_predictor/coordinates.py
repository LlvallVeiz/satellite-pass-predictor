"""Coordinate conversions: ECI <-> ECEF <-> geodetic, topocentric look angles.

Conventions
-----------
* ECI: J2000 mean equator / mean equinox, x towards vernal equinox.
  (Precession/nutation between epoch and J2000 are ignored; the
  Kepler+J2 model works in this frame, consistent with how TLE mean
  elements are used in practice for a first-order predictor.)
* ECEF: WGS84, x through Greenwich meridian, z through north pole.
  ECI -> ECEF is a rotation about z by the GMST (Meeus eq. 12.4).
  Polar motion is ignored (< 1 arcsecond effect).
* Geodetic latitude/longitude/height use the WGS84 ellipsoid
  (a = 6378137.0 m, 1/f = 298.257223563), solved by the standard
  fixed-point iteration (converges in a handful of steps).
* Topocentric look angles use the SEZ (south-east-zenith) frame:
  azimuth measured clockwise from north, elevation above the local
  horizon.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime

import numpy as np

from .constants import (
    DEG2RAD,
    EARTH_ECC2,
    EARTH_RADIUS_KM,
    EARTH_RADIUS_M,
    RAD2DEG,
)
from .time_utils import gmst_radians


def eci_to_ecef(position_eci_km: np.ndarray, time: datetime) -> np.ndarray:
    """Rotate an ECI position (km) into ECEF (km) using the GMST."""
    theta = gmst_radians(time)
    c, s = math.cos(theta), math.sin(theta)
    rot = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    return rot @ np.asarray(position_eci_km, dtype=float)


def ecef_to_eci(position_ecef_km: np.ndarray, time: datetime) -> np.ndarray:
    """Rotate an ECEF position (km) back into ECI (km)."""
    theta = gmst_radians(time)
    c, s = math.cos(theta), math.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return rot @ np.asarray(position_ecef_km, dtype=float)


def geodetic_to_ecef(lat_deg: float, lon_deg: float, height_m: float) -> np.ndarray:
    """Geodetic (lat, lon, height) -> ECEF position in km (WGS84)."""
    lat = lat_deg * DEG2RAD
    lon = lon_deg * DEG2RAD
    sin_lat = math.sin(lat)
    n = EARTH_RADIUS_M / math.sqrt(1.0 - EARTH_ECC2 * sin_lat * sin_lat)
    x = (n + height_m) * math.cos(lat) * math.cos(lon)
    y = (n + height_m) * math.cos(lat) * math.sin(lon)
    z = (n * (1.0 - EARTH_ECC2) + height_m) * sin_lat
    return np.array([x, y, z]) / 1000.0


@dataclass(frozen=True)
class GeodeticPosition:
    """Geodetic position on the WGS84 ellipsoid."""

    latitude_deg: float
    longitude_deg: float
    height_m: float = 0.0


def ecef_to_geodetic(position_ecef_km: np.ndarray) -> GeodeticPosition:
    """ECEF position (km) -> geodetic latitude/longitude/height (WGS84).

    Fixed-point iteration on latitude (standard textbook method);
    converges to < 1e-9 deg in a few iterations for LEO altitudes.
    """
    x, y, z = (float(v) * 1000.0 for v in np.asarray(position_ecef_km, dtype=float))
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - EARTH_ECC2))
    for _ in range(50):
        sin_lat = math.sin(lat)
        n = EARTH_RADIUS_M / math.sqrt(1.0 - EARTH_ECC2 * sin_lat * sin_lat)
        height = p / math.cos(lat) - n
        new_lat = math.atan2(z, p * (1.0 - EARTH_ECC2 * n / (n + height)))
        if abs(new_lat - lat) < 1e-12:
            lat = new_lat
            break
        lat = new_lat
    sin_lat = math.sin(lat)
    n = EARTH_RADIUS_M / math.sqrt(1.0 - EARTH_ECC2 * sin_lat * sin_lat)
    height = p / math.cos(lat) - n
    return GeodeticPosition(
        latitude_deg=lat * RAD2DEG,
        longitude_deg=lon * RAD2DEG,
        height_m=height,
    )


@dataclass(frozen=True)
class LookAngles:
    """Topocentric look angles from a ground station to a satellite."""

    azimuth_deg: float  # clockwise from north, [0, 360)
    elevation_deg: float  # above local horizon, [-90, 90]
    range_km: float
    range_rate_km_s: float = 0.0  # positive = receding


def ecef_to_look_angles(sat_ecef_km: np.ndarray, station_ecef_km: np.ndarray,
                        lat_deg: float, lon_deg: float,
                        range_rate_km_s: float = 0.0) -> LookAngles:
    """Look angles from a station (geodetic lat/lon) to a satellite.

    Both positions are ECEF in km.  Uses the SEZ frame.
    """
    rho = np.asarray(sat_ecef_km, dtype=float) - np.asarray(station_ecef_km, dtype=float)
    rng = float(np.linalg.norm(rho))
    lat = lat_deg * DEG2RAD
    lon = lon_deg * DEG2RAD
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    # SEZ components.
    s = sin_lat * cos_lon * rho[0] + sin_lat * sin_lon * rho[1] - cos_lat * rho[2]
    e = -sin_lon * rho[0] + cos_lon * rho[1]
    z = cos_lat * cos_lon * rho[0] + cos_lat * sin_lon * rho[1] + sin_lat * rho[2]
    elevation = math.asin(max(-1.0, min(1.0, z / rng)))
    azimuth = math.atan2(e, -s)  # from north, clockwise
    azimuth_deg = azimuth * RAD2DEG % 360.0
    if azimuth_deg >= 360.0:  # tiny negative angles round to 360.0
        azimuth_deg = 0.0
    return LookAngles(
        azimuth_deg=azimuth_deg,
        elevation_deg=elevation * RAD2DEG,
        range_km=rng,
        range_rate_km_s=range_rate_km_s,
    )


def slant_range_km(altitude_km: float, elevation_deg: float) -> float:
    """Slant range from a ground station to a satellite at a given elevation.

    Geometry of the triangle Earth-centre / station / satellite
    (spherical Earth, radius 6378.137 km).  By the law of cosines,
    with el the elevation above the horizon:

        rho = -Re sin(el) + sqrt((Re+h)^2 - (Re cos el)^2)
    """
    el = elevation_deg * DEG2RAD
    re = EARTH_RADIUS_KM
    r = re + altitude_km
    # rho^2 + 2 rho Re sin(el) - (r^2 - Re^2) = 0  -> positive root.
    return -re * math.sin(el) + math.sqrt(r * r - (re * math.cos(el)) ** 2)

"""Tests for coordinate conversions (coordinates.py).

Reference values
----------------
* WGS84: equator radius 6378.137 km; polar radius a(1-f) = 6356.752314 km.
* Geodetic (0, 0, 0) -> ECEF (6378.137, 0, 0) km.
* Slant range at 90 deg elevation equals the altitude; at 0 deg it is
  sqrt((Re+h)^2 - Re^2) (Pythagoras).
"""

from __future__ import annotations

import math
from datetime import datetime

import numpy as np
import pytest

from satellite_pass_predictor.constants import EARTH_RADIUS_KM
from satellite_pass_predictor.coordinates import (
    ecef_to_geodetic,
    ecef_to_eci,
    ecef_to_look_angles,
    eci_to_ecef,
    geodetic_to_ecef,
    slant_range_km,
)


def test_geodetic_equator():
    ecef = geodetic_to_ecef(0.0, 0.0, 0.0)
    assert ecef == pytest.approx([EARTH_RADIUS_KM, 0.0, 0.0], abs=1e-9)


def test_geodetic_north_pole():
    ecef = geodetic_to_ecef(90.0, 0.0, 0.0)
    assert ecef == pytest.approx([0.0, 0.0, 6356.752314245], abs=1e-6)


def test_geodetic_longitude_90():
    ecef = geodetic_to_ecef(0.0, 90.0, 0.0)
    assert ecef == pytest.approx([0.0, EARTH_RADIUS_KM, 0.0], abs=1e-9)


def test_geodetic_height_adds_radially_at_equator():
    ecef = geodetic_to_ecef(0.0, 0.0, 1000.0)
    assert ecef[0] == pytest.approx(EARTH_RADIUS_KM + 1.0, abs=1e-9)


def test_geodetic_round_trip():
    for lat, lon, h in ((-6.2, 106.8, 8.0), (51.5, -0.12, 50.0),
                        (-33.9, 151.2, 100.0), (80.0, -40.0, 5000.0),
                        (0.0, 0.0, 400000.0)):
        ecef = geodetic_to_ecef(lat, lon, h)
        back = ecef_to_geodetic(ecef)
        assert back.latitude_deg == pytest.approx(lat, abs=1e-7)
        assert back.longitude_deg == pytest.approx(lon, abs=1e-7)
        assert back.height_m == pytest.approx(h, abs=1e-3)


def test_eci_ecef_round_trip():
    t = datetime(2026, 11, 16, 12, 0, 0)
    pos = np.array([5000.0, 3000.0, 2000.0])
    back = ecef_to_eci(eci_to_ecef(pos, t), t)
    assert back == pytest.approx(pos, abs=1e-9)


def test_eci_ecef_rotation_is_pure_z():
    # A pure z-rotation preserves the z component and the norm.
    t = datetime(2026, 11, 16, 12, 0, 0)
    pos = np.array([5000.0, 3000.0, 2000.0])
    ecef = eci_to_ecef(pos, t)
    assert ecef[2] == pytest.approx(2000.0, abs=1e-9)
    assert np.linalg.norm(ecef) == pytest.approx(np.linalg.norm(pos), abs=1e-9)


def test_eci_ecef_gmst_zero_case():
    # At J2000 (GMST = 280.46061837 deg) the rotation is known exactly.
    t = datetime(2000, 1, 1, 12, 0, 0)
    pos = np.array([EARTH_RADIUS_KM, 0.0, 0.0])
    ecef = eci_to_ecef(pos, t)
    theta = math.radians(280.46061837)
    expected = np.array([EARTH_RADIUS_KM * math.cos(theta),
                         -EARTH_RADIUS_KM * math.sin(theta), 0.0])
    assert ecef == pytest.approx(expected, abs=1e-6)


def test_slant_range_overhead_equals_altitude():
    assert slant_range_km(418.4, 90.0) == pytest.approx(418.4, abs=1e-9)


def test_slant_range_horizon_pythagoras():
    h = 418.4
    expected = math.sqrt((EARTH_RADIUS_KM + h) ** 2 - EARTH_RADIUS_KM ** 2)
    assert slant_range_km(h, 0.0) == pytest.approx(expected, abs=1e-9)


def test_slant_range_decreases_with_elevation():
    ranges = [slant_range_km(500.0, el) for el in (0.0, 10.0, 30.0, 60.0, 90.0)]
    assert ranges == sorted(ranges, reverse=True)


def test_look_angles_overhead():
    # Satellite directly above the station (same lat/lon, 400 km up).
    lat, lon = -6.2, 106.8
    station = geodetic_to_ecef(lat, lon, 0.0)
    sat = geodetic_to_ecef(lat, lon, 400_000.0)
    look = ecef_to_look_angles(sat, station, lat, lon)
    assert look.elevation_deg == pytest.approx(90.0, abs=0.01)
    assert look.range_km == pytest.approx(400.0, abs=0.01)


def test_look_angles_horizon():
    # Satellite on the horizon due east of the station: elevation ~0.
    # Build it from the station's local east direction.
    lat, lon = 0.0, 0.0
    station = geodetic_to_ecef(lat, lon, 0.0)
    east = np.array([0.0, 1.0, 0.0])  # at lon 0, lat 0, east is +y
    sat = station + east * 1000.0
    look = ecef_to_look_angles(sat, station, lat, lon)
    assert look.azimuth_deg == pytest.approx(90.0, abs=1.0)
    assert abs(look.elevation_deg) < 10.0  # geometry: 1000 km east, ~9 deg below


def test_look_angles_azimuth_range():
    lat, lon = -6.2, 106.8
    station = geodetic_to_ecef(lat, lon, 0.0)
    for dlat, dlon in ((5.0, 0.0), (-5.0, 0.0), (0.0, 5.0), (0.0, -5.0)):
        sat = geodetic_to_ecef(lat + dlat, lon + dlon, 400_000.0)
        look = ecef_to_look_angles(sat, station, lat, lon)
        assert 0.0 <= look.azimuth_deg < 360.0
        assert -90.0 <= look.elevation_deg <= 90.0
        assert look.range_km > 0.0


def test_look_angles_north_is_zero_azimuth():
    # Satellite due north of the station at the same longitude.
    lat, lon = -6.2, 106.8
    station = geodetic_to_ecef(lat, lon, 0.0)
    sat = geodetic_to_ecef(lat + 5.0, lon, 400_000.0)
    look = ecef_to_look_angles(sat, station, lat, lon)
    assert look.azimuth_deg == pytest.approx(0.0, abs=2.0)

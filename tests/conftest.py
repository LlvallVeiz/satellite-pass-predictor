"""Shared fixtures for the satellite-pass-predictor test suite."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from satellite_pass_predictor.ground_station import BANDUNG, JAKARTA, SINGAPORE, GroundStation
from satellite_pass_predictor.orbit import Propagator
from satellite_pass_predictor.tle import load_tle_file

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


@pytest.fixture(scope="session")
def tles():
    """The four real TLEs from data/satellites.tle."""
    return load_tle_file(str(DATA_DIR / "satellites.tle"))


@pytest.fixture(scope="session")
def iss_tle(tles):
    return next(t for t in tles if t.sat_number == 25544)


@pytest.fixture(scope="session")
def noaa_tle(tles):
    return next(t for t in tles if t.sat_number == 33591)


@pytest.fixture(scope="session")
def css_tle(tles):
    return next(t for t in tles if t.sat_number == 48274)


@pytest.fixture(scope="session")
def so50_tle(tles):
    return next(t for t in tles if t.sat_number == 27607)


@pytest.fixture(scope="session")
def iss_prop(iss_tle):
    return Propagator(iss_tle)


@pytest.fixture(scope="session")
def jakarta():
    return JAKARTA


@pytest.fixture(scope="session")
def stations():
    return [JAKARTA, BANDUNG, SINGAPORE]


@pytest.fixture(scope="session")
def day_start():
    """2026-11-16 00:00 UTC — the day of the TLE epochs."""
    return datetime(2026, 11, 16, 0, 0, 0)

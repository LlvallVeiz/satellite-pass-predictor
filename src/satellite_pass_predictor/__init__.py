"""satellite-pass-predictor: predict satellite passes over ground stations.

Kepler + secular-J2 propagation from TLEs, topocentric look angles,
pass prediction, Doppler shift, footprint geometry, and Markdown/HTML
reports.  See README.md and docs/RESEARCH.md.
"""

from __future__ import annotations

from .batch import BatchResult, predict_batch
from .coordinates import GeodeticPosition, LookAngles
from .doppler import doppler_shift_hz, range_rate_km_s
from .footprint import Footprint, compute_footprint
from .ground_station import GroundStation
from .orbit import Propagator, StateVector
from .passes import Pass, predict_passes
from .tle import TLE, TLEError, load_tle_file, parse_tle, parse_tle_text
from .time_utils import datetime_to_julian_date, gmst_degrees

__version__ = "1.0.0"

__all__ = [
    "__version__",
    "BatchResult",
    "Footprint",
    "GeodeticPosition",
    "GroundStation",
    "LookAngles",
    "Pass",
    "Propagator",
    "StateVector",
    "TLE",
    "TLEError",
    "compute_footprint",
    "datetime_to_julian_date",
    "doppler_shift_hz",
    "gmst_degrees",
    "load_tle_file",
    "parse_tle",
    "parse_tle_text",
    "predict_batch",
    "predict_passes",
    "range_rate_km_s",
]

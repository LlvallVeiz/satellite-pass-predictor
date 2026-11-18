"""Doppler shift of a satellite downlink as seen from a ground station.

The (non-relativistic) Doppler shift of a carrier at frequency f0 is

    delta_f = -f0 * (range rate) / c

where the range rate is the rate of change of the station-to-satellite
distance.  Positive range rate (satellite receding) lowers the received
frequency.  The range rate is computed analytically in the ECEF frame:

    v_sat_ECEF = R(GMST) v_ECI - omega_earth x r_ECEF
    range rate = rho_hat . v_sat_ECEF

Note the station is fixed in the ECEF frame, so its ECEF velocity
is zero; the (inertial) station velocity omega_earth x r_station
must NOT be subtracted here.  A numerical derivative of the range
is also provided as an independent cross-check.
"""

from __future__ import annotations

import math
from datetime import datetime, timedelta

import numpy as np

from .constants import EARTH_ROTATION_RATE_RAD_S, SPEED_OF_LIGHT_M_S
from .coordinates import eci_to_ecef
from .ground_station import GroundStation
from .orbit import Propagator, StateVector

OMEGA_EARTH = np.array([0.0, 0.0, EARTH_ROTATION_RATE_RAD_S])


def ecef_velocity_km_s(state: StateVector) -> np.ndarray:
    """Satellite velocity in the rotating ECEF frame, km/s."""
    theta = _gmst(state.time)
    c, s = math.cos(theta), math.sin(theta)
    rot = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    v_eci_rotated = rot @ state.velocity_km_s
    r_ecef = eci_to_ecef(state.position_km, state.time)
    return v_eci_rotated - np.cross(OMEGA_EARTH, r_ecef)


def _gmst(time: datetime) -> float:
    from .time_utils import gmst_radians
    return gmst_radians(time)


def range_rate_km_s(state: StateVector, station: GroundStation) -> float:
    """Range rate (km/s) of a satellite state as seen from a station.

    Positive = satellite receding from the station.
    """
    r_ecef = eci_to_ecef(state.position_km, state.time)
    rho = r_ecef - station.ecef_km
    rho_hat = rho / float(np.linalg.norm(rho))
    v_sat = ecef_velocity_km_s(state)
    return float(np.dot(rho_hat, v_sat))


def range_rate_numerical_km_s(propagator: Propagator, time: datetime,
                              station: GroundStation,
                              step_s: float = 1.0) -> float:
    """Range rate by centred finite difference of the slant range.

    Independent cross-check of :func:`range_rate_km_s`.
    """
    def rng(t: datetime) -> float:
        st = propagator.propagate(t)
        r = eci_to_ecef(st.position_km, t)
        return float(np.linalg.norm(r - station.ecef_km))

    return (rng(time + timedelta(seconds=step_s))
            - rng(time - timedelta(seconds=step_s))) / (2.0 * step_s)


def doppler_shift_hz(state: StateVector, station: GroundStation,
                     frequency_hz: float) -> float:
    """Doppler shift (Hz) of a downlink carrier at ``frequency_hz``.

    Positive shift = received frequency higher than transmitted
    (satellite approaching).
    """
    rr = range_rate_km_s(state, station)  # km/s
    return -frequency_hz * (rr * 1000.0) / SPEED_OF_LIGHT_M_S


def doppler_at_pass(propagator: Propagator, station: GroundStation,
                    times: list[datetime], frequency_hz: float) -> list[float]:
    """Doppler shift (Hz) along a sampled pass."""
    return [doppler_shift_hz(propagator.propagate(t), station, frequency_hz)
            for t in times]

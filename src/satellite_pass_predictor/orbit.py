"""Orbit propagation: Kepler solver plus secular J2 corrections.

Model
-----
Mean elements from a TLE are propagated with the classical two-body
(Kepler) solution, with the secular rates caused by Earth's oblateness
(J2) applied to the node, perigee and mean anomaly.  The rates are the
standard first-order expressions (Vallado, "Fundamentals of
Astrodynamics and Applications", 4th ed., ch. 9; also Hoots &
Roehrich, Spacetrack Report No. 3, 1980, SGP secular terms):

    p  = a (1 - e^2)
    dOmega/dt = -3/2 J2 (Re/p)^2 n cos i
    domega/dt =  3/4 J2 (Re/p)^2 n (5 cos^2 i - 1)
    dM/dt     =  n + 3/4 J2 (Re/p)^2 n sqrt(1 - e^2) (3 cos^2 i - 1)

with n the mean motion.  The semi-major axis, eccentricity and
inclination are held constant (first-order secular theory).  Drag
(B*) and all periodic J2 terms are neglected, so this is NOT SGP4:
accuracy is a few kilometres over a day or two near the epoch for
near-circular LEO.  The README states this honestly.

Kepler's equation M = E - e sin E is solved by Newton-Raphson.
Position/velocity use the standard perifocal-to-ECI rotation
R3(Omega) R1(i) R3(omega) (Vallado ch. 2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np

from .constants import (
    EARTH_J2,
    EARTH_RADIUS_KM,
    MU_EARTH_KM3_S2,
)
from .tle import TLE


@dataclass(frozen=True)
class StateVector:
    """Position (km) and velocity (km/s) in ECI at a given time."""

    time: datetime
    position_km: np.ndarray  # shape (3,)
    velocity_km_s: np.ndarray  # shape (3,)

    @property
    def radius_km(self) -> float:
        """Distance from Earth's centre, km."""
        return float(np.linalg.norm(self.position_km))

    @property
    def speed_km_s(self) -> float:
        """Speed, km/s."""
        return float(np.linalg.norm(self.velocity_km_s))

    @property
    def specific_energy_km2_s2(self) -> float:
        """Specific orbital energy, km^2/s^2 (constant in two-body motion)."""
        r = self.radius_km
        v2 = float(np.dot(self.velocity_km_s, self.velocity_km_s))
        return v2 / 2.0 - MU_EARTH_KM3_S2 / r

    @property
    def angular_momentum_km2_s(self) -> np.ndarray:
        """Specific angular momentum vector, km^2/s (constant, two-body)."""
        return np.cross(self.position_km, self.velocity_km_s)


@dataclass(frozen=True)
class J2Rates:
    """Secular J2 rates in rad/s."""

    d_raan: float
    d_arg_perigee: float
    d_mean_anomaly: float  # extra rate on top of mean motion n


def j2_secular_rates(a_km: float, e: float, i_rad: float, n_rad_s: float) -> J2Rates:
    """First-order secular J2 rates for the node, perigee, mean anomaly."""
    p = a_km * (1.0 - e * e)
    factor = 1.5 * EARTH_J2 * (EARTH_RADIUS_KM / p) ** 2 * n_rad_s
    cos_i = math.cos(i_rad)
    d_raan = -factor * cos_i
    d_argp = 0.5 * factor * (5.0 * cos_i * cos_i - 1.0)
    d_mean = 0.5 * factor * math.sqrt(1.0 - e * e) * (3.0 * cos_i * cos_i - 1.0)
    return J2Rates(d_raan=d_raan, d_arg_perigee=d_argp, d_mean_anomaly=d_mean)


def solve_kepler(mean_anomaly_rad: float, eccentricity: float,
                 tol: float = 1e-12, max_iter: int = 100) -> float:
    """Solve Kepler's equation M = E - e sin E for E (radians).

    Newton-Raphson, starting from E = M (good for e < 0.8; all TLE
    eccentricities in our data are < 0.01, and the solver is tested
    to e = 0.9).  Returns the eccentric anomaly wrapped to [0, 2*pi).
    """
    if not 0.0 <= eccentricity < 1.0:
        raise ValueError(f"eccentricity must be in [0, 1), got {eccentricity}")
    m = mean_anomaly_rad % (2.0 * math.pi)
    e_anom = m if eccentricity < 0.8 else math.pi
    for _ in range(max_iter):
        f = e_anom - eccentricity * math.sin(e_anom) - m
        f_prime = 1.0 - eccentricity * math.cos(e_anom)
        step = f / f_prime
        e_anom -= step
        if abs(step) < tol:
            break
    else:
        raise RuntimeError("Kepler solver did not converge")
    return e_anom % (2.0 * math.pi)


def _rotation_matrix(raan: float, inc: float, argp: float) -> np.ndarray:
    """Perifocal -> ECI rotation R3(Omega) R1(i) R3(omega)."""
    co, so = math.cos(raan), math.sin(raan)
    ci, si = math.cos(inc), math.sin(inc)
    cw, sw = math.cos(argp), math.sin(argp)
    return np.array([
        [co * cw - so * sw * ci, -co * sw - so * cw * ci, so * si],
        [so * cw + co * sw * ci, -so * sw + co * cw * ci, -co * si],
        [sw * si, cw * si, ci],
    ])


def elements_to_state(a_km: float, e: float, i_rad: float, raan_rad: float,
                      argp_rad: float, mean_anom_rad: float,
                      time: datetime) -> StateVector:
    """Convert Keplerian elements at one instant to an ECI state vector."""
    ecc_anom = solve_kepler(mean_anom_rad, e)
    cos_e, sin_e = math.cos(ecc_anom), math.sin(ecc_anom)
    sqrt_1me2 = math.sqrt(1.0 - e * e)
    # True anomaly.
    cos_nu = (cos_e - e) / (1.0 - e * cos_e)
    sin_nu = (sqrt_1me2 * sin_e) / (1.0 - e * cos_e)
    nu = math.atan2(sin_nu, cos_nu)
    # Radius and perifocal position.
    r = a_km * (1.0 - e * cos_e)
    pos_pf = np.array([r * math.cos(nu), r * math.sin(nu), 0.0])
    # Perifocal velocity (Vallado ch. 2).
    p = a_km * (1.0 - e * e)
    vel_pf = np.array([
        -math.sqrt(MU_EARTH_KM3_S2 / p) * math.sin(nu),
        math.sqrt(MU_EARTH_KM3_S2 / p) * (e + math.cos(nu)),
        0.0,
    ])
    rot = _rotation_matrix(raan_rad, i_rad, argp_rad)
    return StateVector(
        time=time,
        position_km=rot @ pos_pf,
        velocity_km_s=rot @ vel_pf,
    )


class Propagator:
    """Kepler + secular-J2 propagator for one TLE.

    Elements are advanced linearly from the epoch:
    Omega(t) = Omega0 + dOmega/dt * dt, and likewise for omega and M
    (the latter on top of n * dt).  a, e, i are constant.
    """

    def __init__(self, tle: TLE):
        self.tle = tle
        self.a_km = tle.semi_major_axis_km
        self.e = tle.eccentricity
        self.i_rad = tle.inclination
        self.n_rad_s = tle.mean_motion_rad_s
        self.rates = j2_secular_rates(self.a_km, self.e, self.i_rad, self.n_rad_s)

    def elements_at(self, time: datetime) -> tuple[float, float, float, float]:
        """(a, e, i, raan, argp, M) in (km, -, rad, rad, rad, rad) at time."""
        dt_s = (time - self.tle.epoch).total_seconds()
        raan = self.tle.raan + self.rates.d_raan * dt_s
        argp = self.tle.arg_perigee + self.rates.d_arg_perigee * dt_s
        mean_anom = (self.tle.mean_anomaly
                     + (self.n_rad_s + self.rates.d_mean_anomaly) * dt_s)
        return self.a_km, self.e, self.i_rad, raan, argp, mean_anom

    def propagate(self, time: datetime) -> StateVector:
        """ECI position/velocity at ``time``."""
        a, e, i, raan, argp, mean_anom = self.elements_at(time)
        return elements_to_state(a, e, i, raan, argp, mean_anom, time)

    def propagate_many(self, times: list[datetime]) -> list[StateVector]:
        """Propagate to a list of times."""
        return [self.propagate(t) for t in times]


def integrate_two_body_rk4(state: StateVector, dt_s: float,
                           steps: int) -> list[StateVector]:
    """Numerically integrate the pure two-body problem (RK4, fixed step).

    Independent cross-check for the analytic Kepler propagator: with
    J2 switched off, both must agree to integration tolerance over an
    orbit.  No perturbations are included here.
    """
    def accel(r: np.ndarray) -> np.ndarray:
        return -MU_EARTH_KM3_S2 * r / float(np.linalg.norm(r)) ** 3

    out = [state]
    r = state.position_km.copy()
    v = state.velocity_km_s.copy()
    t = state.time
    for _ in range(steps):
        k1v = accel(r)
        k1r = v
        k2v = accel(r + 0.5 * dt_s * k1r)
        k2r = v + 0.5 * dt_s * k1v
        k3v = accel(r + 0.5 * dt_s * k2r)
        k3r = v + 0.5 * dt_s * k2v
        k4v = accel(r + dt_s * k3r)
        k4r = v + dt_s * k3v
        r = r + (dt_s / 6.0) * (k1r + 2 * k2r + 2 * k3r + k4r)
        v = v + (dt_s / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
        t = t + timedelta(seconds=dt_s)
        out.append(StateVector(time=t, position_km=r.copy(), velocity_km_s=v.copy()))
    return out

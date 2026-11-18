"""Tests for orbit propagation (orbit.py).

Reference values are hand-computed from the TLE mean motion and the
standard J2 secular-rate formulas (Vallado ch. 9):

* ISS mean motion 15.49415609 rev/day -> n = 0.0011267 rad/s,
  a = 6796.5 km, period 92.938 min.
* J2 nodal regression for the ISS: about -4.95 deg/day
  (the textbook value quoted for the ISS is "about 5 deg/day").
* J2 perigee precession for the ISS: about +3.70 deg/day.
"""

from __future__ import annotations

import math
from datetime import datetime, timedelta

import numpy as np
import pytest

from satellite_pass_predictor.constants import MU_EARTH_KM3_S2
from satellite_pass_predictor.orbit import (
    Propagator,
    elements_to_state,
    integrate_two_body_rk4,
    j2_secular_rates,
    solve_kepler,
)


def test_kepler_circular():
    # e = 0 -> E = M exactly.
    for m in (0.0, 1.0, 3.0, 5.5):
        assert solve_kepler(m, 0.0) == pytest.approx(m % (2 * math.pi))


def test_kepler_residual():
    for e in (0.1, 0.5, 0.9):
        for m in (0.3, 1.0, 2.5, 4.0, 5.9):
            ecc = solve_kepler(m, e)
            assert ecc - e * math.sin(ecc) - m % (2 * math.pi) == pytest.approx(0.0, abs=1e-10)


def test_kepler_invalid_eccentricity():
    with pytest.raises(ValueError):
        solve_kepler(1.0, 1.0)
    with pytest.raises(ValueError):
        solve_kepler(1.0, -0.1)


def test_j2_rates_iss(iss_prop):
    rates = iss_prop.rates
    # deg/day: nodal regression ~ -4.95, perigee precession ~ +3.70.
    assert rates.d_raan * 86400 * 180 / math.pi == pytest.approx(-4.95, abs=0.05)
    assert rates.d_arg_perigee * 86400 * 180 / math.pi == pytest.approx(3.70, abs=0.05)
    # Mean-anomaly extra rate is small for the ISS.
    assert abs(rates.d_mean_anomaly) < 1e-6


def test_j2_rates_equatorial_orbit_has_no_regression():
    # cos(i) = 0 at i = 90 deg -> no nodal regression.
    rates = j2_secular_rates(7000.0, 0.0, math.radians(90.0), 0.0011)
    assert rates.d_raan == pytest.approx(0.0, abs=1e-12)


def test_j2_rates_critical_inclination_no_perigee_precession():
    # 5 cos^2(i) - 1 = 0 at i = 63.43 deg -> no perigee precession.
    rates = j2_secular_rates(7000.0, 0.01, math.radians(63.4349), 0.0011)
    assert rates.d_arg_perigee == pytest.approx(0.0, abs=1e-9)


def test_iss_state_at_epoch(iss_prop):
    state = iss_prop.propagate(iss_prop.tle.epoch)
    # Radius within [a(1-e), a(1+e)].
    a = iss_prop.a_km
    e = iss_prop.e
    assert a * (1 - e) - 0.01 <= state.radius_km <= a * (1 + e) + 0.01
    # ISS speed ~7.66 km/s.
    assert state.speed_km_s == pytest.approx(7.66, abs=0.05)


def test_elements_to_state_circular_equatorial():
    # Circular equatorial orbit: r = a on the x-axis at M = 0, v along +y.
    state = elements_to_state(7000.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              datetime(2026, 1, 1))
    assert state.position_km == pytest.approx([7000.0, 0.0, 0.0], abs=1e-9)
    v_circ = math.sqrt(MU_EARTH_KM3_S2 / 7000.0)
    assert state.velocity_km_s == pytest.approx([0.0, v_circ, 0.0], abs=1e-9)


def test_elements_to_state_apogee():
    # e = 0.1, M = pi -> apogee: r = a(1+e), speed = sqrt(mu(1-e)/(a(1+e))).
    a, e = 8000.0, 0.1
    state = elements_to_state(a, e, 0.0, 0.0, 0.0, math.pi, datetime(2026, 1, 1))
    assert state.radius_km == pytest.approx(a * (1 + e), abs=1e-6)
    v_ap = math.sqrt(MU_EARTH_KM3_S2 * (1 - e) / (a * (1 + e)))
    assert state.speed_km_s == pytest.approx(v_ap, abs=1e-6)


def test_specific_energy_conserved(iss_prop):
    # Energy -mu/(2a) is constant in this model (a is fixed).
    t0 = iss_prop.tle.epoch
    e0 = iss_prop.propagate(t0).specific_energy_km2_s2
    e1 = iss_prop.propagate(t0 + timedelta(minutes=92.9383)).specific_energy_km2_s2
    assert e1 == pytest.approx(e0, rel=1e-9)
    assert e0 == pytest.approx(-MU_EARTH_KM3_S2 / (2 * iss_prop.a_km), rel=1e-9)


def test_angular_momentum_magnitude_conserved(iss_prop):
    t0 = iss_prop.tle.epoch
    h0 = np.linalg.norm(iss_prop.propagate(t0).angular_momentum_km2_s)
    h1 = np.linalg.norm(iss_prop.propagate(t0 + timedelta(hours=5)).angular_momentum_km2_s)
    assert h1 == pytest.approx(h0, rel=1e-9)


def test_position_returns_after_one_period_approx(iss_prop):
    # After one Keplerian period the satellite is back near its start
    # (J2 moves the node/perigee slightly, so "near", not exact).
    t0 = iss_prop.tle.epoch
    p0 = iss_prop.propagate(t0).position_km
    p1 = iss_prop.propagate(t0 + timedelta(seconds=iss_prop.tle.period_seconds)).position_km
    assert np.linalg.norm(p1 - p0) < 60.0  # km; J2 drift over one orbit


def test_kepler_vs_rk4_cross_validation():
    # Independent check: analytic Kepler propagation vs numerical RK4
    # integration of the two-body problem (no J2 in either).
    a, e = 7000.0, 0.05
    i, raan, argp, m0 = 0.9, 0.5, 0.3, 0.7
    t0 = datetime(2026, 1, 1)
    n = math.sqrt(MU_EARTH_KM3_S2 / a ** 3)
    state0 = elements_to_state(a, e, i, raan, argp, m0, t0)
    period = 2 * math.pi / n
    quarter = period / 4.0
    dt = 5.0
    steps = int(quarter / dt)
    rk = integrate_two_body_rk4(state0, dt, steps)[-1]
    analytic = elements_to_state(a, e, i, raan, argp, m0 + n * dt * steps,
                                 t0 + timedelta(seconds=dt * steps))
    assert np.linalg.norm(rk.position_km - analytic.position_km) < 1.0  # km
    assert np.linalg.norm(rk.velocity_km_s - analytic.velocity_km_s) < 0.002  # km/s


def test_rk4_conserves_energy():
    a, e = 7000.0, 0.1
    state0 = elements_to_state(a, e, 0.5, 0.2, 0.1, 0.0, datetime(2026, 1, 1))
    traj = integrate_two_body_rk4(state0, 10.0, 60)
    energies = [s.specific_energy_km2_s2 for s in traj]
    assert max(energies) - min(energies) < 1e-6


def test_propagate_many(iss_prop):
    times = [iss_prop.tle.epoch + timedelta(minutes=m) for m in (0, 10, 20)]
    states = iss_prop.propagate_many(times)
    assert len(states) == 3
    assert states[0].time == times[0]
    assert states[2].radius_km > 6000.0


def test_elements_advance_with_j2(iss_prop):
    # RAAN decreases (regression) over a day for a prograde LEO orbit.
    t0 = iss_prop.tle.epoch
    _, _, _, raan0, _, _ = iss_prop.elements_at(t0)
    _, _, _, raan1, _, _ = iss_prop.elements_at(t0 + timedelta(days=1))
    delta_deg = (raan1 - raan0) * 180 / math.pi
    assert delta_deg == pytest.approx(-4.95, abs=0.05)


def test_noaa_period_longer_than_iss(noaa_tle, iss_tle):
    assert noaa_tle.period_minutes > iss_tle.period_minutes

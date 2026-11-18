# satellite-pass-predictor

Predict satellite passes over ground stations from real TLE data
(Python): Kepler propagation with secular J2 corrections, ECI→ECEF→
geodetic conversions, topocentric look angles, pass prediction
(AOS/LOS, maximum elevation, duration), Doppler shift of a
downlink, footprint/coverage geometry, multi-satellite ×
multi-station batch runs, and Markdown/HTML reports.

> **Status:** v1.0.0 — **122/122 tests pass.** Example outputs in
> [`examples/RESULTS.md`](examples/RESULTS.md) are real runs on
> real TLEs (ISS, NOAA 19, CSS/Tianhe, SO-50 — epoch 2026-11-16,
> downloaded from CelesTrak).

## Key results (real runs, 2026-11-16, mask 10°)

| Quantity | Value |
|---|---|
| ISS passes over Jakarta in 24 h | **3** |
| Best ISS pass over Jakarta | max el **50.0°**, duration **6.6 min** (AOS 18:50:25 UTC) |
| Longest pass in the batch | NOAA 19 over Jakarta/Bandung, **11.0 min** |
| Total passes, 4 sats × 3 stations in 24 h | **32** |
| Doppler, ISS 437 MHz downlink | **+9.62 kHz** at AOS, **−9.63 kHz** at LOS, **−27 Hz** at max elevation |
| ISS orbital period (from TLE mean motion) | **92.94 min** |
| ISS footprint radius (10° mask) | **1387 km** |
| J2 nodal regression, ISS | **−4.95°/day** |

## Background / theory

A **TLE** (two-line element set) is two 69-character lines
describing a satellite's orbit at one instant (the epoch): mean
motion, eccentricity, inclination, RAAN, argument of perigee and
mean anomaly, plus drag terms. Each line ends in a modulo-10
checksum, which this package validates on load.

**Propagation** solves Kepler's equation M = E − e·sin E
(Newton–Raphson) and converts the elements to an ECI state vector
with the standard perifocal→ECI rotation. The node, perigee and
mean anomaly are advanced with the first-order **secular J2 rates**
(Vallado ch. 9; Spacetrack Report No. 3):

    dΩ/dt = −(3/2)·J2·(Re/p)²·n·cos i
    dω/dt =  (3/4)·J2·(Re/p)²·n·(5cos²i − 1)

**Coordinates:** ECI→ECEF is a rotation by the GMST (Meeus
eq. 12.4); ECEF→geodetic uses the WGS84 ellipsoid; look angles
(azimuth/elevation/range) use the SEZ frame. **Passes** are found
by a coarse elevation scan plus bisection refinement of AOS/LOS
and golden-section search for maximum elevation. **Doppler** is
Δf = −f₀·(range rate)/c with the range rate computed analytically
in the ECEF frame. **Footprint** uses the spherical-Earth horizon
geometry λ = acos(Re/(Re+h)).

> **This is not SGP4.** TLE mean elements are fitted for SGP4, so
> a Kepler+J2 propagator is an approximation: a few kilometres over
> a day or two near the epoch for near-circular LEO. It is the right
> tool for pass planning and teaching; it is not precise ephemeris.
> See [`docs/RESEARCH.md`](docs/RESEARCH.md) for the full
> discussion, sources, and validation strategy.

## Project structure

```
satellite-pass-predictor/
├── src/satellite_pass_predictor/
│   ├── constants.py      # WGS84 / Vallado constants
│   ├── time_utils.py     # Julian Date, GMST (Meeus 12.4)
│   ├── tle.py            # TLE parser + checksum validation
│   ├── orbit.py          # Kepler solver, J2 rates, propagator, RK4 check
│   ├── coordinates.py    # ECI/ECEF/geodetic, look angles, slant range
│   ├── ground_station.py # GroundStation + example stations
│   ├── doppler.py        # range rate + Doppler shift
│   ├── footprint.py      # coverage geometry
│   ├── passes.py         # pass prediction (scan + bisection + golden section)
│   ├── batch.py          # multi-satellite x multi-station runs
│   ├── report.py         # Markdown + HTML reports
│   └── cli.py            # command-line interface
├── data/satellites.tle   # real TLEs from CelesTrak (epoch 2026-11-16)
├── tests/                # 122 pytest tests
├── examples/             # 5 runnable scenarios + RESULTS.md
└── docs/RESEARCH.md      # sources, formulas, validation
```

## Installation

```bash
pip install -r requirements.txt   # numpy only
pip install -e .                  # optional: installs the CLI
```

Python ≥ 3.10, NumPy ≥ 1.26. No other dependencies.

## Usage

```bash
# Summarise the TLE file
python -m satellite_pass_predictor tles data/satellites.tle

# Position/velocity of the ISS at a time
python -m satellite_pass_predictor propagate data/satellites.tle \
    --sat 25544 --time 2026-11-16T13:27:14

# ISS passes over Jakarta in the next 24 h (mask 10 deg)
python -m satellite_pass_predictor passes data/satellites.tle \
    --sat 25544 --station-name Jakarta --lat -6.2 --lon 106.8 \
    --start 2026-11-16T00:00:00 --hours 24

# Batch: all satellites x stations, Markdown + HTML reports
python -m satellite_pass_predictor batch data/satellites.tle \
    --stations "Jakarta,-6.2,106.8,8" "Bandung,-6.9,107.6,768" "Singapore,1.35,103.82,15" \
    --start 2026-11-16T00:00:00 --hours 24 --mask 10 \
    --frequency 437000000 --output report.md --html report.html

# Doppler shift along the best pass
python -m satellite_pass_predictor doppler data/satellites.tle \
    --sat 25544 --station-name Jakarta --lat -6.2 --lon 106.8 \
    --start 2026-11-16T00:00:00 --hours 24 --frequency 437000000

# Footprint geometry for an altitude
python -m satellite_pass_predictor footprint --altitude 418.4 --mask 10
```

Or run the examples directly:

```bash
python3 examples/run_iss_jakarta.py
python3 examples/run_multi_station.py   # writes examples/output/report.md/.html
python3 examples/run_doppler.py
python3 examples/run_footprint.py
python3 examples/run_ground_track.py
```

## Example output (real run)

```
ISS (ZARYA) over Jakarta, 2026-11-16 UTC, mask 10 deg
Orbital period 92.94 min, mean altitude 418.4 km
Passes in 24 h: 3

ISS (ZARYA) over Jakarta: AOS 2026-11-16 05:37:51 az 137.2 deg, LOS 2026-11-16 05:38:51 az 120.0 deg, max el 10.2 deg at 05:38:21, duration 1.0 min
ISS (ZARYA) over Jakarta: AOS 2026-11-16 07:12:10 az 247.0 deg, LOS 2026-11-16 07:17:47 az 2.8 deg, max el 25.2 deg at 07:14:59, duration 5.6 min
ISS (ZARYA) over Jakarta: AOS 2026-11-16 18:50:25 az 338.1 deg, LOS 2026-11-16 18:56:59 az 131.4 deg, max el 50.0 deg at 18:53:42, duration 6.6 min
```

Full numbers: [`examples/RESULTS.md`](examples/RESULTS.md).

## Design decisions

- **Kepler + secular J2, not SGP4** — a deliberate, documented
  scope lock: full SGP4 is a large, subtle implementation; the
  chosen model is transparent, testable, and accurate enough for
  pass planning near the epoch. B\* is parsed and reported but not
  used; the README says so instead of pretending otherwise.
- **NumPy only** — vectors stay readable; no heavy astrodynamics
  dependency. Everything else is stdlib.
- **Analytic range rate, numerically cross-checked** — the Doppler
  module carries both an analytic ECEF range rate and a finite-
  difference one; the test suite requires them to agree. (This
  cross-check caught a real bug during development: the station's
  inertial velocity was wrongly subtracted in the ECEF frame.)
- **Real TLEs, checksummed** — the example data is genuine CelesTrak
  data with epoch 2026-11-16; the parser rejects any line whose
  modulo-10 checksum fails.

## Limitations & future work

- Not SGP4: positions drift kilometres per day from a true SGP4
  ephemeris; predictions degrade away from the epoch. A full SGP4
  implementation (Vallado 2006 code) is the natural next step.
- No drag model (B\* unused), no lunar/solar or higher-order gravity
  terms, no atmospheric refraction of the elevation angle.
- GMST uses UTC instead of UT1 and ignores nutation (sub-second).
- ECI frame is J2000 mean equator/equinox; precession between
  J2000 and the date is ignored.
- Pass search assumes one elevation maximum per pass (true for
  near-circular LEO; exotic eccentric orbits could violate it).

## Tests

```bash
python3 -m pytest -q        # 122 tests
```

Coverage: TLE parsing/checksums, Julian Date/GMST reference values,
Kepler solver residuals, J2 rates vs hand-computed values, energy
and angular-momentum conservation, Kepler-vs-RK4 cross-validation,
coordinate round-trips, slant-range geometry, Doppler sign/magnitude
and analytic-vs-numeric agreement, footprint closed forms, pass
structure (AOS/LOS at the mask, max elevation inside the pass),
batch totals, reports, and the CLI.

## License

MIT — see [LICENSE](LICENSE).

# Research notes — satellite-pass-predictor

How the physics in this package was chosen, what each formula is,
and where it comes from. Every claim below is tied to a source that
was actually read during this build (2026-11-17); nothing is cited
from memory alone.

## 1. Two-Line Element sets (TLE)

A TLE is two 69-character fixed-width lines describing one object's
orbit at one instant (the *epoch*). Line 1 carries the catalogue
number, epoch, drag terms (ndot, nddot, B\*) and an element-set
number; line 2 carries the six mean orbital elements plus mean
motion in revolutions per day. Column 69 of each line is a
modulo-10 checksum: add all digits in columns 1–68, count each
minus sign as 1, ignore letters, blanks, periods and plus signs;
the checksum is the last digit of the sum. Two-digit years 57–99
map to 1957–1999, 00–56 to 2000–2056.

Crucially, TLE elements are **mean elements fitted for SGP4** —
they are not osculating Keplerian elements, and a TLE only reaches
its stated accuracy when propagated with SGP4/SDP4 (Hoots &
Roehrich 1980; Vallado et al. 2006). This package deliberately
does **not** implement SGP4 (see §2); the consequences are stated
honestly in the README.

Source: Kelso, T. S., "Frequently Asked Questions: Two-Line Element
Set Format", CelesTrak, January 1998 (mirror read 2026-11-17).

## 2. Propagation model: Kepler + secular J2 (not SGP4)

The propagator solves Kepler's equation M = E − e·sin E by
Newton–Raphson, converts to an ECI state vector with the standard
perifocal→ECI rotation R3(Ω)·R1(i)·R3(ω), and advances the node,
perigee and mean anomaly with the first-order secular rates due to
Earth's oblateness (J2):

    p = a(1 − e²)
    dΩ/dt = −(3/2)·J2·(Re/p)²·n·cos i
    dω/dt =  (3/4)·J2·(Re/p)²·n·(5cos²i − 1)
    dM/dt =  n + (3/4)·J2·(Re/p)²·n·√(1−e²)·(3cos²i − 1)

with J2 = 0.0010826269, Re = 6378.137 km, μ = 398600.4418 km³/s².
These are the secular terms of the SGP theory as documented in
Spacetrack Report No. 3 and in Vallado ch. 9. Semi-major axis,
eccentricity and inclination are held constant; drag (B\*) and all
periodic J2 terms are neglected.

Sanity check against the literature: for the ISS the formulas give
dΩ/dt ≈ −4.95°/day — the value quoted in textbooks and tracking
guides as "about 5° per day westward" — and dω/dt ≈ +3.70°/day.
The code reproduces both (see `tests/test_orbit.py`).

Accuracy: for near-circular LEO a few kilometres over a day or two
near the epoch. Do not use it for precise ephemeris, conjunction
assessment, or predictions far from the epoch — use a real SGP4
implementation (e.g. the Vallado 2006 code) for that.

Sources: Hoots, F. R. & Roehrich, R. L., "Spacetrack Report No. 3:
Models for Propagation of NORAD Element Sets", U.S. Air Force,
December 1980; Vallado, D. A., Crawford, P., Hujsak, R. & Kelso,
T. S., "Revisiting Spacetrack Report #3", AIAA 2006-6753 (Rev. 1);
Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
4th ed., Microcosm Press, 2013, ch. 2 and 9.

## 3. Time: Julian Date and GMST

Julian Date conversion follows Meeus ch. 7. Greenwich Mean
Sidereal Time uses Meeus eq. (12.4), the IAU-adopted expression:

    GMST = 280.46061837 + 360.98564736629·(JD − 2451545.0)
           + 0.000387933·T² − T³/38710000        [degrees]

with T = (JD − 2451545.0)/36525. UT1 vs UTC (< 0.9 s) and nutation
(GMST vs GAST) are ignored — sub-second effects, far below this
model's other errors. Verified: GMST at J2000.0 is exactly
280.46061837°, and at 2025-06-21 20:00 UTC the formula gives
210.2666°, matching an independent Meeus-based calculator.

Source: Meeus, J., "Astronomical Algorithms", 2nd ed.,
Willmann-Bell, 1998, ch. 11–12.

## 4. Reference frame and geodetic conversions (WGS84)

ECI is the J2000 mean-equator/mean-equinox frame. ECI→ECEF is a
rotation about z by the GMST; polar motion is ignored. Geodetic
latitude/longitude/height use the WGS84 ellipsoid
(a = 6378137.0 m, 1/f = 298.257223563) with the standard
fixed-point latitude iteration. Topocentric look angles use the
SEZ (south-east-zenith) frame: azimuth clockwise from north,
elevation above the local horizon.

Source: National Geospatial-Intelligence Agency, "Department of
Defense World Geodetic System 1984", NIMA TR8350.2, 3rd ed., 2000
(WGS84 ellipsoid parameters).

## 5. Footprint / coverage geometry

For a satellite at altitude h over a spherical Earth of radius Re,
the horizon subtends the Earth central angle λ = acos(Re/(Re+h));
the footprint radius along the surface is Re·λ and the nadir angle
is asin(Re/(Re+h)). With a ground elevation mask el, the visible
central angle is λ = acos(Re/(Re+h)·cos el) − el. This is the
standard spherical-Earth coverage construction used in satellite
communications textbooks (cf. Vallado ch. 5).

## 6. Doppler shift

Non-relativistic Doppler: Δf = −f₀·(range rate)/c, with the range
rate computed analytically in the ECEF frame (the station is
fixed in ECEF, so its ECEF velocity is zero — an earlier draft
subtracted the station's inertial velocity and was caught by the
analytic-vs-numeric cross-check in `tests/test_doppler.py`). The
numerical finite-difference range rate agrees with the analytic
one to ~0.003 km/s.

## 7. Data

The TLEs in `data/satellites.tle` are real element sets downloaded
from CelesTrak (`https://celestrak.org/NORAD/elements/gp.php`,
query `CATNR=<number>&FORMAT=tle`) on 2026-11-17, all with epoch
2026-11-16: ISS (ZARYA) 25544, NOAA 19 33591, CSS (TIANHE) 48274,
SAUDISAT 1C (SO-50) 27607. TLE data is public orbital data
published by the U.S. Space Force 18th Space Defense Squadron and
redistributed by CelesTrak. The parser validates the modulo-10
checksum of every line on load.

## 8. Validation strategy

* **Hand-computed values:** ISS period 92.938 min from mean motion
  15.49415609 rev/day; J2 rates −4.95/+3.70°/day; GMST reference
  values; footprint geometry from the closed-form expressions.
* **Conservation laws:** specific energy −μ/(2a) and |h| are
  constant along the propagated orbit (a, e, i fixed by the model).
* **Independent re-implementation:** analytic Kepler propagation
  vs RK4 numerical integration of the two-body problem agree to
  < 1 km over a quarter orbit; analytic range rate vs finite
  difference agree to ~0.003 km/s; ECEF↔geodetic round-trips
  recover the input to < 1e-7 deg.
* **Structural pass checks:** elevation at refined AOS/LOS equals
  the mask to 0.1°, maximum elevation lies strictly inside the
  pass, durations are in the LEO range, and a station at latitude
  89° sees zero ISS passes (inclination 51.6°).

## References

1. Kelso, T. S. "Frequently Asked Questions: Two-Line Element Set
   Format." CelesTrak, January 1998.
   https://www.scss.tcd.ie/Stephen.Farrell/ipn/background/computers-and-satellites/v04n03.html
2. Hoots, F. R. & Roehrich, R. L. "Spacetrack Report No. 3: Models
   for Propagation of NORAD Element Sets." U.S. Air Force,
   December 1980. (Hosted by CelesTrak as `spacetrk.pdf`.)
3. Vallado, D. A., Crawford, P., Hujsak, R. & Kelso, T. S.
   "Revisiting Spacetrack Report #3." AIAA 2006-6753, Rev. 1, 2006.
   https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev1.pdf
4. Vallado, D. A. "Fundamentals of Astrodynamics and Applications."
   4th ed., Microcosm Press, 2013.
5. Meeus, J. "Astronomical Algorithms." 2nd ed., Willmann-Bell,
   1998. https://archive.org/stream/astronomicalalgorithmsjeanmeeus1991/Astronomical%20Algorithms-%20Jean%20Meeus%20%281991%29_djvu.txt
6. National Geospatial-Intelligence Agency. "Department of Defense
   World Geodetic System 1984." NIMA TR8350.2, 3rd ed., 2000.
7. CelesTrak. "GP element sets (TLE) query service."
   https://celestrak.org/NORAD/elements/gp.php

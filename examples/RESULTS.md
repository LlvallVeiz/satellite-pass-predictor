# Example results (real runs)

All numbers below are the actual output of the example scripts in
this directory, run on 2026-11-17 with the real TLEs in
`../data/satellites.tle` (epoch 2026-11-16). Window for the pass
predictions: 2026-11-16 00:00–24:00 UTC, elevation mask 10 deg.

## `run_iss_jakarta.py` — ISS over Jakarta, 24 h

| Pass | AOS (UTC) | LOS (UTC) | Max el | Duration |
|---|---|---|---|---|
| 1 | 05:37:51 | 05:38:51 | 10.2° | 1.0 min |
| 2 | 07:12:10 | 07:17:47 | 25.2° | 5.6 min |
| 3 | 18:50:25 | 18:56:59 | 50.0° | 6.6 min |

- **Passes in 24 h: 3.** Best pass: AOS 18:50:25 az 338.1°, max
  elevation 50.0° at 18:53:42, LOS 18:56:59 az 131.4°.
- Doppler at 437 MHz along the best pass: **+9.62 kHz at AOS**,
  **−9.63 kHz at LOS**, **−27 Hz at maximum elevation** (≈ 0, as
  theory requires at closest approach).
- Orbital period 92.94 min, mean altitude 418.4 km.

## `run_multi_station.py` — 4 satellites x 3 stations, 24 h

**Total passes: 32.**

| Satellite | Station | Passes | Best max el | Longest | Total visible |
|---|---|---|---|---|---|
| ISS (ZARYA) | Jakarta | 3 | 50.0° | 6.6 min | 13.2 min |
| ISS (ZARYA) | Bandung | 3 | 52.6° | 6.6 min | 14.6 min |
| ISS (ZARYA) | Singapore | 2 | 88.1° | 6.7 min | 13.4 min |
| NOAA 19 | Jakarta | 3 | 82.7° | 11.0 min | 27.2 min |
| NOAA 19 | Bandung | 3 | 75.1° | 11.0 min | 27.5 min |
| NOAA 19 | Singapore | 3 | 60.8° | 10.9 min | 26.3 min |
| CSS (TIANHE) | Jakarta | 3 | 77.6° | 6.3 min | 14.9 min |
| CSS (TIANHE) | Bandung | 3 | 78.3° | 6.3 min | 15.4 min |
| CSS (TIANHE) | Singapore | 2 | 57.9° | 6.2 min | 12.4 min |
| SAUDISAT 1C (SO-50) | Jakarta | 2 | 56.2° | 8.8 min | 17.2 min |
| SAUDISAT 1C (SO-50) | Bandung | 3 | 60.4° | 8.6 min | 18.2 min |
| SAUDISAT 1C (SO-50) | Singapore | 2 | 85.8° | 9.3 min | 17.8 min |

Reports written by this run: `output/report.md`, `output/report.html`.

## `run_footprint.py` — coverage geometry

| Satellite | Altitude | Ground-track speed | Footprint radius, mask 0° | Footprint radius, mask 10° | Footprint radius, mask 30° |
|---|---|---|---|---|---|
| ISS (ZARYA) | 418.4 km | 7.19 km/s | 2250 km | 1387 km | 628 km |
| NOAA 19 | 847.4 km | 6.56 km/s | 3120 km | 2184 km | 1129 km |
| CSS (TIANHE) | 385.8 km | 7.24 km/s | 2165 km | 1311 km | 585 km |
| SAUDISAT 1C (SO-50) | 618.8 km | 6.88 km/s | 2703 km | 1797 km | 876 km |

ISS footprint area: 15,736,391 km² (horizon) / 6,015,787 km²
(10° mask). NOAA 19, the highest of the four, has the largest
footprint (29,978,054 km² at the horizon).

## `run_ground_track.py` — ISS ground track, one orbit from epoch

From epoch 2026-11-16 13:27:14 UTC (lat 0.08°, lon 8.67°), the ISS
reaches max latitude 51.25° (13:52:14 UTC, lon 103.66°) and min
latitude −51.78° (14:37:14 UTC, lon −96.88°), matching its 51.63°
inclination. Geodetic height varies 414–434 km over the orbit
(mean altitude 418.4 km above a spherical Earth; the spread is the
WGS84 ellipsoid's oblateness). J2 nodal regression: −4.95°/day.

## `run_doppler.py` — Doppler curve, best ISS pass over Jakarta

| Time (UTC) | El | Range | Shift (437 MHz) |
|---|---|---|---|
| 18:50:25 | 10.0° | 1497 km | +9.62 kHz |
| 18:51:25 | 17.8° | 1111 km | +9.04 kHz |
| 18:52:25 | 30.6° | 767 km | +7.36 kHz |
| 18:53:25 | 48.3° | 553 km | +2.20 kHz |
| 18:54:25 | 40.8° | 623 km | −5.17 kHz |
| 18:55:25 | 24.0° | 914 km | −8.36 kHz |
| 18:56:25 | 14.0° | 1283 km | −9.38 kHz |

Classic S-curve: positive while approaching, zero at closest
approach, negative while receding.

# Changelog

## 1.0.0 — 2026-11-17

Initial release.

- TLE parser with modulo-10 checksum validation (real CelesTrak
  TLEs for ISS, NOAA 19, CSS/Tianhe, SO-50; epoch 2026-11-16).
- Kepler propagator with secular J2 rates (nodal regression,
  perigee precession, mean-anomaly rate); RK4 two-body integrator
  as an independent cross-check.
- ECI/ECEF/geodetic (WGS84) conversions, GMST (Meeus 12.4),
  topocentric look angles, slant range.
- Pass prediction: coarse scan + bisection AOS/LOS refinement +
  golden-section maximum elevation; track sampling.
- Doppler shift of a downlink (analytic ECEF range rate,
  cross-checked against finite differences).
- Footprint/coverage geometry (central angle, radius, area,
  nadir angle) and ground-track speed.
- Multi-satellite x multi-station batch runs; Markdown + HTML
  reports with a Doppler section.
- CLI (`tles`, `propagate`, `passes`, `batch`, `footprint`,
  `doppler`); 5 runnable examples with real results.
- 122/122 pytest tests pass.

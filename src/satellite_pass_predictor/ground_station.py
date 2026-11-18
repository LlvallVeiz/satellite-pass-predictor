"""Ground stations: location, ECEF position, and look-angle helper."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .coordinates import LookAngles, ecef_to_look_angles, eci_to_ecef, geodetic_to_ecef
from .orbit import StateVector


@dataclass(frozen=True)
class GroundStation:
    """A ground station: name + geodetic position (WGS84)."""

    name: str
    latitude_deg: float
    longitude_deg: float
    height_m: float = 0.0

    def __post_init__(self) -> None:
        if not -90.0 <= self.latitude_deg <= 90.0:
            raise ValueError(f"latitude out of range: {self.latitude_deg}")
        if not -180.0 <= self.longitude_deg <= 180.0:
            raise ValueError(f"longitude out of range: {self.longitude_deg}")

    @property
    def ecef_km(self) -> np.ndarray:
        """Station position in ECEF, km."""
        return geodetic_to_ecef(self.latitude_deg, self.longitude_deg, self.height_m)

    def look_angles(self, state: StateVector,
                    range_rate_km_s: float = 0.0) -> LookAngles:
        """Azimuth/elevation/range to a satellite state (ECI)."""
        sat_ecef = eci_to_ecef(state.position_km, state.time)
        return ecef_to_look_angles(
            sat_ecef, self.ecef_km,
            self.latitude_deg, self.longitude_deg,
            range_rate_km_s=range_rate_km_s,
        )

    def elevation_deg(self, state: StateVector) -> float:
        """Elevation angle (deg) to a satellite state."""
        return self.look_angles(state).elevation_deg


# --- Example stations used by the examples/ and tests ----------------------
JAKARTA = GroundStation("Jakarta", -6.2, 106.8, 8.0)
BANDUNG = GroundStation("Bandung", -6.9, 107.6, 768.0)
SINGAPORE = GroundStation("Singapore", 1.35, 103.82, 15.0)

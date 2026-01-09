"""Utilities for the Roman Galactic Bulge Time Domain Survey (GBTDS).

This module provides lightweight helpers for working with common GBTDS field
layouts used in community tooling. The primary function is
`gbtds_field_from_radec`, which assigns an (approximate) GBTDS field identifier
to a given sky position by selecting the nearest configured field center.

Notes
-----
GBTDS "fields" are survey pointings on the sky (not SIAF apertures). The field
definitions here are center coordinates in Galactic longitude/latitude (l, b)
drawn from the open-source `gbtds_optimizer` tooling. Since the per-field
footprints/roll angles are not part of those layout files, the mapping here is
based on nearest-center assignment and (optionally) a maximum angular distance.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Sequence, overload

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord


@dataclass(frozen=True)
class GbtdsFieldCenter:
    """A single GBTDS field center definition."""

    field: str
    l_deg: float
    b_deg: float
    fixed: bool = False


# Field center definitions in Galactic lon/lat (degrees).
# Sources: https://github.com/mtpenny/gbtds_optimizer (field_layouts/*.centers)
_LAYOUTS: dict[str, tuple[GbtdsFieldCenter, ...]] = {
    # 9-field baseline-style layout; fields labeled 1..9
    "layout_9f_3": (
        GbtdsFieldCenter("1", 1.1215, -1.14, False),
        GbtdsFieldCenter("2", 0.7134, -1.14, False),
        GbtdsFieldCenter("3", 0.3053, -1.14, False),
        GbtdsFieldCenter("4", 1.5296, -1.14, False),
        GbtdsFieldCenter("5", -0.1028, -1.14, False),
        GbtdsFieldCenter("6", 0.3053, -1.93, False),
        GbtdsFieldCenter("7", 0.7134, -1.93, False),
        GbtdsFieldCenter("8", -0.1028, -1.93, False),
        GbtdsFieldCenter("9", 1.1215, -1.93, False),
    ),
    # 4-field subset; fields labeled 3,5,6,8 (no fixed column in source file)
    "layout_4f_1": (
        GbtdsFieldCenter("3", 0.3053, -1.14, False),
        GbtdsFieldCenter("5", -0.1028, -1.14, False),
        GbtdsFieldCenter("6", 0.3053, -1.93, False),
        GbtdsFieldCenter("8", -0.1028, -1.93, False),
    ),
    # 4-field subset plus an explicit Galactic Center point
    "layout_4f_1gc": (
        GbtdsFieldCenter("3", 0.3053, -1.14, False),
        GbtdsFieldCenter("5", -0.1028, -1.14, False),
        GbtdsFieldCenter("6", 0.3053, -1.93, False),
        GbtdsFieldCenter("8", -0.1028, -1.93, False),
        GbtdsFieldCenter("GC", 0.0, 0.0, True),
    ),
    # 8-field layout with a fixed GC point included in the source file
    "layout_8f_1overlap": (
        GbtdsFieldCenter("1", 1.1215, -1.14, False),
        GbtdsFieldCenter("2", 0.7134, -1.14, False),
        GbtdsFieldCenter("3", 0.3053, -1.14, False),
        GbtdsFieldCenter("4", -0.1028, -1.14, False),
        GbtdsFieldCenter("5", 1.1215, -1.79721, False),
        GbtdsFieldCenter("6", 0.7134, -1.79721, False),
        GbtdsFieldCenter("7", 0.3053, -1.79721, False),
        GbtdsFieldCenter("8", -0.1028, -1.79721, False),
        GbtdsFieldCenter("GC", 0.0, 0.0, True),
    ),
    # Layouts labeled 0..N with one fixed "GC-ish" point at (l,b)=(0,-0.125)
    "layout_40395": (
        GbtdsFieldCenter("0", -0.417948, -1.2, False),
        GbtdsFieldCenter("1", -0.008974, -1.2, False),
        GbtdsFieldCenter("2", 0.4, -1.2, False),
        GbtdsFieldCenter("3", 0.808974, -1.2, False),
        GbtdsFieldCenter("4", 1.217948, -1.2, False),
        GbtdsFieldCenter("5", 0.0, -0.125, True),
    ),
    "layout_88997": (
        GbtdsFieldCenter("0", -0.622435, -1.2, False),
        GbtdsFieldCenter("1", -0.213461, -1.2, False),
        GbtdsFieldCenter("2", 0.195513, -1.2, False),
        GbtdsFieldCenter("3", 0.604487, -1.2, False),
        GbtdsFieldCenter("4", 1.013461, -1.2, False),
        GbtdsFieldCenter("5", 1.422435, -1.2, False),
        GbtdsFieldCenter("6", 0.0, -0.125, True),
    ),
    "layout_163000": (
        GbtdsFieldCenter("0", 0.497757, -1.005558, False),
        GbtdsFieldCenter("1", 0.90673, -1.005558, False),
        GbtdsFieldCenter("2", 1.315704, -1.005558, False),
        GbtdsFieldCenter("3", -0.115704, -1.794442, False),
        GbtdsFieldCenter("4", 0.29327, -1.794442, False),
        GbtdsFieldCenter("5", 0.702243, -1.794442, False),
        GbtdsFieldCenter("6", 1.111217, -1.794442, False),
        GbtdsFieldCenter("7", 0.0, -0.125, True),
    ),
    "layout_wo407196": (
        GbtdsFieldCenter("0", 0.497757, -1.005558, False),
        GbtdsFieldCenter("1", 0.90673, -1.005558, False),
        GbtdsFieldCenter("2", 1.315704, -1.005558, False),
        GbtdsFieldCenter("3", -0.115704, -1.794442, False),
        GbtdsFieldCenter("4", 0.29327, -1.794442, False),
        GbtdsFieldCenter("5", 0.702243, -1.794442, False),
        GbtdsFieldCenter("6", 1.111217, -1.794442, False),
    ),
}

# Convenience aliases (more ergonomic than typing "layout_...").
_ALIASES: dict[str, str] = {
    "9f_3": "layout_9f_3",
    "4f_1": "layout_4f_1",
    "4f_1gc": "layout_4f_1gc",
    "8f_1overlap": "layout_8f_1overlap",
    "40395": "layout_40395",
    "88997": "layout_88997",
    "163000": "layout_163000",
    "wo407196": "layout_wo407196",
}


def get_gbtds_layout_names() -> tuple[str, ...]:
    """Return available GBTDS layout names (including aliases)."""
    return tuple(sorted(set(_LAYOUTS.keys()) | set(_ALIASES.keys())))


def get_gbtds_layout(layout: str = "layout_9f_3") -> tuple[GbtdsFieldCenter, ...]:
    """Return the field-center definitions for a named layout."""
    key = _ALIASES.get(layout, layout)
    try:
        return _LAYOUTS[key]
    except KeyError as exc:
        raise KeyError(
            f"Unknown GBTDS layout {layout!r}. Available: {get_gbtds_layout_names()}"
        ) from exc


def _as_icrs_skycoord(
    ra: object,
    dec: object | None,
    *,
    frame: str,
    unit: u.Unit | tuple[u.Unit, u.Unit] | None,
) -> SkyCoord:
    if isinstance(ra, SkyCoord):
        if dec is not None:
            raise TypeError("If ra is a SkyCoord, dec must be None.")
        return ra
    if dec is None:
        raise TypeError("dec is required unless ra is an astropy.coordinates.SkyCoord.")
    return SkyCoord(ra, dec, frame=frame, unit=unit)


@overload
def gbtds_field_from_radec(
    ra: SkyCoord | float | Sequence[float] | u.Quantity,
    dec: None | float | Sequence[float] | u.Quantity,
    *,
    layout: str = ...,
    include_fixed: bool = ...,
    frame: str = ...,
    unit: u.Unit | tuple[u.Unit, u.Unit] | None = ...,
    max_sep: u.Quantity | None = ...,
    return_sep: Literal[False] = ...,
) -> object: ...


@overload
def gbtds_field_from_radec(
    ra: SkyCoord | float | Sequence[float] | u.Quantity,
    dec: None | float | Sequence[float] | u.Quantity,
    *,
    layout: str = ...,
    include_fixed: bool = ...,
    frame: str = ...,
    unit: u.Unit | tuple[u.Unit, u.Unit] | None = ...,
    max_sep: u.Quantity | None = ...,
    return_sep: Literal[True],
) -> tuple[object, u.Quantity]: ...


def gbtds_field_from_radec(
    ra: SkyCoord | float | Sequence[float] | u.Quantity,
    dec: None | float | Sequence[float] | u.Quantity,
    *,
    layout: str = "layout_9f_3",
    include_fixed: bool = True,
    frame: str = "icrs",
    unit: u.Unit | tuple[u.Unit, u.Unit] | None = u.deg,
    max_sep: u.Quantity | None = None,
    return_sep: bool = False,
):
    """Assign a GBTDS field label to a sky position.

    Parameters
    ----------
    ra, dec : float, array-like, Quantity, or SkyCoord
        Sky position. If `ra` is a `~astropy.coordinates.SkyCoord`, then `dec`
        must be None and the coordinate is used directly.
    layout : str
        Which GBTDS layout to use. See `get_gbtds_layout_names()`.
    include_fixed : bool
        Include any centers marked as "fixed" in the layout (e.g. explicit "GC"
        markers) when selecting the nearest field.
    frame : str
        Astropy coordinate frame name for the input `ra, dec` values (ignored
        if `ra` is a `SkyCoord`).
    unit : Unit or (Unit, Unit) or None
        Units for the input `ra, dec` if not provided as Quantity/SkyCoord.
        Defaults to degrees.
    max_sep : Quantity or None
        If provided, return None for any position whose nearest field center is
        farther away than this separation.
    return_sep : bool
        If True, also return the separation to the selected field center.

    Returns
    -------
    field : object
        Field label(s). For scalar input returns a single string (or None).
        For array input returns a NumPy object array.
    sep : Quantity, optional
        Separation(s) to the selected center, returned if `return_sep=True`.
    """
    centers = get_gbtds_layout(layout)
    if not include_fixed:
        centers = tuple(c for c in centers if not c.fixed)
        if len(centers) == 0:
            raise ValueError(f"Layout {layout!r} has no non-fixed centers to select from.")

    target = _as_icrs_skycoord(ra, dec, frame=frame, unit=unit)
    target_gal = target.galactic

    field_labels = np.array([c.field for c in centers], dtype=object)
    center_coords = SkyCoord(
        l=[c.l_deg for c in centers] * u.deg,
        b=[c.b_deg for c in centers] * u.deg,
        frame="galactic",
    )

    # Astropy does not always broadcast (n,) against (m,) the way NumPy does.
    # Ensure a trailing "field" axis for separation results.
    if target_gal.isscalar:
        sep = target_gal.separation(center_coords)  # (nfields,)
    else:
        sep = target_gal[..., np.newaxis].separation(center_coords[np.newaxis, ...])  # (..., nfields)

    idx = np.argmin(sep, axis=-1)

    # Collect labels + min separation with broadcasting-safe indexing
    chosen_fields = field_labels[idx]
    min_sep = np.take_along_axis(sep, np.expand_dims(idx, axis=-1), axis=-1)[..., 0]

    if max_sep is not None:
        max_sep_q = u.Quantity(max_sep)
        mask = min_sep > max_sep_q
        chosen_fields = np.where(mask, None, chosen_fields).astype(object)

    # Preserve scalar-ness
    if target_gal.isscalar:
        # chosen_fields may be a Python scalar or a 0-d ndarray (after np.where)
        if isinstance(chosen_fields, np.ndarray):
            chosen_fields_out: object = chosen_fields.item()
        else:
            chosen_fields_out = chosen_fields
        min_sep_out = min_sep
    else:
        chosen_fields_out = chosen_fields
        min_sep_out = min_sep

    if return_sep:
        return chosen_fields_out, min_sep_out
    return chosen_fields_out


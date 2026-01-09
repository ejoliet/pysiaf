import astropy.units as u
from astropy.coordinates import SkyCoord

from pysiaf.utils.gbtds import gbtds_field_from_radec, get_gbtds_layout


def test_gbtds_field_from_radec_matches_centers_layout_9f_3():
    centers = get_gbtds_layout("layout_9f_3")

    for c in centers:
        if c.fixed:
            continue
        icrs = SkyCoord(l=c.l_deg * u.deg, b=c.b_deg * u.deg, frame="galactic").icrs
        field = gbtds_field_from_radec(icrs, None, layout="layout_9f_3")
        assert field == c.field


def test_gbtds_field_from_radec_max_sep_masks_far_positions():
    # Pick something far from the Bulge fields
    field = gbtds_field_from_radec(ra=0.0, dec=0.0, layout="layout_9f_3", max_sep=1.0 * u.deg)
    assert field is None


def test_gbtds_field_from_radec_vectorized():
    centers = [c for c in get_gbtds_layout("layout_9f_3") if not c.fixed]
    icrs0 = SkyCoord(l=centers[0].l_deg * u.deg, b=centers[0].b_deg * u.deg, frame="galactic").icrs
    icrs1 = SkyCoord(l=centers[1].l_deg * u.deg, b=centers[1].b_deg * u.deg, frame="galactic").icrs

    fields = gbtds_field_from_radec([icrs0.ra.deg, icrs1.ra.deg], [icrs0.dec.deg, icrs1.dec.deg])
    assert list(fields) == [centers[0].field, centers[1].field]


"""
Tests for the `IdealShock` class.
"""

from rgfrosh import IdealShock

from numpy.testing import assert_allclose
import pytest


class TestRatios:
    r"""
    Test `IdealShock` against values derived from example calculations in Table II.3 of
    Gaydon and Hurle [^1].

    [^1]: Gaydon, A. G. and I. R. Hurle (1963). The shock tube in high-temperature chemical
    physics, Reinhold Publishing Corporation.
    """

    rtol = 0.5e-2

    @staticmethod
    @pytest.mark.parametrize(
        ("gamma", "M1", "P2_1"),
        [(7 / 5, 2.95, 10), (7 / 5, 6.56, 50), (5 / 3, 2.87, 10), (5 / 3, 6.34, 50)],
    )
    def test_incident_pressure_ratio(gamma, M1, P2_1):
        assert_allclose(
            IdealShock.incident_pressure_ratio(M1, gamma), P2_1, rtol=TestRatios.rtol
        )

    @staticmethod
    @pytest.mark.parametrize(
        ("gamma", "M1", "T2_1"),
        [
            (7 / 5, 2.95, 2.62),
            (7 / 5, 6.56, 9.31),
            (5 / 3, 2.87, 3.42),
            (5 / 3, 6.34, 13.4),
        ],
    )
    def test_incident_temperature_ratio(gamma, M1, T2_1):
        assert_allclose(
            IdealShock.incident_temperature_ratio(M1, gamma), T2_1, rtol=TestRatios.rtol
        )

    @staticmethod
    @pytest.mark.parametrize(
        ("gamma", "M1", "rho2_1"),
        [
            (7 / 5, 2.95, 3.82),
            (7 / 5, 6.56, 5.37),
            (5 / 3, 2.87, 2.92),
            (5 / 3, 6.34, 3.72),
        ],
    )
    def test_incident_density_ratio(gamma, M1, rho2_1):
        assert_allclose(
            IdealShock.incident_density_ratio(M1, gamma), rho2_1, rtol=TestRatios.rtol
        )

    @staticmethod
    @pytest.mark.parametrize(
        ("gamma", "M1", "P2_1", "P5_2"),
        [
            (7 / 5, 2.95, 10, 4.95),
            (7 / 5, 6.56, 50, 7.12),
            (5 / 3, 2.87, 10, 4.22),
            (5 / 3, 6.34, 50, 5.54),
        ],
    )
    def test_reflected_pressure_ratio(gamma, M1, P2_1, P5_2):
        assert_allclose(
            IdealShock.reflected_pressure_ratio(M1, gamma),
            P5_2 * P2_1,
            rtol=TestRatios.rtol,
        )

    @staticmethod
    @pytest.mark.parametrize(
        ("gamma", "M1", "T2_1", "T5_2"),
        [
            (7 / 5, 2.95, 2.62, 1.76),
            pytest.param(7 / 5, 6.56, 9.31, 2.28, marks=pytest.mark.xfail),
            (5 / 3, 2.87, 3.42, 1.94),
            pytest.param(5 / 3, 6.34, 13.4, 2.37, marks=pytest.mark.xfail),
        ],
    )
    def test_reflected_temperature_ratio(gamma, M1, T2_1, T5_2):
        assert_allclose(
            IdealShock.reflected_temperature_ratio(M1, gamma),
            T5_2 * T2_1,
            rtol=TestRatios.rtol,
        )

    @staticmethod
    @pytest.mark.parametrize(
        ("gamma", "M1", "VR_S"),
        [
            (7 / 5, 2.95, 0.423),
            (7 / 5, 6.56, 0.351),
            (5 / 3, 2.87, 0.589),
            (5 / 3, 6.34, 0.517),
        ],
    )
    def test_reflected_velocity_ratio(gamma, M1, VR_S):
        assert_allclose(
            IdealShock.reflected_velocity_ratio(M1, gamma), VR_S, rtol=TestRatios.rtol
        )


@pytest.mark.parametrize("gamma,MW", [(7 / 5, 28), (5 / 3, 40)])  # [Nitrogen, Argon]
@pytest.mark.parametrize("M", [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
def test_consistency(gamma, MW, M):
    """
    Tests for consistency between the initialization approaches using the following steps:

    1. Initialize an `IdealShock` object using `M`, `T1`, and `P1`
    2. Initialize another `IdealShock` object using the calculated `T5` and `P5` from step 1
    3. Check that all properties are consistent between the objects from steps 1 and 2
    """

    from_initial = IdealShock(gamma, MW, M=M, T1=300, P1=101325)
    from_target = IdealShock(
        gamma, MW, T5=from_initial.T5, P5=from_initial.P5, T1=from_initial.T1
    )

    assert_allclose(from_initial.u1, from_target.u1)
    assert_allclose(from_initial.T1, from_target.T1)
    assert_allclose(from_initial.P1, from_target.P1)
    assert_allclose(from_initial.rho1, from_target.rho1)

    assert_allclose(from_initial.u2, from_target.u2)
    assert_allclose(from_initial.T2, from_target.T2)
    assert_allclose(from_initial.P2, from_target.P2)
    assert_allclose(from_initial.rho2, from_target.rho2)

    assert_allclose(from_initial.u5, from_target.u5)
    assert_allclose(from_initial.T5, from_target.T5)
    assert_allclose(from_initial.P5, from_target.P5)
    assert_allclose(from_initial.rho5, from_target.rho5)

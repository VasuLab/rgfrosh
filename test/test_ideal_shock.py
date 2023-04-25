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

    !!! Note
        Higher precision values were substituted based on fixed values of $\gamma$ and
        $P_2/P_1$; the original values are included in comments. Therefore, this test
        operates as a regression test.

    [^1]: Gaydon, A. G. and I. R. Hurle (1963). The shock tube in high-temperature chemical
    physics, Reinhold Publishing Corporation.
    """

    atol = 1e-3

    arg_names = "gamma,M1,P2_1,P5_2,T2_1,T5_2,rho2_1,rho5_2,VR_S"
    data = (
        # 2.95, 10, 4.95, 2.62, 1.76, 3.82, 2.82, 0.423
        (7 / 5, 2.95200, 10, 4.9375, 2.6230, 1.7634, 3.8125, 2.82, 0.423),
        # 6.56, 50, 7.12, 9.31, 2.28, 5.37, 3.11, 0.351
        (
            7 / 5,
            6.55744,
            50,
            7.1249,
            9.3022,
            2.1375,
            5.3750,
            3.11,
            0.351,
        ),  # T5/T2 discrepancy
        # 2.87, 10, 4.22, 3.42, 1.94, 2.92, 2.17, 0.589
        (5 / 3, 2.86356, 10, 4.2144, 3.4147, 1.9386, 2.9286, 2.17, 0.589),
        # 6.34, 50, 5.54, 13.4, 2.37, 3.72, 2.31, 0.517
        (
            5 / 3,
            6.34035,
            50,
            5.5369,
            13.4326,
            2.2813,
            3.7222,
            2.31,
            0.517,
        ),  # T5/T2 discrepancy
    )

    @staticmethod
    @pytest.mark.parametrize(arg_names, data)
    def test_incident_pressure_ratio(
        gamma, M1, P2_1, P5_2, T2_1, T5_2, rho2_1, rho5_2, VR_S
    ):
        assert_allclose(
            IdealShock.incident_pressure_ratio(M1, gamma), P2_1, atol=TestRatios.atol
        )

    @staticmethod
    @pytest.mark.parametrize(arg_names, data)
    def test_incident_temperature_ratio(
        gamma, M1, P2_1, P5_2, T2_1, T5_2, rho2_1, rho5_2, VR_S
    ):
        assert_allclose(
            IdealShock.incident_temperature_ratio(M1, gamma), T2_1, atol=TestRatios.atol
        )

    @staticmethod
    @pytest.mark.parametrize(arg_names, data)
    def test_incident_density_ratio(
        gamma, M1, P2_1, P5_2, T2_1, T5_2, rho2_1, rho5_2, VR_S
    ):
        assert_allclose(
            IdealShock.incident_density_ratio(M1, gamma), rho2_1, atol=TestRatios.atol
        )

    @staticmethod
    @pytest.mark.parametrize(arg_names, data)
    def test_reflected_pressure_ratio(
        gamma, M1, P2_1, P5_2, T2_1, T5_2, rho2_1, rho5_2, VR_S
    ):
        assert_allclose(
            IdealShock.reflected_pressure_ratio(M1, gamma)
            / IdealShock.incident_pressure_ratio(M1, gamma),
            P5_2,
            atol=TestRatios.atol,
        )

    @staticmethod
    @pytest.mark.parametrize(arg_names, data)
    def test_reflected_temperature_ratio(
        gamma, M1, P2_1, P5_2, T2_1, T5_2, rho2_1, rho5_2, VR_S
    ):
        assert_allclose(
            IdealShock.reflected_temperature_ratio(M1, gamma)
            / IdealShock.incident_temperature_ratio(M1, gamma),
            T5_2,
            atol=TestRatios.atol,
        )

    @staticmethod
    @pytest.mark.parametrize(arg_names, data)
    def test_reflected_velocity_ratio(
        gamma, M1, P2_1, P5_2, T2_1, T5_2, rho2_1, rho5_2, VR_S
    ):
        assert_allclose(
            IdealShock.reflected_velocity_ratio(M1, gamma), VR_S, atol=TestRatios.atol
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

"""
Tests for the `IdealShock` class.
"""

from rgfrosh import IdealShock

from numpy.testing import assert_allclose
import pytest


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

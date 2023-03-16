"""
Test FrozenShock output against output of the original RGFROSH[^1] for the
Peng-Robinson EOS for nitrogen (from table in the addendum).

[^1]: Davidson, D. F. and R. K. Hanson (1996). "Real Gas Corrections in Shock Tube
Studies at High Pressures." Israel Journal of Chemistry 36(3): 321-326.
"""

from rgfrosh import FrozenShock
from rgfrosh.interface import CPInterface
from numpy.testing import assert_allclose
import CoolProp as CP
import pytest


testdata = [
    (762, 1, 814.9, 19.95, 8.303),
    (762, 10, 809.0, 197.37, 78.143),
    (1000, 1, 1208.4, 45.56, 12.732),
    (1000, 10, 1191.9, 450.47, 116.194),
    (1523, 1, 2338, 145.0, 20.817),
    (1523, 10, 2286, 1432, 181.60),
]


@pytest.fixture
def gas():
    state = CP.AbstractState("PR", "Nitrogen")
    state.specify_phase(CP.iphase_supercritical_gas)
    return CPInterface(state)


@pytest.mark.parametrize("u1,P1,T5,P5,rho5", testdata)
def test_reflected_conditions(gas, u1, P1, T5, P5, rho5):
    shock = FrozenShock(gas, u1, 293, P1 * 101325)
    assert_allclose(
        [shock.T5, shock.P5 / 101325, shock.rho5],
        [T5, P5, rho5],
        rtol=1e-3,
    )

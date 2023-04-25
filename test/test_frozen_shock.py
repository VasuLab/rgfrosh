"""
Tests for the `FrozenShock` class.
"""

from rgfrosh import ThermoInterface, IdealShock, FrozenShock
from rgfrosh.thermo import CPInterface
from rgfrosh.constants import GAS_CONSTANT

import cantera as ct
import CoolProp as CP

from numpy.testing import assert_almost_equal, assert_allclose
import pytest


FrozenShock.rtol = 1e-9


class PerfectGas(ThermoInterface):
    """Thermo interface for a calorically perfect gas."""

    def __init__(self, gamma, MW):
        self._T = 300
        self._P = 101325
        self.gamma = gamma
        self.MW = MW

    @property
    def TP(self):
        return self._T, self._P

    @TP.setter
    def TP(self, value):
        self._T, self._P = value

    @property
    def mean_molecular_weight(self):
        return self.MW

    @property
    def density_mass(self):
        return self._P / (GAS_CONSTANT / self.MW * self._T)

    @property
    def cp_mass(self):
        return self.gamma / (self.gamma - 1) * GAS_CONSTANT / self.MW

    @property
    def enthalpy_mass(self):
        return self.cp_mass * self._T

    @property
    def isothermal_compressibility(self):
        return 1 / self._P

    @property
    def thermal_expansion_coeff(self):
        return 1 / self._T

    @property
    def sound_speed(self):
        return (self.gamma * GAS_CONSTANT / self.MW * self._T) ** 0.5


@pytest.mark.parametrize("gamma,MW", [(7 / 5, 28), (5 / 3, 40)])  # [Nitrogen, Argon]
@pytest.mark.parametrize("M", [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
class TestIdealCase:
    """
    Tests that the `FrozenShock` implementation correctly reduces to the ideal shock equations[^1]
    for the case of a calorically perfect gas `ThermoInterface`.

    [^1]: Gaydon, A. G. and I. R. Hurle (1963). The shock tube in high-temperature chemical physics.
    """

    @pytest.fixture
    def ideal_case(self, M, gamma, MW):
        perfect_gas = PerfectGas(gamma, MW)
        return FrozenShock(
            perfect_gas, u1=M * perfect_gas.sound_speed, T1=300, P1=101325
        )

    def test_incident_temperature_ratio(self, ideal_case, M, gamma, MW):
        assert_almost_equal(
            ideal_case.T2 / ideal_case.T1,
            IdealShock.incident_temperature_ratio(M, gamma),
        )

    def test_incident_pressure_ratio(self, ideal_case, M, gamma, MW):
        assert_almost_equal(
            ideal_case.P2 / ideal_case.P1,
            IdealShock.incident_pressure_ratio(M, gamma),
        )

    def test_incident_density_ratio(self, ideal_case, M, gamma, MW):
        assert_almost_equal(
            ideal_case.rho2 / ideal_case.rho1,
            IdealShock.incident_density_ratio(M, gamma),
        )

    def test_incident_velocity_ratio(self, ideal_case, M, gamma, MW):
        assert_almost_equal(
            ideal_case.u1 / ideal_case.u2,
            IdealShock.incident_density_ratio(M, gamma),
        )

    def test_reflected_temperature_ratio(self, ideal_case, M, gamma, MW):
        assert_almost_equal(
            ideal_case.T5 / ideal_case.T1,
            IdealShock.reflected_temperature_ratio(M, gamma),
        )

    def test_reflected_pressure_ratio(self, ideal_case, M, gamma, MW):
        assert_almost_equal(
            ideal_case.P5 / ideal_case.P1,
            IdealShock.reflected_pressure_ratio(M, gamma),
        )

    @pytest.mark.xfail
    def test_reflected_velocity_ratio(self, ideal_case, M, gamma, MW):
        assert_almost_equal(
            ideal_case.u5 / ideal_case.u1,
            IdealShock.reflected_velocity_ratio(M, gamma),
        )


class TestRealGas:
    """
    Test FrozenShock output against output of the original RGFROSH[^2] for the
    Peng-Robinson EOS for nitrogen (from table in the addendum).

    [^2]: Davidson, D. F. and R. K. Hanson (1996). "Real Gas Corrections in Shock Tube
    Studies at High Pressures." Israel Journal of Chemistry 36(3): 321-326.
    """

    data = [
        (762, 1, 814.9, 19.95, 8.303),
        (762, 10, 809.0, 197.37, 78.143),
        (1000, 1, 1208.4, 45.56, 12.732),
        (1000, 10, 1191.9, 450.47, 116.194),
        (1523, 1, 2338, 145.0, 20.817),
        (1523, 10, 2286, 1432, 181.60),
    ]

    @staticmethod
    @pytest.fixture
    def gas():
        state = CP.AbstractState("PR", "Nitrogen")
        state.specify_phase(CP.iphase_supercritical_gas)
        return CPInterface(state)

    @staticmethod
    @pytest.mark.parametrize("u1,P1,T5,P5,rho5", data)
    def test_reflected_conditions(gas, u1, P1, T5, P5, rho5):
        shock = FrozenShock(gas, u1=u1, T1=293, P1=P1 * 101325)
        assert_allclose(
            [shock.T5, shock.P5 / 101325, shock.rho5],
            [T5, P5, rho5],
            rtol=1e-3,
        )


@pytest.mark.parametrize("gas", [ct.Nitrogen(), ct.CarbonDioxide()])
@pytest.mark.parametrize("u1", [400, 450, 500, 550, 600, 650, 700, 750, 800])
def test_consistency(gas, u1):
    """
    Tests for consistency between the initialization approaches using the following steps:

    1. Initialize a `FrozenShock` object using `M`, `T1`, and `P1`
    2. Initialize another `FrozenShock` object using the calculated `T5` and `P5` from step 1
    3. Check that all properties are consistent between the objects from steps 1 and 2
    """

    gas = ct.Nitrogen()
    from_initial = FrozenShock(gas, u1=u1, T1=300, P1=101325)
    from_target = FrozenShock(
        gas, T5=from_initial.T5, P5=from_initial.P5, T1=from_initial.T1
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

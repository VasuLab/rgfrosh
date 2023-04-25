"""
Tests for the `FrozenShock` class.
"""

from rgfrosh import ThermoInterface, IdealShock, FrozenShock
from rgfrosh.thermo import CPInterface
from rgfrosh.constants import GAS_CONSTANT

import CoolProp as CP

from numpy.testing import assert_almost_equal, assert_allclose
import pytest


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

    test_data = [
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
    @pytest.mark.parametrize("u1,P1,T5,P5,rho5", test_data)
    def test_reflected_conditions(gas, u1, P1, T5, P5, rho5):
        shock = FrozenShock(gas, u1=u1, T1=293, P1=P1 * 101325)
        assert_allclose(
            [shock.T5, shock.P5 / 101325, shock.rho5],
            [T5, P5, rho5],
            rtol=1e-3,
        )

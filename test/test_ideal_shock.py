"""
Tests that the FrozenShock implementation correctly reduces to the ideal shock equations[^1]
for a calorically perfect gas ThermoInterface.

[1]: Gaydon, A. G. and I. R. Hurle (1963). The shock tube in high-temperature chemical physics.
"""

from rgfrosh import ThermoInterface, FrozenShock
from rgfrosh.constants import GAS_CONSTANT
from numpy.testing import assert_almost_equal
import pytest


class PerfectGas(ThermoInterface):
    """Thermo interface for a calorically perfect gas."""

    def __init__(self, gamma, MW):
        self._T = 300
        self._P = 101325
        self._gamma = gamma
        self._MW = MW

    @property
    def TP(self):
        return self._T, self._P

    @TP.setter
    def TP(self, value):
        self._T, self._P = value

    @property
    def mean_molecular_weight(self):
        return self._MW

    @property
    def density_mass(self):
        return self._P / (GAS_CONSTANT / self._MW * self._T)

    @property
    def cp_mass(self):
        return self._gamma / (self._gamma - 1) * GAS_CONSTANT / self._MW

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
        return (self._gamma * GAS_CONSTANT / self._MW * self._T) ** 0.5


def incident_pressure_ratio(M, gamma):
    return (2 * gamma * M**2 - (gamma - 1)) / (gamma + 1)


def incident_density_ratio(M, gamma):
    return (gamma + 1) * M**2 / ((gamma - 1) * M**2 + 2)


@pytest.mark.parametrize("gamma,MW", [(7 / 5, 28), (5 / 3, 40)])  # [Nitrogen, Argon]
@pytest.mark.parametrize("M", [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
class TestClass:
    @pytest.fixture
    def ideal_shock(self, M, gamma, MW):
        perfect_gas = PerfectGas(gamma, MW)
        perfect_gas.TP = 300, 101325
        return FrozenShock(perfect_gas, M * perfect_gas.sound_speed, *perfect_gas.TP)

    def test_incident_temperature_ratio(self, ideal_shock, M, gamma, MW):
        assert_almost_equal(
            ideal_shock.T2 / ideal_shock.T1,
            (gamma * M**2 - (gamma - 1) / 2)
            * ((gamma - 1) / 2 * M**2 + 1)
            / ((gamma + 1) / 2 * M) ** 2,
        )

    def test_incident_pressure_ratio(self, ideal_shock, M, gamma, MW):
        assert_almost_equal(
            ideal_shock.P2 / ideal_shock.P1, incident_pressure_ratio(M, gamma)
        )

    def test_incident_density_ratio(self, ideal_shock, M, gamma, MW):
        assert_almost_equal(
            ideal_shock.rho2 / ideal_shock.rho1,
            incident_density_ratio(M, gamma),
        )

    def test_incident_velocity_ratio(self, ideal_shock, M, gamma, MW):
        assert_almost_equal(
            ideal_shock.u1 / ideal_shock.u2,
            incident_density_ratio(M, gamma),
        )

    def test_reflected_temperature_ratio(self, ideal_shock, M, gamma, MW):
        assert_almost_equal(
            ideal_shock.T5 / ideal_shock.T1,
            (2 * (gamma - 1) * M**2 + 3 - gamma)
            * ((3 * gamma - 1) * M**2 - 2 * (gamma - 1))
            / ((gamma + 1) * M) ** 2,
        )

    def test_reflected_pressure_ratio(self, ideal_shock, M, gamma, MW):
        assert_almost_equal(
            ideal_shock.P5 / ideal_shock.P1,
            incident_pressure_ratio(M, gamma)
            * ((3 * gamma - 1) * M**2 - 2 * (gamma - 1))
            / ((gamma - 1) * M**2 + 2),
        )

    @pytest.mark.xfail
    def test_reflected_velocity_ratio(self, ideal_shock, M, gamma, MW):
        assert_almost_equal(
            ideal_shock.u5 / ideal_shock.u1,
            (2 + 2 / (gamma - 1) / incident_pressure_ratio(M, gamma))
            / ((gamma + 1) / (gamma - 1) - 1 / incident_pressure_ratio(M, gamma)),
        )

from .constants import GAS_CONSTANT

from abc import abstractmethod
from typing import Protocol, Tuple


class ThermoInterface(Protocol):
    """
    Class defining the required interface for calculating mixture thermodynamic
    properties.
    """

    @property
    @abstractmethod
    def TP(self) -> Tuple[float, float]:
        """Get/set temperature [K] and pressure [Pa]."""

    @TP.setter
    @abstractmethod
    def TP(self, value: Tuple[float, float]):
        """Get/set temperature [K] and pressure [Pa]."""

    @property
    @abstractmethod
    def mean_molecular_weight(self) -> float:
        """The mean molecular weight (molar mass) [kg/kmol]."""

    @property
    @abstractmethod
    def density_mass(self) -> float:
        """(Mass) density [kg/m^3^]."""

    @property
    @abstractmethod
    def cp_mass(self) -> float:
        """Specific heat capacity at constant pressure [J/kg/K]."""

    @property
    @abstractmethod
    def enthalpy_mass(self) -> float:
        """Specific enthalpy [J/kg]."""

    @property
    @abstractmethod
    def isothermal_compressibility(self) -> float:
        """Isothermal compressibility [1/Pa]."""

    @property
    @abstractmethod
    def thermal_expansion_coeff(self) -> float:
        """Thermal expansion coefficient [1/K]."""


def compressibility_factor(
    thermo: ThermoInterface, T: float = None, P: float = None
) -> float:
    r"""
    Calculates the compressibility factor of a gas at a specified state:

    $$
    Z = \frac{P}{\rho RT}
    $$

    Parameters:
        thermo: Thermodynamic interface.
        T: Temperature [K].
        P: Pressure [Pa].

    """
    if T and P:
        thermo.TP = T, P
    else:
        T, P = thermo.TP

    return P / (thermo.density_mass * GAS_CONSTANT / thermo.mean_molecular_weight * T)


try:
    import CoolProp as CP

    class CPInterface(ThermoInterface):
        """
        :octicons-tag-24: 0.1.4

        !!! Note
            Only available if `CoolProp` is installed.

        Wrapper for `CoolProp.AbstractState` objects that accounts for
        differing property names and automatically adjusts units for
        `mean_molecular_weight`.
        """

        def __init__(self, state: CP.AbstractState):
            self.state = state

        @property
        def mean_molecular_weight(self):
            return self.state.molar_mass() * 1e3  # [kg/mol] to [kg/kmol]

        @property
        def TP(self):
            return self.state.T(), self.state.p()

        @TP.setter
        def TP(self, value):
            T, P = value
            self.state.update(CP.PT_INPUTS, P, T)

        @property
        def density_mass(self):
            return self.state.rhomass()

        @property
        def cp_mass(self):
            return self.state.cpmass()

        @property
        def enthalpy_mass(self):
            return self.state.hmass()

        @property
        def isothermal_compressibility(self):
            return self.state.isothermal_compressibility()

        @property
        def thermal_expansion_coeff(self):
            return self.state.isobaric_expansion_coefficient()

except ImportError:
    pass

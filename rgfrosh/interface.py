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
        raise NotImplementedError

    @TP.setter
    @abstractmethod
    def TP(self, value: Tuple[float, float]):
        raise NotImplementedError

    @property
    @abstractmethod
    def mean_molecular_weight(self) -> float:
        """The mean molecular weight (molar mass) [kg/kmol]."""
        raise NotImplementedError

    @property
    @abstractmethod
    def density_mass(self) -> float:
        """(Mass) density [kg/m^3^]."""
        raise NotImplementedError

    @property
    @abstractmethod
    def cp_mass(self) -> float:
        """Specific heat capacity at constant pressure [J/kg/K]."""
        raise NotImplementedError

    @property
    @abstractmethod
    def enthalpy_mass(self) -> float:
        """Specific enthalpy [J/kg]."""
        raise NotImplementedError

    @property
    @abstractmethod
    def isothermal_compressibility(self) -> float:
        """Isothermal compressibility [1/Pa]."""
        raise NotImplementedError

    @property
    @abstractmethod
    def thermal_expansion_coeff(self) -> float:
        """Thermal expansion coefficient [1/K]."""
        raise NotImplementedError


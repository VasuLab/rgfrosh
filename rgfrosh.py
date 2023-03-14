__version__ = "0.1.3"

from abc import abstractmethod
from dataclasses import dataclass
from typing import Tuple, Protocol

import numpy as np
from tabulate import tabulate


GAS_CONSTANT: float = 8314.46261815324
"""Universal gas constant [J/kmol/K]."""

max_iter: int = 1000
"""Maximum number of iterations for the solver."""

rtol: float = 1E-6
"""Relative tolerance for the solver convergence criteria."""


class ConvergenceError(Exception):
    """
    Exception raised when the iterative solver fails to converge.
    """


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


def compressibility_factor(thermo: ThermoInterface, T: float = None, P: float = None):
    r"""
    Calculates the compressibility factor of a gas at a specified state:

    $$
    Z = \frac{P}{\rho R_{specific} T}
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


@dataclass(init=False)
class FrozenShock:
    """
    Dataclass with fields for relevant properties - velocity, temperature, pressure, and
    density - for each region associated with a reflecting shock.

    !!! Warning
        The state of the `thermo` argument is modified when the `FrozenShock` object is initialized and
        whenever the thermodynamic interface at a given state ([`state1`][rgfrosh.FrozenShock.state1],
        [`state2`][rgfrosh.FrozenShock.state2], [`state5`][rgfrosh.FrozenShock.state5]) is accessed.

    """

    u1: float
    """Incident shock velocity [m/s]."""
    T1: float
    """Initial temperature [K]."""
    P1: float
    """Initial pressure [Pa]."""
    rho1: float
    """Initial density [kg/m^3^]."""

    u2: float
    """Velocity behind the incident shock (shock-fixed) [m/s]."""
    U2: float
    """Velocity behind the incident shock (lab frame) [m/s]."""
    T2: float
    """Temperature behind the incident shock [K]."""
    P2: float
    """Pressure behind the incident shock [Pa]."""
    rho2: float
    """Density behind the incident shock [kg/m^3^]."""

    u5: float
    """Reflected shock velocity [m/s]."""
    T5: float
    """Temperature behind the reflected shock [K]."""
    P5: float
    """Pressure behind the reflected shock [Pa]."""
    rho5: float
    """Density behind the reflected shock [kg/m^3^]."""

    def __init__(
        self,
        thermo: ThermoInterface,
        u1: float,
        T1: float,
        P1: float,
    ):
        """
        Solves for the post-incident-shock and post-reflected-shock conditions
        given the initial conditions and incident shock velocity.

        Parameters:
            thermo: Thermodynamic interface.
            u1: Incident shock velocity [m/s].
            T1: Initial temperature [K].
            P1: Initial pressure [Pa].

        """

        self.thermo = thermo
        self.thermo.TP = T1, P1

        self.u1 = u1
        self.P1 = P1
        self.T1 = T1
        self.rho1 = thermo.density_mass

        self.u2, self.T2, self.P2, self.rho2 = FrozenShock.incident_conditions(
            self.thermo, u1, T1, P1)
        self.U2 = self.u1 - self.u2

        self.u5, self.T5, self.P5, self.rho5 = FrozenShock.reflected_conditions(
            thermo, u1, P1, self.u2, self.T2, self.P2)

    def __str__(self):
        return tabulate(
            [["1", self.u1, self.T1, self.P1, self.rho1],
             ["2", self.u2, self.T2, self.P2, self.rho2],
             ["5", self.u5, self.T5, self.P5, self.rho5]],
            headers=["State", "Velocity [m/s]", "Temperature [K]", "Pressure [Pa]", "Density [kg/m\u00b3]"],
            tablefmt="rounded_outline",
            floatfmt=("", ".1f", ".1f", ".3e", ".4g"),
        )

    @property
    def state1(self):
        """Thermodynamic interface at initial conditions."""
        self.thermo.TP = self.T1, self.P1
        return self.thermo

    @property
    def state2(self):
        """Thermodynamic interface at post-incident-shock conditions."""
        self.thermo.TP = self.T2, self.P2
        return self.thermo

    @property
    def state5(self):
        """Thermodynamic interface at post-reflected-shock conditions."""
        self.thermo.TP = self.T5, self.P5
        return self.thermo

    @staticmethod
    def incident_conditions(
        thermo: ThermoInterface,
        u1: float,
        T1: float,
        P1: float,
    ) -> Tuple[float, float, float, float]:
        """
        Solves for the conditions behind the incident shock.

        Parameters:
            thermo: Thermodynamic interface.
            u1: Incident shock velocity [m/s].
            T1: Initial temperature [K].
            P1: Initial pressure [Pa].

        Returns:
            u2: Velocity behind the incident shock (shock-fixed) [m/s].
            T2: Temperature behind the incident shock [K].
            P2: Pressure behind the incident shock [Pa].
            rho2: Density behind the incident shock [kg/m^3^].

        Exceptions:
            ConvergenceError: If the relative change in `T5` and `P5` is not below the
                [`rtol`][rgfrosh.rtol] within
                [`max_iter`][rgfrosh.max_iter] iterations.

        """

        thermo.TP = T1, P1
        h1 = thermo.enthalpy_mass
        nu1 = 1 / thermo.density_mass
        cp1 = thermo.cp_mass

        R = GAS_CONSTANT / thermo.mean_molecular_weight
        gamma1 = cp1 / (cp1 - R)

        # Calculate an initial guess of P2 and T2 with the ideal gas assumption
        M1 = u1 / (gamma1 * R * T1) ** 0.5
        T2 = T1 * (gamma1 * M1 ** 2 - (gamma1 - 1) / 2) * \
             ((gamma1 - 1) / 2 * M1 ** 2 + 1) / ((gamma1 + 1) / 2 * M1) ** 2
        P2 = P1 * (2 * gamma1 * M1 ** 2 - (gamma1 - 1)) / (gamma1 + 1)

        for i in range(max_iter):
            # Calculate thermodynamic properties at T2, P2 guess
            thermo.TP = T2, P2
            h2 = thermo.enthalpy_mass
            nu2 = 1 / thermo.density_mass
            cp2 = thermo.cp_mass
            beta2 = thermo.thermal_expansion_coeff
            kappa2 = thermo.isothermal_compressibility

            f1 = (P2 / P1 - 1) + (u1 ** 2 / (P1 * nu1)) * (nu2 / nu1 - 1)
            f2 = ((h2 - h1) / (1 / 2 * u1 ** 2)) + (nu2 ** 2 / nu1 ** 2 - 1)

            df1_dT2 = u1 ** 2 * nu2 * beta2 / (P1 * nu1 ** 2)
            df1_dP2 = 1 / P1 - u1 ** 2 * nu2 * kappa2 / (P1 * nu1 ** 2)

            df2_dT2 = 2 * cp2 / u1 ** 2 + 2 * (nu2 / nu1) ** 2 * beta2
            df2_dP2 = 2 * nu2 * (1 - T2 * beta2) / u1 ** 2 - 2 * nu2 ** 2 * kappa2 / nu1 ** 2

            deltaT2, deltaP2 = np.matmul(
                np.linalg.inv(np.array([
                    [df1_dT2, df1_dP2],
                    [df2_dT2, df2_dP2]
                ])),
                np.array([f1, f2]),
            )

            converged = (
                abs(deltaT2) <= T2 * rtol and
                abs(deltaP2) <= P2 * rtol
            )

            T2 -= deltaT2
            P2 -= deltaP2

            if converged:
                thermo.TP = T2, P2
                nu2 = 1 / thermo.density_mass
                u2 = u1 * nu2 / nu1
                return u2, T2, P2, thermo.density_mass

        raise ConvergenceError

    @staticmethod
    def reflected_conditions(
        thermo: ThermoInterface,
        u1: float,
        P1: float,
        u2: float,
        T2: float,
        P2: float,
    ) -> Tuple[float, float, float, float]:
        """
        Solves for the conditions behind the reflected shock.

        Parameters:
            thermo: Thermodynamic interface.
            u1: Incident shock velocity [m/s].
            P1: Initial pressure [Pa].
            u2: Velocity behind the incident shock (shock-fixed) [m/s].
            T2: Temperature behind the incident shock [K].
            P2: Pressure behind the incident shock [Pa].

        Returns:
            u5: Reflected shock velocity [m/s].
            T5: Temperature behind the reflected shock [K].
            P5: Pressure behind the reflected shock [Pa].
            rho5: Density behind the reflected shock [kg/m^3^].

        Exceptions:
            ConvergenceError: If the relative change in `T5` and `P5` is not below the
                [`rtol`][rgfrosh.rtol] within
                [`max_iter`][rgfrosh.max_iter] iterations.

        """

        thermo.TP = T2, P2
        h2 = thermo.enthalpy_mass
        nu2 = 1 / thermo.density_mass

        # Calculate an initial guess of P5 and T5 with the ideal gas assumption
        cp2 = thermo.cp_mass
        R = GAS_CONSTANT / thermo.mean_molecular_weight
        gamma2 = cp2 / (cp2 - R)

        eta2 = (gamma2 + 1) / (gamma2 - 1)
        P5 = P2 * (eta2 + 2 - P1 / P2) / (1 + eta2 * P1 / P2)
        T5 = T2 * P5 / P2 * (eta2 + P5 / P2) / (1 + eta2 * P5 / P2)

        for i in range(max_iter):
            thermo.TP = T5, P5

            h5 = thermo.enthalpy_mass
            nu5 = 1 / thermo.density_mass
            cp5 = thermo.cp_mass
            beta5 = thermo.thermal_expansion_coeff
            kappa5 = thermo.isothermal_compressibility

            f3 = (P5 / P2 - 1) + (u1 - u2) ** 2 / P2 / (nu5 - nu2)
            f4 = 2 * (h5 - h2) / (u1 - u2) ** 2 + (nu5 + nu2) / (nu5 - nu2)

            df3_dT5 = -nu5 * beta5 * (u1 - u2) ** 2 / (P2 * (nu5 - nu2) ** 2)
            df3_dP5 = 1 / P2 + nu5 * kappa5 * (u1 - u2) ** 2 / (P2 * (nu5 - nu2) ** 2)

            df4_dT5 = 2 * cp5 / (u1 - u2) ** 2 \
                      - 2 * nu2 * nu5 * beta5 / (nu5 - nu2) ** 2
            df4_dP5 = (
                2 * nu5 * (1 - T5 * beta5) / (u1 - u2) ** 2
                + 2 * nu2 * nu5 * kappa5 / (nu5 - nu2) ** 2
            )

            deltaT5, deltaP5 = np.matmul(
                np.linalg.inv(np.array([
                    [df3_dT5, df3_dP5],
                    [df4_dT5, df4_dP5]
                ])),
                np.array([f3, f4]),
            )

            converged = (
                abs(deltaT5) <= T5 * rtol and
                abs(deltaP5) <= P5 * rtol
            )

            T5 -= deltaT5
            P5 -= deltaP5

            if converged:
                thermo.TP = T5, P5
                nu5 = 1 / thermo.density_mass
                u5 = (u1 - u2) / (nu2 / nu5 - 1)
                return u5, T5, P5, thermo.density_mass

        raise ConvergenceError

    @classmethod
    def target_conditions(
        cls,
        thermo: ThermoInterface,
        T5: float,
        P5: float,
        T1: float = 300
    ):
        """
        Solves for the initial pressure and incident shock velocity given the target
        post-reflected-shock conditions.

        Parameters:
            thermo: Thermodynamic interface.
            T5: Temperature behind the reflected shock [K].
            P5: Pressure behind the reflected shock [Pa].
            T1: Initial temperature [K].

        Exceptions:
            ConvergenceError: If the relative change in `u1`, `P1`, `T2`, and `P2`
                is not below the [`rtol`][rgfrosh.rtol]
                within [`max_iter`][rgfrosh.max_iter] iterations.

        """

        thermo.TP = T5, P5
        h5 = thermo.enthalpy_mass
        nu5 = 1 / thermo.density_mass

        thermo.TP = T1, 101325
        cp1 = thermo.cp_mass
        R = GAS_CONSTANT / thermo.mean_molecular_weight
        gamma1 = cp1 / (cp1 - R)
        a1 = (gamma1 * GAS_CONSTANT / thermo.mean_molecular_weight * T1) ** 0.5

        a = 2 * (gamma1 - 1) * (3 * gamma1 - 1)
        b = (3 * gamma1 - 1) * (3 - gamma1) - 4 * (gamma1 - 1) ** 2 - (
                gamma1 + 1) ** 2 * T5 / T1
        c = -2 * (gamma1 - 1) * (3 - gamma1)

        MS = ((-b + (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)) ** 0.5

        u1 = MS * a1
        P1 = P5 * (gamma1 + 1) / (2 * gamma1 * MS ** 2 - (gamma1 - 1)) * (
                (gamma1 - 1) * MS ** 2 + 2) / (
                     (3 * gamma1 - 1) * MS ** 2 - 2 * (gamma1 - 1))

        P2 = P1 * (2 * gamma1 * MS ** 2 - (gamma1 - 1)) / (gamma1 + 1)
        T2 = T1 * ((gamma1 * MS ** 2 - (gamma1 - 1) / 2) *
                   ((gamma1 - 1) / 2 * MS ** 2 + 1) / ((gamma1 + 1) / 2 * MS) ** 2)

        for i in range(max_iter):
            thermo.TP = T1, P1
            h1 = thermo.enthalpy_mass
            nu1 = 1 / thermo.density_mass
            beta1 = thermo.thermal_expansion_coeff
            kappa1 = thermo.isothermal_compressibility

            thermo.TP = T2, P2
            h2 = thermo.enthalpy_mass
            nu2 = 1 / thermo.density_mass
            cp2 = thermo.cp_mass
            beta2 = thermo.thermal_expansion_coeff
            kappa2 = thermo.isothermal_compressibility

            f1 = P2 / P1 - 1 + u1 ** 2 / (P1 * nu1) * (nu2 / nu1 - 1)
            f2 = 2 * (h2 - h1) / u1 ** 2 + (nu2 / nu1) ** 2 - 1
            f3 = (P5 / P2 - 1) + (u1 * (1 - nu2 / nu1)) ** 2 / P2 / (nu5 - nu2)
            f4 = 2 * (h5 - h2) / (u1 * (1 - nu2 / nu1)) ** 2 + (nu5 + nu2) / (nu5 - nu2)

            df1_du1 = 2 * u1 / (P1 * nu1 ** 2) * (nu2 - nu1)
            df1_dP1 = -P2 / P1 ** 2 + u1 ** 2 / (P1 ** 2 * nu1 ** 2) \
                      * (P1 * kappa1 * (2 * nu2 - nu1) - (nu2 - nu1))
            df1_dT2 = u1 ** 2 * nu2 * beta2 / (P1 * nu1 ** 2)
            df1_dP2 = 1 / P1 - u1 ** 2 * nu2 * kappa2 / (P1 * nu1 ** 2)

            df2_du1 = -4 * (h2 - h1) / u1 ** 3
            df2_dP1 = 2 * nu2 ** 2 * kappa1 / nu1 ** 2 - 2 * nu1 * (
                    1 - T1 * beta1) / u1 ** 2
            df2_dT2 = 2 * cp2 / u1 ** 2 + 2 * (nu2 / nu1) ** 2 * beta2
            df2_dP2 = 2 * nu2 * (1 - T2 * beta2) / u1 ** 2 - 2 * nu2 ** 2 * kappa2 / nu1 ** 2

            df3_du1 = 2 * u1 * (1 - nu2 / nu1) ** 2 / (P2 * (nu5 - nu2))
            df3_dP1 = -2 * u1 ** 2 * nu2 * kappa1 * (1 - nu2 / nu1) / (
                    P2 * nu1 * (nu5 - nu2))
            df3_dT2 = u1 ** 2 * nu2 * beta2 / (P2 * nu1 ** 2) * \
                      (nu1 - nu2) / (nu5 - nu2) ** 2 * (nu1 + nu2 - 2 * nu5)
            df3_dP2 = -P5 / P2 ** 2 - u1 ** 2 / (nu1 ** 2 * P2) * (nu1 - nu2) / (nu5 - nu2) \
                      * ((nu2 * kappa2 * (nu1 + nu2 - 2 * nu5)) / (nu5 - nu2) + (nu1 - nu2) / P2)

            df4_du1 = -4 * (h5 - h2) / (u1 ** 3 * (1 - nu2 / nu1) ** 2)
            df4_dP1 = 4 * nu2 * kappa1 * (h5 - h2) / (u1 ** 2 * nu1) / (
                    1 - nu2 / nu1) ** 3
            df4_dT2 = 2 * nu1 ** 2 / u1 ** 2 * (2 * nu2 * beta2 * (h5 - h2) - cp2 * (nu1 - nu2)) \
                      / (nu1 - nu2) ** 3 + 2 * nu2 * nu5 * beta2 / (nu5 - nu2) ** 2
            df4_dP2 = -2 * nu1 ** 2 * nu2 / u1 ** 2 * (2 * kappa2 * (h5 - h2) + (nu1 - nu2) * (1 - T2 * beta2)) \
                      / (nu1 - nu2) ** 3 - 2 * nu2 * nu5 * kappa2 / (nu5 - nu2) ** 2

            delta_u1, delta_P1, delta_T2, delta_P2 = np.matmul(
                np.linalg.inv(np.array([
                    [df1_du1, df1_dP1, df1_dT2, df1_dP2],
                    [df2_du1, df2_dP1, df2_dT2, df2_dP2],
                    [df3_du1, df3_dP1, df3_dT2, df3_dP2],
                    [df4_du1, df4_dP1, df4_dT2, df4_dP2]
                ])),
                np.array([f1, f2, f3, f4]),
            )

            converged = (
                abs(delta_u1) <= u1 * rtol and
                abs(delta_P1) <= P1 * rtol and
                abs(delta_T2) <= T2 * rtol and
                abs(delta_P2) <= P2 * rtol
            )

            u1 -= delta_u1
            P1 -= delta_P1
            T2 -= delta_T2
            P2 -= delta_P2

            if converged:
                return cls(thermo, u1, T1, P1)

        raise ConvergenceError

---
title: 'rgfrosh: A Python frozen shock solver for ideal and real gas equations of state.'
tags:
  - Python
  - shock tube
  - normal shock
  - frozen shock
  - ideal gas
  - real gas
authors:
  - name: Cory Kinney
    orcid: 0000-0001-8774-3090
    affiliation: 1
  - name: Subith Vasu
    affiliation: 1

affiliations:
 - name: Department of Mechanical and Aerospace Engineering, University of Central Florida, Orlando, FL USA
   index: 1
date: 26 April 2024
bibliography: paper.bib
---

# Summary

Combustion is a complex process that is influenced by numerous factors, including temperature, pressure, and 
chemical composition. Combustion modeling relies on chemical kinetics mechanisms that detail how chemical reaction rates
involving many intermediate species vary with temperature and pressure. Experimental data is essential for measuring 
reaction rates, refining mechanisms, and validating predictions from these models for certain fuels and conditions.

Shock tubes are ideal experimental facilities for performing fundamental research in combustion because they use shock 
waves to impart step changes in temperature and pressure to a test gas mixture. The gas behind the reflected shock is
relatively stagnant, and the system is considered adiabatic on the relevant time-scales; therefore, the shock tube 
closely resembles a constant volume zero-dimensional reactor, allowing for comparison between experimental measurements 
and model predictions. Accurate calculation of the reflected shock conditions is essential for correctly interpreting 
experimental results. This is accomplished by solving the reflected shock equations under an assumed equation of state 
using the known initial state of the gas and the measured shock velocity.

Fundamental research into combustion at extremely high-pressure conditions, such as for rocket engines or direct-fire 
supercritical CO2 power cycles [@kinney:2022], requires consideration of real gas effects in experimental measurements 
and model simulations. A shock solver supporting real gas equations of state called RGFROSH was previously developed in 
FORTRAN by @davidson:1996; however, there existed a need for a modernized open-source version, as the original is not 
generally available nor would it be readily compatible with modern tools. Thus, the present `rgfrosh` was developed in 
Python as a solution that uses modern thermodynamic libraries to enable further research into high-pressure combustion.

# Features

`rgfrosh` provides a simple interface for solving the reflected shock equations. Two models are provided - `IdealShock` 
for calorically perfect gases, and `FrozenShock` for arbitrary equations of state. The former is primarily included for 
comparison and validation purposes, while the latter is the primary focus of the package as it allows for accurate 
calculation of the reflected shock conditions for real gases. 

To remain as lightweight and extensible as possible, `rgfrosh` relies on external packages for the key thermodynamic 
functions required by the `FrozenShock` solver. The required interface is defined by the `ThermoInterface` protocol 
class, which was written to provide native support for Cantera [@cantera]. Additionally, an interface is provided to 
wrap CoolProp [@coolprop], which itself has backend support for NIST REFPROP [@refprop], for use with the solver. 
These two compatible packages enable support for a wide range of equations of state which should cover the majority of 
use cases; however, any user-defined class that implements the simple `ThermoInterface` protocol can be used with the 
solver.

The primary use cases for `rgfrosh` are experiment postprocessing and experiment planning. The `solve_incident` and 
`solve_reflected` methods implement the Newton-Raphson solver detailed by @davidson:1996 for calculating the incident 
and reflected shock conditions, respectively, from the initial conditions and the experimentally measured shock
velocity. The `solve_initial` method implements the algorithm derived by the author [see @kinney:2022{A.3}] for 
calculating the initial pressure and incident shock velocity, temperature, and pressure from the initial temperature and
target reflected shock temperature and pressure for an experiment.

# Future work

Future work includes the consideration of vibrational non-equilibrium in the shock solvers. Current
solver functionality would be classified as equilibrium-equilibrium (EE) mode, referring to the incident 
and reflected shock, respectively; frozen-equilibrium and frozen-frozen modes are planned.

# Acknowledgements

We would like to acknowledge D. F. Davidson and R. K. Hanson for authoring the original software this 
work is based on, granting permission to use the package name, and for providing validation data for comparison.

# References

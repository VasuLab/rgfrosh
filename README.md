# PyRGFROSH

PyRGFROSH is a **Py**thon implementation of the **R**eal **G**as **FRO**zen **SH**ock 
equations[^1] solver, which calculates the conditions behind the incident and reflected 
shocks in a shock tube for an arbitrary mixture and equation of state. PyRGFROSH 
requires a thermodynamic interface for calculating mixture properties as a function of 
temperature and pressure and currently supports:

- [Cantera](https://github.com/cantera/cantera)
- [CoolProp](https://github.com/CoolProp/CoolProp) (see [example](https://VasuLab.github.io/PyRGFROSH/guide/#coolprop)) 
- [User-defined interfaces](https://VasuLab.github.io/PyRGFROSH/guide/#user-defined-interfaces)

The original RGFROSH was developed in FORTRAN at Stanford University by D. F. Davidson 
and R. K. Hanson using real gas subroutines for CHEMKIN[^2][^3].

## Documentation

PyRGFROSH's [documentation](https://VasuLab.github.io/PyRGFROSH) provides a detailed 
[API reference](https://VasuLab.github.io/PyRGFROSH/reference) for the package.

## Installation

PyRGFROSH can be installed using

```
pip install rgfrosh
```


[^1]: Davidson, D.F. and Hanson, R.K. (1996), Real Gas Corrections in Shock Tube Studies 
at High Pressures. Isr. J. Chem., 36: 321-326. https://doi.org/10.1002/ijch.199600044
[^2]: P. Barry Butler, "Real Gas Equations of State for Chemkin" Sandia Report No. 
SAND88-3188 (1988). https://doi.org/10.2172/6224858
[^3]: R. G. Schmitt, P. B. Butler, N. B. French "Chemkin real gas: a Fortran package for 
analaysis of thermodynamic properties and chemical kinetics in non-ideal systems," 
U. of Iowa Report UIME PPB 93-006 (1994).

# Python Real Gas FROzen SHock (RGFROSH) Solver

!!! cite ""
    
    This project is a solver for the frozen shock equations[^1] developed in Python at the
    University of Central Florida. The original RGFROSH was developed in FORTRAN at Stanford 
    University by D. F. Davidson and R. K. Hanson using real gas subroutines for 
    CHEMKIN[^2]^,^[^3]. 
    

Python RGFROSH (`rgfrosh`) is a package for calculating conditions behind incident and reflected shock in
a shock tube for an arbitrary equation of state. RGFROSH requires a thermodynamic interface 
for calculating mixture properties as a function of temperature and pressure and currently supports:

- [Cantera](https://github.com/cantera/cantera) (natively)
- [CoolProp](https://github.com/CoolProp/CoolProp) (see [example](guide/#coolprop-example)) 
- [User-defined interfaces](guide/#user-defined-interfaces)



## Installation

RGFROSH can be installed using

```
pip install rgfrosh
```

which also installs required dependencies.

[^1]: Davidson, D.F. and Hanson, R.K. (1996), Real Gas Corrections in Shock Tube Studies 
at High Pressures. Isr. J. Chem., 36: 321-326. 
[https://doi.org/10.1002/ijch.199600044](https://doi.org/10.1002/ijch.199600044)
[^2]: P. Barry Butler, "Real Gas Equations of State for Chemkin" Sandia Report No. 
SAND88-3188 (1988). [https://doi.org/10.2172/6224858](https://doi.org/10.2172/6224858)
[^3]: R. G. Schmitt, P. B. Butler, N. B. French "Chemkin real gas: a Fortran package for 
analaysis of thermodynamic properties and chemical kinetics in non-ideal systems," 
U. of Iowa Report UIME PPB 93-006 (1994).

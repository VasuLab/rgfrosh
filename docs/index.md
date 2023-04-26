# Python Real Gas FROzen SHock (RGFROSH)

!!! cite ""
    
    This project is a solver for the frozen shock equations[^1] developed in Python at the
    University of Central Florida. The original RGFROSH was developed in FORTRAN at Stanford 
    University by D. F. Davidson and R. K. Hanson using real gas subroutines for 
    CHEMKIN[^2]^,^[^3]. 
    

Python RGFROSH (`rgfrosh`) is a package for calculating conditions behind incident and reflected shock in
a shock tube for an arbitrary equation of state. See the [user guide](guide/getting-started) for 
information on how to get started with `rgfrosh`.

## Installation

Python RGFROSH can be installed using

```
pip install rgfrosh
```

which also installs required dependencies. Cantera or CoolProp are optional and must 
be installed separately if desired.

## Contributing

For any bugs or feature requests, create an issue on the 
[issue tracker](https://github.com/VasuLab/RGFROSH/issues). 

After cloning the repository, the development environment can be set up with

```
pip install -r requirements.txt
```

Before creating a pull request, be sure to lint

```
black .
```

and run the automated tests

```
pytest
```

These checks will be performed automatically for all pull requests along
with test coverage comparisons.

## Cite

To cite Python RGFROSH go to the [GitHub repository](https://github.com/VasuLab/RGFROSH) and click 
<kbd>Cite this repository</kbd> on the right side of the page to export the citation for the latest 
release of `rgfrosh`. It is also encouraged to cite the original paper[^1] for the frozen shock 
equations that this work is a derivative of.


[^1]: Davidson, D.F. and Hanson, R.K. (1996), Real Gas Corrections in Shock Tube Studies 
at High Pressures. Isr. J. Chem., 36: 321-326. 
[https://doi.org/10.1002/ijch.199600044](https://doi.org/10.1002/ijch.199600044)
[^2]: P. Barry Butler, "Real Gas Equations of State for Chemkin" Sandia Report No. 
SAND88-3188 (1988). [https://doi.org/10.2172/6224858](https://doi.org/10.2172/6224858)
[^3]: R. G. Schmitt, P. B. Butler, N. B. French "Chemkin real gas: a Fortran package for 
analaysis of thermodynamic properties and chemical kinetics in non-ideal systems," 
U. of Iowa Report UIME PPB 93-006 (1994).

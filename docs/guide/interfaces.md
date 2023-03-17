### Cantera

`rgfrosh` natively supports `cantera.ThermoPhase` objects, as shown in the
[examples](../../guide/#1-experiment-analysis).

!!! Tip "Notice"
    The implementations of the Redlich-Kwong (RK) EOS and Peng-Robinson (PR) EOS
    in Cantera do not yet include definitions for `isothermal_compressibility` and 
    `thermal_expansion_coeff` 
    (see [Cantera/cantera/pull/1421](https://github.com/Cantera/cantera/pull/1421)), thus
    only the ideal gas EOS is currently supported.

### CoolProp

The [`CPInterface`](../../reference/interface/#rgfrosh.interface.CPInterface) wrapper is 
provided to enable straightforward use of `CoolProp.AbstractState` objects with `rgfrosh`: 

```python hl_lines="7"
from rgfrosh import FrozenShock
from rgfrosh.interface import CPInterface
import CoolProp as CP

state = CP.AbstractState("PR", "Argon")
state.specify_phase(CP.iphase_supercritical_gas)  # (1)!
shock = FrozenShock(CPInterface(state), 763, 293, 10 * 101325)
print(shock)
```

1. Specifying phase is necessary if the cubic has three roots

which prints the following results

```
┌─────────┬──────────────────┬───────────────────┬─────────────────┬───────────────────┐
│   State │   Velocity [m/s] │   Temperature [K] │   Pressure [Pa] │   Density [kg/m³] │
├─────────┼──────────────────┼───────────────────┼─────────────────┼───────────────────┤
│       1 │            763.0 │             293.0 │       1.013e+06 │             16.78 │
│       2 │            295.7 │             764.6 │       6.996e+06 │             43.3  │
│       5 │            465.7 │            1372.9 │       2.587e+07 │             86.74 │
└─────────┴──────────────────┴───────────────────┴─────────────────┴───────────────────┘
```


### User-Defined Interfaces

Custom thermodynamic routines are supported using the 
[`ThermoInterface`][rgfrosh.ThermoInterface] class, which utilizes structural 
subtyping (see [PEP 544](https://www.python.org/dev/peps/pep-0544/)) through the 
built-in [`typing.Protocol`](https://docs.python.org/3/library/typing.html#typing.Protocol)
base class. Thus, any class that implements the required methods is considered a 
subtype, even if it is not an explicit subclass of `ThermoInterface`; however, 
explicitly subclassing `ThermoInterface` is recommended for user-defined 
interfaces to ensure all required methods are defined.


*[EOS]: Equation of State

## Cantera

`rgfrosh` has native support for `cantera.ThermoPhase` instances.

!!! Tip "Notice"
    The implementations of `isothermal_compressibility` and `thermal_expansion_coeff` for the
    Redlich-Kwong (RK) EOS and Peng-Robinson (PR) EOS in Cantera were added in 
    [`931de24`](https://github.com/Cantera/cantera/commit/931de2450e34226a431846a07d7db075939cd7ec).
    Thus, to use the real gas EOS, you must either wait for the Cantera 3.0 release or build 
    Cantera from source.


## CoolProp

The [`CPInterface`](../../reference/thermo/#rgfrosh.thermo.CPInterface) wrapper is
provided to enable straightforward use of `CoolProp.AbstractState` instances with `rgfrosh`:

```python hl_lines="8"
from rgfrosh import FrozenShock
from rgfrosh.thermo import CPInterface
import CoolProp as CP

state = CP.AbstractState("PR", "Argon")  # Peng-Robinson
state.specify_phase(CP.iphase_supercritical_gas)  # (1)!

shock = FrozenShock(CPInterface(state), u1=763, T1=293, P1=(10 * 101325))
```

1. Specifying phase is necessary if the cubic has three roots

## User-defined Interfaces

If you would like to use `rgfrosh` with thermodynamic routines other than Cantera or 
CoolProp, all you have to do is pass an object that matches the pattern defined by the 
[`ThermoInterface`][rgfrosh.ThermoInterface] class.

!!! Note
    The [`ThermoInterface`][rgfrosh.ThermoInterface] class utilizes structural 
    subtyping (see [PEP 544](https://www.python.org/dev/peps/pep-0544/)) through the 
    built-in [`typing.Protocol`](https://docs.python.org/3/library/typing.html#typing.Protocol)
    base class. Thus, any class that implements the required methods is considered a 
    subtype, even if it is not an explicit subclass of `ThermoInterface`; however, 
    explicitly subclassing `ThermoInterface` is recommended for user-defined 
    interfaces to ensure all required methods are defined.

Below is a simple example implementation of an interface for a calorically perfect gas:

```python
from rgfrosh import ThermoInterface
from rgfrosh.constants import GAS_CONSTANT


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
```

The interface can then be used with [`FrozenShock`](../../reference/shock/frozen):

```python hl_lines="4"
from rgfrosh import FrozenShock

argon = PerfectGas(5/3, 40)
frozen_shock = FrozenShock(argon, u1=750, P1=101325)

print(f"T5 = {frozen_shock.T5:.1f} K, P5 = {frozen_shock.P5 / 101325:.2f} atm")
```

which yields the reflected shock conditions:

```
T5 = 1353.9 K, P5 = 23.60 atm
```

For the case of a calorically perfect gas, the frozen shock solution should simplify to the
ideal shock solution. We can use the custom interface with the 
[`IdealShock.from_thermo`][rgfrosh.IdealShock.from_thermo] classmethod to verify this:

```python hl_lines="3"
from rgfrosh import IdealShock

ideal_shock = IdealShock.from_thermo(argon, u1=750, P1=101325)

print(f"T5 = {frozen_shock.T5:.1f} K, P5 = {frozen_shock.P5 / 101325:.2f} atm")
```

which yields approximately the same solution:

```
T5 = 1353.9 K, P5 = 23.60 atm
```


*[EOS]: Equation of State

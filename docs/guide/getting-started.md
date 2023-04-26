The main function of `rgfrosh` is to solve the conservation equations across an incident
and reflected shock in a shock tube. `rgfrosh` does not itself perform any calculations to evaluate 
the thermodynamic properties of the gas; therefore, an interface to an external set of thermodynamics 
routines is required. The following are currently supported:

- [Cantera](https://github.com/cantera/cantera) (natively)
- [CoolProp](https://github.com/CoolProp/CoolProp) (through the built-in [wrapper][rgfrosh.thermo.CPInterface])
- User-defined interfaces

For more detailed information and examples see the [thermodynamic interfaces](../thermodynamic-interfaces)
page.

### Experiment Analysis

The primary use case of `rgfrosh` is using the [`FrozenShock`](../../reference/shock/frozen)
class can to calculate shock conditions for an experiment from the initial conditions and 
measured shock velocity. For example, using Cantera's built-in carbon dioxide EOS:

```py hl_lines="4"
from rgfrosh import FrozenShock
import cantera as ct

shock = FrozenShock(ct.CarbonDioxide(), u1=1000, P1=101325)  # (1)!

print(f"T2 = {shock.T2:.1f} K, P2 = {shock.P2 / 101325:.2f} atm")
print(f"T5 = {shock.T5:.1f} K, P5 = {shock.P5 / 101325:.2f} atm")
```

1. If not given, `T1` is `300` by default.

which outputs the following incident and reflected shock conditions:

```
T2 = 774.8 K, P2 = 15.82 atm
T5 = 1225.1 K, P5 = 113.33 atm
```

???+ Note
    The [`IdealShock`](../../reference/shock/ideal) class implements the exact same constructor
    pattern, except it requires the positional arguments `gamma` and `MW` instead of a 
    `ThermoInterface` object. Therefore, the same calculation would look like:

    ```python hl_lines="3"
    from rgfrosh import IdealShock

    shock = IdealShock(1.29, 44, u1=1000, P1=101325)  # Carbon Dioxide

    print(f"T2 = {shock.T2:.1f} K, P2 = {shock.P2 / 101325:.2f} atm")
    print(f"T5 = {shock.T5:.1f} K, P5 = {shock.P5 / 101325:.2f} atm")
    ```

    which outputs:

    ```
    T2 = 873.2 K, P2 = 15.28 atm
    T5 = 1559.5 K, P5 = 99.03 atm
    ```

    Comparison with the frozen shock solution for the same example demonstrates 
    the unsuitability of the ideal shock equations for certain conditions.

### Experiment Planning
    
The other use case of `rgfrosh` is to plan experiments by calculating the required 
initial conditions for target  reflected shock conditions:

```py hl_lines="4"
from rgfrosh import FrozenShock
import cantera as ct

shock = FrozenShock(ct.CarbonDioxide(), T5=1100, P5=200e5, T1=295)

print(f"P1 = {shock.P1 / 133.322:.0f} torr, u1 = {shock.u1:.1f} m/s")
```

which outputs the required fill pressure and incident shock velocity:

```
P1 = 1640 torr, u1 = 920.9 m/s
```

*[EOS]: Equation of State

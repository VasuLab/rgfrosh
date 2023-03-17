# User Guide

??? Info "Notation"

    The standard shock tube region notation is observed in `rgfrosh`:
    
    | Region | Description                           |
    | ------ | ------------------------------------- |
    | 1      | Initial driven gas state              |
    | 2      | Post-incident-shock driven gas state  | 
    | 3      | Expanded driver gas state             | 
    | 4      | Initial driver gas state              |
    | 5      | Post-reflected-shock driven gas state |

RGFROSH has two primary features:

### 1. Experiment Analysis

The [`FrozenShock`][rgfrosh.shock.FrozenShock] class can be used to calculate shock conditions 
for an experiment from the initial conditions and measured shock velocity. For example:

```py hl_lines="4"
from rgfrosh import FrozenShock
import cantera as ct

shock = FrozenShock(ct.CarbonDioxide(), 1000, 295, 101325)
print(shock)
```

which outputs the following table of conditions:

```commandline
┌─────────┬──────────────────┬───────────────────┬─────────────────┬───────────────────┐
│   State │   Velocity [m/s] │   Temperature [K] │   Pressure [Pa] │   Density [kg/m³] │
├─────────┼──────────────────┼───────────────────┼─────────────────┼───────────────────┤
│       1 │           1000.0 │             295.0 │       1.013e+05 │             1.828 │
│       2 │            163.5 │             771.2 │       1.630e+06 │            11.18  │
│       5 │            244.4 │            1222.3 │       1.174e+07 │            49.43  │
└─────────┴──────────────────┴───────────────────┴─────────────────┴───────────────────┘
```

### 2. Experiment Planning
    
The [`FrozenShock.target_conditions`][rgfrosh.shock.FrozenShock.target_conditions] 
class method can be used to plan experiments by calculating the required 
initial conditions for target post-reflected-shock conditions:

```py hl_lines="4"
from rgfrosh import FrozenShock
import cantera as ct

shock = FrozenShock.target_conditions(ct.CarbonDioxide(), 1100, 200e5)
print(f"P1 = {shock.P1 / 133.322:.0f} torr, "
      f"u1 = {shock.u1:.1f} m/s")
```

which outputs the required fill pressure (`P1`) and incident shock velocity (`u1`):

```
P1 = 1687 torr, u1 = 918.9 m/s
```

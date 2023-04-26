# `rgfrosh` Reference

## API

The following classes are available at the package level:

- [`FrozenShock`](shock/#rgfrosh.shock.FrozenShock)
- [`IdealShock`](shock/#rgfrosh.shock.IdealShock)
- [`ThermoInterface`](thermo/#rgfrosh.thermo.ThermoInterface)
- [`ConvergenceError`](errors/#errors.ConvergenceError)

Definitions are contained within the following submodules:

- [`rgfrosh.shock`](shock/)
- [`rgfrosh.thermo`](thermo/)
- [`rgfrosh.errors`](errors/)
- [`rgfrosh.constants`](constants/)

---

## Theory

The standard shock tube region notation is observed:
    
| Region | Description                           |
| ------ | ------------------------------------- |
| 1      | Initial driven gas state              |
| 2      | Post-incident-shock driven gas state  | 
| 3      | Expanded driver gas state             | 
| 4      | Initial driver gas state              |
| 5      | Post-reflected-shock driven gas state |

### Conservation Equations

#### Incident Shock

The conservation of mass, momentum, and energy equations across a normal shock are:

$$
\rho_1 u_1 = \rho_2 u_2
$$

$$
P_1 + \rho_1 {u_1}^2 = P_2 + \rho_2 {u_2}^2
$$

$$
h_1 + \tfrac{1}{2} {u_1}^2 = h_2 + \tfrac{1}{2} {u_2}^2
$$

#### Reflected Shock

Similarly, the conservation equations for the reflected shock are:

$$
\rho_2 u_2' = \rho_5 u_5
$$

$$
P_2 + \rho_2 {u_2'}^2 = P_5 + \rho_5 {u_5}^2
$$

$$
h_2 + \tfrac{1}{2} {u_2'}^2 = h_5 + \tfrac{1}{2} {u_5}^2
$$

where $u_2'$ is the velocity of the gas in the incident shock region relative to the reflected shock

$$
u_2' = u_5 + u_1 - u_2
$$

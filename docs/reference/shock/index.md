# `rgfrosh.shock`

The core functionality of `rgfrosh` is built on the [`Shock`](#rgfrosh.shock.Shock) class, which
defines the state of the different regions of the reflected shock.
[`Shock`](#rgfrosh.shock.Shock) is a frozen data class - its attributes cannot be assigned after 
instantiation - and thus it is not very useful to instantiate it directly. Two different models
are defined to calculate shock conditions: 

- [Ideal shock](ideal/) - the analytical solution of the [conservation equations](../#conservation-equations) 
  for ideal (calorically perfect) gases.
- [Frozen shock](frozen/) - a numerical solver for a gas with an arbitrary equation of state.

---

::: rgfrosh.shock.Shock


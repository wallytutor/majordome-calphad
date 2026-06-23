# To-Do

**Short-term goals:**

- Improve robustness of the optimization routines and document its behavior.

- Create a database of substances that I work on a regular basis for ease access.

**Long-term goals:**

- See *thermodynamics/calphad-devel* in knowledge-base repository.

## 2026-05-17

- [x] Provide a simple data loading interface with possibility of substance selection in Python.

- [x] Add a substance printer method to allow inspecting properties at a given state. Goal of replacing `thermo_evaluation.rs`.

- [x] Create a wrapper API allowing for specification of system composition in terms of compount mass or moles (or respective proportions) and internally normalizing the elemental composition to one mole of matter.

- [x] Export `autodiff` to Python and add example showing that Gibbs energy derivative with temperature matches entropy at reference state. Goal of replacing `autodiff_evaluation.rs`.

- [x] Make species tabulation a default feature of substances. Goal to replace `species_tabulation.rs` and expose method in Python.

- [x] Create module `linalg` for methods not related to the core equilibrium calculations (e.g. `dense_gaussian_solver`).

- [x] Rename local equilibrium solver to `equilibrate_stoichiometric` as it is dedicated to stoichiometric systems.

- [x] Expose `equilibrate_stoichiometric` to Python.

## 2026-05-20

- [x] Simplify `SystemComposition` constructors to take a `HashMap` database.

- [x] Add `into_elemental_fractions` helper to `SystemComposition`.

- [x] Implement the `Equilibrium` class with full thermodynamic details and a gorgeous report formatter.

- [x] Make Python the source of truth: `equilibrate_stoichiometric` now accepts a `SystemComposition` object directly from Python; the Rust-only variant uses `&HashMap` as before.

## 2026-05-21

- [x] Migrate `composition_tabulation` and `hallstedt1990` into a unified `usage.rs` at the crate root; remove `sample/` binaries.

- [x] Rename `data/data.lua` → `data/sample/simple-calcination.lua`; propagate path everywhere.

- [x] Fill WIP sections in `calphad.qmd` for temperature scan and Hallstedt composition scan.

- [x] Rewrite README with authoritative status and capability table.

## Next steps — Ionic Liquid Model (Hallstedt 1990)

The TDB defines an `IONIC_LIQ` phase with two sublattices:

```
PHASE IONIC_LIQ:Y % 2 2 3 !
  CONSTITUENT IONIC_LIQ:Y : AL+3,CA+2 : O-2,VA : !
```

Implementing this requires the following stages:

- [ ] **Ion species registry** – Add `IonSpecies` to the Lua DSL with a charge
  field; extend `load_substances_from_lua` to recognise and load them without
  error.

- [ ] **Sublattice model struct** – Introduce a `SublatticePhase` type holding:
  - site multiplicities (`p`, `q`)
  - a list of cation species with site fractions `y_c`
  - a list of anion/vacancy species with site fractions `y_a`
  - end-member Gibbs parameters `G(IONIC_LIQ, cation:anion)`
  - interaction (Redlich–Kister) parameters `L(…;0)`, `L(…;1)`

- [ ] **Gibbs energy of mixing** – Implement the two-sublattice ionic liquid
  expression:
  - Ideal configurational entropy: `RT [p Σ y_c ln y_c + q Σ y_a ln y_a]`
  - Excess term: `Σ y_c y_c' y_a (L0 + L1(y_c - y_c')) `  (Redlich–Kister)

- [ ] **Internal site-fraction optimizer** – Given temperature and overall
  composition, minimize the sublattice Gibbs energy w.r.t. site fractions
  subject to charge neutrality and site-fraction sum constraints (likely a
  Newton–Raphson inner loop).

- [ ] **Outer equilibrium coupling** – Feed the liquid-phase Gibbs energy
  (now a function of internal site fractions and T) into the outer
  `equilibrate_stoichiometric` as a black-box phase Gibbs value, or extend
  the equilibrium solver to handle non-stoichiometric (solution) phases.

- [ ] **Python exposure** – Wrap `SublatticePhase` in PyO3 so that the
  liquid-phase degree of freedom is accessible from Python, similarly to
  how `Substance` is currently exposed.

- [ ] **Validation** – Reproduce Figure 1 of Hallstedt (1990): liquidus and
  eutectic temperatures in the CaO–Al₂O₃ binary at 1 bar.

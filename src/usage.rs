// | `loading_data`            | CaCO₃/CaO/CO₂ + Diaspore    | Report thermodynamic properties
// | `autodiff_verification`   | Calcite                      | Verify dG/dT = -S via autodiff
// | `species_tabulation_demo` | CaCO₃/CaO/CO₂ + Diaspore    | Property tabulation 298–1200 K
// | `single_equilibrium`      | CaCO₃ + AlOOH          | Single-point Gibbs minimization
// | `composition_tabulation`  | CaCO₃ + AlOOH                | Temperature scan 300–1200 K
// | `hallstedt1990_scan`      | CaO–Al₂O₃ (Hallstedt 1990)  | Composition scan at 1400 K
use crate::calphad::T_REF;
use crate::calphad::core::SystemComposition;
use crate::calphad::data::DatabaseLoader;
use crate::calphad::equil::equilibrate_stoichiometric;
use crate::num::autodiff::{Dual, diff};

fn main() {
    loading_data();
    autodiff_verification();
    species_tabulation_demo();
    single_equilibrium();
    composition_tabulation();
    hallstedt1990_scan();
}

fn loading_data() {
    println!("=== Thermodynamic Properties ===\n");

    let db = DatabaseLoader::new("sample/simple-calcination.lua".to_string(), None).unwrap();

    for (_name, phase) in db.data.iter() {
        println!("{}\n", phase.report(300.0));
    }
}

fn autodiff_verification() {
    println!("=== Automatic Differentiation Verification ===\n");

    let db = DatabaseLoader::new("sample/simple-calcination.lua".to_string(), None).unwrap();
    let calcite = db.data.get("Calcite").unwrap();

    let t_eval = 300.0;
    let g = |t: Dual<f64>| calcite.gibbs(t);

    let dg = diff(g, t_eval);
    let neg_s = -calcite.entropy(t_eval);

    println!("Calcite at {:.2} K:", t_eval);
    println!("  dG/dT (autodiff) = {:.6} J/(mol.K)", dg);
    println!("  -S(T) (computed) = {:.6} J/(mol.K)", neg_s);
    println!("  Difference       = {:.6e} J/(mol.K)", (dg - neg_s).abs());
    println!();
}

fn species_tabulation_demo() {
    println!("=== Species Thermodynamic Tabulation ===\n");

    let db = DatabaseLoader::new("sample/simple-calcination.lua".to_string(), None).unwrap();
    let names = ["Calcite", "Lime", "CO2", "Diaspore", "H2O", "Al2O3"];

    for name in &names {
        if let Some(substance) = db.data.get(*name) {
            println!("{}\n", substance.tabulate(None, None, None));
        }
    }
}

fn single_equilibrium() {
    let db = DatabaseLoader::new("sample/simple-calcination.lua".to_string(), None).unwrap();

    let names = ["Calcite", "Lime", "CO2", "Diaspore", "H2O", "Al2O3"];
    let species: Vec<_> = names.iter().map(|n| db.data.get(*n).unwrap()).collect();

    let t = 1173.15_f64;
    let p = 101325.0_f64;

    let mut proportions = std::collections::HashMap::new();
    proportions.insert("Calcite".to_string(), 1.0);
    proportions.insert("Diaspore".to_string(), 1.0);

    let comp = SystemComposition::from_compound_moles(db.data.clone(), proportions).unwrap();
    let b = comp.into_elemental_fractions();

    println!("\n=== Generic CALPHAD Local Equilibrium ===");
    println!("T = {} K, P = {} Pa", t, p);

    let eq = equilibrate_stoichiometric(&species, &b, t, p);
    println!("{}", eq.report());
}

/// Temperature scan for the simple calcination system (CaCO3 + Diaspore).
///
/// Reports phase amounts and mass-specific enthalpy change from 300 K to 1200 K,
/// using 298.15 K as the reference state.
fn composition_tabulation() {
    let db = DatabaseLoader::new("sample/simple-calcination.lua".to_string(), None).unwrap();

    let names = ["Calcite", "Lime", "CO2", "Diaspore", "H2O", "Al2O3"];
    let species: Vec<_> = names.iter().map(|n| db.data.get(*n).unwrap()).collect();

    let mut proportions = std::collections::HashMap::new();
    proportions.insert("Calcite".to_string(), 1.0);
    proportions.insert("Diaspore".to_string(), 1.0);

    let comp = SystemComposition::from_compound_moles(db.data.clone(), proportions).unwrap();
    let b = comp.into_elemental_fractions();

    // Reference state at T_REF for enthalpy normalization
    let eq_ref = equilibrate_stoichiometric(&species, &b, T_REF, 101325.0);
    let mut h_ref = 0.0;
    let mut m_sys = 0.0;

    for s in &species {
        let n = eq_ref.amounts.get(&s.name).copied().unwrap_or(0.0);
        h_ref += n * s.enthalpy(T_REF);
        m_sys += n * s.molar_mass;
    }

    println!("\n=== Composition Tabulation (300 K – 1200 K) ===");
    println!(
        "{:<10} | {:<12} | {:<10} | {:<10} | {:<14} | {:<10} | {:<10} | {:<14}",
        "T (K)",
        "CaCO3 (mol)",
        "CaO (mol)",
        "CO2 (mol)",
        "Diaspore (mol)",
        "H2O (mol)",
        "Al2O3 (mol)",
        "dH (J/g)"
    );
    println!(
        "{:-<10}-+-{:-<12}-+-{:-<10}-+-{:-<10}-+-{:-<14}-+-{:-<10}-+-{:-<10}-+-{:-<14}-",
        "", "", "", "", "", "", "", ""
    );

    let t_min = 300.0;
    let t_max = 1200.0;
    let t_inc = 100.0;
    let p = 101325.0;

    let mut t = t_min;

    while t <= t_max {
        let eq = equilibrate_stoichiometric(&species, &b, t, p);

        let mut h_sys = 0.0;

        for s in &species {
            let n = eq.amounts.get(&s.name).copied().unwrap_or(0.0);
            h_sys += n * s.enthalpy(t);
        }

        let dh = if m_sys > 0.0 {
            (h_sys - h_ref) / m_sys
        } else {
            0.0
        };

        println!(
            "{:<10.2} | {:<12.6} | {:<10.6} | {:<10.6} | {:<14.6} | {:<10.6} | {:<10.6} | {:<14.2}",
            t,
            eq.amounts.get("Calcite").copied().unwrap_or(0.0),
            eq.amounts.get("Lime").copied().unwrap_or(0.0),
            eq.amounts.get("CO2").copied().unwrap_or(0.0),
            eq.amounts.get("Diaspore").copied().unwrap_or(0.0),
            eq.amounts.get("H2O").copied().unwrap_or(0.0),
            eq.amounts.get("Al2O3").copied().unwrap_or(0.0),
            dh,
        );

        t += t_inc;
    }

    println!();
}

/// Composition scan for the CaO–Al2O3 system from Hallstedt (1990).
///
/// Varies x(Al2O3) from 0 to 1 at 1400 K and reports the normalized
/// mole fraction of each stoichiometric phase at equilibrium.
fn hallstedt1990_scan() {
    let db = DatabaseLoader::new("hallstedt1990.lua".to_string(), None)
        .expect("Failed to load hallstedt1990.lua");

    let names = ["HALITE", "CORUNDUM", "C3A1", "C1A1", "C1A2", "C1A6"];
    let species: Vec<_> = names
        .iter()
        .map(|n| db.data.get(*n).expect(&format!("{} not found", n)))
        .collect();

    let t_eval = 1400.0;
    let p_eval = 101325.0;

    println!(
        "\n=== Hallstedt (1990) CaO–Al2O3 Phase Stability at {:.0} K ===",
        t_eval
    );
    println!(
        "{:<10} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10}",
        "x(Al2O3)", "HALITE", "C3A1", "C1A1", "C1A2", "C1A6", "CORUNDUM"
    );
    println!(
        "{:-<10}-+-{:-<10}-+-{:-<10}-+-{:-<10}-+-{:-<10}-+-{:-<10}-+-{:-<10}-",
        "", "", "", "", "", "", ""
    );

    let mut x = 0.0_f64;

    while x <= 1.0001 {
        let mut proportions = std::collections::HashMap::new();
        proportions.insert("HALITE".to_string(), 1.0 - x);
        proportions.insert("CORUNDUM".to_string(), x);

        let comp = SystemComposition::from_compound_moles(db.data.clone(), proportions).unwrap();
        let b = comp.into_elemental_fractions();

        let eq = equilibrate_stoichiometric(&species, &b, t_eval, p_eval);
        let total: f64 = names
            .iter()
            .map(|n| eq.amounts.get(*n).copied().unwrap_or(0.0))
            .sum();

        let norm = |name: &str| {
            if total > 0.0 {
                eq.amounts.get(name).copied().unwrap_or(0.0) / total
            } else {
                0.0
            }
        };

        println!(
            "{:<10.2} | {:<10.4} | {:<10.4} | {:<10.4} | {:<10.4} | {:<10.4} | {:<10.4}",
            x,
            norm("HALITE"),
            norm("C3A1"),
            norm("C1A1"),
            norm("C1A2"),
            norm("C1A6"),
            norm("CORUNDUM"),
        );

        x += 0.05;
    }

    println!();
}

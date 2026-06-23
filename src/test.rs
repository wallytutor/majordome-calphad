use super::core::AggregationType;
use super::core::SystemComposition;
use super::core::extract_elements;
use super::core::get_atomic_weight;
use super::data::DatabaseLoader;
use super::data::load_substances_from_lua;
use super::equilibrium::equilibrate_stoichiometric;
use super::mixture::FractionConverter;
use super::mixture::InterstitialConverter;

use majordome_numerical::prelude::{Dual, diff};

const SIMPLE_CALCINATION: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../majordome/majordome/data/calphad/sample/simple-calcination.lua"
);

// --------------------------------------------------------------------------

#[test]
fn test_load_from_lua() {
    let result = load_substances_from_lua(SIMPLE_CALCINATION);
    assert!(
        result.is_ok(),
        "Failed to load from lua: {:?}",
        result.err()
    );
    let map = result.unwrap();
    assert!(map.len() >= 6);
    assert!(map.contains_key("Calcite"));
    assert!(map.contains_key("CO2"));
    assert_eq!(
        map.get("Calcite").unwrap().aggregation_type,
        AggregationType::Solid
    );
    assert_eq!(
        map.get("CO2").unwrap().aggregation_type,
        AggregationType::Gas
    );
}

#[test]
fn test_database_loader() {
    let loader = DatabaseLoader::new(SIMPLE_CALCINATION.to_string(), None).unwrap();
    assert_eq!(loader.path, SIMPLE_CALCINATION);
    assert!(loader.phases.len() >= 6);
    assert!(loader.phases.contains(&"Calcite".to_string()));

    let data = &loader.data;
    assert!(data.len() >= 6);
    assert!(data.contains_key("Calcite"));

    // Test loading phase list filtering
    let loader_filtered = DatabaseLoader::new(
        SIMPLE_CALCINATION.to_string(),
        Some(vec!["Calcite".to_string(), "CO2".to_string()]),
    )
    .unwrap();
    assert_eq!(loader_filtered.data.len(), 2);
    assert!(loader_filtered.data.contains_key("Calcite"));
    assert!(loader_filtered.data.contains_key("CO2"));
    assert!(!loader_filtered.data.contains_key("Lime"));

    // Test load_compound retrieving data from cache
    let compound = loader_filtered
        .load_compound("Calcite".to_string())
        .unwrap();
    assert_eq!(compound.name, "Calcite");

    let compound2 = loader.load_compound("Lime".to_string()).unwrap();
    assert_eq!(compound2.name, "Lime");
}

// --------------------------------------------------------------------------

#[test]
fn test_thermo_derivatives() {
    let db = load_substances_from_lua(SIMPLE_CALCINATION).unwrap();
    let calcite = db.get("Calcite").unwrap();
    let g = |t: Dual<f64>| calcite.gibbs(t);

    let expected = -calcite.entropy(300.0);
    let actual = diff(g, 300.0);
    assert!(
        (expected - actual).abs() < 1e-9,
        "expected={}, actual={}",
        expected,
        actual
    );
}

#[test]
fn test_h2o_nasa7() {
    let db = load_substances_from_lua(SIMPLE_CALCINATION).unwrap();
    let h2o = db.get("H2O").unwrap();
    let val = h2o.cp(300.0);
    assert!(val > 0.0);
}

#[test]
fn test_system_composition() {
    let db = load_substances_from_lua(SIMPLE_CALCINATION).unwrap();

    let mut proportions = std::collections::HashMap::new();
    proportions.insert("Calcite".to_string(), 1.0);
    proportions.insert("Diaspore".to_string(), 1.0);

    let comp = SystemComposition::from_compound_moles(db, proportions).unwrap();
    let elements = comp.elements();
    assert_eq!(elements, vec!["Al", "C", "Ca", "H", "O"]);

    let b = comp.fractions();
    assert!((b[0] - 1.0 / 9.0).abs() < 1e-9);
    assert!((b[1] - 1.0 / 9.0).abs() < 1e-9);
    assert!((b[2] - 1.0 / 9.0).abs() < 1e-9);
    assert!((b[3] - 1.0 / 9.0).abs() < 1e-9);
    assert!((b[4] - 5.0 / 9.0).abs() < 1e-9);

    let rep = comp.report();
    assert!(rep.contains("Input Method: Compound Moles"));
}

// --------------------------------------------------------------------------

#[test]
fn test_extract_elements() {
    let db = load_substances_from_lua(SIMPLE_CALCINATION).unwrap();
    let calcite = db.get("Calcite").unwrap();
    let diaspore = db.get("Diaspore").unwrap();
    let elements = extract_elements(&[calcite, diaspore]);
    assert_eq!(elements, vec!["Al", "C", "Ca", "H", "O"]);
}

#[test]
fn test_equilibrate_stoichiometric() {
    let db = load_substances_from_lua(SIMPLE_CALCINATION).unwrap();
    let calcite = db.get("Calcite").unwrap();
    let lime = db.get("Lime").unwrap();
    let co2 = db.get("CO2").unwrap();
    let diaspore = db.get("Diaspore").unwrap();
    let h2o = db.get("H2O").unwrap();
    let al2o3 = db.get("Al2O3").unwrap();

    let species = vec![calcite, lime, co2, diaspore, h2o, al2o3];

    let mut b = std::collections::HashMap::new();
    b.insert("Ca".to_string(), 1.0);
    b.insert("C".to_string(), 1.0);
    b.insert("Al".to_string(), 1.0);
    b.insert("H".to_string(), 1.0);
    b.insert("O".to_string(), 5.0);

    let eq = equilibrate_stoichiometric(&species, &b, 1173.15, 1.0);
    let phi = &eq.amounts;
    assert_eq!(phi.len(), 6);
    assert!((phi.get("Calcite").copied().unwrap_or(0.0)).abs() < 1e-4);
    assert!((phi.get("Lime").copied().unwrap_or(0.0) - 1.0).abs() < 1e-4);
    assert!((phi.get("CO2").copied().unwrap_or(0.0) - 1.0).abs() < 1e-4);
    assert!((phi.get("Diaspore").copied().unwrap_or(0.0)).abs() < 1e-4);
    assert!((phi.get("H2O").copied().unwrap_or(0.0) - 0.5).abs() < 1e-4);
    assert!((phi.get("Al2O3").copied().unwrap_or(0.0) - 0.5).abs() < 1e-4);
}

// --------------------------------------------------------------------------

#[test]
fn test_fraction_converter() {
    let elements = ["Fe".to_string(), "C".to_string()];
    let converter = FractionConverter::new(&elements).unwrap();
    let mass_frac0 = vec![0.99, 0.01];
    let mole_frac1 = converter.mass_to_mole_fraction(&mass_frac0).unwrap();
    let mass_frac2 = converter.mole_to_mass_fraction(&mole_frac1).unwrap();

    assert!((mass_frac0[0] - mass_frac2[0]).abs() < 1.0e-10);
    assert!((mass_frac0[1] - mass_frac2[1]).abs() < 1.0e-10);
}

#[test]
fn test_interstitial_converter() {
    let rho_fe = 7.874;
    let m_fe = get_atomic_weight("Fe").unwrap();
    let converter = InterstitialConverter::new(m_fe, rho_fe);

    let mole_frac0 = vec![0.01, 0.02];
    let concentration1 = converter.mole_fraction_to_concentration(&mole_frac0);
    let mole_frac1 = converter.concentration_to_mole_fraction(&concentration1);

    assert!((mole_frac0[0] - mole_frac1[0]).abs() < 1.0e-10);
    assert!((mole_frac0[1] - mole_frac1[1]).abs() < 1.0e-10);
}

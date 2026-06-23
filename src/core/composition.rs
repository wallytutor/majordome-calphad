use super::get_atomic_weight;
use super::substance::Substance;
use pyo3::prelude::*;
use std::collections::HashMap;

/// Extracts and sorts unique element symbols from a list of compound substances.
///
/// # Arguments
/// * `compounds` - A slice of references to `Substance` structures.
///
/// # Returns
/// A sorted vector of unique element symbols (e.g. `["C", "Fe", "N"]`).
pub fn extract_elements(compounds: &[&Substance]) -> Vec<String> {
    let mut uniq = std::collections::HashSet::new();

    for s in compounds {
        for el in s.elements.keys() {
            uniq.insert(el.clone());
        }
    }

    let mut res: Vec<String> = uniq.into_iter().collect();
    res.sort();
    res
}

/// Represents the overall composition of a thermodynamic system.
///
/// It converts initial input proportions (moles or masses of compounds/elements) into
/// elemental mole fractions normalized to 1 mole of atoms.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct SystemComposition {
    pub(crate) elements: Vec<String>,
    pub(crate) fractions: Vec<f64>,
    pub(crate) input_method: String,
    pub(crate) input_proportions: HashMap<String, f64>,
}

#[pymethods]
impl SystemComposition {
    /// Builds a system composition from compound molar amounts.
    ///
    /// # Arguments
    /// * `species` - A dictionary map of substance name to `Substance` records.
    /// * `proportions` - A map of substance name to its quantity in moles (`mol`).
    ///
    /// # Errors
    /// Returns `PyValueError` if any substance specified in `proportions` is not found in `species`.
    #[staticmethod]
    pub fn from_compound_moles(
        species: HashMap<String, Substance>,
        proportions: HashMap<String, f64>,
    ) -> PyResult<Self> {
        let mut element_moles = HashMap::new();
        let mut input_proportions = HashMap::new();

        for (name, moles) in &proportions {
            let substance = species.get(name).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "Substance '{}' specified in composition not found in species database",
                    name
                ))
            })?;

            input_proportions.insert(substance.name.clone(), *moles);

            for (el, coeff) in &substance.elements {
                *element_moles.entry(el.clone()).or_insert(0.0) += moles * coeff;
            }
        }

        let subs: Vec<&Substance> = species.values().collect();
        let elements = extract_elements(&subs);

        let total_moles: f64 = element_moles.values().sum();
        let fractions = if total_moles > 0.0 {
            elements
                .iter()
                .map(|el| element_moles.get(el).copied().unwrap_or(0.0) / total_moles)
                .collect()
        } else {
            vec![0.0; elements.len()]
        };

        Ok(Self {
            elements,
            fractions,
            input_method: "Compound Moles".to_string(),
            input_proportions,
        })
    }

    /// Builds a system composition from compound mass amounts.
    ///
    /// # Arguments
    /// * `species` - A dictionary map of substance name to `Substance` records.
    /// * `proportions` - A map of substance name to its mass in grams (`g`).
    ///
    /// # Errors
    /// Returns `PyValueError` if any substance specified in `proportions` is not found in `species`.
    #[staticmethod]
    pub fn from_compound_masses(
        species: HashMap<String, Substance>,
        proportions: HashMap<String, f64>,
    ) -> PyResult<Self> {
        let mut element_moles = HashMap::new();
        let mut input_proportions = HashMap::new();

        for (name, mass) in &proportions {
            let substance = species.get(name).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "Substance '{}' specified in composition not found in species database",
                    name
                ))
            })?;

            input_proportions.insert(substance.name.clone(), *mass);
            let moles = mass / substance.molar_mass;

            for (el, coeff) in &substance.elements {
                *element_moles.entry(el.clone()).or_insert(0.0) += moles * coeff;
            }
        }

        let subs: Vec<&Substance> = species.values().collect();
        let elements = extract_elements(&subs);

        let total_moles: f64 = element_moles.values().sum();
        let fractions = if total_moles > 0.0 {
            elements
                .iter()
                .map(|el| element_moles.get(el).copied().unwrap_or(0.0) / total_moles)
                .collect()
        } else {
            vec![0.0; elements.len()]
        };

        Ok(Self {
            elements,
            fractions,
            input_method: "Compound Masses".to_string(),
            input_proportions,
        })
    }

    /// Converts computed elemental fractions into a dictionary map.
    ///
    /// # Returns
    /// A map of element symbols to their computed mole fractions (dimensionless,
    /// normalized to 1.0).
    #[pyo3(name = "into_elemental_fractions")]
    pub fn into_elemental_fractions(&self) -> HashMap<String, f64> {
        let mut map = HashMap::new();

        for (el, frac) in self.elements.iter().zip(self.fractions.iter()) {
            map.insert(el.clone(), *frac);
        }

        map
    }

    /// Builds a system composition directly from elemental moles.
    ///
    /// # Arguments
    /// * `elements_in` - A map of element symbol strings to their quantities in moles (`mol`).
    #[staticmethod]
    pub fn from_elemental_moles(elements_in: HashMap<String, f64>) -> PyResult<Self> {
        let mut element_moles = HashMap::new();
        let mut input_proportions = HashMap::new();

        for (el, moles) in elements_in {
            input_proportions.insert(el.clone(), moles);
            *element_moles.entry(el).or_insert(0.0) += moles;
        }

        let mut elements: Vec<String> = element_moles.keys().cloned().collect();
        elements.sort();

        let total_moles: f64 = element_moles.values().sum();
        let fractions = if total_moles > 0.0 {
            elements
                .iter()
                .map(|el| element_moles.get(el).copied().unwrap_or(0.0) / total_moles)
                .collect()
        } else {
            vec![0.0; elements.len()]
        };

        Ok(Self {
            elements,
            fractions,
            input_method: "Elemental Moles".to_string(),
            input_proportions,
        })
    }

    /// Builds a system composition directly from elemental masses.
    ///
    /// # Arguments
    /// * `elements_in` - A map of element symbol strings to their masses in grams (`g`).
    ///
    /// # Errors
    /// Returns `PyValueError` if any element symbol is unrecognized in the atomic weights database.
    #[staticmethod]
    pub fn from_elemental_masses(elements_in: HashMap<String, f64>) -> PyResult<Self> {
        let mut element_moles = HashMap::new();
        let mut input_proportions = HashMap::new();

        for (el, mass) in elements_in {
            input_proportions.insert(el.clone(), mass);
            let atomic_w = get_atomic_weight(&el).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!("Unknown element symbol: '{}'", el))
            })?;

            let moles = mass / atomic_w;
            *element_moles.entry(el).or_insert(0.0) += moles;
        }

        let mut elements: Vec<String> = element_moles.keys().cloned().collect();
        elements.sort();

        let total_moles: f64 = element_moles.values().sum();
        let fractions = if total_moles > 0.0 {
            elements
                .iter()
                .map(|el| element_moles.get(el).copied().unwrap_or(0.0) / total_moles)
                .collect()
        } else {
            vec![0.0; elements.len()]
        };

        Ok(Self {
            elements,
            fractions,
            input_method: "Elemental Masses".to_string(),
            input_proportions,
        })
    }

    /// Gets the list of elements sorted alphabetically.
    #[getter]
    pub fn elements(&self) -> Vec<String> {
        self.elements.clone()
    }

    /// Gets the list of computed dimensionless mole fractions corresponding to the elements list.
    #[getter]
    pub fn fractions(&self) -> Vec<f64> {
        self.fractions.clone()
    }

    /// Gets the input composition parsing method.
    #[getter]
    pub fn input_method(&self) -> String {
        self.input_method.clone()
    }

    /// Gets the raw input proportions map.
    #[getter]
    pub fn input_proportions(&self) -> HashMap<String, f64> {
        self.input_proportions.clone()
    }

    /// Generates a readable string report summarizing system composition.
    pub fn report(&self) -> String {
        let mut lines = Vec::new();
        lines.push("=== System Composition Report ===".to_string());
        lines.push(format!("Input Method: {}", self.input_method));
        lines.push("\nInput Proportions:".to_string());

        let mut sorted_inputs: Vec<(&String, &f64)> = self.input_proportions.iter().collect();
        sorted_inputs.sort_by(|a, b| a.0.cmp(b.0));

        for (name, prop) in sorted_inputs {
            lines.push(format!("  * {}: {:.6}", name, prop));
        }

        lines.push(
            "\nComputed Elemental Mole Fractions (normalized to 1 mole of atoms):".to_string(),
        );

        for (i, el) in self.elements.iter().enumerate() {
            lines.push(format!("  * {:>2}: {:.6}", el, self.fractions[i]));
        }

        lines.join("\n")
    }
}

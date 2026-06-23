use super::constant::T_REF;
use super::parameterization::{AggregationType, TemperatureRange};
use majordome_numerical::prelude::*;
use majordome_utilities::prelude::*;
use pyo3::prelude::*;
use std::collections::HashMap;

/// Represents a pure chemical substance or compound with standardized thermodynamic parameterizations.
///
/// Provides functions to compute temperature-dependent thermodynamic properties:
/// specific heat ($C_p$), enthalpy ($H$), entropy ($S$), and Gibbs free energy ($G$).
///
/// # Units of Properties
/// - **Temperature ($T$)**: Kelvin ($\text{K}$)
/// - **Molar Mass**: grams per mole ($\text{g/mol}$)
/// - **Molar Volume**: cubic centimeters per mole ($\text{cm}^3\text{/mol}$)
/// - **Specific Heat ($C_p$)**: Joules per mole-Kelvin ($\text{J/(mol·K)}$)
/// - **Entropy ($S$, $S_0$)**: Joules per mole-Kelvin ($\text{J/(mol·K)}$)
/// - **Enthalpy ($H$)**: Joules per mole ($\text{J/mol}$)
/// - **Gibbs Free Energy ($G$)**: Joules per mole ($\text{J/mol}$)
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct Substance {
    pub name: String,
    pub molar_mass: f64,
    pub molar_volume: f64,
    pub s0: f64,
    pub ranges: Vec<TemperatureRange>,
    pub elements: HashMap<String, f64>,
    pub reference: String,
    pub aggregation_type: AggregationType,
}

impl Substance {
    /// Selects the appropriate `TemperatureRange` model for the given temperature.
    ///
    /// # Arguments
    /// * `t` - Temperature in Kelvin ($\text{K}$).
    pub fn get_range(&self, t: f64) -> &TemperatureRange {
        for r in &self.ranges {
            if t >= r.t_min && t <= r.t_max {
                return r;
            }
        }

        if t < self.ranges[0].t_min {
            return &self.ranges[0];
        }

        match self.ranges.last() {
            Some(r) => r,
            None => &self.ranges[0],
        }
    }

    /// Computes the molar isobaric specific heat ($C_p$) at the given temperature.
    ///
    /// # Arguments
    /// * `t` - Temperature in Kelvin ($\text{K}$).
    ///
    /// # Returns
    /// Specific heat in $\text{J/(mol·K)}$.
    pub fn cp<T: Numeric>(&self, t: T) -> T {
        let t_val = t.to_f64();
        let range = self.get_range(t_val);
        range.model.cp(t)
    }

    /// Computes the molar enthalpy ($H$) at the given temperature.
    ///
    /// # Arguments
    /// * `t` - Temperature in Kelvin ($\text{K}$).
    ///
    /// # Returns
    /// Molar enthalpy in $\text{J/mol}$.
    pub fn enthalpy<T: Numeric>(&self, t: T) -> T {
        let t_val = t.to_f64();
        let range = self.get_range(t_val);
        range.model.enthalpy(t)
    }

    /// Computes the absolute molar entropy ($S$) at the given temperature.
    ///
    /// # Arguments
    /// * `t` - Temperature in Kelvin ($\text{K}$).
    ///
    /// # Returns
    /// Absolute entropy in $\text{J/(mol·K)}$.
    pub fn entropy<T: Numeric>(&self, t: T) -> T {
        let t_val = t.to_f64();
        let range = self.get_range(t_val);
        range.model.entropy(t)
    }

    /// Computes the molar Gibbs free energy ($G = H - T \cdot S$) at the given temperature.
    ///
    /// # Arguments
    /// * `t` - Temperature in Kelvin ($\text{K}$).
    ///
    /// # Returns
    /// Molar Gibbs free energy in $\text{J/mol}$.
    pub fn gibbs<T: Numeric>(&self, t: T) -> T {
        self.enthalpy(t) - t * self.entropy(t)
    }

    /// Tabulates thermodynamic properties ($C_p$, $S$, enthalpy relative to reference temperature 298.15 K, etc.)
    /// over a specified temperature range.
    ///
    /// # Arguments
    /// * `t_min` - Minimum temperature in Kelvin ($\text{K}$). Defaults to the substance's minimum range boundary.
    /// * `t_max` - Maximum temperature in Kelvin ($\text{K}$). Defaults to the substance's maximum range boundary.
    /// * `step` - Temperature step increment in Kelvin ($\text{K}$). Defaults to $100\text{ K}$.
    pub fn tabulate(&self, t_min: Option<f64>, t_max: Option<f64>, step: Option<f64>) -> String {
        let t_start =
            t_min.unwrap_or_else(|| self.ranges.first().map(|r| r.t_min).unwrap_or(200.0));
        let t_end = t_max.unwrap_or_else(|| self.ranges.last().map(|r| r.t_max).unwrap_or(3000.0));
        let step_val = step.unwrap_or(100.0);

        let mut temps = Vec::new();
        temps.push(t_start);
        temps.push(t_end);

        if T_REF >= t_start && T_REF <= t_end {
            temps.push(T_REF);
        }

        for r in &self.ranges {
            if r.t_min >= t_start && r.t_min <= t_end {
                temps.push(r.t_min);
            }
            if r.t_max >= t_start && r.t_max <= t_end {
                temps.push(r.t_max);
            }
        }

        let first_multiple = (t_start / step_val).ceil() * step_val;
        let mut current = first_multiple;

        while current <= t_end {
            temps.push(current);
            current += step_val;
        }

        temps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut unique_temps = Vec::new();

        for &t in &temps {
            if unique_temps.is_empty()
                || (t - unique_temps.last().copied().unwrap_or(0.0)).abs() > 1e-6
            {
                unique_temps.push(t);
            }
        }

        let mut lines = Vec::new();
        lines.push(format!(
            "--- {} ({:?}) ---\n",
            self.name, self.aggregation_type
        ));
        lines.push(format!(
            "{:<8} | {:<12} | {:<12} | {:<14} | {:<14}",
            "T (K)", "Cp", "S", "-(G-H298)/T", "H-H298"
        ));
        lines.push(format!(
            "{:-<8}-+-{:-<12}-+-{:-<12}-+-{:-<14}-+-{:-<14}-",
            "", "", "", "", ""
        ));

        let h_ref = self.enthalpy(T_REF);

        for t in unique_temps {
            let cp_val = self.cp(t);
            let h_val = self.enthalpy(t);
            let s_val = self.entropy(t);
            let g_val = self.gibbs(t);

            let free_energy_func = if t > 1e-6 { -(g_val - h_ref) / t } else { 0.0 };
            let h_diff = h_val - h_ref;

            lines.push(format!(
                "{:<8.2} | {:<12.4} | {:<12.4} | {:<14.4} | {:<14.2}",
                t, cp_val, s_val, free_energy_func, h_diff
            ));
        }

        lines.join("\n")
    }
}

#[pymethods]
impl Substance {
    /// Name of the substance.
    #[getter]
    pub fn name(&self) -> String {
        self.name.clone()
    }

    /// Molar mass in grams per mole ($\text{g/mol}$).
    #[getter]
    pub fn molar_mass(&self) -> f64 {
        self.molar_mass
    }

    /// Molar volume in cubic centimeters per mole ($\text{cm}^3\text{/mol}$).
    #[getter]
    pub fn molar_volume(&self) -> f64 {
        self.molar_volume
    }

    /// Reference absolute entropy $S_0$ in $\text{J/(mol·K)}$.
    #[getter]
    pub fn s0(&self) -> f64 {
        self.s0
    }

    /// Map of element symbols to their stoichiometric coefficients (dimensionless ratios).
    #[getter]
    pub fn elements(&self) -> HashMap<String, f64> {
        self.elements.clone()
    }

    /// Bibliographic or database reference string.
    #[getter]
    pub fn reference(&self) -> String {
        self.reference.clone()
    }

    /// State aggregation type description (e.g. "Gas", "Liquid", "Solid").
    #[getter]
    pub fn aggregation_type(&self) -> String {
        format!("{:?}", self.aggregation_type)
    }

    /// Standard molar enthalpy of formation at $298.15\text{ K}$ ($\text{J/mol}$).
    #[getter]
    pub fn formation_enthalpy(&self) -> f64 {
        self.enthalpy(298.15)
    }

    /// Computes specific heat $C_p$ at temperature `t` ($\text{K}$).
    /// Accepts a float or a dual number for automatic differentiation.
    #[pyo3(name = "cp")]
    pub fn cp_py<'py>(
        &self,
        py: Python<'py>,
        t: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(dual_t) = t.extract::<PyDual>() {
            let res = self.cp(dual_t.inner);
            let bound = PyDual { inner: res }.into_pyobject(py)?;
            Ok(bound.into_any())
        } else if let Ok(f64_t) = t.extract::<f64>() {
            let res = self.cp(f64_t);
            let bound = res.into_pyobject(py)?;
            Ok(bound.into_any())
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "Expected f64 or Dual",
            ))
        }
    }

    /// Computes enthalpy $H$ at temperature `t` ($\text{K}$).
    /// Accepts a float or a dual number for automatic differentiation.
    #[pyo3(name = "enthalpy")]
    pub fn enthalpy_py<'py>(
        &self,
        py: Python<'py>,
        t: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(dual_t) = t.extract::<PyDual>() {
            let res = self.enthalpy(dual_t.inner);
            let bound = PyDual { inner: res }.into_pyobject(py)?;
            Ok(bound.into_any())
        } else if let Ok(f64_t) = t.extract::<f64>() {
            let res = self.enthalpy(f64_t);
            let bound = res.into_pyobject(py)?;
            Ok(bound.into_any())
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "Expected f64 or Dual",
            ))
        }
    }

    /// Computes absolute entropy $S$ at temperature `t` ($\text{K}$).
    /// Accepts a float or a dual number for automatic differentiation.
    #[pyo3(name = "entropy")]
    pub fn entropy_py<'py>(
        &self,
        py: Python<'py>,
        t: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(dual_t) = t.extract::<PyDual>() {
            let res = self.entropy(dual_t.inner);
            let bound = PyDual { inner: res }.into_pyobject(py)?;
            Ok(bound.into_any())
        } else if let Ok(f64_t) = t.extract::<f64>() {
            let res = self.entropy(f64_t);
            let bound = res.into_pyobject(py)?;
            Ok(bound.into_any())
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "Expected f64 or Dual",
            ))
        }
    }

    /// Computes Gibbs free energy $G$ at temperature `t` ($\text{K}$).
    /// Accepts a float or a dual number for automatic differentiation.
    #[pyo3(name = "gibbs")]
    pub fn gibbs_py<'py>(
        &self,
        py: Python<'py>,
        t: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(dual_t) = t.extract::<PyDual>() {
            let res = self.gibbs(dual_t.inner);
            let bound = PyDual { inner: res }.into_pyobject(py)?;
            Ok(bound.into_any())
        } else if let Ok(f64_t) = t.extract::<f64>() {
            let res = self.gibbs(f64_t);
            let bound = res.into_pyobject(py)?;
            Ok(bound.into_any())
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "Expected f64 or Dual",
            ))
        }
    }

    /// Generates a complete thermodynamic properties report at temperature `t` ($\text{K}$).
    pub fn report(&self, t: f64) -> String {
        let tabulate = |label: &str, units: &str, value: &str| {
            format!("{:15} | {:10} | {}", label, units, value)
        };

        let mut lines = Vec::new();
        lines.push(format!(
            "--- {} ({:?}) ---\n",
            self.name, self.aggregation_type
        ));

        lines.push(tabulate(
            "Molar mass",
            "g/mol",
            &exponential_fmt(self.molar_mass),
        ));

        lines.push(tabulate(
            "Molar volume",
            "cm³/mol",
            &exponential_fmt(self.molar_volume),
        ));

        lines.push(tabulate(
            "Entropy S0",
            "J/(mol.K)",
            &exponential_fmt(self.s0),
        ));

        lines.push(tabulate(
            "Delta-Hf298",
            "kJ/mol",
            exponential_fmt(self.formation_enthalpy() / 1000.0).as_str(),
        ));

        let mut sorted_elements: Vec<(&String, &f64)> = self.elements.iter().collect();
        sorted_elements.sort_by(|a, b| a.0.cmp(b.0));

        lines.push("\nElements:".to_string());

        for (el, coeff) in sorted_elements {
            lines.push(format!("  * {:>2}: {:.6}", el, coeff));
        }

        lines.push(format!("\nProperties at {:.2} K:", t));
        lines.push(tabulate(
            " Specific heat",
            "J/(mol.K)",
            exponential_fmt(self.cp(t)).as_str(),
        ));
        lines.push(tabulate(
            " Entropy",
            "J/(mol.K)",
            exponential_fmt(self.entropy(t)).as_str(),
        ));
        lines.push(tabulate(
            " Gibbs energy",
            "kJ/mol",
            exponential_fmt(self.gibbs(t) / 1000.0).as_str(),
        ));
        lines.push(tabulate(
            " Enthalpy",
            "kJ/mol",
            exponential_fmt(self.enthalpy(t) / 1000.0).as_str(),
        ));

        lines.push(format!("\nReference: {}", self.reference));
        lines.join("\n")
    }

    /// Tabulates thermodynamic properties for Python interface.
    #[pyo3(name = "tabulate", signature = (t_min = None, t_max = None, step = None))]
    pub fn tabulate_py(&self, t_min: Option<f64>, t_max: Option<f64>, step: Option<f64>) -> String {
        self.tabulate(t_min, t_max, step)
    }
}

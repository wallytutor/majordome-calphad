use super::core::AggregationType;
use super::core::R_GAS;
use super::core::Substance;
use super::core::SystemComposition;
use super::core::extract_elements;
use majordome_numerical::prelude::*;
use pyo3::prelude::*;
use std::collections::HashMap;

/// Represents the resolved thermodynamic state and chemical equilibrium of a system.
///
/// Contains computed molar amounts, individual phase Gibbs energies, and global properties
/// at a specified temperature and pressure.
///
/// # Units of Properties
/// - **amounts**: Moles of each substance present at equilibrium ($\text{mol}$).
/// - **gibbs_energies**: Gibbs free energies of the pure phases at conditions ($\text{J/mol}$).
/// - **total_gibbs**: Global minimized Gibbs free energy of the entire system ($\text{J}$).
/// - **temperature**: System temperature in Kelvin ($\text{K}$).
/// - **pressure**: System pressure in Pascals ($\text{Pa}$).
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct Equilibrium {
    #[pyo3(get)]
    pub amounts: HashMap<String, f64>,
    #[pyo3(get)]
    pub gibbs_energies: HashMap<String, f64>,
    #[pyo3(get)]
    pub total_gibbs: f64,
    #[pyo3(get)]
    pub temperature: f64,
    #[pyo3(get)]
    pub pressure: f64,
    pub substances: Vec<Substance>,
}

#[pymethods]
impl Equilibrium {
    /// Creates a new `Equilibrium` state.
    #[new]
    pub fn new(
        amounts: HashMap<String, f64>,
        gibbs_energies: HashMap<String, f64>,
        total_gibbs: f64,
        temperature: f64,
        pressure: f64,
        substances: Vec<Substance>,
    ) -> Self {
        Self {
            amounts,
            gibbs_energies,
            total_gibbs,
            temperature,
            pressure,
            substances,
        }
    }

    /// Generates a highly detailed text report showing equilibrium conditions,
    /// total thermodynamic Gibbs energy, and phase assemblage yields.
    #[pyo3(name = "report")]
    pub fn report(&self) -> String {
        let mut s = String::new();
        s.push_str(
            "================================================================================\n",
        );
        s.push_str(
            "                           CHEMICAL EQUILIBRIUM REPORT                          \n",
        );
        s.push_str(
            "================================================================================\n",
        );

        s.push_str("  CONDITIONS:\n");
        s.push_str(&format!(
            "    Temperature .....: {:.2} K ({:.2} °C)\n",
            self.temperature,
            self.temperature - 273.15
        ));
        s.push_str(&format!(
            "    Pressure ........: {:.5} bar ({:.1} Pa)\n",
            self.pressure / 100000.0,
            self.pressure
        ));

        s.push_str("\n  GLOBAL THERMODYNAMIC PROPERTIES:\n");
        s.push_str(&format!(
            "    Total Gibbs Energy (G) : {:.4} J\n",
            self.total_gibbs
        ));

        s.push_str("\n  PHASE ASSEMBLAGE DATA:\n");
        s.push_str(&format!(
            "    {:<15} | {:<8} | {:<12} | {:<15} | {:<15}\n",
            "Phase Name", "Status", "Amount (mol)", "Pure Gibbs (J)", "Total Gibbs (J)"
        ));
        s.push_str(
            "    ----------------------------------------------------------------------------\n",
        );

        let mut sorted_substances = self.substances.clone();
        sorted_substances.sort_by(|a, b| a.name.cmp(&b.name));

        for sub in &sorted_substances {
            let amount = self.amounts.get(&sub.name).copied().unwrap_or(0.0);
            let status = if amount > 1e-6 { "Stable" } else { "Inactive" };
            let pure_g = self.gibbs_energies.get(&sub.name).copied().unwrap_or(0.0);
            let contrib_g = amount * pure_g;

            s.push_str(&format!(
                "    {:<15} | {:<8} | {:<12.6E} | {:<15.4} | {:<15.4}\n",
                sub.name, status, amount, pure_g, contrib_g
            ));

            let mut el_strs = Vec::new();
            let mut sorted_els: Vec<_> = sub.elements.keys().collect();
            sorted_els.sort();

            for el in sorted_els {
                let coeff = sub.elements.get(el).copied().unwrap_or(0.0);
                el_strs.push(format!("{}: {}", el, coeff));
            }

            s.push_str(&format!("      Composition: {}\n", el_strs.join(", ")));
        }

        s.push_str(
            "================================================================================\n",
        );

        s
    }

    fn __repr__(&self) -> String {
        self.report()
    }

    fn __str__(&self) -> String {
        self.report()
    }
}

/// Evaluates stoichiometric equilibrium via a Python interface.
///
/// # Arguments
/// * `species` - Map of substance name to `Substance` records.
/// * `composition` - A `SystemComposition` specifying input proportions.
/// * `t` - Temperature in Kelvin ($\text{K}$).
/// * `p` - Pressure in Pascals ($\text{Pa}$, defaults to $101325.0\text{ Pa} = 1\text{ atm}$).
#[pyfunction]
#[pyo3(
    name = "equilibrate_stoichiometric",
    signature = (species, composition, t, p = 101325.0)
)]
pub fn equilibrate_stoichiometric_py(
    species: HashMap<String, Substance>,
    composition: SystemComposition,
    t: f64,
    p: f64,
) -> Equilibrium {
    let refs: Vec<&Substance> = species.values().collect();
    let b: HashMap<String, f64> = composition.into_elemental_fractions();

    equilibrate_stoichiometric(&refs, &b, t, p)
}

/// Find a particular solution to the mass balance equations A * phi = b.
/// This uses a simple gradient descent on the squared error to find ANY solution
/// that satisfies the elemental constraints, regardless of Gibbs energy.
/// Used as a fallback when the support-based minimization fails.
///
/// # Arguments
/// * `a` - Stoichiometric matrix ($A$), elements coefficient ratios.
/// * `b` - Normalized system elemental constraints.
/// * `n_s` - Total species.
/// * `n_e` - Total elements.
pub fn find_particular_solution(a: &[Vec<f64>], b: &[f64], n_s: usize, n_e: usize) -> Vec<f64> {
    let mut phi = vec![0.0; n_s];
    let lr = 0.01;

    for _ in 0..20000 {
        let mut grad = vec![0.0; n_s];

        for i in 0..n_e {
            let mut err = -b[i];
            for k in 0..n_s {
                err += a[i][k] * phi[k];
            }
            for k in 0..n_s {
                grad[k] += err * a[i][k];
            }
        }

        for k in 0..n_s {
            phi[k] -= lr * grad[k];
        }
    }

    phi
}

/// Find the global chemical equilibrium state of stoichiometric phases.
/// Minimizes the total Gibbs free energy subject to elemental conservation.
///
/// # Arguments
/// * `species` - Slices of `Substance` database references.
/// * `b` - Map of element symbol to its computed normalized mole fraction constraints.
/// * `t` - Temperature in Kelvin ($\text{K}$).
/// * `p` - Pressure in Pascals ($\text{Pa}$).
pub fn equilibrate_stoichiometric(
    species: &[&Substance],
    b: &HashMap<String, f64>,
    t: f64,
    p: f64,
) -> Equilibrium {
    let elements = extract_elements(species);
    let n_s = species.len();
    let n_e = elements.len();

    // 1. Precompute molar Gibbs energies
    let mut g_0 = vec![0.0; n_s];

    for i in 0..n_s {
        let s = species[i];
        let mut g = s.gibbs(t);

        if s.aggregation_type == AggregationType::Gas {
            g += R_GAS * t * (p / 101325.0).ln();
        }

        g_0[i] = g;
    }

    // 2. Build stoichiometry matrix A (n_e x n_s)
    let mut a = vec![vec![0.0; n_s]; n_e];

    for i in 0..n_e {
        for j in 0..n_s {
            a[i][j] = species[j]
                .elements
                .get(&elements[i])
                .copied()
                .unwrap_or(0.0);
        }
    }

    // Construct mass balance vector b_vec (n_e) matching the extracted elements
    let mut b_vec = vec![0.0; n_e];

    for i in 0..n_e {
        b_vec[i] = b.get(&elements[i]).copied().unwrap_or(0.0);
    }

    let mut best_phi = vec![0.0; n_s];
    let mut min_g = f64::INFINITY;
    let mut found_solution = false;

    // 3. Support-Based Solver
    for mask in 1..(1 << n_s) {
        let mut active = Vec::new();

        for k in 0..n_s {
            if (mask & (1 << k)) != 0 {
                active.push(k);
            }
        }

        if active.len() > n_e {
            continue;
        }

        let n_a = active.len();
        let mut m = vec![vec![0.0; n_a]; n_a];
        let mut rhs = vec![0.0; n_a];

        for i in 0..n_a {
            for j in 0..n_a {
                let mut sum = 0.0;
                for row in a.iter().take(n_e) {
                    sum += row[active[i]] * row[active[j]];
                }
                m[i][j] = sum;
            }
            m[i][i] += 1e-9;

            let mut sum_rhs = 0.0;
            for k in 0..n_e {
                sum_rhs += a[k][active[i]] * b_vec[k];
            }
            rhs[i] = sum_rhs;
        }

        if let Some(x) = dense_gaussian_solver(m, rhs) {
            let mut phi = vec![0.0; n_s];
            let mut valid = true;

            for (j, &val) in x.iter().enumerate() {
                if val < -1e-6 {
                    valid = false;
                    break;
                }
                phi[active[j]] = val;
            }

            if valid {
                let mut max_err: f64 = 0.0;

                for i in 0..n_e {
                    let mut sum = 0.0;
                    for k in 0..n_s {
                        sum += a[i][k] * phi[k];
                    }
                    max_err = max_err.max((sum - b_vec[i]).abs());
                }

                if max_err < 1e-5 {
                    let mut g = 0.0;
                    for k in 0..n_s {
                        g += phi[k] * g_0[k];
                    }
                    let score = g + (active.len() as f64) * 1e-6;

                    if score < min_g {
                        min_g = score;
                        best_phi = phi;
                        found_solution = true;
                    }
                }
            }
        }
    }

    let final_phi = if found_solution {
        best_phi
    } else {
        find_particular_solution(&a, &b_vec, n_s, n_e)
    };

    let mut amounts = HashMap::new();

    for i in 0..n_s {
        amounts.insert(species[i].name.clone(), final_phi[i]);
    }

    let mut gibbs_energies = HashMap::new();

    for i in 0..n_s {
        gibbs_energies.insert(species[i].name.clone(), g_0[i]);
    }

    let mut total_gibbs = 0.0;

    for i in 0..n_s {
        total_gibbs += final_phi[i] * g_0[i];
    }

    let substances: Vec<Substance> = species.iter().map(|&s| s.clone()).collect();

    Equilibrium {
        amounts,
        gibbs_energies,
        total_gibbs,
        temperature: t,
        pressure: p,
        substances,
    }
}

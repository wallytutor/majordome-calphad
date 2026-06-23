use crate::core::get_atomic_weight;
use std::fmt;


/// Error types encountered during mixture conversion operations.
#[derive(Debug)]
pub enum MixtureError {
    /// Raised when a given composition vector length does not match the converter elements.
    InvalidCompositionLength { expected: usize, found: usize },
    /// Raised when an element name is not recognized or not found in the atomic weights database.
    UnknownElement(String),
}


impl fmt::Display for MixtureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MixtureError::InvalidCompositionLength { expected, found } => {
                write!(
                    f,
                    "Invalid composition length: expected {}, found {}",
                    expected, found
                )
            }
            MixtureError::UnknownElement(el) => write!(f, "Unknown element: {}", el),
        }
    }
}


impl std::error::Error for MixtureError {}


/// A converter utility to translate composition vectors between mass (weight) fraction
/// and mole fraction representations for a mixture of elements.
///
/// Mass fraction (usually denoted as $Y_k$) is dimensionless (representing $\text{wt}\%$ / 100).
/// Mole fraction (usually denoted as $X_k$) is dimensionless (representing $\text{mol}\%$ / 100).
#[derive(Debug, Clone)]
pub struct FractionConverter {
    molar_masses: Vec<f64>,
}


impl FractionConverter {
    /// Creates a new `FractionConverter` for the specified list of elements.
    ///
    /// # Arguments
    /// * `elements` - A slice of element symbol strings (e.g. `&["C", "N", "Fe"]`).
    ///
    /// # Errors
    /// Returns `MixtureError::UnknownElement` if any element symbol is unrecognized in the atomic weights database.
    pub fn new(elements: &[String]) -> Result<Self, MixtureError> {
        let mut molar_masses = Vec::with_capacity(elements.len());

        for el in elements {
            if let Some(mass) = get_atomic_weight(el) {
                molar_masses.push(mass);
            } else {
                return Err(MixtureError::UnknownElement(el.clone()));
            }
        }

        Ok(Self { molar_masses })
    }

    /// Validates that a composition slice matches the expected number of elements.
    fn validate(&self, comp: &[f64]) -> Result<(), MixtureError> {
        if comp.len() != self.molar_masses.len() {
            Err(MixtureError::InvalidCompositionLength {
                expected: self.molar_masses.len(),
                found: comp.len(),
            })
        } else {
            Ok(())
        }
    }

    /// Computes the mean mixture molar mass (in $\text{g/mol}$) from a mass fraction composition.
    fn mean_molar_mass_from_mass(&self, y: &[f64]) -> f64 {
        let sum: f64 = y
            .iter()
            .zip(&self.molar_masses)
            .map(|(yk, wk)| yk / wk)
            .sum();

        1.0 / sum
    }

    /// Computes the mean mixture molar mass (in $\text{g/mol}$) from a mole fraction composition.
    fn mean_molar_mass_from_mole(&self, x: &[f64]) -> f64 {
        let sum: f64 = x
            .iter()
            .zip(&self.molar_masses)
            .map(|(xk, wk)| xk * wk)
            .sum();

        sum
    }

    /// Converts a mass fraction vector into a mole fraction vector.
    ///
    /// # Arguments
    /// * `y` - Dimensionless mass fractions (values between 0.0 and 1.0, summing to 1.0).
    ///
    /// # Errors
    /// Returns `MixtureError::InvalidCompositionLength` if the length of `y` doesn't match elements.
    pub fn mass_to_mole_fraction(&self, y: &[f64]) -> Result<Vec<f64>, MixtureError> {
        self.validate(y)?;
        let m = self.mean_molar_mass_from_mass(y);

        Ok(y.iter()
            .zip(&self.molar_masses)
            .map(|(yk, wk)| m * yk / wk)
            .collect())
    }

    /// Converts a mole fraction vector into a mass fraction vector.
    ///
    /// # Arguments
    /// * `x` - Dimensionless mole fractions (values between 0.0 and 1.0, summing to 1.0).
    ///
    /// # Errors
    /// Returns `MixtureError::InvalidCompositionLength` if the length of `x` doesn't match elements.
    pub fn mole_to_mass_fraction(&self, x: &[f64]) -> Result<Vec<f64>, MixtureError> {
        self.validate(x)?;
        let m = self.mean_molar_mass_from_mole(x);

        Ok(x.iter()
            .zip(&self.molar_masses)
            .map(|(xk, wk)| xk * wk / m)
            .collect())
    }
}


/// A converter utility to translate between mole fractions and volumetric concentrations (in $\text{g/m}^3$)
/// for interstitial species dissolved in a bulk matrix (e.g. Carbon and Nitrogen in Iron).
///
/// Concentration is in units of $\text{g/m}^3$ (grams per cubic meter).
/// Mole fraction is dimensionless.
#[derive(Debug, Clone)]
pub struct InterstitialConverter {
    /// Molar mass of the bulk metal matrix element (e.g. $55.845\text{ g/mol}$ for pure Iron).
    pub matrix_mean_molar_mass: f64,
    /// Density of the pure matrix element (e.g. $7870.0\text{ kg/m}^3 = 7.87\times 10^6\text{ g/m}^3$).
    pub matrix_free_density: f64,
}


impl InterstitialConverter {
    /// Creates a new `InterstitialConverter`.
    ///
    /// # Arguments
    /// * `matrix_mean_molar_mass` - Molar mass of bulk metal matrix element in $\text{g/mol}$.
    /// * `matrix_free_density` - Density of the pure matrix element in $\text{g/m}^3$ (note: $\text{kg/m}^3 \times 10^3$).
    pub fn new(matrix_mean_molar_mass: f64, matrix_free_density: f64) -> Self {
        Self {
            matrix_mean_molar_mass,
            matrix_free_density,
        }
    }

    /// Converts interstitial mole fractions into volumetric concentrations.
    ///
    /// # Arguments
    /// * `x` - Interstitial mole fractions (dimensionless, length $M$).
    ///
    /// # Returns
    /// A vector of concentrations in $\text{g/m}^3$.
    pub fn mole_fraction_to_concentration(&self, x: &[f64]) -> Vec<f64> {
        let sum_x: f64 = x.iter().sum();
        let matrix_mole_fraction = 1.0 - sum_x;
        let matrix_concentration =
            self.matrix_free_density / (matrix_mole_fraction * self.matrix_mean_molar_mass);

        x.iter().map(|xi| xi * matrix_concentration).collect()
    }

    /// Converts interstitial volumetric concentrations back into mole fractions.
    ///
    /// # Arguments
    /// * `c` - Volumetric concentrations in $\text{g/m}^3$, length $M$.
    ///
    /// # Returns
    /// A vector of dimensionless mole fractions.
    pub fn concentration_to_mole_fraction(&self, c: &[f64]) -> Vec<f64> {
        let matrix_concentration = self.matrix_free_density / self.matrix_mean_molar_mass;
        let sum_c: f64 = c.iter().sum();
        let total_concentration = sum_c + matrix_concentration;

        c.iter().map(|ci| ci / total_concentration).collect()
    }
}

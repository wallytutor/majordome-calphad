#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
#![deny(clippy::panic)]
#![deny(unused_must_use)]
#![deny(warnings)]

use pyo3::prelude::*;

#[cfg(test)]
mod test;

mod parameters;

pub mod core;
pub mod data;
pub mod equilibrium;
pub mod mixture;
pub mod prelude;

#[pymodule(name = "calphad")]
pub mod calphad {
    #[pymodule_export]
    use super::core::Substance;

    #[pymodule_export]
    use super::core::SystemComposition;

    #[pymodule_export]
    use super::data::DatabaseLoader;

    #[pymodule_export]
    use super::data::add_data_directory_py as add_data_directory;

    #[pymodule_export]
    use super::data::list_data_directories_py as list_data_directories;

    #[pymodule_export]
    use super::equilibrium::Equilibrium;

    #[pymodule_export]
    use super::equilibrium::equilibrate_stoichiometric_py as equilibrate_stoichiometric;
}

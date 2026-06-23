use pyo3::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::RwLock;

/// Global list of directories searched when loading data files by relative path.
pub static DATA_DIRS: RwLock<Vec<PathBuf>> = RwLock::new(Vec::new());

/// Resolve a (possibly relative) file path against the global data-directory
/// registry.  Returns the first match that exists on disk; if none matches,
/// returns the original path unchanged so callers get a meaningful OS error.
pub fn resolve_data_path(path: &str) -> PathBuf {
    let p = Path::new(path);
    if p.is_absolute() {
        return p.to_path_buf();
    }
    let dirs = DATA_DIRS.read().unwrap_or_else(|p| p.into_inner());
    for dir in dirs.iter() {
        let candidate = dir.join(p);
        if candidate.exists() {
            return candidate;
        }
    }
    // Fallback: let the caller get a natural "file not found" error.
    p.to_path_buf()
}

/// Register an additional directory to search when resolving data-file paths.
///
/// Call this before `DatabaseLoader::new` (or `load_substances_from_lua`) so
/// that relative paths such as `"data/sample/simple-calcination.lua"` are
/// looked up inside every registered directory in the order they were added.
pub fn add_data_directory_rs(dir: impl AsRef<Path>) {
    let path = dir.as_ref().to_path_buf();
    DATA_DIRS
        .write()
        .unwrap_or_else(|p| p.into_inner())
        .push(path);
}

/// List the directories that are searched when loading data files by relative path.
pub fn list_data_directories_rs() {
    let dirs = DATA_DIRS.read().unwrap_or_else(|p| p.into_inner());
    if dirs.is_empty() {
        println!("\nNo data directories registered.");
        return;
    }

    println!("\nRegistered data directories:");

    for dir in dirs.iter() {
        println!("  - {}", dir.display());
    }
}

/// Python-facing wrapper: register a search directory for data files.
///
/// # Example (Python)
///
/// ```text
/// import majordome.calphad as calphad
/// calphad.add_data_directory("/path/to/my/data")
/// db = calphad.DatabaseLoader("my-phases.lua")
/// ```
#[pyfunction(name = "add_data_directory")]
pub fn add_data_directory_py(path: String) {
    add_data_directory_rs(path);
}

/// Python-facing wrapper: list search directories for data files.
#[pyfunction(name = "list_data_directories")]
pub fn list_data_directories_py() {
    list_data_directories_rs();
}

use super::core::AggregationType;
use super::core::Parameterization;
use super::core::Substance;
use super::core::TemperatureRange;
use super::core::get_atomic_weight;
use super::core::{P_REF, R_GAS, T_REF};
use mlua::prelude::*;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::fs::read_to_string;

mod path;
pub use path::add_data_directory_py;
pub use path::list_data_directories_py;
use path::resolve_data_path;

impl FromLua for AggregationType {
    fn from_lua(value: LuaValue, _lua: &Lua) -> LuaResult<Self> {
        let s = value
            .as_string()
            .ok_or_else(|| LuaError::FromLuaConversionError {
                from: "Value",
                to: "AggregationType".to_string(),
                message: Some("Expected string".to_string()),
            })?
            .to_str()?
            .to_string();

        match s.as_str() {
            "Solid" => Ok(AggregationType::Solid),
            "Liquid" => Ok(AggregationType::Liquid),
            "Gas" => Ok(AggregationType::Gas),
            _ => Err(LuaError::FromLuaConversionError {
                from: "String",
                to: "AggregationType".to_string(),
                message: Some(format!("Unknown aggregation type: {}", s)),
            }),
        }
    }
}

impl FromLua for Parameterization {
    fn from_lua(value: LuaValue, _lua: &Lua) -> LuaResult<Self> {
        let table = value
            .as_table()
            .ok_or_else(|| LuaError::FromLuaConversionError {
                from: "Value",
                to: "Parameterization".to_string(),
                message: Some("Expected table".to_string()),
            })?;

        let model_type: String = table.get("type")?;

        match model_type.as_str() {
            "MaierKelley" => Ok(Parameterization::MaierKelley {
                a: table.get("a")?,
                b: table.get("b")?,
                c: table.get("c")?,
                h_ref: table.get("h_ref")?,
                s_ref: table.get("s_ref").unwrap_or(0.0),
            }),

            "NASA7" => Ok(Parameterization::NASA7 {
                a1: table.get("a1")?,
                a2: table.get("a2")?,
                a3: table.get("a3")?,
                a4: table.get("a4")?,
                a5: table.get("a5")?,
                a6: table.get("a6")?,
                a7: table.get("a7")?,
            }),
            "NASA9" => Ok(Parameterization::NASA9 {
                a1: table.get("a1")?,
                a2: table.get("a2")?,
                a3: table.get("a3")?,
                a4: table.get("a4")?,
                a5: table.get("a5")?,
                a6: table.get("a6")?,
                a7: table.get("a7")?,
                a8: table.get("a8")?,
                a9: table.get("a9")?,
            }),
            "Shomate" => Ok(Parameterization::Shomate {
                a: table.get("a")?,
                b: table.get("b")?,
                c: table.get("c")?,
                d: table.get("d")?,
                e: table.get("e")?,
                f: table.get("f")?,
                g: table.get("g")?,
                h: table.get("h")?,
            }),
            "GibbsPolynomial" => Ok(Parameterization::GibbsPolynomial {
                a: table.get("a")?,
                b: table.get("b")?,
                c: table.get("c")?,
                d: table.get("d")?,
                e: table.get("e")?,
                f: table.get("f")?,
                g: table.get("g")?,
            }),
            "Compound" => {
                let components_table: LuaTable = table.get("components")?;
                let mut components = Vec::new();
                for pair in components_table.sequence_values::<LuaTable>() {
                    let pair = pair?;
                    let substance: Substance = pair.get(1)?;
                    let coeff: f64 = pair.get(2)?;
                    components.push((Box::new(substance), coeff));
                }
                let deviation: Parameterization = table.get("deviation")?;
                Ok(Parameterization::Compound {
                    components,
                    deviation: Box::new(deviation),
                })
            }

            _ => Err(LuaError::FromLuaConversionError {
                from: "Table",
                to: "Parameterization".to_string(),
                message: Some(format!("Unknown parameterization type: {}", model_type)),
            }),
        }
    }
}

impl FromLua for TemperatureRange {
    fn from_lua(value: LuaValue, _lua: &Lua) -> LuaResult<Self> {
        let table = value
            .as_table()
            .ok_or_else(|| LuaError::FromLuaConversionError {
                from: "Value",
                to: "TemperatureRange".to_string(),
                message: Some("Expected table".to_string()),
            })?;

        Ok(TemperatureRange {
            t_min: table.get("t_min")?,
            t_max: table.get("t_max")?,
            model: table.get("model")?,
        })
    }
}

impl FromLua for Substance {
    fn from_lua(value: LuaValue, _lua: &Lua) -> LuaResult<Self> {
        let table = value
            .as_table()
            .ok_or_else(|| LuaError::FromLuaConversionError {
                from: "Value",
                to: "Substance".to_string(),
                message: Some("Expected table".to_string()),
            })?;

        let name: String = table.get("name")?;
        let s0: f64 = table.get("s0").unwrap_or(0.0);
        let mut ranges: Vec<TemperatureRange> = table.get("ranges")?;

        for range in &mut ranges {
            if let Parameterization::MaierKelley { ref mut s_ref, .. } = range.model
                && *s_ref == 0.0
            {
                *s_ref = s0;
            }
        }

        let elements: HashMap<String, f64> = table.get("elements")?;
        let aggregation_type: AggregationType = table
            .get("aggregation_type")
            .unwrap_or(AggregationType::Solid);

        let mut molar_mass = 0.0;
        for (symbol, count) in &elements {
            let weight = get_atomic_weight(symbol).ok_or_else(|| {
                LuaError::RuntimeError(format!("Unknown element symbol: {}", symbol))
            })?;
            molar_mass += weight * count;
        }

        let molar_volume = if aggregation_type == AggregationType::Gas {
            (R_GAS * T_REF / P_REF) * 1e6
        } else {
            table.get("molar_volume").unwrap_or(0.0)
        };

        Ok(Substance {
            name,
            molar_mass,
            molar_volume,
            s0: table.get("s0").unwrap_or(0.0),
            ranges,
            elements,
            reference: table.get("reference").unwrap_or_else(|_| "".to_string()),
            aggregation_type,
        })
    }
}

fn get_coeff(t: &LuaTable, i: usize, name: &str) -> LuaResult<f64> {
    match t.get::<f64>(i) {
        Ok(v) => Ok(v),
        Err(_) => t.get::<f64>(name),
    }
}

fn create_lua_maier_kelley(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|lua, args: LuaMultiValue| {
        let table = lua.create_table()?;
        table.set("type", "MaierKelley")?;
        if args.len() == 1
            && let Some(t) = args[0].as_table()
        {
            table.set("a", get_coeff(t, 1, "a")?)?;
            table.set("b", get_coeff(t, 2, "b")?)?;
            table.set("c", get_coeff(t, 3, "c")?)?;
            table.set("h_ref", get_coeff(t, 4, "h_ref")?)?;
            return Ok(table);
        }
        let mut iter = args.into_iter();
        table.set(
            "a",
            f64::from_lua(iter.next().unwrap_or(LuaValue::Nil), lua)?,
        )?;
        table.set(
            "b",
            f64::from_lua(iter.next().unwrap_or(LuaValue::Nil), lua)?,
        )?;
        table.set(
            "c",
            f64::from_lua(iter.next().unwrap_or(LuaValue::Nil), lua)?,
        )?;
        table.set(
            "h_ref",
            f64::from_lua(iter.next().unwrap_or(LuaValue::Nil), lua)?,
        )?;
        Ok(table)
    })
}

fn create_lua_nasa7(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|lua, args: LuaMultiValue| {
        let table = lua.create_table()?;
        table.set("type", "NASA7")?;
        if args.len() == 1
            && let Some(t) = args[0].as_table()
        {
            for i in 1..=7 {
                let name = format!("a{}", i);
                table.set(name.clone(), get_coeff(t, i, &name)?)?;
            }
            return Ok(table);
        }
        for (i, arg) in args.into_iter().enumerate() {
            if i < 7 {
                table.set(format!("a{}", i + 1), f64::from_lua(arg, lua)?)?;
            }
        }
        Ok(table)
    })
}

fn create_lua_nasa9(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|lua, args: LuaMultiValue| {
        let table = lua.create_table()?;
        table.set("type", "NASA9")?;
        if args.len() == 1
            && let Some(t) = args[0].as_table()
        {
            for i in 1..=9 {
                let name = format!("a{}", i);
                table.set(name.clone(), get_coeff(t, i, &name)?)?;
            }
            return Ok(table);
        }
        for (i, arg) in args.into_iter().enumerate() {
            if i < 9 {
                table.set(format!("a{}", i + 1), f64::from_lua(arg, lua)?)?;
            }
        }
        Ok(table)
    })
}

fn create_lua_shomate(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|lua, args: LuaMultiValue| {
        let table = lua.create_table()?;
        table.set("type", "Shomate")?;
        let keys = ["a", "b", "c", "d", "e", "f", "g", "h"];
        if args.len() == 1
            && let Some(t) = args[0].as_table()
        {
            for (i, key) in keys.iter().enumerate() {
                table.set(*key, get_coeff(t, i + 1, key)?)?;
            }
            return Ok(table);
        }
        for (i, arg) in args.into_iter().enumerate() {
            if i < 8 {
                table.set(keys[i], f64::from_lua(arg, lua)?)?;
            }
        }
        Ok(table)
    })
}

fn create_lua_gibbs_polynomial(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|lua, args: LuaMultiValue| {
        let table = lua.create_table()?;
        table.set("type", "GibbsPolynomial")?;
        let keys = ["a", "b", "c", "d", "e", "f", "g"];
        if args.len() == 1
            && let Some(t) = args[0].as_table()
        {
            for (i, key) in keys.iter().enumerate() {
                table.set(*key, get_coeff(t, i + 1, key)?)?;
            }
            return Ok(table);
        }
        for (i, arg) in args.into_iter().enumerate() {
            if i < 7 {
                table.set(keys[i], f64::from_lua(arg, lua)?)?;
            }
        }
        Ok(table)
    })
}

fn create_lua_compound(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|_lua, table: LuaTable| {
        if !table.contains_key("type")? {
            table.set("type", "Compound")?;
        }
        if !table.contains_key("deviation")? {
            let dev = _lua.create_table()?;
            dev.set("type", "GibbsPolynomial")?;
            dev.set("a", 0.0)?;
            dev.set("b", 0.0)?;
            dev.set("c", 0.0)?;
            dev.set("d", 0.0)?;
            dev.set("e", 0.0)?;
            dev.set("f", 0.0)?;
            dev.set("g", 0.0)?;
            table.set("deviation", dev)?;
        }
        Ok(table)
    })
}

fn create_lua_temperature_range(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|lua, (t_min, t_max, model): (f64, f64, LuaValue)| {
        let table = lua.create_table()?;
        table.set("t_min", t_min)?;
        table.set("t_max", t_max)?;
        table.set("model", model)?;
        Ok(table)
    })
}

fn create_lua_substance(lua: &Lua) -> LuaResult<LuaFunction> {
    lua.create_function(|_lua, table: LuaTable| {
        if !table.contains_key("reference")? {
            table.set("reference", "")?;
        }
        if !table.contains_key("aggregation_type")? {
            table.set("aggregation_type", "Solid")?;
        }
        Ok(table)
    })
}

/// Loads a database of substances from a specialized Lua data file.
///
/// # Arguments
/// * `path` - Path string to the Lua script file containing substance tables.
pub fn load_substances_from_lua(path: &str) -> LuaResult<HashMap<String, Substance>> {
    let resolved = resolve_data_path(path);

    let lua = Lua::new();
    let globals = lua.globals();

    let mk_fn = create_lua_maier_kelley(&lua)?;
    globals.set("MaierKelley", mk_fn)?;

    let nasa7_fn = create_lua_nasa7(&lua)?;
    globals.set("NASA7", nasa7_fn)?;

    let nasa9_fn = create_lua_nasa9(&lua)?;
    globals.set("NASA9", nasa9_fn)?;

    let shomate_fn = create_lua_shomate(&lua)?;
    globals.set("Shomate", shomate_fn)?;

    let gibbs_polynomial_fn = create_lua_gibbs_polynomial(&lua)?;
    globals.set("GibbsPolynomial", gibbs_polynomial_fn)?;

    let compound_fn = create_lua_compound(&lua)?;
    globals.set("Compound", compound_fn)?;

    let temperature_range_fn = create_lua_temperature_range(&lua)?;

    globals.set("TemperatureRange", temperature_range_fn.clone())?;
    globals.set("Range", temperature_range_fn)?;

    let substance_fn = create_lua_substance(&lua)?;
    globals.set("Substance", substance_fn)?;

    let content =
        read_to_string(&resolved).map_err(|e| LuaError::ExternalError(std::sync::Arc::new(e)))?;

    let map: HashMap<String, Substance> = lua.load(&content).eval()?;
    Ok(map)
}

/// Coordinates loading and retrieving substance thermodynamic database files (Lua formatted).
/// Supports filtering active thermodynamic phases.
#[pyclass]
pub struct DatabaseLoader {
    pub path: String,
    pub phases: Vec<String>,
    pub data: HashMap<String, Substance>,
}

#[pymethods]
impl DatabaseLoader {
    /// Loads a database from the given file path.
    ///
    /// # Arguments
    /// * `path` - The database filename or absolute file path.
    /// * `phases` - An optional filter list of phase names to load (e.g. `Some(vec!["Fe(s)", "Fe(l)"])`).
    ///
    /// # Errors
    /// Returns a `PyValueError` if database loading or parsing fails.
    #[new]
    #[pyo3(signature = (path, phases = None))]
    pub fn new(path: String, phases: Option<Vec<String>>) -> PyResult<Self> {
        let resolved_path = resolve_data_path(&path).to_string_lossy().into_owned();

        let mut raw_data = load_substances_from_lua(&resolved_path)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let phases = if let Some(ref selected_phases) = phases {
            raw_data.retain(|k, _| selected_phases.contains(k));
            selected_phases.clone()
        } else {
            let mut p = raw_data.keys().cloned().collect::<Vec<String>>();
            p.sort();
            p
        };

        Ok(Self {
            path,
            phases,
            data: raw_data,
        })
    }

    /// Gets the loaded database file path.
    #[getter]
    pub fn path(&self) -> String {
        self.path.clone()
    }

    /// Gets the list of loaded/filtered active phases.
    #[getter]
    pub fn phases(&self) -> Vec<String> {
        self.phases.clone()
    }

    /// Returns the loaded database as a Python dictionary map of `Substance` records.
    pub fn get_data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        for phase_name in &self.phases {
            if let Some(v) = self.data.get(phase_name) {
                dict.set_item(phase_name, v.clone())?;
            }
        }

        Ok(dict)
    }

    /// Loads a single pure substance by name from the database map.
    ///
    /// # Errors
    /// Returns `PyKeyError` if the compound name is not found in the database.
    pub fn load_compound(&self, name: String) -> PyResult<Substance> {
        self.data.get(&name).cloned().ok_or_else(|| {
            pyo3::exceptions::PyKeyError::new_err(format!(
                "Compound '{}' not found in database",
                name
            ))
        })
    }
}

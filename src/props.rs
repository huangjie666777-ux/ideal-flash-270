//! Component data, system validation and pure-component property evaluation.

use crate::FlashError;
use serde::{Deserialize, Serialize};

/// Reference temperature for the enthalpy correlations [K].
pub const T_REF: f64 = 298.15;

/// Antoine parameters with validity range.
/// log10(Psat/bar) = a - b / (T + c), T in [t_min, t_max] (K).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Antoine {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub t_min: f64,
    pub t_max: f64,
}

/// One feed component.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Component {
    pub name: String,
    /// Feed mole fraction, must be > 0; all z must sum to 1.
    pub z: f64,
    pub antoine: Antoine,
    /// Constant liquid heat capacity [J/(mol*K)], must be > 0.
    pub cp_l: f64,
    /// Constant vapor heat capacity [J/(mol*K)], must be > 0.
    pub cp_v: f64,
    /// Reference heat of vaporization at T_REF [J/mol], must be > 0.
    pub l_ref: f64,
}

/// Validated flash system: components plus the common validity range.
#[derive(Debug, Clone)]
pub struct FlashSystem {
    pub components: Vec<Component>,
    /// Intersection of all Antoine validity ranges [K].
    pub t_common_min: f64,
    pub t_common_max: f64,
}

fn require_finite(value: f64, what: &str) -> Result<(), FlashError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(FlashError::InvalidInput(format!("{what} is not finite")))
    }
}

/// Validate raw component data and build a FlashSystem.
pub fn validate_system(components: Vec<Component>) -> Result<FlashSystem, FlashError> {
    let n = components.len();
    if !(2..=16).contains(&n) {
        return Err(FlashError::InvalidInput(format!(
            "component count {n} out of range 2..=16"
        )));
    }
    for (i, c) in components.iter().enumerate() {
        if c.name.trim().is_empty() {
            return Err(FlashError::InvalidInput(format!(
                "component {i} has an empty name"
            )));
        }
        if components[..i].iter().any(|o| o.name == c.name) {
            return Err(FlashError::InvalidInput(format!(
                "duplicate component name {:?}",
                c.name
            )));
        }
        require_finite(c.z, "z")?;
        if c.z <= 0.0 {
            return Err(FlashError::InvalidInput(format!(
                "component {:?}: z must be positive",
                c.name
            )));
        }
        let a = &c.antoine;
        for (v, what) in [
            (a.a, "antoine.a"),
            (a.b, "antoine.b"),
            (a.c, "antoine.c"),
            (a.t_min, "antoine.t_min"),
            (a.t_max, "antoine.t_max"),
            (c.cp_l, "cp_l"),
            (c.cp_v, "cp_v"),
            (c.l_ref, "l_ref"),
        ] {
            require_finite(v, what)?;
        }
        if a.t_min >= a.t_max {
            return Err(FlashError::InvalidInput(format!(
                "component {:?}: antoine t_min must be < t_max",
                c.name
            )));
        }
        if c.cp_l <= 0.0 || c.cp_v <= 0.0 {
            return Err(FlashError::InvalidInput(format!(
                "component {:?}: cp_l and cp_v must be positive",
                c.name
            )));
        }
        if c.l_ref <= 0.0 {
            return Err(FlashError::InvalidInput(format!(
                "component {:?}: l_ref must be positive",
                c.name
            )));
        }
    }
    let z_sum: f64 = components.iter().map(|c| c.z).sum();
    if (z_sum - 1.0).abs() > 1e-9 {
        return Err(FlashError::InvalidInput(format!(
            "mole fractions sum to {z_sum}, expected 1"
        )));
    }

    let t_common_min = components
        .iter()
        .map(|c| c.antoine.t_min)
        .fold(f64::NEG_INFINITY, f64::max);
    let t_common_max = components
        .iter()
        .map(|c| c.antoine.t_max)
        .fold(f64::INFINITY, f64::min);
    if t_common_min >= t_common_max {
        return Err(FlashError::InvalidInput(
            "no common Antoine validity range across components".into(),
        ));
    }

    // T + C must stay positive over the common range (no extrapolation).
    for c in &components {
        if t_common_min + c.antoine.c <= 0.0 {
            return Err(FlashError::InvalidInput(format!(
                "component {:?}: T + C <= 0 within common range",
                c.name
            )));
        }
    }
    // hV - hL = L + (CpV - CpL)(T - T_REF) must stay positive over the
    // common range; its minimum over the range is at one of the endpoints.
    for c in &components {
        let dcp = c.cp_v - c.cp_l;
        for &t in &[t_common_min, t_common_max] {
            if c.l_ref + dcp * (t - T_REF) <= 0.0 {
                return Err(FlashError::InvalidInput(format!(
                    "component {:?}: hV - hL <= 0 within common range",
                    c.name
                )));
            }
        }
    }

    Ok(FlashSystem {
        components,
        t_common_min,
        t_common_max,
    })
}

impl FlashSystem {
    /// Ensure T lies inside the common validity range.
    pub fn check_temperature(&self, t: f64) -> Result<(), FlashError> {
        if !t.is_finite() {
            return Err(FlashError::InvalidInput("temperature is not finite".into()));
        }
        if t < self.t_common_min || t > self.t_common_max {
            return Err(FlashError::OutOfRange(format!(
                "T = {t} K outside common range [{}, {}] K",
                self.t_common_min, self.t_common_max
            )));
        }
        Ok(())
    }

    /// Saturation pressure [bar] of component i at T (validity checked).
    pub fn psat(&self, i: usize, t: f64) -> Result<f64, FlashError> {
        self.check_temperature(t)?;
        let a = &self.components[i].antoine;
        Ok(10f64.powf(a.a - a.b / (t + a.c)))
    }

    /// Liquid molar enthalpy of component i [J/mol].
    pub fn h_l(&self, i: usize, t: f64) -> f64 {
        self.components[i].cp_l * (t - T_REF)
    }

    /// Vapor molar enthalpy of component i [J/mol].
    pub fn h_v(&self, i: usize, t: f64) -> f64 {
        self.components[i].l_ref + self.components[i].cp_v * (t - T_REF)
    }
}

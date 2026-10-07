//! ideal_flash270: ideal vapor-liquid flash (TP and PH) for 2..=16 components.
//!
//! Units: temperature K, pressure bar, enthalpy J/mol.
//!   log10(Psat/bar) = A - B / (T + C)
//!   hL = CpL * (T - 298.15)
//!   hV = L + CpV * (T - 298.15)

pub mod energy;
pub mod equilibrium;
pub mod props;
pub mod serde_io;

pub use energy::ph_flash;
pub use equilibrium::tp_flash;
pub use props::{validate_system, Antoine, Component, FlashSystem, T_REF};
pub use serde_io::{flash_ph_json, flash_tp_json};

use serde::Serialize;
use std::fmt;

/// Errors returned by validation and flash calculations.
#[derive(Debug, Clone, PartialEq)]
pub enum FlashError {
    InvalidInput(String),
    OutOfRange(String),
    DegenerateEquilibrium(String),
    NotBracketed(String),
    NotConverged(String),
    Json(String),
}

impl fmt::Display for FlashError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FlashError::InvalidInput(m) => write!(f, "invalid input: {m}"),
            FlashError::OutOfRange(m) => write!(f, "property out of range: {m}"),
            FlashError::DegenerateEquilibrium(m) => write!(f, "degenerate equilibrium: {m}"),
            FlashError::NotBracketed(m) => write!(f, "target not bracketed: {m}"),
            FlashError::NotConverged(m) => write!(f, "iteration did not converge: {m}"),
            FlashError::Json(m) => write!(f, "json error: {m}"),
        }
    }
}

impl std::error::Error for FlashError {}

impl Serialize for FlashError {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&self.to_string())
    }
}

/// Phase state of a flash result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Phase {
    Liquid,
    Vapor,
    TwoPhase,
}

/// Result of a flash calculation. Compositions follow input component order;
/// a missing phase has an empty composition vector.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FlashResult {
    /// Equilibrium temperature [K].
    pub temperature: f64,
    /// Pressure [bar].
    pub pressure: f64,
    /// Phase state.
    pub phase: Phase,
    /// Vapor fraction beta (molar), 0.0 for liquid, 1.0 for vapor.
    pub beta: f64,
    /// Liquid-phase mole fractions (empty if no liquid phase).
    pub x: Vec<f64>,
    /// Vapor-phase mole fractions (empty if no vapor phase).
    pub y: Vec<f64>,
    /// Mixture molar enthalpy [J/mol].
    pub enthalpy: f64,
}

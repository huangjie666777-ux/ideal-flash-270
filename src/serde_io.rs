//! JSON input/output for TP and PH flashes.

use crate::energy::ph_flash;
use crate::equilibrium::tp_flash;
use crate::props::{validate_system, Component};
use crate::{FlashError, FlashResult};
use serde::{Deserialize, Serialize};

/// JSON request for a TP flash.
#[derive(Debug, Deserialize)]
pub struct TpRequest {
    pub components: Vec<Component>,
    /// Temperature [K].
    pub temperature: f64,
    /// Pressure [bar].
    pub pressure: f64,
}

/// JSON request for a PH flash.
#[derive(Debug, Deserialize)]
pub struct PhRequest {
    pub components: Vec<Component>,
    /// Pressure [bar].
    pub pressure: f64,
    /// Target mixture molar enthalpy [J/mol].
    pub enthalpy: f64,
    /// Temperature bracket lower bound [K].
    pub t_lo: f64,
    /// Temperature bracket upper bound [K].
    pub t_hi: f64,
}

/// JSON response: either a result or an error message.
#[derive(Debug, Serialize)]
#[serde(untagged)]
pub enum FlashResponse {
    Ok { result: FlashResult },
    Err { error: FlashError },
}

fn respond(r: Result<FlashResult, FlashError>) -> String {
    let resp = match r {
        Ok(result) => FlashResponse::Ok { result },
        Err(error) => FlashResponse::Err { error },
    };
    serde_json::to_string_pretty(&resp).unwrap_or_else(|e| {
        format!(r#"{{"error":"json error: serialization failed: {e}"}}"#)
    })
}

/// Run a TP flash from a JSON string, returning a JSON string.
pub fn flash_tp_json(input: &str) -> String {
    respond(flash_tp(input))
}

/// Run a PH flash from a JSON string, returning a JSON string.
pub fn flash_ph_json(input: &str) -> String {
    respond(flash_ph(input))
}

/// Typed TP flash from JSON input.
pub fn flash_tp(input: &str) -> Result<FlashResult, FlashError> {
    let req: TpRequest =
        serde_json::from_str(input).map_err(|e| FlashError::Json(e.to_string()))?;
    let sys = validate_system(req.components)?;
    tp_flash(&sys, req.temperature, req.pressure)
}

/// Typed PH flash from JSON input.
pub fn flash_ph(input: &str) -> Result<FlashResult, FlashError> {
    let req: PhRequest =
        serde_json::from_str(input).map_err(|e| FlashError::Json(e.to_string()))?;
    let sys = validate_system(req.components)?;
    ph_flash(&sys, req.pressure, req.enthalpy, req.t_lo, req.t_hi)
}

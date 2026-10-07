use crate::energy::ph_flash;
use crate::equilibrium::{tp_flash, FlashState, Phase};
use crate::error::{FlashError, FlashResult};
use crate::properties::{Component, Mixture};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize)]
pub struct ComponentInput {
    pub name: String,
    #[serde(alias = "Z")]
    pub z: f64,
    #[serde(alias = "A")]
    pub a: f64,
    #[serde(alias = "B")]
    pub b: f64,
    #[serde(alias = "C")]
    pub c: f64,
    #[serde(alias = "tmin")]
    pub t_min: f64,
    #[serde(alias = "tmax")]
    pub t_max: f64,
    #[serde(alias = "CpL")]
    pub cp_l: f64,
    #[serde(alias = "CpV")]
    pub cp_v: f64,
    #[serde(alias = "L", alias = "l")]
    pub latent_heat: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct TpRequest {
    pub components: Vec<ComponentInput>,
    #[serde(alias = "T")]
    pub temperature: f64,
    #[serde(alias = "P")]
    pub pressure: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct PhRequest {
    pub components: Vec<ComponentInput>,
    #[serde(alias = "P")]
    pub pressure: f64,
    #[serde(alias = "H")]
    pub enthalpy: f64,
    #[serde(alias = "T_lo", alias = "Tlo")]
    pub t_lo: f64,
    #[serde(alias = "T_hi", alias = "Thi")]
    pub t_hi: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct FlashResponse {
    pub temperature: f64,
    pub pressure: f64,
    pub phase: &'static str,
    pub beta: f64,
    pub component_names: Vec<String>,
    pub liquid_composition: Vec<f64>,
    pub vapor_composition: Vec<f64>,
    pub enthalpy: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ErrorResponse {
    pub error: String,
    pub kind: &'static str,
}

pub fn build_mixture(inputs: Vec<ComponentInput>) -> FlashResult<Mixture> {
    Mixture::new(inputs.into_iter().map(into_component).collect())
}

fn into_component(input: ComponentInput) -> Component {
    Component {
        name: input.name,
        z: input.z,
        a: input.a,
        b: input.b,
        c: input.c,
        t_min: input.t_min,
        t_max: input.t_max,
        cp_l: input.cp_l,
        cp_v: input.cp_v,
        latent_heat: input.latent_heat,
    }
}

fn phase_name(phase: Phase) -> &'static str {
    match phase {
        Phase::Liquid => "liquid",
        Phase::Vapor => "vapor",
        Phase::TwoPhase => "two_phase",
    }
}

fn response(mixture: &Mixture, state: FlashState) -> FlashResponse {
    FlashResponse {
        temperature: state.temperature,
        pressure: state.pressure,
        phase: phase_name(state.phase),
        beta: state.beta,
        component_names: mixture.components.iter().map(|c| c.name.clone()).collect(),
        liquid_composition: state.liquid_composition,
        vapor_composition: state.vapor_composition,
        enthalpy: state.enthalpy,
    }
}

pub fn tp_from_request(request: TpRequest) -> FlashResult<FlashResponse> {
    let mixture = build_mixture(request.components)?;
    let state = tp_flash(&mixture, request.temperature, request.pressure)?;
    Ok(response(&mixture, state))
}

pub fn ph_from_request(request: PhRequest) -> FlashResult<FlashResponse> {
    let mixture = build_mixture(request.components)?;
    let state = ph_flash(
        &mixture,
        request.pressure,
        request.enthalpy,
        request.t_lo,
        request.t_hi,
    )?;
    Ok(response(&mixture, state))
}

pub fn run_json(input: &str) -> String {
    match run_json_value(input) {
        Ok(value) => serde_json::to_string_pretty(&value).unwrap_or_else(|error| {
            serde_json::to_string(&ErrorResponse {
                error: error.to_string(),
                kind: "json",
            })
            .unwrap()
        }),
        Err(error) => {
            let kind = match error {
                FlashError::InvalidInput(_) => "invalid_input",
                FlashError::OutOfRange(_) => "out_of_range",
                FlashError::Degenerate(_) => "degenerate",
                FlashError::NoBracket(_) => "no_bracket",
                FlashError::NotConverged(_) => "not_converged",
                FlashError::Json(_) => "json",
            };
            serde_json::to_string_pretty(&ErrorResponse {
                error: error.to_string(),
                kind,
            })
            .unwrap()
        }
    }
}

fn run_json_value(input: &str) -> FlashResult<FlashResponse> {
    let value: serde_json::Value = serde_json::from_str(input)?;
    match value.get("kind").and_then(|v| v.as_str()) {
        Some("tp") => tp_from_request(serde_json::from_value(value)?),
        Some("ph") => ph_from_request(serde_json::from_value(value)?),
        Some(other) => Err(FlashError::InvalidInput(format!(
            "unknown calculation kind '{other}'"
        ))),
        None => Err(FlashError::InvalidInput("missing calculation kind".into())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TP_JSON: &str = r#"{
        "kind": "tp",
        "temperature": 278.0,
        "pressure": 1.0,
        "components": [
            {"name":"light","z":0.5,"a":2.5,"b":900.0,"c":100.0,"t_min":250.0,"t_max":450.0,"cp_l":100.0,"cp_v":80.0,"L":30000.0},
            {"name":"heavy","z":0.5,"a":2.5,"b":1000.0,"c":100.0,"t_min":250.0,"t_max":450.0,"cp_l":100.0,"cp_v":80.0,"L":30000.0}
        ]
    }"#;

    #[test]
    fn json_tp_round_trip() {
        let output = run_json(TP_JSON);
        let value: serde_json::Value = serde_json::from_str(&output).unwrap();
        assert_eq!(value["phase"], "two_phase");
        assert_eq!(value["component_names"][0], "light");
        assert_eq!(value["liquid_composition"].as_array().unwrap().len(), 2);
        assert_eq!(value["vapor_composition"].as_array().unwrap().len(), 2);
    }

    #[test]
    fn json_ph_accepts_title_case_aliases() {
        let input = r#"{
            "kind": "ph",
            "P": 1.0,
            "H": 10000.0,
            "T_lo": 250.0,
            "T_hi": 300.0,
            "components": [
                {"name":"light","Z":0.5,"A":2.5,"B":900.0,"C":100.0,"t_min":250.0,"t_max":450.0,"CpL":100.0,"CpV":80.0,"L":30000.0},
                {"name":"heavy","Z":0.5,"A":2.5,"B":1000.0,"C":100.0,"t_min":250.0,"t_max":450.0,"CpL":100.0,"CpV":80.0,"L":30000.0}
            ]
        }"#;
        let value: serde_json::Value = serde_json::from_str(&run_json(input)).unwrap();
        assert!((value["temperature"].as_f64().unwrap() - 279.3348924952).abs() < 1.0e-8);
        assert_eq!(value["phase"], "two_phase");
    }
}

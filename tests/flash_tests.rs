use ideal_flash270::props::{validate_system, Antoine, Component};
use ideal_flash270::{flash_ph_json, flash_tp_json, ph_flash, tp_flash, FlashError, Phase};

fn comp(name: &str, z: f64, a: f64, b: f64, c: f64, tmin: f64, tmax: f64,
        cpl: f64, cpv: f64, l: f64) -> Component {
    Component {
        name: name.into(),
        z,
        antoine: Antoine { a, b, c, t_min: tmin, t_max: tmax },
        cp_l: cpl,
        cp_v: cpv,
        l_ref: l,
    }
}

fn benz_tol() -> Vec<Component> {
    vec![
        comp("benzene", 0.4, 4.01814, 1203.835, -53.226, 280.0, 380.0, 136.0, 82.0, 30760.0),
        comp("toluene", 0.6, 4.07827, 1343.943, -53.773, 280.0, 410.0, 157.0, 103.0, 33180.0),
    ]
}

fn sys() -> ideal_flash270::FlashSystem {
    validate_system(benz_tol()).unwrap()
}

#[test]
fn validation_rejects_bad_input() {
    // too few components
    assert!(validate_system(vec![benz_tol().remove(0)]).is_err());
    // duplicate names
    let mut v = benz_tol();
    v[1].name = "benzene".into();
    assert!(validate_system(v).is_err());
    // z not summing to 1
    let mut v = benz_tol();
    v[0].z = 0.5;
    assert!(validate_system(v).is_err());
    // non-positive z
    let mut v = benz_tol();
    v[0].z = -0.1;
    v[1].z = 1.1;
    assert!(validate_system(v).is_err());
    // non-finite value
    let mut v = benz_tol();
    v[0].cp_l = f64::NAN;
    assert!(validate_system(v).is_err());
    // non-positive Cp
    let mut v = benz_tol();
    v[1].cp_v = 0.0;
    assert!(validate_system(v).is_err());
    // T + C <= 0 in common range
    let mut v = benz_tol();
    v[0].antoine.t_min = 50.0;
    v[0].antoine.c = -53.226;
    v[1].antoine.t_min = 50.0;
    assert!(validate_system(v).is_err());
    // hV - hL <= 0 in common range
    let mut v = benz_tol();
    v[0].l_ref = 100.0;
    v[0].cp_v = 1.0;
    assert!(validate_system(v).is_err());
}

#[test]
fn tp_two_phase_conserves_and_normalizes() {
    let s = sys();
    let r = tp_flash(&s, 370.0, 1.01325).unwrap();
    assert_eq!(r.phase, Phase::TwoPhase);
    assert!(r.beta > 0.0 && r.beta < 1.0);
    let zs = [0.4, 0.6];
    for i in 0..2 {
        let mix = (1.0 - r.beta) * r.x[i] + r.beta * r.y[i];
        assert!((mix - zs[i]).abs() < 1e-9, "component balance failed");
    }
    let sx: f64 = r.x.iter().sum();
    let sy: f64 = r.y.iter().sum();
    assert!((sx - 1.0).abs() < 1e-9);
    assert!((sy - 1.0).abs() < 1e-9);
}

#[test]
fn tp_single_phases() {
    let s = sys();
    // Low T: subcooled liquid.
    let r = tp_flash(&s, 290.0, 5.0).unwrap();
    assert_eq!(r.phase, Phase::Liquid);
    assert_eq!(r.beta, 0.0);
    assert!(r.y.is_empty());
    assert_eq!(r.x, vec![0.4, 0.6]);
    // High T, low P: superheated vapor.
    let r = tp_flash(&s, 375.0, 0.5).unwrap();
    assert_eq!(r.phase, Phase::Vapor);
    assert_eq!(r.beta, 1.0);
    assert!(r.x.is_empty());
    assert_eq!(r.y, vec![0.4, 0.6]);
}

#[test]
fn tp_degenerate_all_k_one_errors() {
    // Two identical components => all K = 1 at any T.
    let mut v = benz_tol();
    v[1] = v[0].clone();
    v[1].name = "benzene2".into();
    v[1].z = 0.6;
    let s = validate_system(v).unwrap();
    // With identical components, K_i = 1 for all i when P = Psat(T).
    let p = s.psat(0, 370.0).unwrap();
    let err = tp_flash(&s, 370.0, p).unwrap_err();
    assert!(matches!(err, FlashError::DegenerateEquilibrium(_)));
}

#[test]
fn tp_out_of_range_errors() {
    let s = sys();
    let err = tp_flash(&s, 500.0, 1.0).unwrap_err();
    assert!(matches!(err, FlashError::OutOfRange(_)));
}

#[test]
fn ph_recovers_tp_temperature() {
    let s = sys();
    let tp = tp_flash(&s, 370.0, 1.01325).unwrap();
    let r = ph_flash(&s, 1.01325, tp.enthalpy, 280.0, 410.0).unwrap();
    assert!((r.temperature - 370.0).abs() < 1e-6);
    assert_eq!(r.phase, Phase::TwoPhase);
    assert!((r.enthalpy - tp.enthalpy).abs() < 1e-3);
}

#[test]
fn ph_crosses_phase_boundaries() {
    let s = sys();
    // Target a single-phase liquid state reached from a bracket that
    // spans the two-phase region.
    let tp = tp_flash(&s, 300.0, 5.0).unwrap();
    assert_eq!(tp.phase, Phase::Liquid);
    let r = ph_flash(&s, 5.0, tp.enthalpy, 280.0, 410.0).unwrap();
    assert!((r.temperature - 300.0).abs() < 1e-6);
    assert_eq!(r.phase, Phase::Liquid);
}

#[test]
fn ph_unbracketed_target_errors() {
    let s = sys();
    // Huge target enthalpy, above anything in the bracket.
    let err = ph_flash(&s, 1.0, 1.0e9, 280.0, 410.0).unwrap_err();
    assert!(matches!(err, FlashError::NotBracketed(_)));
}

#[test]
fn ph_endpoint_not_accepted() {
    let s = sys();
    // Target exactly the enthalpy at the bracket endpoint.
    let tp = tp_flash(&s, 280.0, 1.01325).unwrap();
    let err = ph_flash(&s, 1.01325, tp.enthalpy, 280.0, 410.0).unwrap_err();
    assert!(matches!(err, FlashError::NotBracketed(_)));
}

#[test]
fn json_roundtrip_tp_and_ph() {
    let input = r#"{
        "components": [
            {"name":"benzene","z":0.4,
             "antoine":{"a":4.01814,"b":1203.835,"c":-53.226,"t_min":280.0,"t_max":380.0},
             "cp_l":136.0,"cp_v":82.0,"l_ref":30760.0},
            {"name":"toluene","z":0.6,
             "antoine":{"a":4.07827,"b":1343.943,"c":-53.773,"t_min":280.0,"t_max":410.0},
             "cp_l":157.0,"cp_v":103.0,"l_ref":33180.0}
        ],
        "temperature": 370.0, "pressure": 1.01325
    }"#;
    let out = flash_tp_json(input);
    let v: serde_json::Value = serde_json::from_str(&out).unwrap();
    assert_eq!(v["result"]["phase"], "two_phase");
    let h = v["result"]["enthalpy"].as_f64().unwrap();

    let ph_input = format!(
        r#"{{"components":{},"pressure":1.01325,"enthalpy":{h},"t_lo":280.0,"t_hi":410.0}}"#,
        serde_json::to_string(&serde_json::from_str::<serde_json::Value>(input).unwrap()["components"]).unwrap()
    );
    let out = flash_ph_json(&ph_input);
    let v: serde_json::Value = serde_json::from_str(&out).unwrap();
    assert!((v["result"]["temperature"].as_f64().unwrap() - 370.0).abs() < 1e-6);
}

#[test]
fn json_error_reported() {
    let out = flash_tp_json("{not json");
    assert!(out.contains("error"));
    let out = flash_tp_json(r#"{"components":[],"temperature":300.0,"pressure":1.0}"#);
    assert!(out.contains("error"));
}

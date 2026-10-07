//! Demo: TP and PH flashes for a benzene/toluene mixture, via JSON.

use ideal_flash270::{flash_ph_json, flash_tp_json};

fn main() {
    // Benzene / toluene, Antoine constants (Psat in bar, T in K).
    let components = r#"[
        {"name":"benzene","z":0.4,
         "antoine":{"a":4.01814,"b":1203.835,"c":-53.226,"t_min":280.0,"t_max":380.0},
         "cp_l":136.0,"cp_v":82.0,"l_ref":30760.0},
        {"name":"toluene","z":0.6,
         "antoine":{"a":4.07827,"b":1343.943,"c":-53.773,"t_min":280.0,"t_max":410.0},
         "cp_l":157.0,"cp_v":103.0,"l_ref":33180.0}
    ]"#;

    let tp_input = format!(
        r#"{{"components":{components},"temperature":370.0,"pressure":1.01325}}"#
    );
    println!("=== TP flash (T=370 K, P=1.01325 bar) ===");
    let tp_out = flash_tp_json(&tp_input);
    println!("{tp_out}");

    // Use the TP mixture enthalpy as the PH target: should recover 360 K.
    let v: serde_json::Value = serde_json::from_str(&tp_out).unwrap();
    let h = v["result"]["enthalpy"].as_f64().unwrap();
    let ph_input = format!(
        r#"{{"components":{components},"pressure":1.01325,"enthalpy":{h},"t_lo":280.0,"t_hi":410.0}}"#
    );
    println!("=== PH flash (P=1.01325 bar, h={h:.3} J/mol) ===");
    println!("{}", flash_ph_json(&ph_input));
}

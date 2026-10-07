use ideal_flash270::{ph_flash, run_json, tp_flash, Component, Mixture};

fn main() {
    let components = vec![
        Component {
            name: "light".into(),
            z: 0.5,
            a: 2.5,
            b: 900.0,
            c: 100.0,
            t_min: 250.0,
            t_max: 450.0,
            cp_l: 100.0,
            cp_v: 80.0,
            latent_heat: 30_000.0,
        },
        Component {
            name: "heavy".into(),
            z: 0.5,
            a: 2.5,
            b: 1000.0,
            c: 100.0,
            t_min: 250.0,
            t_max: 450.0,
            cp_l: 100.0,
            cp_v: 80.0,
            latent_heat: 30_000.0,
        },
    ];
    let mixture = Mixture::new(components).expect("valid mixture");

    let tp = tp_flash(&mixture, 278.0, 1.0).expect("TP flash");
    println!(
        "TP: phase={:?}, beta={:.10}, h={:.6} J/mol",
        tp.phase, tp.beta, tp.enthalpy
    );

    let ph = ph_flash(&mixture, 1.0, 10_000.0, 250.0, 300.0).expect("PH flash");
    println!(
        "PH: T={:.10} K, phase={:?}, beta={:.10}",
        ph.temperature, ph.phase, ph.beta
    );

    let request = r#"{
      "kind": "tp",
      "temperature": 278.0,
      "pressure": 1.0,
      "components": [
        {"name":"light","z":0.5,"a":2.5,"b":900.0,"c":100.0,"t_min":250.0,"t_max":450.0,"cp_l":100.0,"cp_v":80.0,"L":30000.0},
        {"name":"heavy","z":0.5,"a":2.5,"b":1000.0,"c":100.0,"t_min":250.0,"t_max":450.0,"cp_l":100.0,"cp_v":80.0,"L":30000.0}
      ]
    }"#;
    println!("\nJSON result:\n{}", run_json(request));
}

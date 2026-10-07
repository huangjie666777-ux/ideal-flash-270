use crate::error::{FlashError, FlashResult};
use crate::properties::Mixture;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    Liquid,
    Vapor,
    TwoPhase,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FlashState {
    pub temperature: f64,
    pub pressure: f64,
    pub phase: Phase,
    pub beta: f64,
    pub liquid_composition: Vec<f64>,
    pub vapor_composition: Vec<f64>,
    pub enthalpy: f64,
}

fn rachford_rice(z: &[f64], k: &[f64], beta: f64) -> f64 {
    z.iter()
        .zip(k)
        .map(|(&zi, &ki)| zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)))
        .sum()
}

fn mixture_enthalpy(
    mixture: &Mixture,
    phase: Phase,
    beta: f64,
    x: &[f64],
    y: &[f64],
    temperature: f64,
) -> f64 {
    let liquid: f64 = mixture
        .components
        .iter()
        .zip(x)
        .map(|(component, &xi)| xi * component.liquid_enthalpy(temperature))
        .sum();
    let vapor: f64 = mixture
        .components
        .iter()
        .zip(y)
        .map(|(component, &yi)| yi * component.vapor_enthalpy(temperature))
        .sum();
    match phase {
        Phase::Liquid => liquid,
        Phase::Vapor => vapor,
        Phase::TwoPhase => (1.0 - beta) * liquid + beta * vapor,
    }
}

pub fn tp_flash(mixture: &Mixture, temperature: f64, pressure: f64) -> FlashResult<FlashState> {
    mixture.validate_temperature(temperature)?;
    if !pressure.is_finite() || pressure <= 0.0 {
        return Err(FlashError::InvalidInput(
            "pressure must be finite and positive".into(),
        ));
    }

    let z: Vec<f64> = mixture.components.iter().map(|c| c.z).collect();
    let k = mixture.k_values(temperature, pressure)?;
    if k.iter().all(|ki| (ki - 1.0).abs() <= 1.0e-14) {
        return Err(FlashError::Degenerate(
            "all K-values are 1; phase split is underdetermined".into(),
        ));
    }

    let r0 = rachford_rice(&z, &k, 0.0);
    let r1 = rachford_rice(&z, &k, 1.0);

    let (phase, beta, x, y) = if r0 <= 0.0 && r1 < 0.0 {
        (Phase::Liquid, 0.0, z.clone(), Vec::new())
    } else {
        // Vapor: g(0)>0 and g(1)>=0; otherwise g(0)<=0, g(1)>=0 is two-phase.
        if r0 > 0.0 && r1 >= 0.0 {
            (Phase::Vapor, 1.0, Vec::new(), z.clone())
        } else {
            let mut low = 0.0_f64;
            let mut high = 1.0_f64;
            let mut f_low = r0;
            let mut f_high = r1;
            let mut root = None;

            for _ in 0..200 {
                let mid = (low + high) * 0.5;
                let f_mid = rachford_rice(&z, &k, mid);
                if !f_mid.is_finite() {
                    return Err(FlashError::NotConverged(
                        "non-finite Rachford-Rice residual".into(),
                    ));
                }
                if f_mid == 0.0 || high - low <= 1.0e-15 {
                    root = Some(mid);
                    break;
                }
                if f_mid * f_low < 0.0 {
                    high = mid;
                    f_high = f_mid;
                } else if f_mid * f_high < 0.0 {
                    low = mid;
                    f_low = f_mid;
                } else {
                    break;
                }
            }

            let beta = root.ok_or_else(|| {
                FlashError::NotConverged(
                    "Rachford-Rice root did not remain bracketed in [0, 1]".into(),
                )
            })?;
            if beta <= 0.0 || beta >= 1.0 {
                return Err(FlashError::NotConverged(
                    "two-phase root collapsed to a phase boundary".into(),
                ));
            }

            let denominators: Vec<f64> = k.iter().map(|ki| 1.0 + beta * (ki - 1.0)).collect();
            if denominators
                .iter()
                .any(|den| !den.is_finite() || *den <= 0.0)
            {
                return Err(FlashError::NotConverged(
                    "invalid Rachford-Rice composition denominator".into(),
                ));
            }
            let x: Vec<f64> = z
                .iter()
                .zip(&denominators)
                .map(|(&zi, den)| zi / den)
                .collect();
            let y: Vec<f64> = x.iter().zip(&k).map(|(&xi, ki)| ki * xi).collect();
            (Phase::TwoPhase, beta, x, y)
        }
    };

    let enthalpy = mixture_enthalpy(mixture, phase, beta, &x, &y, temperature);
    if !enthalpy.is_finite() {
        return Err(FlashError::NotConverged(
            "non-finite mixture enthalpy".into(),
        ));
    }

    Ok(FlashState {
        temperature,
        pressure,
        phase,
        beta,
        liquid_composition: x,
        vapor_composition: y,
        enthalpy,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::Component;

    fn component(name: &str, z: f64, b: f64) -> Component {
        Component {
            name: name.into(),
            z,
            a: 2.5,
            b,
            c: 100.0,
            t_min: 250.0,
            t_max: 450.0,
            cp_l: 100.0,
            cp_v: 80.0,
            latent_heat: 30_000.0,
        }
    }

    #[test]
    fn detects_single_and_two_phase_states() {
        let mixture = Mixture::new(vec![
            component("light", 0.5, 900.0),
            component("heavy", 0.5, 1000.0),
        ])
        .unwrap();
        assert_eq!(
            tp_flash(&mixture, 278.0, 1.0).unwrap().phase,
            Phase::TwoPhase
        );
        assert_eq!(tp_flash(&mixture, 260.0, 1.0).unwrap().phase, Phase::Liquid);
        assert_eq!(tp_flash(&mixture, 440.0, 1.0).unwrap().phase, Phase::Vapor);
    }

    #[test]
    fn two_phase_conserves_components() {
        let mixture = Mixture::new(vec![
            component("light", 0.5, 900.0),
            component("heavy", 0.5, 1000.0),
        ])
        .unwrap();
        let state = tp_flash(&mixture, 278.0, 1.0).unwrap();
        assert!((state.liquid_composition.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        assert!((state.vapor_composition.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        for i in 0..2 {
            let rebuilt = (1.0 - state.beta) * state.liquid_composition[i]
                + state.beta * state.vapor_composition[i];
            assert!((rebuilt - 0.5).abs() < 1.0e-12);
        }
    }
}

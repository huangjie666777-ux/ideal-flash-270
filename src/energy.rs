use crate::equilibrium::{tp_flash, FlashState};
use crate::error::{FlashError, FlashResult};
use crate::properties::Mixture;

pub fn ph_flash(
    mixture: &Mixture,
    pressure: f64,
    target_enthalpy: f64,
    t_lo: f64,
    t_hi: f64,
) -> FlashResult<FlashState> {
    if !target_enthalpy.is_finite() {
        return Err(FlashError::InvalidInput(
            "target enthalpy must be finite".into(),
        ));
    }
    if !pressure.is_finite() || pressure <= 0.0 {
        return Err(FlashError::InvalidInput(
            "pressure must be finite and positive".into(),
        ));
    }
    if !t_lo.is_finite() || !t_hi.is_finite() || !(t_lo < t_hi) {
        return Err(FlashError::InvalidInput(
            "temperature bracket must satisfy t_lo < t_hi".into(),
        ));
    }
    mixture.validate_temperature(t_lo)?;
    mixture.validate_temperature(t_hi)?;

    let mut low = t_lo;
    let mut high = t_hi;
    let mut f_low = tp_flash(mixture, low, pressure)?.enthalpy - target_enthalpy;
    let mut f_high = tp_flash(mixture, high, pressure)?.enthalpy - target_enthalpy;

    if f_low == 0.0 || f_high == 0.0 {
        return Err(FlashError::NoBracket(
            "target enthalpy is reached at a bracket endpoint; endpoints are not accepted".into(),
        ));
    }
    if f_low * f_high > 0.0 {
        return Err(FlashError::NoBracket(
            "target enthalpy is not enclosed by the supplied temperature bracket".into(),
        ));
    }

    let mut midpoint;
    for _ in 0..100 {
        midpoint = (low + high) * 0.5;
        if midpoint <= low || midpoint >= high {
            break;
        }
        let f_mid = tp_flash(mixture, midpoint, pressure)?.enthalpy - target_enthalpy;
        if f_mid == 0.0 {
            return tp_flash(mixture, midpoint, pressure);
        }
        if f_mid * f_low < 0.0 {
            high = midpoint;
            f_high = f_mid;
        } else if f_mid * f_high < 0.0 {
            low = midpoint;
            f_low = f_mid;
        } else {
            return Err(FlashError::NotConverged(
                "enthalpy residual lost its sign bracket during iteration".into(),
            ));
        }
        if high - low <= 1.0e-10 {
            midpoint = (low + high) * 0.5;
            return tp_flash(mixture, midpoint, pressure);
        }
    }

    Err(FlashError::NotConverged(
        "PH iteration reached its limit without converging".into(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::equilibrium::Phase;
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
    fn solves_enclosed_enthalpy_and_rejects_endpoint() {
        let mixture = Mixture::new(vec![
            component("light", 0.5, 900.0),
            component("heavy", 0.5, 1000.0),
        ])
        .unwrap();
        let result = ph_flash(&mixture, 1.0, 10_000.0, 250.0, 300.0).unwrap();
        assert!((result.enthalpy - 10_000.0).abs() < 1.0e-5);
        assert!(matches!(result.phase, Phase::TwoPhase | Phase::Vapor));

        let endpoint_h = tp_flash(&mixture, 250.0, 1.0).unwrap().enthalpy;
        assert!(ph_flash(&mixture, 1.0, endpoint_h, 250.0, 300.0).is_err());
        assert!(ph_flash(&mixture, 1.0, 1.0e12, 250.0, 300.0).is_err());
    }
}

//! Isothermal (TP) flash via Rachford-Rice.

use crate::props::FlashSystem;
use crate::{FlashError, FlashResult, Phase};

const RR_TOL: f64 = 1e-12;
const RR_MAX_ITER: usize = 200;

/// Rachford-Rice residual: sum_i z_i (K_i - 1) / (1 + beta (K_i - 1)).
fn rr_residual(z: &[f64], k: &[f64], beta: f64) -> f64 {
    z.iter()
        .zip(k)
        .map(|(zi, ki)| zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0)))
        .sum()
}

/// Solve for vapor fraction beta in [0, 1] by bisection.
/// Assumes f(0) > 0 and f(1) < 0 (two-phase).
fn solve_beta(z: &[f64], k: &[f64]) -> Result<f64, FlashError> {
    let mut lo = 0.0;
    let mut hi = 1.0;
    for _ in 0..RR_MAX_ITER {
        let mid = 0.5 * (lo + hi);
        if rr_residual(z, k, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < RR_TOL {
            return Ok(0.5 * (lo + hi));
        }
    }
    Err(FlashError::NotConverged(
        "Rachford-Rice beta solve did not converge".into(),
    ))
}

/// Mixture molar enthalpy [J/mol] given phase compositions and beta.
pub fn mixture_enthalpy(sys: &FlashSystem, t: f64, beta: f64, x: &[f64], y: &[f64]) -> f64 {
    let mut h = 0.0;
    for i in 0..sys.components.len() {
        h += (1.0 - beta) * x[i] * sys.h_l(i, t) + beta * y[i] * sys.h_v(i, t);
    }
    h
}

/// TP flash: equilibrium at given temperature [K] and pressure [bar].
pub fn tp_flash(sys: &FlashSystem, t: f64, p: f64) -> Result<FlashResult, FlashError> {
    if !p.is_finite() || p <= 0.0 {
        return Err(FlashError::InvalidInput(format!(
            "pressure must be positive and finite, got {p}"
        )));
    }
    sys.check_temperature(t)?;

    let n = sys.components.len();
    let z: Vec<f64> = sys.components.iter().map(|c| c.z).collect();
    let k: Vec<f64> = (0..n).map(|i| sys.psat(i, t).map(|ps| ps / p)).collect::<Result<_, _>>()?;

    if k.iter().all(|&ki| (ki - 1.0).abs() <= 1e-14) {
        return Err(FlashError::DegenerateEquilibrium(
            "all K-factors equal 1; phase split is indeterminate".into(),
        ));
    }

    // Phase tests: f(0) = sum z*K - 1, f(1) = 1 - sum z/K (with opposite sign).
    let f0 = rr_residual(&z, &k, 0.0); // sum z_i (K_i - 1)
    let f1 = rr_residual(&z, &k, 1.0); // 1 - sum z_i / K_i ... evaluated directly

    let (phase, beta, x, y) = if f0 <= 0.0 {
        // Undersaturated bubble point: single liquid phase.
        (Phase::Liquid, 0.0, z.clone(), Vec::new())
    } else if f1 >= 0.0 {
        // Oversaturated dew point: single vapor phase.
        (Phase::Vapor, 1.0, Vec::new(), z.clone())
    } else {
        let beta = solve_beta(&z, &k)?;
        let x: Vec<f64> = (0..n)
            .map(|i| z[i] / (1.0 + beta * (k[i] - 1.0)))
            .collect();
        let y: Vec<f64> = (0..n).map(|i| k[i] * x[i]).collect();
        (Phase::TwoPhase, beta, x, y)
    };

    let enthalpy = match phase {
        Phase::Liquid => (0..n).map(|i| z[i] * sys.h_l(i, t)).sum(),
        Phase::Vapor => (0..n).map(|i| z[i] * sys.h_v(i, t)).sum(),
        Phase::TwoPhase => mixture_enthalpy(sys, t, beta, &x, &y),
    };

    Ok(FlashResult {
        temperature: t,
        pressure: p,
        phase,
        beta,
        x,
        y,
        enthalpy,
    })
}

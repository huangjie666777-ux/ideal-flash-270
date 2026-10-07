//! PH flash: find equilibrium temperature from pressure and target enthalpy.

use crate::equilibrium::tp_flash;
use crate::props::FlashSystem;
use crate::{FlashError, FlashResult};

const T_TOL: f64 = 1e-9;
const MAX_ITER: usize = 200;

/// PH flash: equilibrium at pressure p [bar] with mixture molar enthalpy
/// h_target [J/mol]. The temperature is bracketed by [t_lo, t_hi] (K),
/// which is intersected with the common validity range. The solution may
/// cross single-phase and two-phase regions. Endpoints are never accepted
/// as the answer.
pub fn ph_flash(
    sys: &FlashSystem,
    p: f64,
    h_target: f64,
    t_lo: f64,
    t_hi: f64,
) -> Result<FlashResult, FlashError> {
    if !p.is_finite() || p <= 0.0 {
        return Err(FlashError::InvalidInput(format!(
            "pressure must be positive and finite, got {p}"
        )));
    }
    if !h_target.is_finite() {
        return Err(FlashError::InvalidInput(
            "target enthalpy is not finite".into(),
        ));
    }
    if !t_lo.is_finite() || !t_hi.is_finite() || t_lo >= t_hi {
        return Err(FlashError::InvalidInput(format!(
            "invalid temperature bracket [{t_lo}, {t_hi}]"
        )));
    }
    let lo = t_lo.max(sys.t_common_min);
    let hi = t_hi.min(sys.t_common_max);
    if lo >= hi {
        return Err(FlashError::OutOfRange(
            "temperature bracket does not overlap the common validity range".into(),
        ));
    }

    let residual = |t: f64| -> Result<f64, FlashError> {
        Ok(tp_flash(sys, t, p)?.enthalpy - h_target)
    };

    let f_lo = residual(lo)?;
    let f_hi = residual(hi)?;
    // Strict bracketing: endpoints are not accepted as the answer.
    if f_lo >= 0.0 || f_hi <= 0.0 {
        return Err(FlashError::NotBracketed(format!(
            "h(T_lo) - h_target = {f_lo}, h(T_hi) - h_target = {f_hi}; \
             target enthalpy not strictly bracketed"
        )));
    }

    // Bisection; h(T) is non-decreasing for positive Cp and L.
    let mut a = lo;
    let mut b = hi;
    let mut t_mid = 0.5 * (a + b);
    let mut converged = false;
    for _ in 0..MAX_ITER {
        t_mid = 0.5 * (a + b);
        let f_mid = residual(t_mid)?;
        if f_mid > 0.0 {
            b = t_mid;
        } else {
            a = t_mid;
        }
        if b - a < T_TOL {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(FlashError::NotConverged(
            "PH temperature solve did not converge".into(),
        ));
    }
    tp_flash(sys, t_mid, p)
}

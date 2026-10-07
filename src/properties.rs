use crate::error::{FlashError, FlashResult};

pub const REFERENCE_TEMPERATURE: f64 = 298.15;

#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    pub name: String,
    pub z: f64,
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub t_min: f64,
    pub t_max: f64,
    pub cp_l: f64,
    pub cp_v: f64,
    pub latent_heat: f64,
}

impl Component {
    pub fn saturation_pressure(&self, temperature: f64) -> FlashResult<f64> {
        if !temperature.is_finite() || temperature + self.c <= 0.0 {
            return Err(FlashError::OutOfRange(format!(
                "component {} cannot be evaluated at T={temperature} K",
                self.name
            )));
        }
        let log_p = self.a - self.b / (temperature + self.c);
        if !log_p.is_finite() {
            return Err(FlashError::OutOfRange(format!(
                "non-finite saturation pressure for component {} at T={temperature} K",
                self.name
            )));
        }
        Ok(10.0_f64.powf(log_p))
    }

    pub fn liquid_enthalpy(&self, temperature: f64) -> f64 {
        self.cp_l * (temperature - REFERENCE_TEMPERATURE)
    }

    pub fn vapor_enthalpy(&self, temperature: f64) -> f64 {
        self.latent_heat + self.cp_v * (temperature - REFERENCE_TEMPERATURE)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Mixture {
    pub components: Vec<Component>,
    pub t_min: f64,
    pub t_max: f64,
}

impl Mixture {
    pub fn new(components: Vec<Component>) -> FlashResult<Self> {
        if !(2..=16).contains(&components.len()) {
            return Err(FlashError::InvalidInput(format!(
                "expected 2 to 16 components, got {}",
                components.len()
            )));
        }

        for component in &components {
            let finite_fields = [
                component.z,
                component.a,
                component.b,
                component.c,
                component.t_min,
                component.t_max,
                component.cp_l,
                component.cp_v,
                component.latent_heat,
            ];
            if finite_fields.iter().any(|value| !value.is_finite()) {
                return Err(FlashError::InvalidInput(format!(
                    "component {} contains a non-finite property",
                    component.name
                )));
            }
            if component.name.trim().is_empty() {
                return Err(FlashError::InvalidInput(
                    "component names must be non-empty".into(),
                ));
            }
            if component.z <= 0.0 {
                return Err(FlashError::InvalidInput(format!(
                    "mole fraction for {} must be positive",
                    component.name
                )));
            }
            if component.t_min > component.t_max {
                return Err(FlashError::InvalidInput(format!(
                    "invalid Antoine temperature range for {}",
                    component.name
                )));
            }
            if component.cp_l <= 0.0 || component.cp_v <= 0.0 {
                return Err(FlashError::InvalidInput(format!(
                    "CpL and CpV must be positive for {}",
                    component.name
                )));
            }
        }

        let mut names: Vec<&str> = components.iter().map(|c| c.name.as_str()).collect();
        names.sort_unstable();
        if names.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(FlashError::InvalidInput(
                "component names must be unique".into(),
            ));
        }

        let z_sum: f64 = components.iter().map(|c| c.z).sum();
        if (z_sum - 1.0).abs() > 1.0e-12 {
            return Err(FlashError::InvalidInput(format!(
                "mole fractions must sum to 1, got {z_sum}"
            )));
        }

        let t_min = components
            .iter()
            .map(|c| c.t_min)
            .fold(f64::NEG_INFINITY, f64::max);
        let t_max = components
            .iter()
            .map(|c| c.t_max)
            .fold(f64::INFINITY, f64::min);
        if !(t_min < t_max) {
            return Err(FlashError::OutOfRange(
                "components have no common interior temperature range".into(),
            ));
        }

        for component in &components {
            if component.c + t_min <= 0.0 || component.c + t_max <= 0.0 {
                return Err(FlashError::OutOfRange(format!(
                    "T+C must be positive throughout the common range for {}",
                    component.name
                )));
            }
            let latent_low = component.vapor_enthalpy(t_min) - component.liquid_enthalpy(t_min);
            let latent_high = component.vapor_enthalpy(t_max) - component.liquid_enthalpy(t_max);
            if !(latent_low > 0.0 && latent_high > 0.0) {
                return Err(FlashError::OutOfRange(format!(
                    "hV-hL must be positive throughout the common range for {}",
                    component.name
                )));
            }
            component.saturation_pressure(t_min)?;
            component.saturation_pressure(t_max)?;
        }

        Ok(Self {
            components,
            t_min,
            t_max,
        })
    }

    pub fn validate_temperature(&self, temperature: f64) -> FlashResult<()> {
        if !temperature.is_finite() || temperature < self.t_min || temperature > self.t_max {
            return Err(FlashError::OutOfRange(format!(
                "T={temperature} K is outside [{}, {}] K",
                self.t_min, self.t_max
            )));
        }
        Ok(())
    }

    pub fn validate_pressure(&self, pressure: f64) -> FlashResult<()> {
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(FlashError::InvalidInput(format!(
                "pressure must be finite and positive, got {pressure} bar"
            )));
        }
        for component in &self.components {
            let psat = component.saturation_pressure(self.t_min)?;
            if !psat.is_finite() || psat <= 0.0 {
                return Err(FlashError::OutOfRange(format!(
                    "invalid saturation pressure for {}",
                    component.name
                )));
            }
        }
        Ok(())
    }

    pub fn k_values(&self, temperature: f64, pressure: f64) -> FlashResult<Vec<f64>> {
        self.validate_temperature(temperature)?;
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(FlashError::InvalidInput(
                "pressure must be finite and positive".into(),
            ));
        }
        self.components
            .iter()
            .map(|component| {
                let psat = component.saturation_pressure(temperature)?;
                Ok(psat / pressure)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    pub fn test_component(name: &str, z: f64, b: f64) -> Component {
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
    fn validates_mixture() {
        let mixture = Mixture::new(vec![
            test_component("a", 0.5, 900.0),
            test_component("b", 0.5, 1000.0),
        ])
        .unwrap();
        assert_eq!(mixture.t_min, 250.0);
        assert_eq!(mixture.t_max, 450.0);
    }

    #[test]
    fn rejects_nonfinite_and_bad_mole_fractions() {
        let mut c1 = test_component("a", 0.5, 900.0);
        let mut c2 = test_component("b", 0.5, 1000.0);
        c2.z = f64::NAN;
        assert!(Mixture::new(vec![c1.clone(), c2]).is_err());
        c2 = test_component("b", 0.6, 1000.0);
        c1.z = 0.5;
        assert!(Mixture::new(vec![c1, c2]).is_err());
    }
}

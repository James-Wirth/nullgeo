use crate::metric::{Mat4, Metric, Vec4};
use crate::metrics::kerr_schild;
use crate::spacetime::{cartesian_circular_four_velocity, CircularOrbits, Spacetime};
use crate::{Error, Result};

#[derive(Clone, Copy, Debug)]
pub struct Schwarzschild {
    m: f64,
}

impl Schwarzschild {
    pub fn new(m: f64) -> Result<Self> {
        if !(m >= 0.0 && m.is_finite()) {
            return Err(Error::InvalidArg(format!(
                "Schwarzschild requires M >= 0 and finite, got M = {m}"
            )));
        }
        Ok(Self { m })
    }

    pub fn mass(&self) -> f64 {
        self.m
    }
}

impl Spacetime for Schwarzschild {
    fn is_captured(&self, x: &Vec4) -> bool {
        self.radius(x) <= 2.0 * self.m * (1.0 + 1e-3)
    }

    fn circular_orbits(&self) -> Option<&dyn CircularOrbits> {
        Some(self)
    }
}

impl CircularOrbits for Schwarzschild {
    fn isco_radius(&self) -> f64 {
        6.0 * self.m
    }

    fn four_velocity(&self, x: &Vec4) -> Option<Vec4> {
        let r = self.radius(x);
        cartesian_circular_four_velocity(self, x, (self.m / (r * r * r)).sqrt())
    }
}

impl Metric for Schwarzschild {
    fn g(&self, x: &Vec4) -> Mat4 {
        let (r, n) = kerr_schild::radial_unit(x);
        kerr_schild::g(self.m / r, n)
    }

    fn g_inv(&self, x: &Vec4) -> Mat4 {
        let (r, n) = kerr_schild::radial_unit(x);
        kerr_schild::g_inv(self.m / r, n)
    }

    fn dg_inv(&self, x: &Vec4) -> [Mat4; 4] {
        let (r, n) = kerr_schild::radial_unit(x);
        kerr_schild::dg_inv(self.m / r, -self.m / (r * r), r, n)
    }
}

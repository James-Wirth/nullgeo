use crate::metric::{Mat4, Metric, Vec4};
use crate::metrics::kerr_schild;
use crate::spacetime::Spacetime;

#[derive(Clone, Copy, Debug)]
pub struct Schwarzschild {
    pub m: f64,
}

impl Spacetime for Schwarzschild {
    fn is_captured(&self, x: &Vec4) -> bool {
        self.radius(x) <= 2.0 * self.m * (1.0 + 1e-3)
    }

    fn isco_radius(&self) -> Option<f64> {
        Some(6.0 * self.m)
    }

    fn disk_emitter(&self, x: &Vec4) -> Option<Vec4> {
        let r = self.radius(x);
        self.circular_emitter(x, (self.m / (r * r * r)).sqrt())
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

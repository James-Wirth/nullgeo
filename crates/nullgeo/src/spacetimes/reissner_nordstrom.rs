use crate::geometry::{Chart, Mat4, Metric, Vec4};
use crate::spacetimes::{kerr_schild, Spacetime};
use crate::{Error, Result};

#[derive(Clone, Copy, Debug)]
pub struct ReissnerNordstrom {
    m: f64,
    q: f64,
}

impl ReissnerNordstrom {
    pub fn new(m: f64, q: f64) -> Result<Self> {
        if q.abs() > m {
            return Err(Error::InvalidArg(format!(
                "Reissner-Nordstrom requires |Q| <= M, got Q = {q}, M = {m}"
            )));
        }
        Ok(Self { m, q })
    }

    pub fn mass(&self) -> f64 {
        self.m
    }

    pub fn charge(&self) -> f64 {
        self.q
    }

    pub fn outer_horizon(&self) -> f64 {
        self.m + (self.m * self.m - self.q * self.q).sqrt()
    }

    fn h(&self, r: f64) -> f64 {
        (self.m - 0.5 * self.q * self.q / r) / r
    }

    fn h_prime(&self, r: f64) -> f64 {
        (self.q * self.q / r - self.m) / (r * r)
    }
}

impl Chart for ReissnerNordstrom {}

impl Spacetime for ReissnerNordstrom {
    fn is_captured(&self, x: &Vec4) -> bool {
        self.radius(x) <= self.outer_horizon() * (1.0 + 1e-3)
    }
}

impl Metric for ReissnerNordstrom {
    fn g(&self, x: &Vec4) -> Mat4 {
        let (r, n) = kerr_schild::radial_unit(x);
        kerr_schild::g(self.h(r), n)
    }

    fn g_inv(&self, x: &Vec4) -> Mat4 {
        let (r, n) = kerr_schild::radial_unit(x);
        kerr_schild::g_inv(self.h(r), n)
    }

    fn dg_inv(&self, x: &Vec4) -> [Mat4; 4] {
        let (r, n) = kerr_schild::radial_unit(x);
        kerr_schild::dg_inv(self.h(r), self.h_prime(r), r, n)
    }
}

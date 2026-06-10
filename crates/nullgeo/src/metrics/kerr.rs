use crate::metric::{Mat4, Metric, Vec4};
use crate::metrics::kerr_schild;
use crate::spacetime::Spacetime;
use crate::{Error, Result};

#[derive(Clone, Copy, Debug)]
pub struct Kerr {
    m: f64,
    a: f64,
}

impl Kerr {
    pub fn new(m: f64, a: f64) -> Result<Self> {
        if a.abs() > m {
            return Err(Error::InvalidArg(format!(
                "Kerr requires |a| <= M, got a = {a}, M = {m}"
            )));
        }
        Ok(Self { m, a })
    }

    pub fn mass(&self) -> f64 {
        self.m
    }

    pub fn spin(&self) -> f64 {
        self.a
    }

    pub fn outer_horizon(&self) -> f64 {
        self.m + (self.m * self.m - self.a * self.a).sqrt()
    }

    fn ks_radius_sq(&self, x: &Vec4) -> f64 {
        let (px, py, pz) = (x[1], x[2], x[3]);
        let rr = px * px + py * py + pz * pz;
        let half = 0.5 * (rr - self.a * self.a);
        (half + (half * half + self.a * self.a * pz * pz).sqrt()).max(0.0)
    }

    fn h_and_l(&self, x: &Vec4) -> (f64, Vec4) {
        let a = self.a;
        let (px, py, pz) = (x[1], x[2], x[3]);
        let r2 = self.ks_radius_sq(x);
        let r = r2.sqrt().max(1e-20);
        let h = self.m * r2 * r / (r2 * r2 + a * a * pz * pz).max(1e-300);
        let denom = r2 + a * a;
        let l = Vec4::new(
            1.0,
            (r * px + a * py) / denom,
            (r * py - a * px) / denom,
            pz / r,
        );
        (h, l)
    }
}

impl Spacetime for Kerr {
    fn radius(&self, x: &Vec4) -> f64 {
        self.ks_radius_sq(x).sqrt()
    }

    fn is_captured(&self, x: &Vec4) -> bool {
        self.radius(x) <= self.outer_horizon() * (1.0 + 1e-3)
    }

    fn isco_radius(&self) -> Option<f64> {
        let chi = self.a / self.m;
        let z1 = 1.0 + (1.0 - chi * chi).cbrt() * ((1.0 + chi).cbrt() + (1.0 - chi).cbrt());
        let z2 = (3.0 * chi * chi + z1 * z1).sqrt();
        let branch = ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).max(0.0).sqrt();
        Some(self.m * (3.0 + z2 - branch))
    }

    fn disk_emitter(&self, x: &Vec4) -> Option<Vec4> {
        let r = self.radius(x);
        let sqrt_m = self.m.sqrt();
        let sense = if self.a >= 0.0 { 1.0 } else { -1.0 };
        let omega = sense * sqrt_m / (r.powf(1.5) + self.a.abs() * sqrt_m);
        self.circular_emitter(x, omega)
    }
}

impl Metric for Kerr {
    fn g(&self, x: &Vec4) -> Mat4 {
        let (h, l) = self.h_and_l(x);
        kerr_schild::eta() + 2.0 * h * l * l.transpose()
    }

    fn g_inv(&self, x: &Vec4) -> Mat4 {
        let (h, mut l) = self.h_and_l(x);
        l[0] = -1.0;
        kerr_schild::eta() - 2.0 * h * l * l.transpose()
    }
}

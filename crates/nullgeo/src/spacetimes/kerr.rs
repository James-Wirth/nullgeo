use crate::geometry::{Chart, Mat4, Metric, Vec4};
use crate::spacetimes::{cartesian_circular_four_velocity, kerr_schild, CircularOrbits, Spacetime};
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

impl Chart for Kerr {
    fn radius(&self, x: &Vec4) -> f64 {
        self.ks_radius_sq(x).sqrt()
    }
}

impl Spacetime for Kerr {
    fn is_captured(&self, x: &Vec4) -> bool {
        self.radius(x) <= self.outer_horizon() * (1.0 + 1e-3)
    }

    fn circular_orbits(&self) -> Option<&dyn CircularOrbits> {
        Some(self)
    }
}

impl CircularOrbits for Kerr {
    fn isco_radius(&self) -> f64 {
        let chi = self.a / self.m;
        let z1 = 1.0 + (1.0 - chi * chi).cbrt() * ((1.0 + chi).cbrt() + (1.0 - chi).cbrt());
        let z2 = (3.0 * chi * chi + z1 * z1).sqrt();
        let branch = ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).max(0.0).sqrt();
        self.m * (3.0 + z2 - branch)
    }

    fn four_velocity(&self, x: &Vec4) -> Option<Vec4> {
        let r = self.radius(x);
        let sqrt_m = self.m.sqrt();
        let sense = if self.a >= 0.0 { 1.0 } else { -1.0 };
        let omega = sense * sqrt_m / (r.powf(1.5) + self.a.abs() * sqrt_m);
        cartesian_circular_four_velocity(self, x, omega)
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

    fn dg_inv(&self, x: &Vec4) -> [Mat4; 4] {
        let a = self.a;
        let (px, py, pz) = (x[1], x[2], x[3]);
        let r2 = self.ks_radius_sq(x);
        let r = r2.sqrt().max(1e-20);
        let sigma = (r2 * r2 + a * a * pz * pz).max(1e-300);
        let h = self.m * r2 * r / sigma;
        let q = r2 + a * a;
        let l = Vec4::new(-1.0, (r * px + a * py) / q, (r * py - a * px) / q, pz / r);
        let ll = l * l.transpose();

        let dr = [r2 * r * px / sigma, r2 * r * py / sigma, r * pz * q / sigma];

        let mut out = [Mat4::zeros(); 4];
        for i in 0..3 {
            let delta = |axis: usize| ((i == axis) as i32) as f64;
            let dsigma = 4.0 * r2 * r * dr[i] + 2.0 * a * a * pz * delta(2);
            let dh = self.m * (3.0 * r2 * dr[i] - r2 * r * dsigma / sigma) / sigma;
            let dl = Vec4::new(
                0.0,
                (dr[i] * px + r * delta(0) + a * delta(1)) / q - l[1] * 2.0 * r * dr[i] / q,
                (dr[i] * py + r * delta(1) - a * delta(0)) / q - l[2] * 2.0 * r * dr[i] / q,
                delta(2) / r - pz * dr[i] / r2,
            );
            let dll = dl * l.transpose() + l * dl.transpose();
            out[i + 1] = -2.0 * (dh * ll + h * dll);
        }
        out
    }
}

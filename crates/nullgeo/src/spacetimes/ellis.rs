use std::f64::consts::FRAC_PI_2;

use crate::geometry::{Chart, Mat3, Mat4, Metric, PhasePoint, RayAlignment, Vec4};
use crate::spacetimes::{SkySide, Spacetime};
use crate::{Error, Result};
use nalgebra::Vector3;

#[derive(Clone, Copy, Debug)]
pub struct Ellis {
    b0: f64,
}

struct SphereFrame {
    radial: Vector3<f64>,
    polar: Vector3<f64>,
    azimuthal: Vector3<f64>,
}

fn sphere_frame(theta: f64, phi: f64) -> SphereFrame {
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi.sin_cos();
    SphereFrame {
        radial: Vector3::new(st * cp, st * sp, ct),
        polar: Vector3::new(ct * cp, ct * sp, -st),
        azimuthal: Vector3::new(-sp, cp, 0.0),
    }
}

impl Ellis {
    pub fn new(b0: f64) -> Result<Self> {
        if !(b0 > 0.0 && b0.is_finite()) {
            return Err(Error::InvalidArg(format!(
                "Ellis throat radius must be positive and finite, got b0 = {b0}"
            )));
        }
        Ok(Self { b0 })
    }

    pub fn throat_radius(&self) -> f64 {
        self.b0
    }

    fn rho_sq(&self, l: f64) -> f64 {
        self.b0 * self.b0 + l * l
    }
}

impl Metric for Ellis {
    fn g(&self, x: &Vec4) -> Mat4 {
        let rho2 = self.rho_sq(x[1]);
        let s2 = x[2].sin().powi(2).max(1e-24);
        Mat4::from_diagonal(&[-1.0, 1.0, rho2, rho2 * s2].into())
    }

    fn g_inv(&self, x: &Vec4) -> Mat4 {
        let rho2 = self.rho_sq(x[1]);
        let s2 = x[2].sin().powi(2).max(1e-24);
        Mat4::from_diagonal(&[-1.0, 1.0, 1.0 / rho2, 1.0 / (rho2 * s2)].into())
    }

    fn dg_inv(&self, x: &Vec4) -> [Mat4; 4] {
        let l = x[1];
        let rho2 = self.rho_sq(l);
        let (st, ct) = x[2].sin_cos();
        let s2 = (st * st).max(1e-24);

        let mut out = [Mat4::zeros(); 4];
        out[1][(2, 2)] = -2.0 * l / (rho2 * rho2);
        out[1][(3, 3)] = -2.0 * l / (rho2 * rho2 * s2);
        out[2][(3, 3)] = -2.0 * ct / (rho2 * s2 * st);
        out
    }
}

impl Chart for Ellis {
    fn embed(&self, x: &Vec4) -> [f64; 3] {
        let rho = self.rho_sq(x[1]).sqrt();
        let p = rho * sphere_frame(x[2], x[3]).radial;
        [p[0], p[1], p[2]]
    }

    fn embed_direction(&self, x: &Vec4, v: &Vec4) -> [f64; 3] {
        let l = x[1];
        let rho = self.rho_sq(l).sqrt();
        let st = x[2].sin();
        let frame = sphere_frame(x[2], x[3]);
        let d = (l / rho) * v[1] * frame.radial
            + rho * v[2] * frame.polar
            + rho * st * v[3] * frame.azimuthal;
        let len = d.dot(&d).sqrt().max(1e-300);
        [d[0] / len, d[1] / len, d[2] / len]
    }

    fn lift_direction(&self, x: &Vec4, d: [f64; 3]) -> Vec4 {
        let rho = self.rho_sq(x[1]).sqrt();
        let st = x[2].sin().max(1e-12);
        let frame = sphere_frame(x[2], x[3]);
        let d = Vector3::new(d[0], d[1], d[2]);
        Vec4::new(
            0.0,
            d.dot(&frame.radial),
            d.dot(&frame.polar) / rho,
            d.dot(&frame.azimuthal) / (rho * st),
        )
    }

    fn radius(&self, x: &Vec4) -> f64 {
        x[1].abs()
    }

    fn align_ray(&self, s: PhasePoint) -> (PhasePoint, RayAlignment) {
        let st = s.x[2].sin().max(1e-12);
        let frame = sphere_frame(s.x[2], s.x[3]);
        let e1 = frame.radial;

        let n = s.p[2] * frame.azimuthal - (s.p[3] / st) * frame.polar;
        let l_total = n.dot(&n).sqrt();

        let e3 = if l_total > 1e-12 * s.p.amax() {
            n / l_total
        } else {
            let pick = if e1[0].abs() < 0.9 {
                Vector3::x()
            } else {
                Vector3::y()
            };
            let w = e1.cross(&pick);
            w / w.dot(&w).sqrt()
        };
        let e2 = e3.cross(&e1);

        let rotation = Mat3::new(
            e1[0], e2[0], e3[0], e1[1], e2[1], e3[1], e1[2], e2[2], e3[2],
        );

        let aligned = PhasePoint {
            x: Vec4::new(s.x[0], s.x[1], FRAC_PI_2, 0.0),
            p: Vec4::new(s.p[0], s.p[1], 0.0, l_total),
        };
        (aligned, RayAlignment::from_rotation(rotation))
    }
}

impl Spacetime for Ellis {
    fn sky_side(&self, x: &Vec4) -> SkySide {
        if x[1] >= 0.0 {
            SkySide::Primary
        } else {
            SkySide::Secondary
        }
    }
}

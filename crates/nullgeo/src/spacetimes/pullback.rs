use crate::geometry::{fd_partials, Chart, Mat4, Metric, Vec4};
use crate::spacetimes::{CircularOrbits, SkySide, Spacetime};
use crate::{Error, Result};

pub trait Transition {
    fn map(&self, y: &Vec4) -> Vec4;

    fn jacobian(&self, y: &Vec4) -> Mat4;

    fn inv_jacobian(&self, y: &Vec4) -> Mat4;

    fn d_inv_jacobian(&self, y: &Vec4) -> [Mat4; 4] {
        fd_partials(|z| self.inv_jacobian(z), y)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Pullback<S, T> {
    base: S,
    transition: T,
}

impl<S, T: Transition> Pullback<S, T> {
    pub fn new(base: S, transition: T) -> Self {
        Self { base, transition }
    }
}

impl<S: Metric, T: Transition> Metric for Pullback<S, T> {
    fn g(&self, y: &Vec4) -> Mat4 {
        let j = self.transition.jacobian(y);
        j.transpose() * self.base.g(&self.transition.map(y)) * j
    }

    fn g_inv(&self, y: &Vec4) -> Mat4 {
        let a = self.transition.inv_jacobian(y);
        a * self.base.g_inv(&self.transition.map(y)) * a.transpose()
    }

    fn dg_inv(&self, y: &Vec4) -> [Mat4; 4] {
        let x = self.transition.map(y);
        let j = self.transition.jacobian(y);
        let a = self.transition.inv_jacobian(y);
        let da = self.transition.d_inv_jacobian(y);
        let g_inv = self.base.g_inv(&x);
        let dg_inv = self.base.dg_inv(&x);

        let mut out = [Mat4::zeros(); 4];
        for (alpha, out_alpha) in out.iter_mut().enumerate() {
            let mut pushed = Mat4::zeros();
            for (mu, dg_mu) in dg_inv.iter().enumerate() {
                pushed += j[(mu, alpha)] * dg_mu;
            }
            *out_alpha = da[alpha] * g_inv * a.transpose()
                + a * g_inv * da[alpha].transpose()
                + a * pushed * a.transpose();
        }
        out
    }
}

impl<S: Chart, T: Transition> Chart for Pullback<S, T> {
    fn embed(&self, y: &Vec4) -> [f64; 3] {
        self.base.embed(&self.transition.map(y))
    }

    fn embed_direction(&self, y: &Vec4, v: &Vec4) -> [f64; 3] {
        let pushed = self.transition.jacobian(y) * v;
        self.base.embed_direction(&self.transition.map(y), &pushed)
    }

    fn lift_direction(&self, y: &Vec4, d: [f64; 3]) -> Vec4 {
        let lifted = self.base.lift_direction(&self.transition.map(y), d);
        self.transition.inv_jacobian(y) * lifted
    }

    fn radius(&self, y: &Vec4) -> f64 {
        self.base.radius(&self.transition.map(y))
    }

    fn equator_distance(&self, y: &Vec4) -> f64 {
        self.base.equator_distance(&self.transition.map(y))
    }

    fn equator_distance_rate(&self, y: &Vec4, v: &Vec4) -> f64 {
        let pushed = self.transition.jacobian(y) * v;
        self.base
            .equator_distance_rate(&self.transition.map(y), &pushed)
    }
}

impl<S: Spacetime, T: Transition> Spacetime for Pullback<S, T> {
    fn is_captured(&self, y: &Vec4) -> bool {
        self.base.is_captured(&self.transition.map(y))
    }

    fn sky_side(&self, y: &Vec4) -> SkySide {
        self.base.sky_side(&self.transition.map(y))
    }

    fn circular_orbits(&self) -> Option<&dyn CircularOrbits> {
        self.base
            .circular_orbits()
            .map(|_| self as &dyn CircularOrbits)
    }
}

impl<S: Spacetime, T: Transition> CircularOrbits for Pullback<S, T> {
    fn isco_radius(&self) -> f64 {
        self.base
            .circular_orbits()
            .expect("exposed only when the base spacetime has circular orbits")
            .isco_radius()
    }

    fn four_velocity(&self, y: &Vec4) -> Option<Vec4> {
        let u = self
            .base
            .circular_orbits()?
            .four_velocity(&self.transition.map(y))?;
        Some(self.transition.inv_jacobian(y) * u)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Rotation {
    forward: Mat4,
    inverse: Mat4,
}

impl Rotation {
    pub fn about_axis(axis: [f64; 3], angle: f64) -> Result<Self> {
        let len_sq = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
        if !(len_sq > 0.0 && len_sq.is_finite() && angle.is_finite()) {
            return Err(Error::InvalidArg(format!(
                "rotation needs a finite nonzero axis and a finite angle, got axis = {axis:?}, angle = {angle}"
            )));
        }
        let len = len_sq.sqrt();
        let k = [axis[0] / len, axis[1] / len, axis[2] / len];
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;

        let r = [
            [
                c + t * k[0] * k[0],
                t * k[0] * k[1] - s * k[2],
                t * k[0] * k[2] + s * k[1],
            ],
            [
                t * k[1] * k[0] + s * k[2],
                c + t * k[1] * k[1],
                t * k[1] * k[2] - s * k[0],
            ],
            [
                t * k[2] * k[0] - s * k[1],
                t * k[2] * k[1] + s * k[0],
                c + t * k[2] * k[2],
            ],
        ];

        let mut forward = Mat4::identity();
        for (i, row) in r.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                forward[(i + 1, j + 1)] = *value;
            }
        }
        Ok(Self {
            forward,
            inverse: forward.transpose(),
        })
    }
}

impl Transition for Rotation {
    fn map(&self, y: &Vec4) -> Vec4 {
        self.forward * y
    }

    fn jacobian(&self, _y: &Vec4) -> Mat4 {
        self.forward
    }

    fn inv_jacobian(&self, _y: &Vec4) -> Mat4 {
        self.inverse
    }

    fn d_inv_jacobian(&self, _y: &Vec4) -> [Mat4; 4] {
        [Mat4::zeros(); 4]
    }
}

use crate::frame::metric_dot;
use crate::metric::{Metric, State4, Vec4};
use nalgebra::{Matrix3, Vector3};

pub type Mat3 = Matrix3<f64>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkySide {
    Primary,
    Secondary,
}

#[derive(Debug, Clone, Copy)]
pub struct RayAlignment(Mat3);

impl RayAlignment {
    pub fn identity() -> Self {
        Self(Mat3::identity())
    }

    pub fn from_rotation(rotation: Mat3) -> Self {
        Self(rotation)
    }

    pub fn apply(&self, d: [f64; 3]) -> [f64; 3] {
        let v = self.0 * Vector3::new(d[0], d[1], d[2]);
        [v[0], v[1], v[2]]
    }
}

pub trait Spacetime: Metric {
    fn cartesian_position(&self, x: &Vec4) -> [f64; 3] {
        [x[1], x[2], x[3]]
    }

    fn radius(&self, x: &Vec4) -> f64 {
        let [px, py, pz] = self.cartesian_position(x);
        (px * px + py * py + pz * pz).sqrt()
    }

    fn equator_distance(&self, x: &Vec4) -> f64 {
        self.cartesian_position(x)[2]
    }

    fn chart_direction(&self, _x: &Vec4, d: [f64; 3]) -> Vec4 {
        Vec4::new(0.0, d[0], d[1], d[2])
    }

    fn cartesian_direction(&self, s: &State4) -> [f64; 3] {
        let v = self.g_inv(&s.x) * s.p;
        let len = (v[1] * v[1] + v[2] * v[2] + v[3] * v[3]).sqrt().max(1e-300);
        [v[1] / len, v[2] / len, v[3] / len]
    }

    fn is_captured(&self, _x: &Vec4) -> bool {
        false
    }

    fn sky_side(&self, _x: &Vec4) -> SkySide {
        SkySide::Primary
    }

    fn align_ray(&self, s: State4) -> (State4, RayAlignment) {
        (s, RayAlignment::identity())
    }

    fn circular_orbits(&self) -> Option<&dyn CircularOrbits> {
        None
    }
}

pub trait CircularOrbits {
    fn isco_radius(&self) -> f64;

    fn four_velocity(&self, x: &Vec4) -> Option<Vec4>;
}

pub fn cartesian_circular_four_velocity<M: Metric + ?Sized>(
    metric: &M,
    x: &Vec4,
    omega: f64,
) -> Option<Vec4> {
    let u = Vec4::new(1.0, -omega * x[2], omega * x[1], 0.0);
    let norm_sq = metric_dot(&metric.g(x), &u, &u);
    (norm_sq < 0.0).then(|| u / (-norm_sq).sqrt())
}

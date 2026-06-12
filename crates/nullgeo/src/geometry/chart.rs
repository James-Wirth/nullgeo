use super::{Mat3, PhasePoint, Vec4};
use nalgebra::Vector3;

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

pub trait Chart {
    fn embed(&self, x: &Vec4) -> [f64; 3] {
        [x[1], x[2], x[3]]
    }

    fn embed_direction(&self, _x: &Vec4, v: &Vec4) -> [f64; 3] {
        let len = (v[1] * v[1] + v[2] * v[2] + v[3] * v[3]).sqrt().max(1e-300);
        [v[1] / len, v[2] / len, v[3] / len]
    }

    fn lift_direction(&self, _x: &Vec4, d: [f64; 3]) -> Vec4 {
        Vec4::new(0.0, d[0], d[1], d[2])
    }

    fn radius(&self, x: &Vec4) -> f64 {
        let [px, py, pz] = self.embed(x);
        (px * px + py * py + pz * pz).sqrt()
    }

    fn equator_distance(&self, x: &Vec4) -> f64 {
        self.embed(x)[2]
    }

    fn equator_distance_rate(&self, _x: &Vec4, v: &Vec4) -> f64 {
        v[3]
    }

    fn align_ray(&self, s: PhasePoint) -> (PhasePoint, RayAlignment) {
        (s, RayAlignment::identity())
    }
}

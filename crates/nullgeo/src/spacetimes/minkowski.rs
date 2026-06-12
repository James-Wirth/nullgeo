use crate::geometry::{Chart, Mat4, Metric, Vec4};
use crate::spacetimes::Spacetime;

#[derive(Clone, Copy, Debug)]
pub struct Minkowski;

impl Chart for Minkowski {}

impl Spacetime for Minkowski {}

impl Metric for Minkowski {
    fn g(&self, _x: &Vec4) -> Mat4 {
        Mat4::from_diagonal(&[-1.0, 1.0, 1.0, 1.0].into())
    }
    fn g_inv(&self, x: &Vec4) -> Mat4 {
        self.g(x)
    }
    fn dg_inv(&self, _x: &Vec4) -> [Mat4; 4] {
        [Mat4::zeros(); 4]
    }
}

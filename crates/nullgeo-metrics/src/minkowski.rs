use nullgeo_core::metric::{Mat4, Metric, Vec4};

pub struct Minkowski;

impl Metric for Minkowski {
    fn g(&self, _x: &Vec4) -> Mat4 {
        Mat4::from_diagonal(&[-1.0, 1.0, 1.0, 1.0].into())
    }
    fn g_inv(&self, x: &Vec4) -> Mat4 { self.g(x) }
    fn dg_inv(&self, _x: &Vec4) -> Option<[Mat4; 4]> { None } 
}

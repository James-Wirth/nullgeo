use nalgebra::{Matrix4, Vector4};

pub type Vec4 = Vector4<f64>;
pub type Mat4 = Matrix4<f64>;

#[derive(Clone, Copy, Debug)]
pub struct State4 {
    pub x: Vec4,
    pub p: Vec4,
}

pub trait Metric {
    fn g(&self, x: &Vec4) -> Mat4;
    fn g_inv(&self, x: &Vec4) -> Mat4;
    fn gamma(&self, _x: &Vec4) -> Option<[[[f64; 4]; 4]; 4]> { None }
    fn dg_inv(&self, _x: &Vec4) -> Option<[Mat4; 4]> { None }
}

pub mod chart;
pub mod frame;
pub mod metric;

pub use chart::{Chart, RayAlignment};
pub use frame::{build_coframe, build_coframe_for, build_coframe_seeded, make_null_covector};
pub use metric::Metric;

use nalgebra::{Matrix3, Matrix4, Vector4};

pub type Vec4 = Vector4<f64>;
pub type Mat4 = Matrix4<f64>;
pub type Mat3 = Matrix3<f64>;

#[derive(Clone, Copy, Debug)]
pub struct PhasePoint {
    pub x: Vec4,
    pub p: Vec4,
}

pub fn inner(g: &Mat4, a: &Vec4, b: &Vec4) -> f64 {
    a.dot(&(g * b))
}

pub fn raise(g_inv: &Mat4, p: &Vec4) -> Vec4 {
    g_inv * p
}

pub fn lower(g: &Mat4, v: &Vec4) -> Vec4 {
    g * v
}

pub fn fd_partials<F: Fn(&Vec4) -> Mat4>(f: F, x: &Vec4) -> [Mat4; 4] {
    const H0: f64 = 7e-4;
    let mut out = [Mat4::zeros(); 4];
    for mu in 0..4 {
        let h = H0 * (1.0 + x[mu].abs());
        let at = |s: f64| {
            let mut y = *x;
            y[mu] += s * h;
            f(&y)
        };
        out[mu] = (at(-2.0) - at(2.0) + 8.0 * (at(1.0) - at(-1.0))) / (12.0 * h);
    }
    out
}

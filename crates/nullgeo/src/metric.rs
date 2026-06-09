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

    fn dg_inv(&self, x: &Vec4) -> [Mat4; 4] {
        fd_dg_inv(|y| self.g_inv(y), x)
    }
}

pub fn fd_dg_inv<F: Fn(&Vec4) -> Mat4>(f: F, x: &Vec4) -> [Mat4; 4] {
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

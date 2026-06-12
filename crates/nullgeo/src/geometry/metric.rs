use super::{fd_partials, Mat4, Vec4};

pub trait Metric {
    fn g(&self, x: &Vec4) -> Mat4;

    fn g_inv(&self, x: &Vec4) -> Mat4;

    fn dg_inv(&self, x: &Vec4) -> [Mat4; 4] {
        fd_partials(|y| self.g_inv(y), x)
    }
}

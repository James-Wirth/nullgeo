use crate::metric::{Mat4, Vec4};

pub fn eta() -> Mat4 {
    Mat4::from_diagonal(&[-1.0, 1.0, 1.0, 1.0].into())
}

pub fn radial_unit(x: &Vec4) -> (f64, [f64; 3]) {
    let (px, py, pz) = (x[1], x[2], x[3]);
    let r = (px * px + py * py + pz * pz).sqrt().max(1e-20);
    (r, [px / r, py / r, pz / r])
}

fn l_dn(n: [f64; 3]) -> Vec4 {
    Vec4::new(1.0, n[0], n[1], n[2])
}

fn l_up(n: [f64; 3]) -> Vec4 {
    Vec4::new(-1.0, n[0], n[1], n[2])
}

pub fn g(h: f64, n: [f64; 3]) -> Mat4 {
    let l = l_dn(n);
    eta() + 2.0 * h * l * l.transpose()
}

pub fn g_inv(h: f64, n: [f64; 3]) -> Mat4 {
    let l = l_up(n);
    eta() - 2.0 * h * l * l.transpose()
}

pub fn dg_inv(h: f64, h_prime: f64, r: f64, n: [f64; 3]) -> [Mat4; 4] {
    let l = l_up(n);
    let ll = l * l.transpose();
    let mut out = [Mat4::zeros(); 4];
    for i in 0..3 {
        let mut dl = Vec4::zeros();
        for a in 0..3 {
            dl[a + 1] = ((i == a) as i32 as f64 - n[i] * n[a]) / r;
        }
        let dll = dl * l.transpose() + l * dl.transpose();
        out[i + 1] = -2.0 * (h_prime * n[i] * ll + h * dll);
    }
    out
}

use super::metric::{Metric, State4, Vec4, Mat4};

#[derive(Debug, Clone, Copy)]
pub struct Tolerances { pub rtol: f64, pub atol: f64 }
impl Default for Tolerances { fn default() -> Self { Self { rtol: 1e-9, atol: 1e-9 } } }

/// dx^mu/ds, dp_mu/ds.
pub fn rhs_hamiltonian<M: Metric>(m: &M, s: &State4) -> (Vec4, Vec4) {
    let ginv = m.g_inv(&s.x);
    let dx = ginv * s.p;

    let mut dp = Vec4::zeros();
    if let Some(dginv) = m.dg_inv(&s.x) {
        for mu in 0..4 {
            let v = dginv[mu] * s.p;    
            let q = s.p.dot(&v);         
            dp[mu] = -0.5 * q;
        }
    } else {
        dp.fill(0.0);
    }

    (dx, dp)
}

pub fn rk4_step<M: Metric>(m: &M, s: &State4, dl: f64) -> State4 {
    let (k1x, k1p) = rhs_hamiltonian(m, s);

    let s2 = State4 { x: s.x + 0.5 * dl * k1x, p: s.p + 0.5 * dl * k1p };
    let (k2x, k2p) = rhs_hamiltonian(m, &s2);

    let s3 = State4 { x: s.x + 0.5 * dl * k2x, p: s.p + 0.5 * dl * k2p };
    let (k3x, k3p) = rhs_hamiltonian(m, &s3);

    let s4 = State4 { x: s.x + dl * k3x, p: s.p + dl * k3p };
    let (k4x, k4p) = rhs_hamiltonian(m, &s4);

    let x = s.x + (dl/6.0) * (k1x + 2.0*k2x + 2.0*k3x + k4x);
    let p = s.p + (dl/6.0) * (k1p + 2.0*k2p + 2.0*k3p + k4p);

    State4 { x, p }
}

pub fn hamiltonian<M: Metric>(m: &M, s: &State4) -> f64 {
    let g_inv = m.g_inv(&s.x);
    0.5 * quad_form(&g_inv, &s.p)
}

#[inline]
fn quad_form(m: &Mat4, p: &Vec4) -> f64 {
    let mp = m * p;
    p.dot(&mp)
}

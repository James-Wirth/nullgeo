use super::metric::{Mat4, Metric, State4, Vec4};

#[derive(Debug, Clone, Copy)]
pub struct Tolerances {
    pub rtol: f64,
    pub atol: f64,
}
impl Default for Tolerances {
    fn default() -> Self {
        Self {
            rtol: 1e-9,
            atol: 1e-9,
        }
    }
}

/// dx^mu/ds, dp_mu/ds.
pub fn rhs_hamiltonian<M: Metric + ?Sized>(m: &M, s: &State4) -> (Vec4, Vec4) {
    let ginv = m.g_inv(&s.x);
    let dx = ginv * s.p;

    let dginv = m.dg_inv(&s.x);
    let mut dp = Vec4::zeros();
    for mu in 0..4 {
        let v = dginv[mu] * s.p;
        dp[mu] = -0.5 * s.p.dot(&v);
    }

    (dx, dp)
}

pub fn rk4_step<M: Metric + ?Sized>(m: &M, s: &State4, dl: f64) -> State4 {
    let (k1x, k1p) = rhs_hamiltonian(m, s);

    let s2 = State4 {
        x: s.x + 0.5 * dl * k1x,
        p: s.p + 0.5 * dl * k1p,
    };
    let (k2x, k2p) = rhs_hamiltonian(m, &s2);

    let s3 = State4 {
        x: s.x + 0.5 * dl * k2x,
        p: s.p + 0.5 * dl * k2p,
    };
    let (k3x, k3p) = rhs_hamiltonian(m, &s3);

    let s4 = State4 {
        x: s.x + dl * k3x,
        p: s.p + dl * k3p,
    };
    let (k4x, k4p) = rhs_hamiltonian(m, &s4);

    let x = s.x + (dl / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    let p = s.p + (dl / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);

    State4 { x, p }
}

#[derive(Debug, Clone, Copy)]
pub struct StepResult {
    pub state: State4,
    pub dx_start: Vec4,
    pub dx_end: Vec4,
    pub err: f64,
    pub dl_used: f64,
    pub dl_next: f64,
    pub accepted: bool,
}

const DP_A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];

const DP_B5: [f64; 7] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
    0.0,
];

const DP_B4: [f64; 7] = [
    5179.0 / 57600.0,
    0.0,
    7571.0 / 16695.0,
    393.0 / 640.0,
    -92097.0 / 339200.0,
    187.0 / 2100.0,
    1.0 / 40.0,
];

const STEP_SAFETY: f64 = 0.9;
const STEP_SHRINK_MIN: f64 = 0.2;
const STEP_GROW_MAX: f64 = 5.0;

pub fn rk45_step<M: Metric + ?Sized>(m: &M, s: &State4, dl: f64, tol: &Tolerances) -> StepResult {
    let mut kx = [Vec4::zeros(); 7];
    let mut kp = [Vec4::zeros(); 7];
    (kx[0], kp[0]) = rhs_hamiltonian(m, s);

    for i in 1..7 {
        let mut stage = *s;
        for j in 0..i {
            let a = DP_A[i - 1][j];
            stage.x += dl * a * kx[j];
            stage.p += dl * a * kp[j];
        }
        (kx[i], kp[i]) = rhs_hamiltonian(m, &stage);
    }

    let mut s5 = *s;
    let mut s4 = *s;
    for i in 0..7 {
        s5.x += dl * DP_B5[i] * kx[i];
        s5.p += dl * DP_B5[i] * kp[i];
        s4.x += dl * DP_B4[i] * kx[i];
        s4.p += dl * DP_B4[i] * kp[i];
    }

    let err = error_norm(s, &s5, &s4, tol);
    let accepted = err <= 1.0;
    let scale = if err > 0.0 {
        (STEP_SAFETY * err.powf(-0.2)).clamp(STEP_SHRINK_MIN, STEP_GROW_MAX)
    } else {
        STEP_GROW_MAX
    };

    StepResult {
        state: if accepted { s5 } else { *s },
        dx_start: kx[0],
        dx_end: kx[6],
        err,
        dl_used: dl,
        dl_next: dl * scale,
        accepted,
    }
}

fn error_norm(s0: &State4, s5: &State4, s4: &State4, tol: &Tolerances) -> f64 {
    let mut sum = 0.0;
    for i in 0..4 {
        let sc_x = tol.atol + tol.rtol * s0.x[i].abs().max(s5.x[i].abs());
        sum += ((s5.x[i] - s4.x[i]) / sc_x).powi(2);
        let sc_p = tol.atol + tol.rtol * s0.p[i].abs().max(s5.p[i].abs());
        sum += ((s5.p[i] - s4.p[i]) / sc_p).powi(2);
    }
    (sum / 8.0).sqrt()
}

pub fn hamiltonian<M: Metric + ?Sized>(m: &M, s: &State4) -> f64 {
    let g_inv = m.g_inv(&s.x);
    0.5 * quad_form(&g_inv, &s.p)
}

#[inline]
fn quad_form(m: &Mat4, p: &Vec4) -> f64 {
    let mp = m * p;
    p.dot(&mp)
}

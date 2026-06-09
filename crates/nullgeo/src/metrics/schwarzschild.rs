use crate::metric::{Mat4, Metric, Vec4};

#[derive(Clone, Copy, Debug)]
pub struct Schwarzschild {
    pub m: f64,
}

impl Schwarzschild {
    #[inline]
    fn eta_cov() -> Mat4 {
        Mat4::from_diagonal(&[-1.0, 1.0, 1.0, 1.0].into())
    }
    #[inline]
    fn eta_con() -> Mat4 {
        Self::eta_cov()
    }

    #[inline]
    fn r_n(x: &Vec4) -> (f64, [f64; 3]) {
        let (px, py, pz) = (x[1], x[2], x[3]);
        let r = (px * px + py * py + pz * pz).sqrt().max(1e-20);
        (r, [px / r, py / r, pz / r])
    }

    #[inline]
    fn l_up(n: [f64; 3]) -> Vec4 {
        Vec4::new(-1.0, n[0], n[1], n[2])
    }
    #[inline]
    fn l_dn(n: [f64; 3]) -> Vec4 {
        Vec4::new(1.0, n[0], n[1], n[2])
    }

    #[inline]
    fn h(self, r: f64) -> f64 {
        self.m / r
    }
}

impl Metric for Schwarzschild {
    fn g(&self, x: &Vec4) -> Mat4 {
        let (r, n) = Self::r_n(x);
        let h = self.h(r);
        let l = Self::l_dn(n);
        let mut g = Self::eta_cov();
        for mu in 0..4 {
            for nu in 0..4 {
                g[(mu, nu)] += 2.0 * h * l[mu] * l[nu];
            }
        }
        g
    }

    fn g_inv(&self, x: &Vec4) -> Mat4 {
        let (r, n) = Self::r_n(x);
        let h = self.h(r);
        let l = Self::l_up(n);
        let mut ginv = Self::eta_con();
        for mu in 0..4 {
            for nu in 0..4 {
                ginv[(mu, nu)] -= 2.0 * h * l[mu] * l[nu];
            }
        }
        ginv
    }

    fn dg_inv(&self, x: &Vec4) -> [Mat4; 4] {
        let (r, n) = Self::r_n(x);
        let h = self.h(r);
        let l = Self::l_up(n);

        // dH/dx^i and dl^a/dx^i (time components of l are constant).
        let dh = |i: usize| -> f64 { -h * n[i] / r };
        let dl = |i: usize, alpha: usize| -> f64 {
            if alpha == 0 {
                return 0.0;
            }
            let a = alpha - 1;
            ((i == a) as i32 as f64 - n[i] * n[a]) / r
        };

        let mut out = [Mat4::zeros(); 4];
        for (mu, m) in out.iter_mut().enumerate().skip(1) {
            let i = mu - 1;
            for a in 0..4 {
                for b in 0..4 {
                    let term = dh(i) * l[a] * l[b] + h * dl(i, a) * l[b] + h * l[a] * dl(i, b);
                    m[(a, b)] = -2.0 * term;
                }
            }
        }
        out
    }
}

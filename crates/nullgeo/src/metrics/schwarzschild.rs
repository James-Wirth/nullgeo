use crate::metric::{Mat4, Metric, Vec4};

#[derive(Clone, Copy, Debug)]
pub struct Schwarzschild {
    pub m: f64, 
}

impl Schwarzschild {
    #[inline]
    fn eta_cov() -> Mat4 { Mat4::from_diagonal(&[-1.0, 1.0, 1.0, 1.0].into()) }
    #[inline]
    fn eta_con() -> Mat4 { Self::eta_cov() } 

    #[inline]
    fn r_xyz(x: &Vec4) -> (f64, [f64; 3], [f64; 3]) {
        let px = x[1]; let py = x[2]; let pz = x[3];
        let r = (px*px + py*py + pz*pz).sqrt().max(1e-20);
        let n = [px/r, py/r, pz/r];
        (r, [px, py, pz], n)
    }

    #[inline]
    fn l_up(n: [f64;3]) -> Vec4 { Vec4::new(1.0, n[0], n[1], n[2]) }     
    #[inline]
    fn l_dn(n: [f64;3]) -> Vec4 { Vec4::new(-1.0, n[0], n[1], n[2]) }    

    #[inline]
    fn h(self, r: f64) -> f64 { self.m / r }
}

impl Metric for Schwarzschild {
    fn g(&self, x: &Vec4) -> Mat4 {
        let (r, _xyz, n) = Self::r_xyz(x);
        let h = self.h(r);
        let l = Self::l_dn(n);
        let eta = Self::eta_cov();
        let mut g = eta;
        for mu in 0..4 {
            for nu in 0..4 {
                g[(mu,nu)] += 2.0 * h * l[mu] * l[nu];
            }
        }
        g
    }

    fn g_inv(&self, x: &Vec4) -> Mat4 {
        let (r, _xyz, n) = Self::r_xyz(x);
        let h = self.h(r);
        let l = Self::l_up(n);
        let eta = Self::eta_con();
        let mut ginv = eta;
        for mu in 0..4 {
            for nu in 0..4 {
                ginv[(mu,nu)] -= 2.0 * h * l[mu] * l[nu];
            }
        }
        ginv
    }

    fn dg_inv(&self, x: &Vec4) -> Option<[Mat4; 4]> {
        let (r, _xyz, n) = Self::r_xyz(x);
        let h = self.h(r);
        let l = Self::l_up(n);

        let dh = |i: usize| -> f64 { -h * n[i] / r }; 
        let dl = |i: usize, alpha: usize| -> f64 {
            if alpha == 0 { return 0.0; }
            let a = alpha - 1; 
            ((i == a) as i32 as f64 - n[i]*n[a]) / r
        };

        let mut out = [Mat4::zeros(), Mat4::zeros(), Mat4::zeros(), Mat4::zeros()];

        for mu in 1..4 {
            let i = mu - 1; 
            let mut m = Mat4::zeros();
            for a in 0..4 {
                for b in 0..4 {
                    let term =
                        dh(i) * l[a] * l[b]
                        + h * dl(i, a) * l[b]
                        + h * l[a] * dl(i, b);
                    m[(a,b)] = -2.0 * term;
                }
            }
            out[mu] = m;
        }
        out[0] = Mat4::zeros();
        Some(out)
    }
}

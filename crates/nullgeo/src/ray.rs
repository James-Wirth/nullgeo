use super::metric::{State4, Vec4, Metric};
use super::integrator::rk4_step;

#[derive(Clone, Debug)]
pub struct RayBundle {
    pub x: Vec<Vec4>,
    pub p: Vec<Vec4>,
}

impl RayBundle {
    pub fn new(n: usize) -> Self {
        Self { x: vec![Vec4::zeros(); n], p: vec![Vec4::zeros(); n] }
    }
    pub fn len(&self) -> usize { self.x.len() }
    pub fn is_empty(&self) -> bool { self.x.is_empty() }

    pub fn state(&self, i: usize) -> State4 { State4 { x: self.x[i], p: self.p[i] } }

    pub fn set_state(&mut self, i: usize, s: State4) {
        self.x[i] = s.x; self.p[i] = s.p;
    }

    pub fn step<M: Metric>(&mut self, m: &M, dl: f64) {
        for i in 0..self.len() {
            let s = self.state(i);
            self.set_state(i, rk4_step(m, &s, dl));
        }
    }
}

#[cfg(feature = "parallel")]
impl RayBundle {
    pub fn step_par<M: Metric + Sync>(&mut self, m: &M, dl: f64)
    where M: Send, {
        use rayon::prelude::*;
        let next: Vec<_> = (0..self.len())
            .into_par_iter()
            .map(|i| rk4_step(m, &State4 { x: self.x[i], p: self.p[i] }, dl))
            .collect();
        for (i, s) in next.into_iter().enumerate() { self.set_state(i, s); }
    }
}

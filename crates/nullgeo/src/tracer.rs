use crate::integrator::{rhs_hamiltonian, rk4_step, rk45_step, StepResult, Tolerances};
use crate::metric::State4;
use crate::spacetime::{SkySide, Spacetime};

#[derive(Debug, Clone, Copy)]
pub struct EquatorialAnnulus {
    pub r_in: f64,
    pub r_out: f64,
}

impl EquatorialAnnulus {
    fn band(&self) -> (f64, f64) {
        (0.9 * self.r_in, 1.1 * self.r_out)
    }

    fn dl_clamp(&self) -> f64 {
        self.r_in / 6.0
    }
}

#[derive(Debug, Clone, Copy)]
pub struct TraceConfig {
    pub tol: Tolerances,
    pub dl_init: f64,
    pub dl_min: f64,
    pub dl_max: f64,
    pub max_steps: usize,
    pub escape_radius: f64,
    pub disk: Option<EquatorialAnnulus>,
}

impl Default for TraceConfig {
    fn default() -> Self {
        Self {
            tol: Tolerances::default(),
            dl_init: 0.1,
            dl_min: 1e-9,
            dl_max: 5.0,
            max_steps: 100_000,
            escape_radius: 1.0e3,
            disk: None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Termination {
    Captured {
        state: State4,
    },
    Escaped {
        side: SkySide,
        dir: [f64; 3],
        state: State4,
    },
    HitSurface {
        state: State4,
    },
    MaxSteps {
        state: State4,
    },
    Stalled {
        state: State4,
    },
}

pub fn trace<S: Spacetime + ?Sized>(
    spacetime: &S,
    start: State4,
    cfg: &TraceConfig,
) -> Termination {
    let (mut s, alignment) = spacetime.align_ray(start);
    let mut dl = cfg.dl_init.clamp(cfg.dl_min, cfg.dl_max);

    for _ in 0..cfg.max_steps {
        if spacetime.is_captured(&s.x) {
            return Termination::Captured { state: s };
        }
        if spacetime.radius(&s.x) > cfg.escape_radius {
            return Termination::Escaped {
                side: spacetime.sky_side(&s.x),
                dir: alignment.apply(spacetime.cartesian_direction(&s)),
                state: s,
            };
        }

        let step = rk45_step(spacetime, &s, dl, &cfg.tol);
        if step.accepted {
            if let Some(annulus) = &cfg.disk {
                if let Some(hit) = annulus_crossing(spacetime, &s, &step, annulus) {
                    return Termination::HitSurface { state: hit };
                }
            }
            s = step.state;
        } else if step.dl_used <= cfg.dl_min {
            return Termination::Stalled { state: s };
        }
        dl = step.dl_next.clamp(cfg.dl_min, cfg.dl_max);
        if let Some(annulus) = &cfg.disk {
            let (lo, hi) = annulus.band();
            let r = spacetime.radius(&s.x);
            if r >= lo && r <= hi {
                dl = dl.min(annulus.dl_clamp());
            }
        }
    }

    Termination::MaxSteps { state: s }
}

fn annulus_crossing<S: Spacetime + ?Sized>(
    spacetime: &S,
    s0: &State4,
    step: &StepResult,
    annulus: &EquatorialAnnulus,
) -> Option<State4> {
    let z0 = spacetime.equator_distance(&s0.x);
    let z1 = spacetime.equator_distance(&step.state.x);
    if z0 * z1 >= 0.0 {
        return None;
    }

    let (lo, hi) = annulus.band();
    let in_band = |r: f64| r >= lo && r <= hi;
    if !in_band(spacetime.radius(&s0.x)) && !in_band(spacetime.radius(&step.state.x)) {
        return None;
    }

    let h = step.dl_used;
    let dz0 = rhs_hamiltonian(spacetime, s0).0[3];
    let dz1 = rhs_hamiltonian(spacetime, &step.state).0[3];
    let sigma = hermite_root(z0, z1, h * dz0, h * dz1);

    let hit = rk4_step(spacetime, s0, sigma * h);
    let r_hit = spacetime.radius(&hit.x);
    (annulus.r_in..=annulus.r_out)
        .contains(&r_hit)
        .then_some(hit)
}

fn hermite_root(z0: f64, z1: f64, m0: f64, m1: f64) -> f64 {
    let value = |sig: f64| {
        let s2 = sig * sig;
        let s3 = s2 * sig;
        (2.0 * s3 - 3.0 * s2 + 1.0) * z0
            + (s3 - 2.0 * s2 + sig) * m0
            + (3.0 * s2 - 2.0 * s3) * z1
            + (s3 - s2) * m1
    };

    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if value(mid) * z0 > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

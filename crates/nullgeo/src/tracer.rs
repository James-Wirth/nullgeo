use crate::integrator::{rk45_step, Tolerances};
use crate::metric::State4;
use crate::spacetime::{SkySide, Spacetime};

#[derive(Debug, Clone, Copy)]
pub struct TraceConfig {
    pub tol: Tolerances,
    pub dl_init: f64,
    pub dl_min: f64,
    pub dl_max: f64,
    pub max_steps: usize,
    pub escape_radius: f64,
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
            s = step.state;
        } else if step.dl_used <= cfg.dl_min {
            return Termination::Stalled { state: s };
        }
        dl = step.dl_next.clamp(cfg.dl_min, cfg.dl_max);
    }

    Termination::MaxSteps { state: s }
}

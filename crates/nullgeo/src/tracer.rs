use crate::geometry::{raise, PhasePoint};
use crate::integrator::{rk45_step, rk4_step, StepResult, Tolerances};
use crate::spacetimes::{SkySide, Spacetime};

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
        state: PhasePoint,
    },
    Escaped {
        side: SkySide,
        dir: [f64; 3],
        state: PhasePoint,
    },
    HitSurface {
        state: PhasePoint,
    },
    MaxSteps {
        state: PhasePoint,
    },
    Stalled {
        state: PhasePoint,
    },
}

#[derive(Debug, Clone, Copy)]
pub struct TraceStats {
    pub affine_length: f64,
    pub coord_time: f64,
    pub equatorial_crossings: usize,
    pub min_radius: f64,
    pub steps_accepted: usize,
    pub steps_rejected: usize,
}

pub fn trace<S: Spacetime + ?Sized>(
    spacetime: &S,
    start: PhasePoint,
    cfg: &TraceConfig,
) -> Termination {
    trace_with_stats(spacetime, start, cfg).0
}

pub fn trace_with_stats<S: Spacetime + ?Sized>(
    spacetime: &S,
    start: PhasePoint,
    cfg: &TraceConfig,
) -> (Termination, TraceStats) {
    let (mut s, alignment) = spacetime.align_ray(start);
    let mut dl = cfg.dl_init.clamp(cfg.dl_min, cfg.dl_max);
    let t_start = s.x[0];
    let mut stats = TraceStats {
        affine_length: 0.0,
        coord_time: 0.0,
        equatorial_crossings: 0,
        min_radius: spacetime.radius(&s.x),
        steps_accepted: 0,
        steps_rejected: 0,
    };

    for _ in 0..cfg.max_steps {
        stats.coord_time = (t_start - s.x[0]).abs();
        if spacetime.is_captured(&s.x) {
            return (Termination::Captured { state: s }, stats);
        }
        if spacetime.radius(&s.x) > cfg.escape_radius {
            let v = raise(&spacetime.g_inv(&s.x), &s.p);
            let escaped = Termination::Escaped {
                side: spacetime.sky_side(&s.x),
                dir: alignment.apply(spacetime.embed_direction(&s.x, &v)),
                state: s,
            };
            return (escaped, stats);
        }

        let step = rk45_step(spacetime, &s, dl, &cfg.tol);
        if step.accepted {
            let crossing =
                spacetime.equator_distance(&s.x) * spacetime.equator_distance(&step.state.x) < 0.0;
            stats.steps_accepted += 1;
            if crossing {
                if let Some(annulus) = &cfg.disk {
                    if let Some((hit, dl_to_hit)) = annulus_crossing(spacetime, &s, &step, annulus)
                    {
                        stats.affine_length += dl_to_hit;
                        stats.coord_time = (t_start - hit.x[0]).abs();
                        stats.min_radius = stats.min_radius.min(spacetime.radius(&hit.x));
                        return (Termination::HitSurface { state: hit }, stats);
                    }
                }
                stats.equatorial_crossings += 1;
            }
            stats.affine_length += step.dl_used;
            stats.min_radius = stats.min_radius.min(spacetime.radius(&step.state.x));
            s = step.state;
        } else {
            stats.steps_rejected += 1;
            if step.dl_used <= cfg.dl_min {
                return (Termination::Stalled { state: s }, stats);
            }
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

    stats.coord_time = (t_start - s.x[0]).abs();
    (Termination::MaxSteps { state: s }, stats)
}

fn annulus_crossing<S: Spacetime + ?Sized>(
    spacetime: &S,
    s0: &PhasePoint,
    step: &StepResult,
    annulus: &EquatorialAnnulus,
) -> Option<(PhasePoint, f64)> {
    let z0 = spacetime.equator_distance(&s0.x);
    let z1 = spacetime.equator_distance(&step.state.x);

    let (lo, hi) = annulus.band();
    let in_band = |r: f64| r >= lo && r <= hi;
    if !in_band(spacetime.radius(&s0.x)) && !in_band(spacetime.radius(&step.state.x)) {
        return None;
    }

    let h = step.dl_used;
    let sigma = hermite_root(
        z0,
        z1,
        h * spacetime.equator_distance_rate(&s0.x, &step.dx_start),
        h * spacetime.equator_distance_rate(&step.state.x, &step.dx_end),
    );

    let hit = rk4_step(spacetime, s0, sigma * h);
    let r_hit = spacetime.radius(&hit.x);
    (annulus.r_in..=annulus.r_out)
        .contains(&r_hit)
        .then_some((hit, sigma * h))
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
    for _ in 0..30 {
        let mid = 0.5 * (lo + hi);
        if value(mid) * z0 > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

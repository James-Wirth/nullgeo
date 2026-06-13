use super::camera::Camera;
use super::scene::Scene;
use crate::geometry::PhasePoint;
use crate::spacetimes::{SkySide, Spacetime};
use crate::tracer::{trace_with_stats, EquatorialAnnulus, Termination, TraceConfig, TraceStats};
use crate::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RayClass {
    Captured,
    EscapedPrimary,
    EscapedSecondary,
    Disk,
    MaxSteps,
    Stalled,
}

#[derive(Debug, Clone, Copy)]
pub enum RayOutcome {
    Captured,
    Escaped {
        side: SkySide,
        dir: [f64; 3],
        g: Option<f64>,
    },
    DiskHit {
        radius: f64,
        g: Option<f64>,
    },
    MaxSteps,
    Stalled,
}

#[derive(Debug, Clone, Copy)]
pub struct RayInfo {
    pub outcome: RayOutcome,
    pub stats: TraceStats,
}

impl RayInfo {
    pub fn class(&self) -> RayClass {
        match self.outcome {
            RayOutcome::Captured => RayClass::Captured,
            RayOutcome::Escaped {
                side: SkySide::Primary,
                ..
            } => RayClass::EscapedPrimary,
            RayOutcome::Escaped {
                side: SkySide::Secondary,
                ..
            } => RayClass::EscapedSecondary,
            RayOutcome::DiskHit { .. } => RayClass::Disk,
            RayOutcome::MaxSteps => RayClass::MaxSteps,
            RayOutcome::Stalled => RayClass::Stalled,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GeometryBuffer {
    pub width: usize,
    pub height: usize,
    pub samples: usize,
    pub annulus: Option<EquatorialAnnulus>,
    pub rays: Vec<RayInfo>,
}

impl GeometryBuffer {
    pub fn ray(&self, sample: usize, pixel: usize) -> &RayInfo {
        &self.rays[sample * self.width * self.height + pixel]
    }

    pub fn primary(&self, pixel: usize) -> &RayInfo {
        self.ray(0, pixel)
    }
}

pub fn trace_geometry<S: Spacetime + Sync + ?Sized>(
    spacetime: &S,
    camera: &Camera,
    scene: &Scene,
    cfg: &TraceConfig,
) -> Result<GeometryBuffer> {
    let mut cfg = *cfg;
    if let Some(disk) = &scene.disk {
        let orbits = spacetime.circular_orbits().ok_or_else(|| {
            Error::InvalidArg("this spacetime does not support an equatorial disk".into())
        })?;
        let r_in = disk.r_in.max(orbits.isco_radius());
        if disk.r_out <= r_in {
            return Err(Error::InvalidArg(format!(
                "disk r_out = {} must exceed inner radius {} (after ISCO clamp)",
                disk.r_out, r_in
            )));
        }
        if let crate::render::DiskModel::Blackbody { t_in, .. } = disk.model {
            if !(t_in > 0.0 && t_in.is_finite()) {
                return Err(Error::InvalidArg(format!(
                    "blackbody disk needs a positive finite t_in, got {t_in}"
                )));
            }
        }
        cfg.disk = Some(EquatorialAnnulus {
            r_in,
            r_out: disk.r_out,
        });
    }

    let u_obs = camera.observer_four_velocity(spacetime)?;

    let probe = |ray: &PhasePoint| -> RayInfo {
        let (termination, stats) = trace_with_stats(spacetime, *ray, &cfg);
        let outcome = match termination {
            Termination::Captured { .. } => RayOutcome::Captured,
            Termination::Escaped { side, dir, state } => {
                let killing_energy = state.p[0];
                let g = (killing_energy > 0.0 && killing_energy.is_finite())
                    .then(|| ray.p.dot(&u_obs) / killing_energy);
                RayOutcome::Escaped { side, dir, g }
            }
            Termination::HitSurface { state } => {
                let g = spacetime
                    .circular_orbits()
                    .and_then(|orbits| orbits.four_velocity(&state.x))
                    .map(|u_em| ray.p.dot(&u_obs) / state.p.dot(&u_em));
                RayOutcome::DiskHit {
                    radius: spacetime.radius(&state.x),
                    g,
                }
            }
            Termination::MaxSteps { .. } => RayOutcome::MaxSteps,
            Termination::Stalled { .. } => RayOutcome::Stalled,
        };
        RayInfo { outcome, stats }
    };

    let (width, height) = camera.spec.res;
    let offsets = camera.subpixel_offsets();
    let samples = offsets.len();
    let mut rays_info = Vec::with_capacity(width * height * samples);

    for offset in offsets {
        let rays = camera.pixel_rays_at(spacetime, offset)?;

        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            rays_info.par_extend(rays.par_iter().map(probe));
        }
        #[cfg(not(feature = "parallel"))]
        rays_info.extend(rays.iter().map(probe));
    }

    Ok(GeometryBuffer {
        width,
        height,
        samples,
        annulus: cfg.disk,
        rays: rays_info,
    })
}

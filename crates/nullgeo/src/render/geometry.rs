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
    pub refined: Vec<Vec<RayInfo>>,
}

impl GeometryBuffer {
    pub fn ray(&self, sample: usize, pixel: usize) -> &RayInfo {
        &self.rays[sample * self.width * self.height + pixel]
    }

    pub fn primary(&self, pixel: usize) -> &RayInfo {
        self.ray(0, pixel)
    }

    pub fn pixel_count(&self) -> usize {
        self.width * self.height
    }

    pub fn for_each_sample<F: FnMut(&RayInfo)>(&self, pixel: usize, mut f: F) {
        let refined = &self.refined[pixel];
        if refined.is_empty() {
            for sample in 0..self.samples {
                f(self.ray(sample, pixel));
            }
        } else {
            for info in refined {
                f(info);
            }
        }
    }

    pub fn sample_count(&self, pixel: usize) -> usize {
        let refined = self.refined[pixel].len();
        if refined == 0 {
            self.samples
        } else {
            refined
        }
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
    let pixels = width * height;
    let generator = camera.ray_generator(spacetime)?;

    let base_offsets = camera.subpixel_offsets_for(camera.spec.supersample);
    let samples = base_offsets.len();
    let mut rays = Vec::with_capacity(pixels * samples);
    for offset in &base_offsets {
        let trace_pixel = |pixel: usize| probe(&generator(pixel, *offset));
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            rays.par_extend((0..pixels).into_par_iter().map(trace_pixel));
        }
        #[cfg(not(feature = "parallel"))]
        rays.extend((0..pixels).map(trace_pixel));
    }

    let mut refined = vec![Vec::new(); pixels];
    if camera.spec.supersample_max > camera.spec.supersample {
        let pixel_angle = camera.spec.fov_deg.to_radians() / width.max(1) as f64;
        let mask = refinement_mask(
            &rays[..pixels],
            width,
            height,
            REFINE_CURVATURE * pixel_angle,
        );
        let flagged: Vec<usize> = (0..pixels).filter(|&p| mask[p]).collect();
        let refine_offsets = camera.subpixel_offsets_for(camera.spec.supersample_max);
        let refine_pixel = |&pixel: &usize| -> (usize, Vec<RayInfo>) {
            let samples = refine_offsets
                .iter()
                .map(|offset| probe(&generator(pixel, *offset)))
                .collect();
            (pixel, samples)
        };
        #[cfg(feature = "parallel")]
        let traced: Vec<(usize, Vec<RayInfo>)> = {
            use rayon::prelude::*;
            flagged.par_iter().map(refine_pixel).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let traced: Vec<(usize, Vec<RayInfo>)> = flagged.iter().map(refine_pixel).collect();
        for (pixel, samples) in traced {
            refined[pixel] = samples;
        }
    }

    Ok(GeometryBuffer {
        width,
        height,
        samples,
        annulus: cfg.disk,
        rays,
        refined,
    })
}

const REFINE_CURVATURE: f64 = 1.0;

const NEIGHBORS: [(isize, isize); 8] = [
    (-1, -1),
    (0, -1),
    (1, -1),
    (-1, 0),
    (1, 0),
    (-1, 1),
    (0, 1),
    (1, 1),
];

const AXES: [[(isize, isize); 2]; 2] = [[(-1, 0), (1, 0)], [(0, -1), (0, 1)]];

fn refinement_mask(
    primary: &[RayInfo],
    width: usize,
    height: usize,
    curvature_threshold: f64,
) -> Vec<bool> {
    let at = |i: isize, j: isize| -> Option<&RayInfo> {
        (i >= 0 && j >= 0 && i < width as isize && j < height as isize)
            .then(|| &primary[j as usize * width + i as usize])
    };
    (0..width * height)
        .map(|pixel| {
            let (i, j) = ((pixel % width) as isize, (pixel / width) as isize);
            let here = &primary[pixel];

            let on_boundary = NEIGHBORS
                .iter()
                .any(|&(di, dj)| at(i + di, j + dj).is_some_and(|n| n.class() != here.class()));
            if on_boundary {
                return true;
            }

            let RayOutcome::Escaped { dir: center, .. } = here.outcome else {
                return false;
            };
            let curvature: f64 = AXES
                .iter()
                .filter_map(|[(adi, adj), (bdi, bdj)]| {
                    let a = escaped_dir(at(i + adi, j + adj)?)?;
                    let b = escaped_dir(at(i + bdi, j + bdj)?)?;
                    Some(second_difference(a, center, b))
                })
                .sum();
            curvature > curvature_threshold
        })
        .collect()
}

fn escaped_dir(info: &RayInfo) -> Option<[f64; 3]> {
    match info.outcome {
        RayOutcome::Escaped { dir, .. } => Some(dir),
        _ => None,
    }
}

fn second_difference(a: [f64; 3], center: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [
        a[0] + b[0] - 2.0 * center[0],
        a[1] + b[1] - 2.0 * center[1],
        a[2] + b[2] - 2.0 * center[2],
    ];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

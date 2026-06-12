pub mod camera;
pub mod disk;
pub mod scene;
pub mod sky;

pub use camera::{Camera, CameraPose, CameraSpec};
pub use disk::Disk;
pub use scene::Scene;
pub use sky::{EquirectImage, SkyMap};

use crate::geometry::PhasePoint;
use crate::spacetimes::Spacetime;
use crate::tracer::{trace, EquatorialAnnulus, Termination, TraceConfig};
use crate::{Error, Result};

#[derive(Debug, Clone)]
pub struct ImageF32 {
    pub width: usize,
    pub height: usize,
    pub data: Vec<[f32; 3]>,
}

pub fn render<S: Spacetime + Sync + ?Sized>(
    spacetime: &S,
    camera: &Camera,
    scene: &Scene,
    cfg: &TraceConfig,
) -> Result<ImageF32> {
    let mut cfg = *cfg;
    let mut r_in = 0.0;
    if let Some(disk) = &scene.disk {
        let orbits = spacetime.circular_orbits().ok_or_else(|| {
            Error::InvalidArg("this spacetime does not support an equatorial disk".into())
        })?;
        r_in = disk.r_in.max(orbits.isco_radius());
        if disk.r_out <= r_in {
            return Err(Error::InvalidArg(format!(
                "disk r_out = {} must exceed inner radius {} (after ISCO clamp)",
                disk.r_out, r_in
            )));
        }
        cfg.disk = Some(EquatorialAnnulus {
            r_in,
            r_out: disk.r_out,
        });
    }

    let u_obs = camera.observer_four_velocity(spacetime)?;

    let shade = |ray: &PhasePoint| -> [f32; 3] {
        match trace(spacetime, *ray, &cfg) {
            Termination::Escaped { side, dir, .. } => scene.sky_for(side).sample(dir),
            Termination::HitSurface { state } => {
                let Some(disk) = &scene.disk else {
                    return [0.0; 3];
                };
                let Some(u_em) = spacetime
                    .circular_orbits()
                    .and_then(|o| o.four_velocity(&state.x))
                else {
                    return [0.0; 3];
                };
                let g_factor = ray.p.dot(&u_obs) / state.p.dot(&u_em);
                let r = spacetime.radius(&state.x);
                let brightness =
                    (g_factor.powf(disk.g_power) * (r / r_in).powf(-disk.emissivity_index)) as f32;
                [brightness; 3]
            }
            _ => [0.0; 3],
        }
    };

    let (width, height) = camera.spec.res;
    let mut data = vec![[0.0f32; 3]; width * height];
    let offsets = camera.subpixel_offsets();
    let weight = 1.0 / offsets.len() as f32;

    for offset in offsets {
        let rays = camera.pixel_rays_at(spacetime, offset)?;

        #[cfg(feature = "parallel")]
        let pass: Vec<[f32; 3]> = {
            use rayon::prelude::*;
            rays.par_iter().map(shade).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let pass: Vec<[f32; 3]> = rays.iter().map(shade).collect();

        for (pixel, sample) in data.iter_mut().zip(&pass) {
            for c in 0..3 {
                pixel[c] += weight * sample[c];
            }
        }
    }

    Ok(ImageF32 {
        width,
        height,
        data,
    })
}

pub fn tone_map(image: &ImageF32, exposure: f32) -> Vec<[u8; 3]> {
    image
        .data
        .iter()
        .map(|c| {
            let mut out = [0u8; 3];
            for (byte, &channel) in out.iter_mut().zip(c) {
                let v = (channel * exposure).max(0.0);
                let v = v / (1.0 + v);
                *byte = (v.powf(1.0 / 2.2) * 255.0 + 0.5).clamp(0.0, 255.0) as u8;
            }
            out
        })
        .collect()
}

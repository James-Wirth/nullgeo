use crate::camera::Camera;
use crate::metric::State4;
use crate::scene::Scene;
use crate::spacetime::Spacetime;
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
        let isco = spacetime.isco_radius().ok_or_else(|| {
            Error::InvalidArg("this spacetime does not support an equatorial disk".into())
        })?;
        r_in = disk.r_in.max(isco);
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

    let rays = camera.pixel_rays(spacetime)?;
    let g_tt = spacetime.g(&camera.pose.position)[(0, 0)];
    let u_obs_t = 1.0 / (-g_tt).sqrt();

    let shade = |ray: &State4| -> [f32; 3] {
        match trace(spacetime, *ray, &cfg) {
            Termination::Escaped { side, dir, .. } => scene.sky_for(side).sample(dir),
            Termination::HitSurface { state } => {
                let Some(disk) = &scene.disk else {
                    return [0.0; 3];
                };
                let Some(u_em) = spacetime.disk_emitter(&state.x) else {
                    return [0.0; 3];
                };
                let g_factor = (ray.p[0] * u_obs_t) / state.p.dot(&u_em);
                let r = spacetime.radius(&state.x);
                let brightness = (g_factor.powf(disk.g_power)
                    * (r / r_in).powf(-disk.emissivity_index))
                    as f32;
                [brightness; 3]
            }
            _ => [0.0; 3],
        }
    };

    #[cfg(feature = "parallel")]
    let data = {
        use rayon::prelude::*;
        rays.par_iter().map(shade).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let data = rays.iter().map(shade).collect();

    let (width, height) = camera.spec.res;
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

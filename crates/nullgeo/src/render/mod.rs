mod camera;
mod disk;
mod geometry;
mod maps;
mod scene;
mod sky;

pub use camera::{Camera, CameraPose, CameraSpec};
pub use disk::Disk;
pub use geometry::{trace_geometry, GeometryBuffer, RayClass, RayInfo, RayOutcome};
pub use maps::{class_color, colorize, shade_map, Colormap, MapField, MapQuantity};
pub use scene::Scene;
pub use sky::{EquirectImage, SkyMap};

use crate::spacetimes::Spacetime;
use crate::tracer::TraceConfig;
use crate::Result;

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
    let buffer = trace_geometry(spacetime, camera, scene, cfg)?;
    Ok(shade_beauty(&buffer, scene))
}

pub fn shade_beauty(buffer: &GeometryBuffer, scene: &Scene) -> ImageF32 {
    let pixels = buffer.width * buffer.height;
    let weight = 1.0 / buffer.samples as f32;
    let r_in = buffer.annulus.map(|annulus| annulus.r_in);

    let mut data = vec![[0.0f32; 3]; pixels];
    for sample in 0..buffer.samples {
        for (pixel, info) in data
            .iter_mut()
            .zip(&buffer.rays[sample * pixels..(sample + 1) * pixels])
        {
            let color = shade_ray(info, scene, r_in);
            for c in 0..3 {
                pixel[c] += weight * color[c];
            }
        }
    }

    ImageF32 {
        width: buffer.width,
        height: buffer.height,
        data,
    }
}

fn shade_ray(info: &RayInfo, scene: &Scene, r_in: Option<f64>) -> [f32; 3] {
    match info.outcome {
        RayOutcome::Escaped { side, dir } => scene.sky_for(side).sample(dir),
        RayOutcome::DiskHit { radius, g } => match (&scene.disk, g, r_in) {
            (Some(disk), Some(g), Some(r_in)) => {
                let brightness =
                    (g.powf(disk.g_power) * (radius / r_in).powf(-disk.emissivity_index)) as f32;
                [brightness; 3]
            }
            _ => [0.0; 3],
        },
        _ => [0.0; 3],
    }
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

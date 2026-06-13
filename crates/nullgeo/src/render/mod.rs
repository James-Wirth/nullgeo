mod camera;
mod color;
mod disk;
mod geometry;
mod maps;
mod scene;
mod sky;

pub use camera::{Camera, CameraPose, CameraSpec};
pub use color::{planck_xyz, quantize16, quantize8, tone_map_curve, xyz_to_linear_srgb, ToneCurve};
pub use disk::{shakura_sunyaev_peak_radius, shakura_sunyaev_temperature, Disk, DiskModel};
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
    let disk_shader = DiskShader::new(scene, buffer.annulus.map(|annulus| annulus.r_in));

    let mut data = vec![[0.0f32; 3]; pixels];
    for sample in 0..buffer.samples {
        for (pixel, info) in data
            .iter_mut()
            .zip(&buffer.rays[sample * pixels..(sample + 1) * pixels])
        {
            let color = shade_ray(info, scene, &disk_shader);
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

enum DiskShader {
    None,
    Stylized {
        r_in: f64,
        emissivity_index: f64,
        g_power: f64,
    },
    Blackbody {
        r_in: f64,
        t_in: f64,
        doppler_beaming: bool,
        redshift_color: bool,
        y_ref: f64,
    },
}

impl DiskShader {
    fn new(scene: &Scene, r_in: Option<f64>) -> Self {
        let (Some(disk), Some(r_in)) = (&scene.disk, r_in) else {
            return DiskShader::None;
        };
        match disk.model {
            DiskModel::Stylized {
                emissivity_index,
                g_power,
            } => DiskShader::Stylized {
                r_in,
                emissivity_index,
                g_power,
            },
            DiskModel::Blackbody {
                t_in,
                doppler_beaming,
                redshift_color,
            } => {
                let t_peak =
                    shakura_sunyaev_temperature(t_in, r_in, shakura_sunyaev_peak_radius(r_in));
                DiskShader::Blackbody {
                    r_in,
                    t_in,
                    doppler_beaming,
                    redshift_color,
                    y_ref: planck_xyz(t_peak)[1],
                }
            }
        }
    }

    fn shade(&self, radius: f64, g: f64) -> [f32; 3] {
        match *self {
            DiskShader::None => [0.0; 3],
            DiskShader::Stylized {
                r_in,
                emissivity_index,
                g_power,
            } => {
                let brightness = (g.powf(g_power) * (radius / r_in).powf(-emissivity_index)) as f32;
                [brightness; 3]
            }
            DiskShader::Blackbody {
                r_in,
                t_in,
                doppler_beaming,
                redshift_color,
                y_ref,
            } => {
                let t_emitted = shakura_sunyaev_temperature(t_in, r_in, radius);
                let t_color = if redshift_color {
                    g * t_emitted
                } else {
                    t_emitted
                };
                let t_bright = if doppler_beaming {
                    g * t_emitted
                } else {
                    t_emitted
                };

                let xyz = planck_xyz(t_color);
                if xyz[1] <= 0.0 || y_ref <= 0.0 {
                    return [0.0; 3];
                }
                let scale = if t_bright == t_color {
                    1.0 / y_ref
                } else {
                    planck_xyz(t_bright)[1] / (xyz[1] * y_ref)
                };
                xyz_to_linear_srgb(xyz.map(|c| c * scale)).map(|c| c as f32)
            }
        }
    }
}

fn shade_ray(info: &RayInfo, scene: &Scene, disk_shader: &DiskShader) -> [f32; 3] {
    match info.outcome {
        RayOutcome::Escaped { side, dir, g } => {
            let sample = scene.sky_for(side).sample(dir);
            match g {
                Some(g) => {
                    let boost = (g * g * g * g) as f32;
                    sample.map(|c| c * boost)
                }
                None => sample,
            }
        }
        RayOutcome::DiskHit { radius, g: Some(g) } => disk_shader.shade(radius, g),
        _ => [0.0; 3],
    }
}

pub fn tone_map(image: &ImageF32, exposure: f32) -> Vec<[u8; 3]> {
    quantize8(&tone_map_curve(image, exposure, ToneCurve::Reinhard))
}

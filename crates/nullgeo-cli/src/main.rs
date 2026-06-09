mod io;

use clap::{Parser, Subcommand, ValueEnum};
use nullgeo::metric::Vec4;
use nullgeo::metrics::minkowski::Minkowski;
use nullgeo::metrics::schwarzschild::Schwarzschild;
use nullgeo::{trace, Camera, CameraPose, CameraSpec, Spacetime, Termination, TraceConfig};
use rayon::prelude::*;

#[derive(Copy, Clone, Debug, ValueEnum)]
enum MetricKind {
    Minkowski,
    Schwarzschild,
}

#[derive(Parser, Debug)]
#[command(name = "nullgeo", about = "nullgeo - general relativistic ray tracing")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Propagate {
        #[arg(long, default_value_t = 0.1)]
        dl: f64,
        #[arg(long, default_value_t = 10)]
        steps: usize,
    },

    Render {
        #[arg(long, default_value_t = 32)]
        width: usize,
        #[arg(long, default_value_t = 32)]
        height: usize,
        #[arg(long, default_value_t = 20.0)]
        fov_deg: f64,
        #[arg(long, default_value_t = 1.0)]
        energy: f64,
        #[arg(long, default_value_t = 0)]
        steps: usize,
        #[arg(long, default_value_t = 0.05)]
        dl: f64,
    },

    Shadow {
        #[arg(long, value_enum, default_value_t=MetricKind::Schwarzschild)]
        metric: MetricKind,
        #[arg(long, default_value_t = 1.0)]
        mass: f64,
        #[arg(long, default_value_t = 256)]
        width: usize,
        #[arg(long, default_value_t = 256)]
        height: usize,
        #[arg(long, default_value_t = 20.0)]
        fov_deg: f64,
        #[arg(long, default_value_t=-15.0)]
        cam_x: f64,
        #[arg(long, default_value_t = 1.0)]
        energy: f64,
        #[arg(long, default_value_t = 0.01)]
        dl: f64,
        #[arg(long, default_value_t = 5000)]
        max_steps: usize,
        #[arg(long, default_value = "shadow.ppm")]
        out: String,
    },
}

fn main() {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Command::Propagate { .. } => {
            eprintln!("'propagate' not yet written");
            std::process::exit(1);
        }
        Command::Render { .. } => {
            eprintln!("'render' not yet written");
            std::process::exit(1);
        }

        Command::Shadow {
            metric,
            mass,
            width,
            height,
            fov_deg,
            cam_x,
            energy,
            dl,
            max_steps,
            out,
        } => {
            let camera = match Camera::new(
                CameraSpec {
                    fov_deg,
                    res: (width, height),
                    energy,
                },
                CameraPose {
                    position: Vec4::new(0.0, cam_x, 0.0, 0.0),
                    look_at: [0.0, 0.0, 0.0],
                    up: [0.0, 0.0, 1.0],
                },
            ) {
                Ok(camera) => camera,
                Err(e) => {
                    eprintln!("invalid camera: {e}");
                    std::process::exit(1);
                }
            };

            let cfg = TraceConfig {
                dl_init: dl,
                max_steps,
                escape_radius: 4.0 * cam_x.abs().max(10.0 * mass.abs()),
                ..TraceConfig::default()
            };

            let result = match metric {
                MetricKind::Minkowski => shadow_image(&Minkowski, &camera, &cfg),
                MetricKind::Schwarzschild => {
                    shadow_image(&Schwarzschild { m: mass }, &camera, &cfg)
                }
            };

            let img = match result {
                Ok(img) => img,
                Err(e) => {
                    eprintln!("trace failed: {e}");
                    std::process::exit(1);
                }
            };

            if let Err(e) = io::write_ppm_gray(&out, width, height, &img) {
                eprintln!("Failed to write {}: {}", out, e);
            } else {
                println!("Wrote {}", out);
            }
        }
    }
}

fn shadow_image<S: Spacetime + Sync>(
    spacetime: &S,
    camera: &Camera,
    cfg: &TraceConfig,
) -> nullgeo::Result<Vec<u8>> {
    let rays = camera.pixel_rays(spacetime)?;
    let img = rays
        .into_par_iter()
        .map(|ray| match trace(spacetime, ray, cfg) {
            Termination::Escaped { .. } => 255u8,
            _ => 0u8,
        })
        .collect();
    Ok(img)
}

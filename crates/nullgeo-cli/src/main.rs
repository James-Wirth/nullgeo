mod io;

use clap::{Parser, Subcommand, ValueEnum};
use nullgeo::integrator::rk4_step;
use nullgeo::metric::{Mat4, Metric, State4, Vec4};
use nullgeo::metrics::minkowski::Minkowski;
use nullgeo::metrics::schwarzschild::Schwarzschild;
use nullgeo::{Camera, CameraSpec};
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
            let cam = Camera::from_spec(CameraSpec {
                fov_deg,
                res: (width, height),
                energy,
            });
            let x0 = Vec4::new(0.0, cam_x, 0.0, 0.0);
            let r_h = 2.0 * mass;
            let x_escape = cam_x + (-cam_x + 5.0 * r_h.max(1.0));

            let img = match metric {
                MetricKind::Minkowski => {
                    shadow_image(&Minkowski, &cam, x0, energy, dl, max_steps, None, x_escape)
                }
                MetricKind::Schwarzschild => shadow_image(
                    &Schwarzschild { m: mass },
                    &cam,
                    x0,
                    energy,
                    dl,
                    max_steps,
                    Some(1.05 * r_h),
                    x_escape,
                ),
            };

            if let Err(e) = io::write_ppm_gray(&out, width, height, &img) {
                eprintln!("Failed to write {}: {}", out, e);
            } else {
                println!("Wrote {}", out);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn shadow_image<M: Metric + Sync>(
    m: &M,
    cam: &Camera,
    x0: Vec4,
    energy: f64,
    dl: f64,
    max_steps: usize,
    r_cap: Option<f64>,
    x_escape: f64,
) -> Vec<u8> {
    cam.generate_directions()
        .par_iter()
        .map(|&n| {
            let mut s = State4 {
                x: x0,
                p: make_null_covector(m, &x0, n, energy),
            };
            let mut captured = false;
            for _ in 0..max_steps {
                s = rk4_step(m, &s, dl);
                let r = (s.x[1] * s.x[1] + s.x[2] * s.x[2] + s.x[3] * s.x[3]).sqrt();
                if let Some(rc) = r_cap {
                    if r < rc {
                        captured = true;
                        break;
                    }
                }
                if s.x[1] > x_escape || r > 1.0e6 {
                    break;
                }
            }
            if captured {
                0u8
            } else {
                255u8
            }
        })
        .collect()
}

#[inline]
fn lower(g: &Mat4, v: &Vec4) -> Vec4 {
    g * v
}
#[inline]
fn dot(g: &Mat4, a: &Vec4, b: &Vec4) -> f64 {
    (a.transpose() * (*g) * (*b))[0]
}
#[inline]
fn norm_timelike(g: &Mat4, v: &Vec4) -> f64 {
    (-dot(g, v, v)).sqrt()
}
#[inline]
fn norm_spacelike(g: &Mat4, v: &Vec4) -> f64 {
    (dot(g, v, v)).sqrt()
}
#[inline]
fn normalize_timelike(g: &Mat4, v: &Vec4) -> Vec4 {
    v / norm_timelike(g, v).max(1e-300)
}
#[inline]
fn normalize_spacelike(g: &Mat4, v: &Vec4) -> Vec4 {
    v / norm_spacelike(g, v).max(1e-300)
}
#[inline]
fn project_out(g: &Mat4, v: &Vec4, onto: &Vec4) -> Vec4 {
    let num = dot(g, v, onto);
    let den = dot(g, onto, onto);
    if den.abs() < 1e-300 {
        *v
    } else {
        v - onto * (num / den)
    }
}

fn build_coframe<M: Metric>(metric: &M, x: &Vec4) -> [Vec4; 4] {
    let g = metric.g(x);

    let mut t = Vec4::new(1.0, 0.0, 0.0, 0.0);
    if g[(0, 0)].abs() < 1e-14 {
        t = Vec4::new(1.0, 1e-8, 0.0, 0.0);
    }
    let mut e0 = normalize_timelike(&g, &t);
    if e0[0] < 0.0 {
        e0 = -e0;
    }

    let ex = project_out(&g, &Vec4::new(0.0, 1.0, 0.0, 0.0), &e0);
    let e1 = normalize_spacelike(&g, &ex);

    let mut ey = project_out(&g, &Vec4::new(0.0, 0.0, 1.0, 0.0), &e0);
    ey = project_out(&g, &ey, &e1);
    let e2 = normalize_spacelike(&g, &ey);

    let mut ez = project_out(&g, &Vec4::new(0.0, 0.0, 0.0, 1.0), &e0);
    ez = project_out(&g, &ez, &e1);
    ez = project_out(&g, &ez, &e2);
    let e3 = normalize_spacelike(&g, &ez);

    [
        lower(&g, &e0),
        lower(&g, &e1),
        lower(&g, &e2),
        lower(&g, &e3),
    ]
}

pub fn make_null_covector<M: Metric>(metric: &M, x: &Vec4, n: [f64; 3], e: f64) -> Vec4 {
    let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-300);
    let nx = n[0] / norm;
    let ny = n[1] / norm;
    let nz = n[2] / norm;

    let th = build_coframe(metric, x);
    let mut p = th[0] * e;
    p += th[1] * (e * nx);
    p += th[2] * (e * ny);
    p += th[3] * (e * nz);
    p
}

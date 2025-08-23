mod io; 

use clap::{Parser, Subcommand, ValueEnum};
use rayon::prelude::*;
use nullgeo_core::metric::{State4, Vec4, Metric};
use nullgeo_core::integrator::{rk4_step};
use nullgeo_core::{Camera, CameraSpec};
use nullgeo_metrics::minkowski::Minkowski;
use nullgeo_metrics::schwarzschild::Schwarzschild;

#[derive(Copy, Clone, Debug, ValueEnum)]
enum MetricKind { Minkowski, Schwarzschild }

#[derive(Parser, Debug)]
#[command(name="nullgeo", about="GR Ray Tracing")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Propagate {
        #[arg(long, default_value_t=0.1)] dl: f64,
        #[arg(long, default_value_t=10)] steps: usize
    },

    Render {
        #[arg(long, default_value_t=32)] width: usize,
        #[arg(long, default_value_t=32)] height: usize,
        #[arg(long, default_value_t=20.0)] fov_deg: f64,
        #[arg(long, default_value_t=1.0)] energy: f64,
        #[arg(long, default_value_t=0)] steps: usize,
        #[arg(long, default_value_t=0.05)] dl: f64
    },

    Shadow {
        #[arg(long, value_enum, default_value_t=MetricKind::Schwarzschild)]
        metric: MetricKind,
        #[arg(long, default_value_t=1.0)]
        mass: f64,
        #[arg(long, default_value_t=256)]
        width: usize,
        #[arg(long, default_value_t=256)]
        height: usize,
        #[arg(long, default_value_t=20.0)]
        fov_deg: f64,
        #[arg(long, default_value_t=-15.0)]
        cam_x: f64,   
        #[arg(long, default_value_t=1.0)]
        energy: f64,
        #[arg(long, default_value_t=0.01)]
        dl: f64,
        #[arg(long, default_value_t=5000)]
        max_steps: usize,
        #[arg(long, default_value = "shadow.ppm")]
        out: String,
    }
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

        Command::Shadow { metric, mass, width, height, fov_deg, cam_x, energy, dl, max_steps, out } => {

        enum AnyMetric { M(Minkowski), S(Schwarzschild) }
        let metric = match metric {
            MetricKind::Minkowski => AnyMetric::M(Minkowski),
            MetricKind::Schwarzschild => AnyMetric::S(Schwarzschild { m: mass }),
        };

        let cam = Camera::from_spec(CameraSpec { fov_deg, res: (width, height), energy });
        let rays_p_cov = cam.generate_rays();

        let dirs_opt: Option<Vec<[f64; 3]>> = match &metric {
            AnyMetric::S(_) => Some(cam.generate_directions()),
            _ => None,
        };

        let x0 = Vec4::new(0.0, cam_x, 0.0, 0.0);
        let r_h = 2.0 * mass;
        let r_cap = 1.05 * r_h;
        let x_far = -cam_x + 5.0 * r_h.max(1.0);

        let mut img = vec![255u8; width * height];

        img.par_iter_mut()
            .enumerate()
            .for_each(|(idx, out_px)| {
                let mut s = State4 { x: x0, p: Vec4::zeros() };

                match &metric {
                    AnyMetric::M(m) => {
                        s.p = rays_p_cov[idx];
                        for _ in 0..max_steps {
                            s = rk4_step(m, &s, dl);
                            if s.x[1] > cam_x + x_far { break; }
                        }
                        *out_px = 255;
                    }
                    AnyMetric::S(m) => {
                        let n = dirs_opt.as_ref().unwrap()[idx];
                        s.p = make_null_covector(m, &s.x, n, energy);

                        let mut captured = false;
                        for _ in 0..max_steps {
                            s = rk4_step(m, &s, dl);

                            let xi = s.x[1]; let yi = s.x[2]; let zi = s.x[3];
                            let r = (xi*xi + yi*yi + zi*zi).sqrt();

                            if r < r_cap { captured = true; break; }
                            if s.x[1] > cam_x + x_far { break; }
                            if r > 1.0e6 { break; }
                        }
                        *out_px = if captured { 0 } else { 255 };
                    }
                }
            });

            if let Err(e) = io::write_ppm_gray(&out, width, height, &img) {
                eprintln!("Failed to write {}: {}", out, e);
            } else {
                println!("Wrote {}", out);
            }
        }
    }
}

use nullgeo_core::metric::{Mat4};
fn make_null_covector<M: Metric>(metric: &M, x: &Vec4, n: [f64;3], e: f64) -> Vec4 {
    let ginv: Mat4 = metric.g_inv(x);
    let qi = [e*n[0], e*n[1], e*n[2]];

    let mut a = ginv[(0,0)];
    if a.abs() < 1e-15 {
        let s = if a >= 0.0 { 1.0 } else { -1.0 }; 
        a = s * 1e-15;
    }

    let b = 2.0 * (ginv[(0,1)]*qi[0] + ginv[(0,2)]*qi[1] + ginv[(0,3)]*qi[2]);
    let c = ginv[(1,1)]*qi[0]*qi[0]
          + ginv[(2,2)]*qi[1]*qi[1]
          + ginv[(3,3)]*qi[2]*qi[2]
          + 2.0*( ginv[(1,2)]*qi[0]*qi[1] + ginv[(1,3)]*qi[0]*qi[2] + ginv[(2,3)]*qi[1]*qi[2] );
    let disc = (b*b - 4.0*a*c).max(0.0);
    let sqrt_disc = disc.sqrt();
    
    let r1 = (-b - sqrt_disc)/(2.0*a);
    let r2 = (-b + sqrt_disc)/(2.0*a);
    let (rneg, rpos) = if r1 <= r2 { (r1, r2) } else { (r2, r1) };
    let mut p0 = if rneg.is_finite() { rneg } else { rpos };
    if p0 > 0.0 && rneg.is_finite() { p0 = rneg; }

    Vec4::new(p0, qi[0], qi[1], qi[2])
}

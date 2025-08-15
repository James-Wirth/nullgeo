use clap::{Parser, Subcommand};
use nullgeo_core::metric::{State4, Vec4};
use nullgeo_core::integrator::{rk4_step, hamiltonian};
use nullgeo_core::{Camera, CameraSpec, RayBundle};
use nullgeo_metrics::minkowski::Minkowski;
use log::{info, warn};

#[derive(Parser, Debug)]
#[command(name="nullgeo", about="nullgeo")]
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
    }
}

fn main() {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Command::Propagate { dl, steps } => {
            let metric = Minkowski;

            let e = 1.0;
            let mut state = State4 { x: Vec4::zeros(), p: Vec4::new(-e, e, 0.0, 0.0) };

            let h0 = hamiltonian(&metric, &state);
            for _ in 0..steps { state = rk4_step(&metric, &state, dl); }
            let h1 = hamiltonian(&metric, &state);

            println!("Final x: {:?}", state.x);
            println!("H0 = {:.15e}, H1 = {:.15e}, ΔH = {:+.3e}", h0, h1, h1 - h0);
        }

        Command::Render { width, height, fov_deg, energy, steps, dl } => {
            let cam = Camera::from_spec(CameraSpec { fov_deg, res: (width, height), energy });
            let ps = cam.generate_rays();
            info!("Generated {} primary rays", ps.len());

            let mut bundle = RayBundle::new(ps.len());
            for (i, p) in ps.into_iter().enumerate() {
                bundle.x[i] = Vec4::zeros();
                bundle.p[i] = p;
            }

            if steps > 0 {
                let metric = Minkowski;
                for _ in 0..steps {
                    #[cfg(feature = "parallel")]
                    { bundle.step_par(&metric, dl); }
                    #[cfg(not(feature = "parallel"))]
                    { bundle.step(&metric, dl); }
                }

                let mut hmean = 0.0;
                for i in 0..bundle.len() {
                    let s = nullgeo_core::metric::State4 { x: bundle.x[i], p: bundle.p[i] };
                    hmean += hamiltonian(&Minkowski, &s);
                }
                hmean /= bundle.len() as f64;
                println!("Bundle rays: {}, steps: {}, ⟨H⟩ ≈ {:.3e}", bundle.len(), steps, hmean);
            } else {
                println!("Bundle rays: {}", bundle.len());
                warn!("No propagation carried out (steps=0).");
            }
        }
    }
}

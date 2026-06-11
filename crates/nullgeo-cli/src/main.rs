mod io;
mod scene_file;

use std::path::Path;

use clap::{Args, Parser, Subcommand};
use nullgeo::metric::Vec4;
use nullgeo::{render, tone_map, Camera, CameraPose, CameraSpec, Scene, SkyMap, TraceConfig};
use scene_file::{
    build_camera, build_disk, build_sky, build_spacetime, build_trace_config, MetricKind,
    MetricSection, OutputFormat, SceneFile,
};

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
        scene: String,
    },

    Shadow(ShadowArgs),
}

#[derive(Args, Debug)]
struct ShadowArgs {
    #[arg(long, value_enum, default_value_t = MetricKind::Schwarzschild)]
    metric: MetricKind,
    #[arg(long, default_value_t = 1.0)]
    mass: f64,
    #[arg(long, default_value_t = 0.0)]
    spin: f64,
    #[arg(long, default_value_t = 0.0)]
    charge: f64,
    #[arg(long, default_value_t = 1.0)]
    b0: f64,
    #[arg(long, default_value_t = 256)]
    width: usize,
    #[arg(long, default_value_t = 256)]
    height: usize,
    #[arg(long, default_value_t = 60.0)]
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
}

fn main() {
    env_logger::init();
    let cli = Cli::parse();

    let result = match cli.command {
        Command::Propagate { .. } => Err("'propagate' not yet written".to_string()),
        Command::Render { scene } => run_render(&scene),
        Command::Shadow(args) => run_shadow(&args),
    };

    if let Err(e) = result {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run_render(path: &str) -> Result<(), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let file: SceneFile =
        toml::from_str(&text).map_err(|e| format!("invalid scene file {path}: {e}"))?;
    let base = Path::new(path).parent().unwrap_or(Path::new("."));

    let spacetime = build_spacetime(&file.metric)?;
    let camera = build_camera(&file.camera)?;
    let scene = Scene {
        sky: build_sky(&file.sky, base)?,
        sky_secondary: file
            .sky_secondary
            .as_ref()
            .map(|s| build_sky(s, base))
            .transpose()?,
        disk: file.disk.as_ref().map(build_disk),
    };
    let cfg = build_trace_config(&file.integrator, file.camera.position);

    let img = render(spacetime.as_ref(), &camera, &scene, &cfg).map_err(|e| e.to_string())?;
    let pixels = tone_map(&img, file.output.exposure);

    let out = &file.output.path;
    match file.output.resolved_format() {
        OutputFormat::Png => {
            let flat: Vec<u8> = pixels.iter().flatten().copied().collect();
            image::save_buffer(
                out,
                &flat,
                img.width as u32,
                img.height as u32,
                image::ExtendedColorType::Rgb8,
            )
            .map_err(|e| format!("failed to write {}: {e}", out.display()))?;
        }
        OutputFormat::Ppm => {
            let path_str = out
                .to_str()
                .ok_or_else(|| format!("non-utf8 output path {}", out.display()))?;
            io::write_ppm_rgb(path_str, img.width, img.height, &pixels)
                .map_err(|e| format!("failed to write {}: {e}", out.display()))?;
        }
    }
    println!("Wrote {}", out.display());
    Ok(())
}

fn run_shadow(args: &ShadowArgs) -> Result<(), String> {
    let camera = Camera::new(
        CameraSpec {
            fov_deg: args.fov_deg,
            res: (args.width, args.height),
            energy: args.energy,
        },
        CameraPose {
            position: Vec4::new(0.0, args.cam_x, 0.0, 0.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .map_err(|e| format!("invalid camera: {e}"))?;

    let cfg = TraceConfig {
        dl_init: args.dl,
        max_steps: args.max_steps,
        escape_radius: 4.0 * args.cam_x.abs().max(10.0 * args.mass.abs()),
        ..TraceConfig::default()
    };

    let spacetime = build_spacetime(&MetricSection {
        kind: args.metric,
        mass: args.mass,
        spin: args.spin,
        charge: args.charge,
        b0: args.b0,
    })?;

    let scene = Scene {
        sky: SkyMap::Uniform([1.0; 3]),
        sky_secondary: None,
        disk: None,
    };
    let img = render(spacetime.as_ref(), &camera, &scene, &cfg)
        .map_err(|e| format!("trace failed: {e}"))?;
    let gray: Vec<u8> = img
        .data
        .iter()
        .map(|c| if c[0] > 0.5 { 255 } else { 0 })
        .collect();

    io::write_ppm_gray(&args.out, args.width, args.height, &gray)
        .map_err(|e| format!("failed to write {}: {e}", args.out))?;
    println!("Wrote {}", args.out);
    Ok(())
}

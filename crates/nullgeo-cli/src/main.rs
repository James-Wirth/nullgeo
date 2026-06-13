mod io;
mod scene_file;

use std::path::Path;

use clap::{Args, Parser, Subcommand};
use nullgeo::geometry::Vec4;
use nullgeo::integrator::{hamiltonian, rk45_step, Tolerances};
use nullgeo::{
    colorize, quantize16, quantize8, render, shade_beauty, shade_map, tone_map_curve,
    trace_geometry, Camera, CameraPose, CameraSpec, Colormap, Scene, SkyMap, TraceConfig,
};
use scene_file::{
    build_camera, build_disk, build_sky, build_spacetime, build_trace_config, MetricKind,
    MetricSection, OutputFormat, OutputSection, SceneFile,
};

#[derive(Parser, Debug)]
#[command(name = "nullgeo", about = "nullgeo - general relativistic ray tracing")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Propagate(PropagateArgs),

    Render { scene: String },

    Shadow(ShadowArgs),
}

#[derive(Args, Debug)]
struct PropagateArgs {
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
    #[arg(
        long,
        required = true,
        value_delimiter = ',',
        allow_hyphen_values = true
    )]
    pos: Vec<f64>,
    #[arg(
        long,
        required = true,
        value_delimiter = ',',
        allow_hyphen_values = true
    )]
    dir: Vec<f64>,
    #[arg(long, default_value_t = 1.0)]
    energy: f64,
    #[arg(long, default_value_t = 1e-9)]
    tol: f64,
    #[arg(long, default_value_t = 100_000)]
    max_steps: usize,
    #[arg(long)]
    escape_radius: Option<f64>,
    #[arg(long)]
    out: Option<String>,
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
    let cli = Cli::parse();

    let result = match cli.command {
        Command::Propagate(args) => run_propagate(&args),
        Command::Render { scene } => run_render(&scene),
        Command::Shadow(args) => run_shadow(&args),
    };

    if let Err(e) = result {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run_propagate(args: &PropagateArgs) -> Result<(), String> {
    let spacetime = build_spacetime(&MetricSection {
        kind: args.metric,
        mass: args.mass,
        spin: args.spin,
        charge: args.charge,
        b0: args.b0,
    })?;

    let (pos, dir) = (&args.pos, &args.dir);
    if pos.len() != 3 || dir.len() != 3 {
        return Err("--pos and --dir each need three components, e.g. --pos=-20,0,0".into());
    }
    let x = Vec4::new(0.0, pos[0], pos[1], pos[2]);
    let embedded = spacetime.embed(&x);
    let look_at = [
        embedded[0] + dir[0],
        embedded[1] + dir[1],
        embedded[2] + dir[2],
    ];
    let dir_len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
    let up = if dir[2].abs() < 0.9 * dir_len {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res: (1, 1),
            energy: args.energy,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: x,
            look_at,
            up,
            velocity: [0.0; 3],
        },
    )
    .map_err(|e| format!("invalid ray: {e}"))?;
    let ray = camera
        .pixel_rays(spacetime.as_ref())
        .map_err(|e| e.to_string())?[0];

    let defaults = TraceConfig::default();
    let escape_radius = args.escape_radius.unwrap_or_else(|| {
        let r = (embedded[0] * embedded[0] + embedded[1] * embedded[1] + embedded[2] * embedded[2])
            .sqrt();
        (4.0 * r).max(100.0)
    });
    let tol = Tolerances {
        rtol: args.tol,
        atol: args.tol,
    };

    let mut out: Box<dyn std::io::Write> = match &args.out {
        Some(path) => {
            Box::new(std::fs::File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?)
        }
        None => Box::new(std::io::stdout().lock()),
    };
    let write_err = |e: std::io::Error| format!("write failed: {e}");

    writeln!(out, "lambda,t,x,y,z,H").map_err(write_err)?;
    let mut emit = |lambda: f64, s: &nullgeo::PhasePoint| -> Result<(), String> {
        let [px, py, pz] = spacetime.embed(&s.x);
        let h = hamiltonian(spacetime.as_ref(), s);
        writeln!(out, "{lambda},{},{px},{py},{pz},{h}", s.x[0]).map_err(write_err)
    };

    let mut s = ray;
    let mut lambda = 0.0;
    let mut dl = defaults.dl_init;
    let mut status = "max steps reached";
    emit(lambda, &s)?;
    for _ in 0..args.max_steps {
        if spacetime.is_captured(&s.x) {
            status = "captured";
            break;
        }
        if spacetime.radius(&s.x) > escape_radius {
            status = "escaped";
            break;
        }
        let step = rk45_step(spacetime.as_ref(), &s, dl, &tol);
        if step.accepted {
            s = step.state;
            lambda += step.dl_used;
            emit(lambda, &s)?;
        } else if step.dl_used <= defaults.dl_min {
            status = "stalled";
            break;
        }
        dl = step.dl_next.clamp(defaults.dl_min, defaults.dl_max);
    }
    eprintln!(
        "{status} at lambda = {lambda}, r = {}",
        spacetime.radius(&s.x)
    );
    Ok(())
}

fn run_render(path: &str) -> Result<(), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let file: SceneFile =
        toml::from_str(&text).map_err(|e| format!("invalid scene file {path}: {e}"))?;
    if file.outputs.is_empty() {
        return Err("scene file needs at least one [[output]] entry".into());
    }
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
        disk: file.disk.as_ref().map(build_disk).transpose()?,
    };
    let cfg = build_trace_config(&file.integrator, file.camera.position);

    let buffer =
        trace_geometry(spacetime.as_ref(), &camera, &scene, &cfg).map_err(|e| e.to_string())?;

    for output in &file.outputs {
        write_output(&buffer, &scene, output)?;
        println!("Wrote {}", output.path.display());
    }
    Ok(())
}

fn write_output(
    buffer: &nullgeo::GeometryBuffer,
    scene: &Scene,
    output: &OutputSection,
) -> Result<(), String> {
    let out = &output.path;
    let path_str = out
        .to_str()
        .ok_or_else(|| format!("non-utf8 output path {}", out.display()))?;
    let format = output.resolved_format();

    let result = match output.kind.map_quantity() {
        None => {
            let img = shade_beauty(buffer, scene);
            let curve = output.tone.map(|t| t.to_curve()).unwrap_or_default();
            match format {
                OutputFormat::Png | OutputFormat::Ppm => {
                    let display = tone_map_curve(&img, output.exposure, curve);
                    match (format, output.bit_depth.unwrap_or(8)) {
                        (_, 8) => write_rgb(
                            &quantize8(&display),
                            img.width,
                            img.height,
                            format,
                            path_str,
                        ),
                        (OutputFormat::Png, 16) => {
                            write_png16(&quantize16(&display), img.width, img.height, path_str)
                        }
                        (_, 16) => Err("bit_depth = 16 needs png output".into()),
                        (_, other) => Err(format!("bit_depth must be 8 or 16, got {other}")),
                    }
                }
                OutputFormat::Pfm => {
                    if output.tone.is_some() || output.bit_depth.is_some() {
                        return Err(format!(
                            "{}: pfm export is raw linear radiance; tone and bit_depth do not \
                             apply",
                            out.display()
                        ));
                    }
                    io::write_pfm_rgb(path_str, img.width, img.height, &img.data)
                        .map_err(|e| e.to_string())
                }
                OutputFormat::Csv => Err("a beauty render cannot be exported as csv".into()),
            }
        }
        Some(quantity) => {
            if output.tone.is_some() || output.bit_depth.is_some() {
                return Err(format!(
                    "{}: tone and bit_depth apply to beauty outputs only",
                    out.display()
                ));
            }
            let field = shade_map(buffer, quantity);
            match format {
                OutputFormat::Png | OutputFormat::Ppm => {
                    let colormap = output
                        .colormap
                        .map(|c| c.to_colormap())
                        .unwrap_or_else(|| Colormap::default_for(quantity));
                    let pixels = colorize(&field, colormap);
                    write_rgb(&pixels, field.width, field.height, format, path_str)
                }
                OutputFormat::Pfm => {
                    io::write_pfm_gray(path_str, field.width, field.height, &field.values)
                        .map_err(|e| e.to_string())
                }
                OutputFormat::Csv => {
                    io::write_csv_matrix(path_str, field.width, field.height, &field.values)
                        .map_err(|e| e.to_string())
                }
            }
        }
    };
    result.map_err(|e| format!("failed to write {}: {e}", out.display()))
}

fn write_rgb(
    pixels: &[[u8; 3]],
    width: usize,
    height: usize,
    format: OutputFormat,
    path: &str,
) -> Result<(), String> {
    match format {
        OutputFormat::Png => {
            let flat: Vec<u8> = pixels.iter().flatten().copied().collect();
            image::save_buffer(
                path,
                &flat,
                width as u32,
                height as u32,
                image::ExtendedColorType::Rgb8,
            )
            .map_err(|e| e.to_string())
        }
        OutputFormat::Ppm => {
            io::write_ppm_rgb(path, width, height, pixels).map_err(|e| e.to_string())
        }
        _ => unreachable!(),
    }
}

fn write_png16(pixels: &[[u16; 3]], width: usize, height: usize, path: &str) -> Result<(), String> {
    let raw: Vec<u16> = pixels.iter().flatten().copied().collect();
    let buffer =
        image::ImageBuffer::<image::Rgb<u16>, _>::from_raw(width as u32, height as u32, raw)
            .ok_or("pixel buffer size mismatch")?;
    buffer
        .save_with_format(path, image::ImageFormat::Png)
        .map_err(|e| e.to_string())
}

fn run_shadow(args: &ShadowArgs) -> Result<(), String> {
    let camera = Camera::new(
        CameraSpec {
            fov_deg: args.fov_deg,
            res: (args.width, args.height),
            energy: args.energy,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
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

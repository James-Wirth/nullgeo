use std::time::Instant;

use nullgeo::geometry::Vec4;
use nullgeo::spacetimes::Kerr;
use nullgeo::{
    render, trace_geometry, Camera, CameraPose, CameraSpec, Disk, Scene, SkyMap, TraceConfig,
};

const WIDTH: usize = 96;
const HEIGHT: usize = 54;

fn camera(supersample: usize, supersample_max: usize) -> Camera {
    Camera::new(
        CameraSpec {
            fov_deg: 24.0,
            res: (WIDTH, HEIGHT),
            energy: 1.0,
            supersample,
            supersample_max,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, -45.0, 0.0, 6.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap()
}

fn scene() -> Scene {
    Scene {
        sky: SkyMap::Uniform([0.0; 3]),
        sky_secondary: None,
        disk: Some(Disk::blackbody(14.0, 10_000.0)),
    }
}

fn config() -> TraceConfig {
    TraceConfig {
        escape_radius: 200.0,
        max_steps: 50_000,
        ..TraceConfig::default()
    }
}

fn time<T>(label: &str, run: impl FnOnce() -> T) -> T {
    let start = Instant::now();
    let value = run();
    println!("{label:<28} {:>8.3} s", start.elapsed().as_secs_f64());
    value
}

fn main() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let scene = scene();
    let cfg = config();
    let pixels = WIDTH * HEIGHT;

    println!("Kerr a=0.9 blackbody disk, {WIDTH}x{HEIGHT}");

    time("uniform 1x1", || {
        render(&kerr, &camera(1, 1), &scene, &cfg).unwrap()
    });
    time("uniform 3x3", || {
        render(&kerr, &camera(3, 3), &scene, &cfg).unwrap()
    });

    let buffer = time("adaptive 1->3", || {
        trace_geometry(&kerr, &camera(1, 3), &scene, &cfg).unwrap()
    });

    let refined = buffer.refined.iter().filter(|r| !r.is_empty()).count();
    let effective = pixels + refined * (9 - 1);
    println!(
        "adaptive refined {refined}/{pixels} pixels, {effective} rays vs {} uniform 3x3",
        pixels * 9
    );
}

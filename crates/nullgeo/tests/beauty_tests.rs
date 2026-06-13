use nullgeo::geometry::Vec4;
use nullgeo::integrator::Tolerances;
use nullgeo::spacetimes::{Kerr, Minkowski, Schwarzschild};
use nullgeo::{
    planck_xyz, shade_beauty, shakura_sunyaev_peak_radius, shakura_sunyaev_temperature,
    trace_geometry, xyz_to_linear_srgb, Camera, CameraPose, CameraSpec, Disk, DiskModel, ImageF32,
    RayOutcome, Scene, SkyMap, TraceConfig,
};

fn camera(
    position: [f64; 3],
    look_at: [f64; 3],
    velocity: [f64; 3],
    fov_deg: f64,
    res: (usize, usize),
) -> Camera {
    Camera::new(
        CameraSpec {
            fov_deg,
            res,
            energy: 1.0,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, position[0], position[1], position[2]),
            look_at,
            up: [0.0, 0.0, 1.0],
            velocity,
        },
    )
    .unwrap()
}

fn black_sky_scene(disk: Disk) -> Scene {
    Scene {
        sky: SkyMap::Uniform([0.0; 3]),
        sky_secondary: None,
        disk: Some(disk),
    }
}

#[test]
fn face_on_blackbody_disk_shades_planck_of_g_times_t() {
    let m = Schwarzschild::new(1.0).unwrap();
    let z_cam = 1000.0;
    let t_in = 1.0e4;
    let cam = camera([0.0, 0.0, z_cam], [12.0, 0.0, 0.0], [0.0; 3], 2.0, (3, 3));
    let scene = black_sky_scene(Disk::blackbody(40.0, t_in));
    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-10,
            atol: 1e-10,
        },
        escape_radius: 4000.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    };

    let buffer = trace_geometry(&m, &cam, &scene, &cfg).unwrap();
    let image = shade_beauty(&buffer, &scene);

    let center = 4;
    let RayOutcome::DiskHit { radius, g: Some(g) } = buffer.primary(center).outcome else {
        panic!("center pixel should hit the disk");
    };

    let g_predicted = (1.0 - 3.0 / radius).sqrt() / (1.0 - 2.0 / z_cam).sqrt();
    assert!((g - g_predicted).abs() < 1e-5 * g_predicted);

    let r_in = buffer.annulus.unwrap().r_in;
    assert_eq!(r_in, 6.0);
    let t_observed = g * shakura_sunyaev_temperature(t_in, r_in, radius);
    let y_ref = planck_xyz(shakura_sunyaev_temperature(
        t_in,
        r_in,
        shakura_sunyaev_peak_radius(r_in),
    ))[1];
    let expected = xyz_to_linear_srgb(planck_xyz(t_observed).map(|c| c / y_ref)).map(|c| c as f32);

    assert_eq!(image.data[center], expected);
    assert!(expected[0] > 0.0 && expected[1] > 0.0);
    assert!(
        expected[0] != expected[1] || expected[1] != expected[2],
        "a blackbody at {t_observed:.0} K should not be pure gray"
    );
}

#[test]
fn suppressing_redshift_color_keeps_the_emitted_chromaticity() {
    let m = Schwarzschild::new(1.0).unwrap();
    let t_in = 1.0e4;
    let cam = camera([0.0, 0.0, 1000.0], [12.0, 0.0, 0.0], [0.0; 3], 2.0, (3, 3));
    let suppressed = Disk {
        r_in: 0.0,
        r_out: 40.0,
        model: DiskModel::Blackbody {
            t_in,
            doppler_beaming: true,
            redshift_color: false,
        },
    };
    let scene = black_sky_scene(suppressed);
    let cfg = TraceConfig {
        escape_radius: 4000.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    };

    let buffer = trace_geometry(&m, &cam, &scene, &cfg).unwrap();
    let image = shade_beauty(&buffer, &scene);

    let center = 4;
    let RayOutcome::DiskHit { radius, .. } = buffer.primary(center).outcome else {
        panic!("center pixel should hit the disk");
    };

    let t_emitted = shakura_sunyaev_temperature(t_in, 6.0, radius);
    let unshifted = xyz_to_linear_srgb(planck_xyz(t_emitted));
    let pixel = image.data[center];
    let chroma_expected = unshifted[2] / unshifted[0];
    let chroma_actual = (pixel[2] / pixel[0]) as f64;
    assert!(
        (chroma_actual - chroma_expected).abs() < 1e-6 * chroma_expected,
        "blue/red ratio {chroma_actual} should match the emitted blackbody {chroma_expected}"
    );
}

fn edge_on_kerr_asymmetries() -> (f64, f64) {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let cam = camera([-100.0, 0.0, 10.0], [0.0; 3], [0.0; 3], 16.0, (24, 24));
    let disk_with = |doppler_beaming| Disk {
        r_in: 0.0,
        r_out: 12.0,
        model: DiskModel::Blackbody {
            t_in: 1.0e4,
            doppler_beaming,
            redshift_color: true,
        },
    };
    let cfg = TraceConfig {
        escape_radius: 500.0,
        max_steps: 200_000,
        ..TraceConfig::default()
    };

    let scene_on = black_sky_scene(disk_with(true));
    let scene_off = black_sky_scene(disk_with(false));
    let buffer = trace_geometry(&kerr, &cam, &scene_on, &cfg).unwrap();
    (
        half_ratio(&shade_beauty(&buffer, &scene_on)),
        half_ratio(&shade_beauty(&buffer, &scene_off)),
    )
}

fn half_ratio(img: &ImageF32) -> f64 {
    let mut left = 0.0;
    let mut right = 0.0;
    for (i, px) in img.data.iter().enumerate() {
        let brightness = (px[0] + px[1] + px[2]) as f64;
        if i % img.width < img.width / 2 {
            left += brightness;
        } else {
            right += brightness;
        }
    }
    left / right
}

#[test]
fn doppler_beaming_toggle_controls_the_brightness_asymmetry() {
    let (with_beaming, without_beaming) = edge_on_kerr_asymmetries();
    assert!(
        with_beaming > 1.5,
        "beaming should brighten the approaching side: {with_beaming}"
    );
    assert!(
        without_beaming < 0.75 * with_beaming,
        "suppressing beaming should flatten the asymmetry: {without_beaming} vs {with_beaming}"
    );
}

#[test]
fn boosted_camera_sees_the_forward_sky_beamed_by_g_fourth() {
    let beta = 0.5;
    let cam = camera(
        [-10.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [beta, 0.0, 0.0],
        40.0,
        (3, 3),
    );
    let scene = Scene {
        sky: SkyMap::Uniform([1.0; 3]),
        sky_secondary: None,
        disk: None,
    };
    let cfg = TraceConfig {
        escape_radius: 50.0,
        ..TraceConfig::default()
    };
    let buffer = trace_geometry(&Minkowski, &cam, &scene, &cfg).unwrap();
    let image = shade_beauty(&buffer, &scene);

    let gamma = 1.0 / (1.0 - beta * beta).sqrt();
    let g_forward = gamma * (1.0 + beta);
    let expected = (g_forward.powi(4)) as f32;
    let center = 4;
    assert!(
        (image.data[center][0] - expected).abs() < 1e-5 * expected,
        "forward pixel {} vs gamma(1+beta))^4 = {expected}",
        image.data[center][0]
    );

    let corner = 0;
    assert!(
        image.data[corner][0] < image.data[center][0],
        "off-axis pixels should be less beamed than the forward direction"
    );
}

#[test]
fn static_camera_sky_carries_the_gravitational_blueshift() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cam = camera([-15.0, 0.0, 0.0], [-30.0, 0.0, 0.0], [0.0; 3], 30.0, (3, 3));
    let scene = Scene {
        sky: SkyMap::Uniform([1.0; 3]),
        sky_secondary: None,
        disk: None,
    };
    let cfg = TraceConfig {
        escape_radius: 2000.0,
        ..TraceConfig::default()
    };
    let buffer = trace_geometry(&m, &cam, &scene, &cfg).unwrap();

    let g_predicted = 1.0 / (1.0 - 2.0 / 15.0_f64).sqrt();
    for pixel in 0..9 {
        let RayOutcome::Escaped { g: Some(g), .. } = buffer.primary(pixel).outcome else {
            panic!("pixel {pixel} looking away from the hole should escape");
        };
        assert!(
            (g - g_predicted).abs() < 1e-6 * g_predicted,
            "pixel {pixel}: g = {g} vs 1/sqrt(1 - 2M/r) = {g_predicted}"
        );
    }
}

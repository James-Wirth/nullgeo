use nullgeo::geometry::{build_coframe, make_null_covector};
use nullgeo::geometry::{Metric, PhasePoint, Vec4};
use nullgeo::integrator::Tolerances;
use nullgeo::spacetimes::kerr::Kerr;
use nullgeo::spacetimes::minkowski::Minkowski;
use nullgeo::spacetimes::schwarzschild::Schwarzschild;
use nullgeo::{
    render, tone_map, trace, Camera, CameraPose, CameraSpec, Chart, Disk, EquatorialAnnulus,
    EquirectImage, ImageF32, Scene, SkyMap, SkySide, Spacetime, Termination, TraceConfig,
};

fn backward_ray<M: Metric>(m: &M, x: Vec4, dir: [f64; 3], energy: f64) -> PhasePoint {
    let coframe = build_coframe(m, &x).unwrap();
    let arriving = make_null_covector(&coframe, [-dir[0], -dir[1], -dir[2]], energy);
    PhasePoint { x, p: -arriving }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalized(a: [f64; 3]) -> [f64; 3] {
    let len = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / len, a[1] / len, a[2] / len]
}

#[test]
fn minkowski_render_matches_direct_sky_lookup() {
    let position = [2.0, -3.0, 1.0];
    let look_at = [3.0, 1.0, 2.0];
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 40.0,
            res: (8, 6),
            energy: 1.0,
            supersample: 1,
        },
        CameraPose {
            position: Vec4::new(0.0, position[0], position[1], position[2]),
            look_at,
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap();
    let scene = Scene {
        sky: SkyMap::Checker {
            angular_size_deg: 23.0,
        },
        sky_secondary: None,
        disk: None,
    };
    let cfg = TraceConfig {
        escape_radius: 50.0,
        ..TraceConfig::default()
    };
    let img = render(&Minkowski, &camera, &scene, &cfg).unwrap();

    let forward = normalized([
        look_at[0] - position[0],
        look_at[1] - position[1],
        look_at[2] - position[2],
    ]);
    let right = normalized(cross(forward, [0.0, 0.0, 1.0]));
    let up = cross(right, forward);

    for (idx, [f, r, u]) in camera.pixel_directions().into_iter().enumerate() {
        let world = [
            f * forward[0] + r * right[0] + u * up[0],
            f * forward[1] + r * right[1] + u * up[1],
            f * forward[2] + r * right[2] + u * up[2],
        ];
        let expected = scene.sky_for(SkySide::Primary).sample(world);
        assert_eq!(img.data[idx], expected, "pixel {idx} mismatch");
    }
}

#[test]
fn supersampled_render_averages_subpixel_sky_samples() {
    let position = [2.0, -3.0, 1.0];
    let look_at = [3.0, 1.0, 2.0];
    let make_camera = |supersample| {
        Camera::new(
            CameraSpec {
                fov_deg: 40.0,
                res: (8, 6),
                energy: 1.0,
                supersample,
            },
            CameraPose {
                position: Vec4::new(0.0, position[0], position[1], position[2]),
                look_at,
                up: [0.0, 0.0, 1.0],
                velocity: [0.0; 3],
            },
        )
        .unwrap()
    };
    let scene = Scene {
        sky: SkyMap::Checker {
            angular_size_deg: 23.0,
        },
        sky_secondary: None,
        disk: None,
    };
    let cfg = TraceConfig {
        escape_radius: 50.0,
        ..TraceConfig::default()
    };
    let camera = make_camera(2);
    let img = render(&Minkowski, &camera, &scene, &cfg).unwrap();
    let single = render(&Minkowski, &make_camera(1), &scene, &cfg).unwrap();

    let forward = normalized([
        look_at[0] - position[0],
        look_at[1] - position[1],
        look_at[2] - position[2],
    ]);
    let right = normalized(cross(forward, [0.0, 0.0, 1.0]));
    let up = cross(right, forward);

    let offsets = camera.subpixel_offsets();
    assert_eq!(offsets.len(), 4);
    let mut expected = vec![[0.0f32; 3]; 8 * 6];
    for offset in offsets {
        for (idx, [f, r, u]) in camera.pixel_directions_at(offset).into_iter().enumerate() {
            let world = [
                f * forward[0] + r * right[0] + u * up[0],
                f * forward[1] + r * right[1] + u * up[1],
                f * forward[2] + r * right[2] + u * up[2],
            ];
            let sample = scene.sky_for(SkySide::Primary).sample(world);
            for (c, &value) in sample.iter().enumerate() {
                expected[idx][c] += 0.25 * value;
            }
        }
    }

    for (idx, (got, want)) in img.data.iter().zip(&expected).enumerate() {
        for c in 0..3 {
            assert!(
                (got[c] - want[c]).abs() < 1e-6,
                "pixel {idx} channel {c}: {} vs {}",
                got[c],
                want[c]
            );
        }
    }
    assert!(
        img.data != single.data,
        "supersampling should smooth at least one checker-edge pixel"
    );
}

fn face_on_ray(m: &Schwarzschild, z0: f64, b_aim: f64) -> PhasePoint {
    let alpha = (b_aim / z0).atan();
    backward_ray(
        m,
        Vec4::new(0.0, 0.0, 0.0, z0),
        [alpha.sin(), 0.0, -alpha.cos()],
        1.0,
    )
}

fn face_on_config() -> TraceConfig {
    TraceConfig {
        tol: Tolerances {
            rtol: 1e-10,
            atol: 1e-10,
        },
        escape_radius: 4000.0,
        max_steps: 1_000_000,
        disk: Some(EquatorialAnnulus {
            r_in: 6.0,
            r_out: 40.0,
        }),
        ..TraceConfig::default()
    }
}

#[test]
fn face_on_disk_redshift_matches_formula() {
    let m = Schwarzschild::new(1.0).unwrap();
    let z_cam = 1000.0;
    let cfg = face_on_config();

    for b_aim in [8.0, 12.0, 20.0] {
        let ray = face_on_ray(&m, z_cam, b_aim);
        let Termination::HitSurface { state } = trace(&m, ray, &cfg) else {
            panic!("ray aimed at b = {b_aim} should hit the disk");
        };
        let r = m.radius(&state.x);
        assert!((6.0..=40.0).contains(&r), "hit radius {r}");
        assert!(
            state.x[3].abs() < 1e-5,
            "hit off the plane: z = {}",
            state.x[3]
        );

        let u_em = m
            .circular_orbits()
            .unwrap()
            .four_velocity(&state.x)
            .unwrap();
        let g_tt_cam = m.g(&ray.x)[(0, 0)];
        let e_obs = ray.p[0] / (-g_tt_cam).sqrt();
        let g_factor = e_obs / state.p.dot(&u_em);
        let predicted = (1.0 - 3.0 / r).sqrt() / (1.0 - 2.0 / z_cam).sqrt();
        assert!(
            (g_factor - predicted).abs() < 2e-6 * predicted,
            "g = {g_factor:.9} vs predicted {predicted:.9} at r = {r}"
        );
    }
}

#[test]
fn plane_crossings_outside_annulus_do_not_hit() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cfg = face_on_config();

    match trace(&m, face_on_ray(&m, 1000.0, 2.0), &cfg) {
        Termination::Captured { .. } => {}
        other => panic!("ray through the disk hole should be captured, got {other:?}"),
    }
    match trace(&m, face_on_ray(&m, 1000.0, 60.0), &cfg) {
        Termination::Escaped { .. } => {}
        other => panic!("ray crossing beyond r_out should escape, got {other:?}"),
    }
}

fn edge_on_disk_image(spin: f64) -> ImageF32 {
    let kerr = Kerr::new(1.0, spin).unwrap();
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 16.0,
            res: (24, 24),
            energy: 1.0,
            supersample: 1,
        },
        CameraPose {
            position: Vec4::new(0.0, -100.0, 0.0, 10.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap();
    let scene = Scene {
        sky: SkyMap::Uniform([0.0; 3]),
        sky_secondary: None,
        disk: Some(Disk::new(12.0)),
    };
    let cfg = TraceConfig {
        escape_radius: 500.0,
        max_steps: 200_000,
        ..TraceConfig::default()
    };
    render(&kerr, &camera, &scene, &cfg).unwrap()
}

fn half_sums(img: &ImageF32) -> (f64, f64) {
    let mut left = 0.0;
    let mut right = 0.0;
    for (i, px) in img.data.iter().enumerate() {
        if i % img.width < img.width / 2 {
            left += px[0] as f64;
        } else {
            right += px[0] as f64;
        }
    }
    (left, right)
}

#[test]
fn edge_on_doppler_asymmetry_flips_with_spin_sign() {
    let (left_pro, right_pro) = half_sums(&edge_on_disk_image(0.9));
    assert!(
        left_pro > 1.5 * right_pro,
        "approaching side should beam brighter: left {left_pro:.3e}, right {right_pro:.3e}"
    );

    let (left_retro, right_retro) = half_sums(&edge_on_disk_image(-0.9));
    assert!(
        right_retro > 1.5 * left_retro,
        "asymmetry should flip with spin: left {left_retro:.3e}, right {right_retro:.3e}"
    );
}

#[test]
fn equirect_sampling_is_exact_at_texel_centers() {
    let (w, h) = (8, 4);
    let mut data = Vec::with_capacity(w * h);
    for j in 0..h {
        for i in 0..w {
            data.push([i as f32, j as f32, (i * j) as f32 + 0.5]);
        }
    }
    let image = EquirectImage::new(w, h, data.clone()).unwrap();
    let sky = SkyMap::Equirect(image.clone());

    for j in 0..h {
        for i in 0..w {
            let u = (i as f64 + 0.5) / w as f64;
            let v = (j as f64 + 0.5) / h as f64;
            assert_eq!(image.sample(u, v), data[j * w + i], "texel ({i},{j})");

            let theta = v * std::f64::consts::PI;
            let phi = u * 2.0 * std::f64::consts::PI - std::f64::consts::PI;
            let dir = [
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ];
            let sampled = sky.sample(dir);
            for (c, &value) in sampled.iter().enumerate() {
                assert!(
                    (value - data[j * w + i][c]).abs() < 1e-3,
                    "sky lookup at texel ({i},{j}) channel {c}: {} vs {}",
                    value,
                    data[j * w + i][c]
                );
            }
        }
    }
}

#[test]
fn tone_map_is_monotone_and_fixes_black() {
    let img = ImageF32 {
        width: 6,
        height: 1,
        data: vec![
            [0.0; 3], [0.05; 3], [0.3; 3], [1.0; 3], [5.0; 3], [100.0; 3],
        ],
    };
    let mapped = tone_map(&img, 1.0);
    assert_eq!(mapped[0], [0, 0, 0]);
    for pair in mapped.windows(2) {
        assert!(pair[0][0] < pair[1][0], "tone map not monotone: {mapped:?}");
    }
    let brighter = tone_map(&img, 4.0);
    for (a, b) in mapped.iter().zip(&brighter).skip(1) {
        assert!(b[0] > a[0], "exposure should brighten nonzero pixels");
    }
}

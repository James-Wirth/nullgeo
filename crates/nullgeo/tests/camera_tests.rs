use nullgeo::geometry::{build_coframe_for, inner};
use nullgeo::geometry::{Metric, Vec4};
use nullgeo::integrator::hamiltonian;
use nullgeo::spacetimes::{Minkowski, Schwarzschild};
use nullgeo::{Camera, CameraPose, CameraSpec, Error};

fn test_camera(res: (usize, usize)) -> Camera {
    Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res,
            energy: 1.0,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, -12.0, 5.0, 3.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap()
}

#[test]
fn coframe_is_orthonormal() {
    let m = Schwarzschild::new(1.0).unwrap();
    let x = Vec4::new(0.0, -12.0, 5.0, 3.0);
    let seeds = [
        Vec4::new(0.0, 0.8, -0.5, 0.1),
        Vec4::new(0.0, 0.2, 0.9, -0.4),
        Vec4::new(0.0, -0.1, 0.3, 1.0),
    ];
    let coframe = build_coframe_for(&m, &x, &Vec4::new(1.0, 0.0, 0.0, 0.0), seeds).unwrap();

    let g_inv = m.g_inv(&x);
    let eta = [-1.0, 1.0, 1.0, 1.0];
    for a in 0..4 {
        for b in 0..4 {
            let expected = if a == b { eta[a] } else { 0.0 };
            let got = inner(&g_inv, &coframe[a], &coframe[b]);
            assert!(
                (got - expected).abs() < 1e-12,
                "coframe[{a}].coframe[{b}] = {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn pixel_rays_are_null() {
    let m = Schwarzschild::new(1.0).unwrap();
    let camera = test_camera((7, 5));
    for ray in camera.pixel_rays(&m).unwrap() {
        let h = hamiltonian(&m, &ray);
        assert!(h.abs() < 1e-12, "H = {h:.3e}");
    }
}

fn boosted_camera(velocity: [f64; 3]) -> Camera {
    Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res: (3, 3),
            energy: 2.5,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, -12.0, 5.0, 3.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity,
        },
    )
    .unwrap()
}

#[test]
fn boosted_camera_rays_are_null_with_observer_frame_energy() {
    let m = Schwarzschild::new(1.0).unwrap();
    let camera = boosted_camera([0.05, -0.1, 0.2]);
    let u = camera.observer_four_velocity(&m).unwrap();

    let g = m.g(&camera.pose.position);
    let u_norm = u.dot(&(g * u));
    assert!((u_norm + 1.0).abs() < 1e-12, "u.u = {u_norm}");

    for ray in camera.pixel_rays(&m).unwrap() {
        let h = hamiltonian(&m, &ray);
        assert!(h.abs() < 1e-12, "H = {h:.3e}");
        let energy = ray.p.dot(&u);
        assert!(
            (energy - 2.5).abs() < 1e-12,
            "observer-frame energy {energy}"
        );
    }
}

#[test]
fn superluminal_camera_rejected_at_ray_construction() {
    let camera = boosted_camera([1.5, 0.0, 0.0]);
    assert!(matches!(
        camera.pixel_rays(&Minkowski),
        Err(Error::NonTimelikeObserver(_))
    ));
}

#[test]
fn pixel_rays_run_backward_in_time() {
    let m = Schwarzschild::new(1.0).unwrap();
    let camera = test_camera((5, 5));
    for ray in camera.pixel_rays(&m).unwrap() {
        let velocity = m.g_inv(&ray.x) * ray.p;
        assert!(velocity[0] < 0.0, "dt/dlambda = {}", velocity[0]);
    }
}

#[test]
fn pixel_rays_carry_killing_energy_of_static_observer() {
    let mass = 1.0;
    let m = Schwarzschild::new(mass).unwrap();
    let camera = test_camera((5, 5));
    let r = (12.0_f64 * 12.0 + 5.0 * 5.0 + 3.0 * 3.0).sqrt();
    let expected = (1.0 - 2.0 * mass / r).sqrt();
    for ray in camera.pixel_rays(&m).unwrap() {
        assert!(
            (ray.p[0] - expected).abs() < 1e-12,
            "q_t = {}, expected {expected}",
            ray.p[0]
        );
    }
}

#[test]
fn minkowski_center_pixel_points_at_look_at() {
    let m = Minkowski;
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 45.0,
            res: (3, 3),
            energy: 2.0,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, -10.0, 4.0, -6.0),
            look_at: [2.0, -1.0, 3.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap();

    let rays = camera.pixel_rays(&m).unwrap();
    let center = rays[4];

    let to_target = [2.0_f64 - (-10.0), -1.0 - 4.0, 3.0 - (-6.0)];
    let len =
        (to_target[0] * to_target[0] + to_target[1] * to_target[1] + to_target[2] * to_target[2])
            .sqrt();

    assert!((center.p[0] - 2.0).abs() < 1e-12);
    for (i, &target) in to_target.iter().enumerate() {
        let expected = 2.0 * target / len;
        assert!(
            (center.p[i + 1] - expected).abs() < 1e-12,
            "spatial component {i}: {} vs {expected}",
            center.p[i + 1]
        );
    }
}

fn sampling_camera(supersample: usize, supersample_max: usize, jitter: bool) -> Camera {
    Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res: (4, 4),
            energy: 1.0,
            supersample,
            supersample_max,
            jitter,
        },
        CameraPose {
            position: Vec4::new(0.0, -12.0, 0.0, 0.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap()
}

#[test]
fn stratified_offsets_are_centered_without_jitter() {
    let camera = sampling_camera(1, 3, false);
    let offsets = camera.subpixel_offsets_for(3);
    assert_eq!(offsets.len(), 9);
    let expected: Vec<(f64, f64)> = (0..9)
        .map(|k| ((k % 3) as f64 + 0.5) / 3.0)
        .zip((0..9).map(|k| ((k / 3) as f64 + 0.5) / 3.0))
        .collect();
    assert_eq!(offsets, expected);
    assert_eq!(camera.subpixel_offsets(), camera.subpixel_offsets_for(1));
}

#[test]
fn jittered_offsets_stay_inside_their_cells() {
    let camera = sampling_camera(3, 3, true);
    let offsets = camera.subpixel_offsets_for(3);
    assert_eq!(offsets.len(), 9);
    let mut moved = false;
    for (k, &(x, y)) in offsets.iter().enumerate() {
        let (cell_x, cell_y) = (k % 3, k / 3);
        assert!((cell_x as f64 / 3.0..=(cell_x as f64 + 1.0) / 3.0).contains(&x));
        assert!((cell_y as f64 / 3.0..=(cell_y as f64 + 1.0) / 3.0).contains(&y));
        if (x - (cell_x as f64 + 0.5) / 3.0).abs() > 1e-12 {
            moved = true;
        }
    }
    assert!(moved, "jitter should perturb at least one sample");
    assert_eq!(offsets, sampling_camera(3, 3, true).subpixel_offsets_for(3));
}

#[test]
fn supersample_max_below_supersample_is_rejected() {
    assert!(matches!(
        Camera::new(
            CameraSpec {
                fov_deg: 60.0,
                res: (4, 4),
                energy: 1.0,
                supersample: 3,
                supersample_max: 2,
                jitter: false,
            },
            CameraPose {
                position: Vec4::new(0.0, -12.0, 0.0, 0.0),
                look_at: [0.0, 0.0, 0.0],
                up: [0.0, 0.0, 1.0],
                velocity: [0.0; 3],
            },
        ),
        Err(Error::InvalidArg(_))
    ));
}

#[test]
fn camera_rejects_bad_configuration() {
    let pose = CameraPose {
        position: Vec4::new(0.0, -10.0, 0.0, 0.0),
        look_at: [0.0, 0.0, 0.0],
        up: [0.0, 0.0, 1.0],
        velocity: [0.0; 3],
    };
    let spec = |fov_deg, res, energy| CameraSpec {
        fov_deg,
        res,
        energy,
        supersample: 1,
        supersample_max: 1,
        jitter: false,
    };

    assert!(matches!(
        Camera::new(spec(0.0, (4, 4), 1.0), pose),
        Err(Error::InvalidArg(_))
    ));
    assert!(matches!(
        Camera::new(spec(60.0, (0, 4), 1.0), pose),
        Err(Error::InvalidArg(_))
    ));
    assert!(matches!(
        Camera::new(spec(60.0, (4, 4), 0.0), pose),
        Err(Error::InvalidArg(_))
    ));
    assert!(matches!(
        Camera::new(
            spec(60.0, (4, 4), 1.0),
            CameraPose {
                up: [-1.0, 0.0, 0.0],
                look_at: [0.0, 0.0, 0.0],
                position: Vec4::new(0.0, -10.0, 0.0, 0.0),
                velocity: [0.0; 3],
            }
        ),
        Err(Error::InvalidArg(_))
    ));
}

#[test]
fn observer_inside_horizon_is_rejected() {
    let m = Schwarzschild::new(1.0).unwrap();
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res: (2, 2),
            energy: 1.0,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, 1.0, 0.0, 0.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap();
    assert!(matches!(
        camera.pixel_rays(&m),
        Err(Error::NonTimelikeObserver(_))
    ));
}

use nullgeo::frame::{build_coframe_seeded, metric_dot};
use nullgeo::integrator::hamiltonian;
use nullgeo::metric::{Metric, Vec4};
use nullgeo::metrics::minkowski::Minkowski;
use nullgeo::metrics::schwarzschild::Schwarzschild;
use nullgeo::{Camera, CameraPose, CameraSpec, Error};

fn test_camera(res: (usize, usize)) -> Camera {
    Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res,
            energy: 1.0,
        },
        CameraPose {
            position: Vec4::new(0.0, -12.0, 5.0, 3.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
        },
    )
    .unwrap()
}

#[test]
fn coframe_is_orthonormal() {
    let m = Schwarzschild { m: 1.0 };
    let x = Vec4::new(0.0, -12.0, 5.0, 3.0);
    let seeds = [
        Vec4::new(0.0, 0.8, -0.5, 0.1),
        Vec4::new(0.0, 0.2, 0.9, -0.4),
        Vec4::new(0.0, -0.1, 0.3, 1.0),
    ];
    let coframe = build_coframe_seeded(&m, &x, seeds).unwrap();

    let g_inv = m.g_inv(&x);
    let eta = [-1.0, 1.0, 1.0, 1.0];
    for a in 0..4 {
        for b in 0..4 {
            let expected = if a == b { eta[a] } else { 0.0 };
            let got = metric_dot(&g_inv, &coframe[a], &coframe[b]);
            assert!(
                (got - expected).abs() < 1e-12,
                "coframe[{a}].coframe[{b}] = {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn pixel_rays_are_null() {
    let m = Schwarzschild { m: 1.0 };
    let camera = test_camera((7, 5));
    for ray in camera.pixel_rays(&m).unwrap() {
        let h = hamiltonian(&m, &ray);
        assert!(h.abs() < 1e-12, "H = {h:.3e}");
    }
}

#[test]
fn pixel_rays_run_backward_in_time() {
    let m = Schwarzschild { m: 1.0 };
    let camera = test_camera((5, 5));
    for ray in camera.pixel_rays(&m).unwrap() {
        let velocity = m.g_inv(&ray.x) * ray.p;
        assert!(velocity[0] < 0.0, "dt/dlambda = {}", velocity[0]);
    }
}

#[test]
fn pixel_rays_carry_killing_energy_of_static_observer() {
    let mass = 1.0;
    let m = Schwarzschild { m: mass };
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
        },
        CameraPose {
            position: Vec4::new(0.0, -10.0, 4.0, -6.0),
            look_at: [2.0, -1.0, 3.0],
            up: [0.0, 0.0, 1.0],
        },
    )
    .unwrap();

    let rays = camera.pixel_rays(&m).unwrap();
    let center = rays[4];

    let to_target = [2.0_f64 - (-10.0), -1.0 - 4.0, 3.0 - (-6.0)];
    let len = (to_target[0] * to_target[0]
        + to_target[1] * to_target[1]
        + to_target[2] * to_target[2])
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

#[test]
fn camera_rejects_bad_configuration() {
    let pose = CameraPose {
        position: Vec4::new(0.0, -10.0, 0.0, 0.0),
        look_at: [0.0, 0.0, 0.0],
        up: [0.0, 0.0, 1.0],
    };
    let spec = |fov_deg, res, energy| CameraSpec {
        fov_deg,
        res,
        energy,
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
            }
        ),
        Err(Error::InvalidArg(_))
    ));
}

#[test]
fn observer_inside_horizon_is_rejected() {
    let m = Schwarzschild { m: 1.0 };
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res: (2, 2),
            energy: 1.0,
        },
        CameraPose {
            position: Vec4::new(0.0, 1.0, 0.0, 0.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
        },
    )
    .unwrap();
    assert!(matches!(
        camera.pixel_rays(&m),
        Err(Error::NonTimelikeObserver(_))
    ));
}

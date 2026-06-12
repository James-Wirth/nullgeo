use std::f64::consts::{FRAC_PI_2, PI};

use nullgeo::geometry::{fd_partials, raise, Mat4, Metric, PhasePoint, Vec4};
use nullgeo::integrator::{hamiltonian, rk45_step, Tolerances};
use nullgeo::spacetimes::ellis::Ellis;
use nullgeo::{trace, Camera, CameraPose, CameraSpec, Chart, SkySide, Termination, TraceConfig};

fn sample_points() -> Vec<Vec4> {
    vec![
        Vec4::new(0.0, 5.0, 1.0, 0.3),
        Vec4::new(1.0, -3.0, 2.0, -1.0),
        Vec4::new(-2.0, 0.0, 0.7, 2.0),
        Vec4::new(0.5, 12.0, 2.5, 4.0),
        Vec4::new(0.0, -0.4, 1.2, 0.0),
    ]
}

fn assert_mat_close(a: &Mat4, b: &Mat4, tol: f64, ctx: &str) {
    for i in 0..4 {
        for j in 0..4 {
            let d = (a[(i, j)] - b[(i, j)]).abs();
            assert!(
                d <= tol,
                "{ctx}: ({i},{j}) differs by {d:.3e} (a={}, b={})",
                a[(i, j)],
                b[(i, j)]
            );
        }
    }
}

fn equatorial_ray(ellis: &Ellis, l0: f64, b: f64) -> PhasePoint {
    let rho2 = ellis.throat_radius().powi(2) + l0 * l0;
    PhasePoint {
        x: Vec4::new(0.0, l0, FRAC_PI_2, 0.0),
        p: Vec4::new(1.0, -(1.0 - b * b / rho2).sqrt(), 0.0, b),
    }
}

#[test]
fn metric_times_inverse_is_identity() {
    let ellis = Ellis::new(1.0).unwrap();
    for x in sample_points() {
        let prod = ellis.g(&x) * ellis.g_inv(&x);
        assert_mat_close(
            &prod,
            &Mat4::identity(),
            1e-12,
            &format!("ellis g*g_inv at {x:?}"),
        );
    }
}

#[test]
fn fd_matches_analytic_dg_inv() {
    let ellis = Ellis::new(1.0).unwrap();
    for x in sample_points() {
        let analytic = ellis.dg_inv(&x);
        let fd = fd_partials(|y| ellis.g_inv(y), &x);
        let scale = analytic.iter().map(|m| m.amax()).fold(1.0_f64, f64::max);
        for (mu, fd_mu) in fd.iter().enumerate() {
            assert_mat_close(
                fd_mu,
                &analytic[mu],
                1e-6 * scale,
                &format!("ellis dg_inv[{mu}] at {x:?}"),
            );
        }
    }
}

#[test]
fn ellis_rejects_nonpositive_throat() {
    assert!(Ellis::new(0.0).is_err());
    assert!(Ellis::new(-1.0).is_err());
    assert!(Ellis::new(f64::INFINITY).is_err());
}

#[test]
fn alignment_maps_state_to_equator_and_rotates_directions_back() {
    let ellis = Ellis::new(1.0).unwrap();
    let states = [
        PhasePoint {
            x: Vec4::new(0.0, 8.0, 1.1, 0.6),
            p: Vec4::new(1.3, 0.4, 2.0, -1.5),
        },
        PhasePoint {
            x: Vec4::new(2.0, -5.0, 2.3, -1.9),
            p: Vec4::new(0.8, -0.9, -3.1, 0.7),
        },
        PhasePoint {
            x: Vec4::new(0.0, 30.0, 0.4, 1.0),
            p: Vec4::new(1.0, -1.0, 0.0, 0.0),
        },
    ];

    for s0 in states {
        let (s, alignment) = ellis.align_ray(s0);

        assert!((s.x[2] - FRAC_PI_2).abs() < 1e-15);
        assert_eq!(s.x[3], 0.0);
        assert_eq!(s.p[2], 0.0);
        let st = s0.x[2].sin();
        let l_total = (s0.p[2].powi(2) + (s0.p[3] / st).powi(2)).sqrt();
        assert!((s.p[3] - l_total).abs() < 1e-12);

        let h0 = hamiltonian(&ellis, &s0);
        let h = hamiltonian(&ellis, &s);
        assert!((h - h0).abs() < 1e-12, "H changed: {h0} -> {h}");

        let d0 = ellis.embed_direction(&s0.x, &raise(&ellis.g_inv(&s0.x), &s0.p));
        let d = alignment.apply(ellis.embed_direction(&s.x, &raise(&ellis.g_inv(&s.x), &s.p)));
        for i in 0..3 {
            assert!(
                (d[i] - d0[i]).abs() < 1e-12,
                "direction component {i}: {} vs {}",
                d[i],
                d0[i]
            );
        }

        let p0 = ellis.embed(&s0.x);
        let p = alignment.apply(ellis.embed(&s.x));
        for i in 0..3 {
            assert!(
                (p[i] - p0[i]).abs() < 1e-12,
                "position component {i}: {} vs {}",
                p[i],
                p0[i]
            );
        }
    }
}

#[test]
fn polar_ray_integrates_pole_free_and_conserves_energy_and_angular_momentum() {
    let ellis = Ellis::new(1.0).unwrap();
    let l0 = 30.0;
    let rho2 = 1.0 + l0 * l0;
    let s0 = PhasePoint {
        x: Vec4::new(0.0, l0, FRAC_PI_2, 0.0),
        p: Vec4::new(1.0, -(1.0 - 9.0 / rho2).sqrt(), -3.0, 0.0),
    };

    let (mut s, _) = ellis.align_ray(s0);
    let tol = Tolerances {
        rtol: 1e-11,
        atol: 1e-11,
    };
    let mut dl = 0.1_f64;
    let mut lambda = 0.0;
    let mut max_drift = 0.0_f64;
    while lambda < 80.0 {
        let step = rk45_step(&ellis, &s, dl, &tol);
        if step.accepted {
            s = step.state;
            lambda += step.dl_used;
            max_drift = max_drift
                .max((s.p[0] - 1.0).abs())
                .max((s.p[3] - 3.0).abs())
                .max(s.p[2].abs())
                .max((s.x[2] - FRAC_PI_2).abs())
                .max(hamiltonian(&ellis, &s).abs());
        }
        dl = step.dl_next.max(1e-9);
    }
    assert!(
        max_drift < 1e-9,
        "max conserved-quantity drift {max_drift:.3e}"
    );
}

#[test]
fn polar_ray_escapes_in_its_original_plane() {
    let ellis = Ellis::new(1.0).unwrap();
    let l0 = 30.0;
    let rho2 = 1.0 + l0 * l0;
    let s0 = PhasePoint {
        x: Vec4::new(0.0, l0, FRAC_PI_2, 0.0),
        p: Vec4::new(1.0, -(1.0 - 9.0 / rho2).sqrt(), -3.0, 0.0),
    };

    let cfg = TraceConfig {
        escape_radius: 100.0,
        ..TraceConfig::default()
    };
    let Termination::Escaped { side, dir, .. } = trace(&ellis, s0, &cfg) else {
        panic!("polar ray with b = 3 should escape");
    };
    assert_eq!(side, SkySide::Primary);
    assert!(dir[1].abs() < 1e-9, "escape left the xz-plane: {dir:?}");
    assert!(dir[0] < 0.0 && dir[2] > 0.0, "unexpected escape {dir:?}");
}

#[test]
fn small_impact_parameter_ray_traverses_the_throat() {
    let ellis = Ellis::new(1.0).unwrap();
    let ray = equatorial_ray(&ellis, 30.0, 0.3);
    let cfg = TraceConfig {
        escape_radius: 100.0,
        ..TraceConfig::default()
    };
    let Termination::Escaped { side, state, .. } = trace(&ellis, ray, &cfg) else {
        panic!("ray aimed at the throat should traverse and escape");
    };
    assert_eq!(side, SkySide::Secondary);
    assert!(state.x[1] < -100.0);
}

#[test]
fn large_impact_parameter_ray_stays_on_primary_side() {
    let ellis = Ellis::new(1.0).unwrap();
    let ray = equatorial_ray(&ellis, 30.0, 5.0);
    let cfg = TraceConfig {
        escape_radius: 100.0,
        ..TraceConfig::default()
    };
    let Termination::Escaped { side, state, .. } = trace(&ellis, ray, &cfg) else {
        panic!("ray with b = 5 b0 should escape");
    };
    assert_eq!(side, SkySide::Primary);
    assert!(state.x[1] > 100.0);
}

#[test]
fn weak_deflection_matches_quadratic_formula() {
    let ellis = Ellis::new(1.0).unwrap();
    let b = 20.0;
    let ray = equatorial_ray(&ellis, 2000.0, b);

    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-12,
            atol: 1e-12,
        },
        escape_radius: 4000.0,
        dl_max: 100.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    };
    let d0 = ellis.embed_direction(&ray.x, &raise(&ellis.g_inv(&ray.x), &ray.p));
    let Termination::Escaped { side, dir, .. } = trace(&ellis, ray, &cfg) else {
        panic!("weak-field ray should escape");
    };
    assert_eq!(side, SkySide::Primary);

    let dot = (d0[0] * dir[0] + d0[1] * dir[1] + d0[2] * dir[2]).clamp(-1.0, 1.0);
    let deflection = dot.acos();
    let predicted = (PI / 4.0) * (1.0 / b).powi(2);
    assert!(
        (deflection - predicted).abs() < 0.02 * predicted,
        "deflection {deflection:.6e} vs predicted {predicted:.6e}"
    );
}

#[test]
fn camera_rays_are_null_and_see_both_sides() {
    let ellis = Ellis::new(1.0).unwrap();
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 30.0,
            res: (3, 3),
            energy: 1.0,
            supersample: 1,
        },
        CameraPose {
            position: Vec4::new(0.0, 20.0, FRAC_PI_2, 0.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap();

    let rays = camera.pixel_rays(&ellis).unwrap();
    assert_eq!(rays.len(), 9);
    for ray in &rays {
        let h = hamiltonian(&ellis, ray);
        assert!(h.abs() < 1e-12, "H = {h:.3e}");
        let v = ellis.g_inv(&ray.x) * ray.p;
        assert!(v[0] < 0.0, "dt/dlambda = {}", v[0]);
    }

    let cfg = TraceConfig {
        escape_radius: 100.0,
        ..TraceConfig::default()
    };
    let Termination::Escaped { side, .. } = trace(&ellis, rays[4], &cfg) else {
        panic!("center pixel should traverse the throat");
    };
    assert_eq!(side, SkySide::Secondary);

    let Termination::Escaped { side, .. } = trace(&ellis, rays[0], &cfg) else {
        panic!("corner pixel should escape");
    };
    assert_eq!(side, SkySide::Primary);
}

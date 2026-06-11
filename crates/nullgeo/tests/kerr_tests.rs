use nullgeo::frame::{build_coframe, make_null_covector};
use nullgeo::integrator::{hamiltonian, rk45_step, Tolerances};
use nullgeo::metric::{Mat4, Metric, State4, Vec4};
use nullgeo::metrics::kerr::Kerr;
use nullgeo::metrics::schwarzschild::Schwarzschild;
use nullgeo::{trace, Spacetime, Termination, TraceConfig};

fn sample_points() -> Vec<Vec4> {
    vec![
        Vec4::new(0.0, 10.0, 0.0, 0.0),
        Vec4::new(1.0, 3.0, 4.0, -2.0),
        Vec4::new(-2.0, -50.0, 30.0, 12.0),
        Vec4::new(0.0, 0.1, 0.05, 3.0),
        Vec4::new(0.0, 2.0, -1.0, 0.5),
        Vec4::new(0.0, 1.5, 0.3, 0.0),
        Vec4::new(5.0, 0.5, -0.2, 0.1),
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

fn backward_ray<M: Metric>(m: &M, x: Vec4, dir: [f64; 3], energy: f64) -> State4 {
    let coframe = build_coframe(m, &x).unwrap();
    let arriving = make_null_covector(&coframe, [-dir[0], -dir[1], -dir[2]], energy);
    State4 { x, p: -arriving }
}

fn xi(s: &State4) -> f64 {
    -(s.x[1] * s.p[2] - s.x[2] * s.p[1]) / s.p[0]
}

fn equatorial_ray_with_xi(kerr: &Kerr, r0: f64, xi_target: f64) -> State4 {
    let x0 = Vec4::new(0.0, -r0, 0.0, 0.0);
    let ray_at = |alpha: f64| backward_ray(kerr, x0, [alpha.cos(), alpha.sin(), 0.0], 1.0);
    let sense = if xi(&ray_at(0.1)) * xi_target > 0.0 {
        1.0
    } else {
        -1.0
    };
    let (mut lo, mut hi) = (1e-4_f64, 1.4_f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if xi(&ray_at(sense * mid)).abs() < xi_target.abs() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let ray = ray_at(sense * 0.5 * (lo + hi));
    assert!(
        (xi(&ray) - xi_target).abs() < 1e-4,
        "aim bisection failed: xi = {}, target = {xi_target}",
        xi(&ray)
    );
    ray
}

fn photon_orbit_radius(m: f64, a: f64, prograde: bool) -> f64 {
    let s = if prograde { -1.0 } else { 1.0 };
    2.0 * m * (1.0 + ((2.0 / 3.0) * (s * a / m).acos()).cos())
}

fn critical_xi(m: f64, a: f64, r_ph: f64) -> f64 {
    -(r_ph.powi(3) - 3.0 * m * r_ph * r_ph + a * a * r_ph + a * a * m) / (a * (r_ph - m))
}

#[test]
fn bardeen_formulas_match_extremal_anchor() {
    let r_ph = photon_orbit_radius(1.0, 1.0, false);
    assert!((r_ph - 4.0).abs() < 1e-12);
    assert!((critical_xi(1.0, 1.0, r_ph) + 7.0).abs() < 1e-12);
    assert!((photon_orbit_radius(1.0, 1.0, true) - 1.0).abs() < 1e-7);
    assert!((photon_orbit_radius(1.0, 0.0, true) - 3.0).abs() < 1e-12);
    assert!((photon_orbit_radius(1.0, 0.0, false) - 3.0).abs() < 1e-12);
}

#[test]
fn kerr_schild_radius_satisfies_quartic() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let a = 0.9_f64;
    for x in sample_points() {
        let r = kerr.radius(&x);
        let rr = x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
        let residual = r.powi(4) - (rr - a * a) * r * r - a * a * x[3] * x[3];
        let scale = r.powi(4) + (rr - a * a).abs() * r * r + a * a * x[3] * x[3] + 1e-30;
        assert!(
            residual.abs() <= 1e-11 * scale,
            "quartic residual {residual:.3e} at {x:?}"
        );
    }
}

#[test]
fn kerr_with_zero_spin_matches_schwarzschild() {
    let kerr = Kerr::new(1.0, 0.0).unwrap();
    let schw = Schwarzschild::new(1.0).unwrap();
    for x in sample_points() {
        assert_mat_close(&kerr.g(&x), &schw.g(&x), 1e-12, &format!("g at {x:?}"));
        assert_mat_close(
            &kerr.g_inv(&x),
            &schw.g_inv(&x),
            1e-12,
            &format!("g_inv at {x:?}"),
        );
    }
}

#[test]
fn kerr_metric_times_inverse_is_identity() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let mut points = sample_points();
    points.push(Vec4::new(0.0, 1.8, 0.3, 0.05));
    points.push(Vec4::new(0.0, -1.9, 0.2, -0.1));
    for x in points {
        let prod = kerr.g(&x) * kerr.g_inv(&x);
        assert_mat_close(
            &prod,
            &Mat4::identity(),
            1e-12,
            &format!("kerr g*g_inv at {x:?} (r = {})", kerr.radius(&x)),
        );
    }
}

#[test]
fn kerr_inside_ergosphere_points_are_inside_ergosphere() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    for x in [
        Vec4::new(0.0, 1.8, 0.3, 0.05),
        Vec4::new(0.0, -1.9, 0.2, -0.1),
    ] {
        let r = kerr.radius(&x);
        assert!(r > kerr.outer_horizon(), "point at r = {r} inside horizon");
        assert!(
            kerr.g(&x)[(0, 0)] > 0.0,
            "point at r = {r} outside ergosphere (g_tt = {})",
            kerr.g(&x)[(0, 0)]
        );
    }
}

#[test]
fn kerr_rejects_overspun_black_hole() {
    assert!(Kerr::new(1.0, 1.1).is_err());
    assert!(Kerr::new(1.0, -1.1).is_err());
    assert!(Kerr::new(1.0, 1.0).is_ok());
}

#[test]
fn kerr_geodesic_conserves_energy_angular_momentum_and_hamiltonian() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let s0 = equatorial_ray_with_xi(&kerr, 15.0, 3.0);
    let energy = |s: &State4| s.p[0];
    let l_z = |s: &State4| s.x[1] * s.p[2] - s.x[2] * s.p[1];

    let tol = Tolerances {
        rtol: 1e-11,
        atol: 1e-11,
    };
    let mut s = s0;
    let mut dl = 0.1_f64;
    let mut lambda = 0.0;
    let mut max_drift = 0.0_f64;
    while lambda < 80.0 {
        let step = rk45_step(&kerr, &s, dl, &tol);
        if step.accepted {
            s = step.state;
            lambda += step.dl_used;
            max_drift = max_drift
                .max((energy(&s) - energy(&s0)).abs())
                .max((l_z(&s) - l_z(&s0)).abs())
                .max(hamiltonian(&kerr, &s).abs());
        }
        dl = step.dl_next.max(1e-9);
    }
    assert!(
        kerr.radius(&s.x) > kerr.outer_horizon() * 1.01,
        "ray fell in: r = {}",
        kerr.radius(&s.x)
    );
    assert!(
        max_drift < 1e-9,
        "max conserved-quantity drift {max_drift:.3e}"
    );
}

#[test]
fn equatorial_photon_orbits_separate_capture_from_escape() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-11,
            atol: 1e-11,
        },
        escape_radius: 100.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    };

    for prograde in [true, false] {
        let r_ph = photon_orbit_radius(1.0, 0.9, prograde);
        let xi_c = critical_xi(1.0, 0.9, r_ph);
        for (factor, expect_captured) in [(1.0 - 1e-3, true), (1.0 + 1e-3, false)] {
            let ray = equatorial_ray_with_xi(&kerr, 30.0, xi_c * factor);
            match trace(&kerr, ray, &cfg) {
                Termination::Captured { .. } => assert!(
                    expect_captured,
                    "prograde={prograde}: xi = xi_c * {factor} should escape"
                ),
                Termination::Escaped { .. } => assert!(
                    !expect_captured,
                    "prograde={prograde}: xi = xi_c * {factor} should be captured"
                ),
                other => {
                    panic!("unexpected termination {other:?} (prograde={prograde}, {factor})")
                }
            }
        }
    }
}

#[test]
fn prograde_rays_deflect_less_than_retrograde() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-11,
            atol: 1e-11,
        },
        escape_radius: 200.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    };

    let deflection = |xi_target: f64| -> f64 {
        let ray = equatorial_ray_with_xi(&kerr, 30.0, xi_target);
        let d0 = kerr.cartesian_direction(&ray);
        let Termination::Escaped { dir, .. } = trace(&kerr, ray, &cfg) else {
            panic!("ray with xi = {xi_target} should escape");
        };
        let dot = d0[0] * dir[0] + d0[1] * dir[1] + d0[2] * dir[2];
        dot.clamp(-1.0, 1.0).acos()
    };

    let prograde = deflection(10.0);
    let retrograde = deflection(-10.0);
    assert!(
        prograde + 0.05 < retrograde,
        "frame dragging sign: prograde deflection {prograde:.4} vs retrograde {retrograde:.4}"
    );
}

use std::f64::consts::PI;

use nullgeo::geometry::{build_coframe, make_null_covector};
use nullgeo::geometry::{Metric, PhasePoint, Vec4};
use nullgeo::integrator::{rk45_step, Tolerances};
use nullgeo::spacetimes::minkowski::Minkowski;
use nullgeo::spacetimes::schwarzschild::Schwarzschild;
use nullgeo::{trace, Termination, TraceConfig};

fn backward_ray<M: Metric>(m: &M, x: Vec4, dir: [f64; 3], energy: f64) -> PhasePoint {
    let coframe = build_coframe(m, &x).unwrap();
    let arriving = make_null_covector(&coframe, [-dir[0], -dir[1], -dir[2]], energy);
    PhasePoint { x, p: -arriving }
}

fn static_frame_direction<M: Metric>(m: &M, s: &PhasePoint) -> [f64; 3] {
    let coframe = build_coframe(m, &s.x).unwrap();
    let v = m.g_inv(&s.x) * s.p;
    let speed = (-coframe[0].dot(&v)).abs();
    [
        coframe[1].dot(&v) / speed,
        coframe[2].dot(&v) / speed,
        coframe[3].dot(&v) / speed,
    ]
}

fn aimed_ray(m: &Schwarzschild, r0: f64, b: f64) -> PhasePoint {
    let alpha = (b * (1.0 - 2.0 * m.mass() / r0).sqrt() / r0).asin();
    backward_ray(
        m,
        Vec4::new(0.0, -r0, 0.0, 0.0),
        [alpha.cos(), alpha.sin(), 0.0],
        1.0,
    )
}

#[test]
fn critical_impact_parameter_separates_capture_from_escape() {
    let m = Schwarzschild::new(1.0).unwrap();
    let b_crit = 27.0_f64.sqrt();
    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-11,
            atol: 1e-11,
        },
        escape_radius: 100.0,
        max_steps: 500_000,
        ..TraceConfig::default()
    };

    for (factor, expect_captured) in [(1.0 - 1e-3, true), (1.0 + 1e-3, false)] {
        let ray = aimed_ray(&m, 30.0, b_crit * factor);
        match trace(&m, ray, &cfg) {
            Termination::Captured { .. } => {
                assert!(expect_captured, "b = b_c * {factor} should escape")
            }
            Termination::Escaped { .. } => {
                assert!(!expect_captured, "b = b_c * {factor} should be captured")
            }
            other => panic!("unexpected termination {other:?} at b = b_c * {factor}"),
        }
    }
}

#[test]
fn weak_field_deflection_matches_second_order_formula() {
    let m = Schwarzschild::new(1.0).unwrap();
    let x0 = Vec4::new(0.0, -2000.0, 50.0, 0.0);
    let ray = backward_ray(&m, x0, [1.0, 0.0, 0.0], 1.0);

    let b = ((x0[1] * ray.p[2] - x0[2] * ray.p[1]) / ray.p[0]).abs();

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
    let Termination::Escaped { state, .. } = trace(&m, ray, &cfg) else {
        panic!("ray with b = {b} should escape");
    };

    let deflection = static_frame_direction(&m, &state)[0].acos();
    let predicted = 4.0 / b + 15.0 * PI / (4.0 * b * b);
    assert!(
        (deflection - predicted).abs() < 0.005 * predicted,
        "deflection {deflection:.6e} vs predicted {predicted:.6e} at b = {b}"
    );
}

#[test]
fn energy_and_angular_momentum_conserved_along_bent_ray() {
    let m = Schwarzschild::new(1.0).unwrap();
    let s0 = aimed_ray(&m, 15.0, 5.25);
    let energy = |s: &PhasePoint| s.p[0];
    let l_z = |s: &PhasePoint| s.x[1] * s.p[2] - s.x[2] * s.p[1];

    let tol = Tolerances {
        rtol: 1e-11,
        atol: 1e-11,
    };
    let mut s = s0;
    let mut dl = 0.1_f64;
    let mut lambda = 0.0;
    let mut max_drift = 0.0_f64;
    while lambda < 60.0 {
        let step = rk45_step(&m, &s, dl, &tol);
        if step.accepted {
            s = step.state;
            lambda += step.dl_used;
            max_drift = max_drift
                .max((energy(&s) - energy(&s0)).abs())
                .max((l_z(&s) - l_z(&s0)).abs());
        }
        dl = step.dl_next.max(1e-9);
    }
    assert!(
        max_drift < 1e-9,
        "max conserved-quantity drift {max_drift:.3e}"
    );
}

#[test]
fn minkowski_rays_travel_straight() {
    let m = Minkowski;
    let dir = [2.0 / 7.0, 3.0 / 7.0, 6.0 / 7.0];
    let x0 = Vec4::new(0.0, 1.0, -2.0, 0.5);
    let ray = backward_ray(&m, x0, dir, 1.0);

    let cfg = TraceConfig {
        escape_radius: 100.0,
        ..TraceConfig::default()
    };
    let Termination::Escaped {
        dir: out_dir,
        state,
        ..
    } = trace(&m, ray, &cfg)
    else {
        panic!("minkowski ray must escape");
    };

    for i in 0..3 {
        assert!(
            (out_dir[i] - dir[i]).abs() < 1e-12,
            "direction component {i} drifted: {} vs {}",
            out_dir[i],
            dir[i]
        );
    }

    let displacement = [state.x[1] - x0[1], state.x[2] - x0[2], state.x[3] - x0[3]];
    let cross = [
        displacement[1] * dir[2] - displacement[2] * dir[1],
        displacement[2] * dir[0] - displacement[0] * dir[2],
        displacement[0] * dir[1] - displacement[1] * dir[0],
    ];
    for (i, c) in cross.iter().enumerate() {
        assert!(c.abs() < 1e-9, "path not straight: cross[{i}] = {c:.3e}");
    }
}

#[test]
fn near_critical_ray_escapes_after_orbiting() {
    let m = Schwarzschild::new(1.0).unwrap();
    let b_crit = 27.0_f64.sqrt();
    let ray = aimed_ray(&m, 30.0, b_crit * (1.0 + 1e-4));

    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-12,
            atol: 1e-12,
        },
        escape_radius: 100.0,
        max_steps: 2_000_000,
        ..TraceConfig::default()
    };
    let Termination::Escaped { state, .. } = trace(&m, ray, &cfg) else {
        panic!("ray just outside critical b must eventually escape");
    };

    let out = static_frame_direction(&m, &state);
    let turning = ray.p[1].signum() != state.p[1].signum() || out[0] < 0.0;
    assert!(turning, "near-critical ray should be strongly deflected");
}

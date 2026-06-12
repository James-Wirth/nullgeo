use std::f64::consts::FRAC_PI_2;

use nullgeo::geometry::{
    build_coframe, fd_partials, make_null_covector, Mat4, Metric, PhasePoint, Vec4,
};
use nullgeo::spacetimes::kerr::Kerr;
use nullgeo::spacetimes::schwarzschild::Schwarzschild;
use nullgeo::spacetimes::{Pullback, Rotation, Transition};
use nullgeo::{trace, Chart, EquatorialAnnulus, Spacetime, Termination, TraceConfig};

fn sample_points() -> Vec<Vec4> {
    vec![
        Vec4::new(0.0, 10.0, 0.0, 0.0),
        Vec4::new(1.0, 3.0, 4.0, -2.0),
        Vec4::new(0.0, -8.0, 6.0, 3.0),
        Vec4::new(2.0, 5.0, -7.0, 4.0),
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

fn tilted_kerr() -> Pullback<Kerr, Rotation> {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let rotation = Rotation::about_axis([1.0, 2.0, 0.5], 0.83).unwrap();
    Pullback::new(kerr, rotation)
}

fn backward_ray<M: Metric>(m: &M, x: Vec4, dir: [f64; 3], energy: f64) -> PhasePoint {
    let coframe = build_coframe(m, &x).unwrap();
    let arriving = make_null_covector(&coframe, [-dir[0], -dir[1], -dir[2]], energy);
    PhasePoint { x, p: -arriving }
}

fn rotate_state(rotation: &Rotation, s: &PhasePoint) -> PhasePoint {
    let to_new = rotation.inv_jacobian(&s.x);
    PhasePoint {
        x: to_new * s.x,
        p: to_new * s.p,
    }
}

#[test]
fn rotation_rejects_degenerate_input() {
    assert!(Rotation::about_axis([0.0, 0.0, 0.0], 1.0).is_err());
    assert!(Rotation::about_axis([1.0, 0.0, 0.0], f64::NAN).is_err());
}

#[test]
fn pullback_metric_inverse_and_derivative_are_consistent() {
    let tilted = tilted_kerr();
    for y in sample_points() {
        let product = tilted.g(&y) * tilted.g_inv(&y);
        assert_mat_close(
            &product,
            &Mat4::identity(),
            1e-12,
            &format!("tilted kerr g*g_inv at {y:?}"),
        );

        let analytic = tilted.dg_inv(&y);
        let fd = fd_partials(|z| tilted.g_inv(z), &y);
        let scale = analytic.iter().map(|m| m.amax()).fold(1.0_f64, f64::max);
        for (mu, fd_mu) in fd.iter().enumerate() {
            assert_mat_close(
                fd_mu,
                &analytic[mu],
                1e-6 * scale,
                &format!("tilted kerr dg_inv[{mu}] at {y:?}"),
            );
        }
    }
}

struct Spherical;

impl Transition for Spherical {
    fn map(&self, y: &Vec4) -> Vec4 {
        let (r, theta, phi) = (y[1], y[2], y[3]);
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        Vec4::new(y[0], r * st * cp, r * st * sp, r * ct)
    }

    fn jacobian(&self, y: &Vec4) -> Mat4 {
        let (r, theta, phi) = (y[1], y[2], y[3]);
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let mut j = Mat4::identity();
        j[(1, 1)] = st * cp;
        j[(1, 2)] = r * ct * cp;
        j[(1, 3)] = -r * st * sp;
        j[(2, 1)] = st * sp;
        j[(2, 2)] = r * ct * sp;
        j[(2, 3)] = r * st * cp;
        j[(3, 1)] = ct;
        j[(3, 2)] = -r * st;
        j[(3, 3)] = 0.0;
        j
    }

    fn inv_jacobian(&self, y: &Vec4) -> Mat4 {
        let (r, theta, phi) = (y[1], y[2], y[3]);
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let mut a = Mat4::identity();
        a[(1, 1)] = st * cp;
        a[(1, 2)] = st * sp;
        a[(1, 3)] = ct;
        a[(2, 1)] = ct * cp / r;
        a[(2, 2)] = ct * sp / r;
        a[(2, 3)] = -st / r;
        a[(3, 1)] = -sp / (r * st);
        a[(3, 2)] = cp / (r * st);
        a[(3, 3)] = 0.0;
        a
    }
}

#[test]
fn schwarzschild_pulled_back_to_spherical_matches_closed_form() {
    let spherical = Pullback::new(Schwarzschild::new(1.0).unwrap(), Spherical);
    for (r, theta, phi) in [(8.0, 1.1, 0.4), (30.0, FRAC_PI_2, -1.3), (3.0, 2.0, 2.9)] {
        let y = Vec4::new(0.0, r, theta, phi);
        let g = spherical.g(&y);

        let mut expected = Mat4::zeros();
        expected[(0, 0)] = -(1.0 - 2.0 / r);
        expected[(0, 1)] = 2.0 / r;
        expected[(1, 0)] = 2.0 / r;
        expected[(1, 1)] = 1.0 + 2.0 / r;
        expected[(2, 2)] = r * r;
        expected[(3, 3)] = r * r * theta.sin().powi(2);

        assert_mat_close(
            &g,
            &expected,
            1e-12 * r * r,
            &format!("ingoing KS spherical form at r = {r}"),
        );
    }
}

#[test]
fn rotated_chart_leaves_geodesics_invariant() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let rotation = Rotation::about_axis([0.3, -1.0, 2.0], 1.21).unwrap();
    let tilted = Pullback::new(kerr, rotation);

    let cfg = TraceConfig {
        escape_radius: 200.0,
        ..TraceConfig::default()
    };

    for b in [7.0, -7.0, 12.0] {
        let ray = backward_ray(&kerr, Vec4::new(0.0, -30.0, b, 0.0), [1.0, 0.0, 0.0], 1.0);
        let Termination::Escaped { dir: d_base, .. } = trace(&kerr, ray, &cfg) else {
            panic!("base ray with b = {b} should escape");
        };
        let Termination::Escaped { dir: d_tilted, .. } =
            trace(&tilted, rotate_state(&rotation, &ray), &cfg)
        else {
            panic!("rotated-chart ray with b = {b} should escape");
        };
        for i in 0..3 {
            assert!(
                (d_base[i] - d_tilted[i]).abs() < 1e-6,
                "escape direction component {i} differs: {} vs {} at b = {b}",
                d_base[i],
                d_tilted[i]
            );
        }
    }

    let plunging = backward_ray(&kerr, Vec4::new(0.0, -30.0, 0.5, 0.0), [1.0, 0.0, 0.0], 1.0);
    assert!(matches!(
        trace(&kerr, plunging, &cfg),
        Termination::Captured { .. }
    ));
    assert!(matches!(
        trace(&tilted, rotate_state(&rotation, &plunging), &cfg),
        Termination::Captured { .. }
    ));
}

#[test]
fn rotated_chart_leaves_disk_hits_and_redshift_invariant() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let rotation = Rotation::about_axis([1.0, 1.0, 1.0], -0.62).unwrap();
    let tilted = Pullback::new(kerr, rotation);

    let cfg = TraceConfig {
        escape_radius: 200.0,
        disk: Some(EquatorialAnnulus {
            r_in: 6.0,
            r_out: 20.0,
        }),
        ..TraceConfig::default()
    };

    let start = Vec4::new(0.0, -30.0, 0.0, 5.0);
    let aim = [20.0_f64, 0.0, -5.0];
    let norm = (aim[0] * aim[0] + aim[2] * aim[2]).sqrt();
    let ray = backward_ray(&kerr, start, [aim[0] / norm, 0.0, aim[2] / norm], 1.0);

    let Termination::HitSurface { state: hit_base } = trace(&kerr, ray, &cfg) else {
        panic!("base ray should hit the disk");
    };
    let Termination::HitSurface { state: hit_tilted } =
        trace(&tilted, rotate_state(&rotation, &ray), &cfg)
    else {
        panic!("rotated-chart ray should hit the disk");
    };

    let r_base = kerr.radius(&hit_base.x);
    let r_tilted = tilted.radius(&hit_tilted.x);
    assert!(
        (r_base - r_tilted).abs() < 1e-6,
        "hit radius differs: {r_base} vs {r_tilted}"
    );

    let u_base = kerr
        .circular_orbits()
        .unwrap()
        .four_velocity(&hit_base.x)
        .unwrap();
    let u_tilted = tilted
        .circular_orbits()
        .unwrap()
        .four_velocity(&hit_tilted.x)
        .unwrap();
    let pu_base = hit_base.p.dot(&u_base);
    let pu_tilted = hit_tilted.p.dot(&u_tilted);
    assert!(
        (pu_base - pu_tilted).abs() < 1e-6 * pu_base.abs(),
        "emitter projection differs: {pu_base} vs {pu_tilted}"
    );
}

#[test]
fn spherical_chart_traces_match_cartesian_traces() {
    let schwarzschild = Schwarzschild::new(1.0).unwrap();
    let spherical = Pullback::new(schwarzschild, Spherical);

    let cfg = TraceConfig {
        escape_radius: 200.0,
        ..TraceConfig::default()
    };

    let x0 = Vec4::new(0.0, -30.0, 8.0, 0.0);
    let ray = backward_ray(&schwarzschild, x0, [1.0, 0.0, 0.0], 1.0);

    let r0 = (x0[1] * x0[1] + x0[2] * x0[2]).sqrt();
    let y0 = Vec4::new(0.0, r0, FRAC_PI_2, x0[2].atan2(x0[1]));
    let pulled = PhasePoint {
        x: y0,
        p: Spherical.jacobian(&y0).transpose() * ray.p,
    };

    let Termination::Escaped { dir: d_base, .. } = trace(&schwarzschild, ray, &cfg) else {
        panic!("Cartesian-chart ray should escape");
    };
    let Termination::Escaped {
        dir: d_spherical, ..
    } = trace(&spherical, pulled, &cfg)
    else {
        panic!("spherical-chart ray should escape");
    };
    for i in 0..3 {
        assert!(
            (d_base[i] - d_spherical[i]).abs() < 1e-5,
            "escape direction component {i} differs: {} vs {}",
            d_base[i],
            d_spherical[i]
        );
    }
}

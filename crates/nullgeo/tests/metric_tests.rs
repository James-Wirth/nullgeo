use nullgeo::geometry::{fd_partials, Mat4, Metric, Vec4};
use nullgeo::spacetimes::kerr::Kerr;
use nullgeo::spacetimes::minkowski::Minkowski;
use nullgeo::spacetimes::reissner_nordstrom::ReissnerNordstrom;
use nullgeo::spacetimes::schwarzschild::Schwarzschild;

fn sample_points() -> Vec<Vec4> {
    vec![
        Vec4::new(0.0, 10.0, 0.0, 0.0),
        Vec4::new(1.0, 3.0, 4.0, -2.0),
        Vec4::new(-2.0, -50.0, 30.0, 12.0),
        Vec4::new(0.0, 1.0, 0.5, -0.3),
        Vec4::new(5.0, 0.7, -0.9, 0.4),
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

#[test]
fn metric_times_inverse_is_identity() {
    let mink = Minkowski;
    let schw = Schwarzschild::new(1.0).unwrap();
    for x in sample_points() {
        let prod_m = mink.g(&x) * mink.g_inv(&x);
        assert_mat_close(&prod_m, &Mat4::identity(), 1e-12, "minkowski g*g_inv");

        let prod_s = schw.g(&x) * schw.g_inv(&x);
        assert_mat_close(
            &prod_s,
            &Mat4::identity(),
            1e-12,
            &format!("schwarzschild g*g_inv at {x:?}"),
        );
    }
}

#[test]
fn reissner_nordstrom_metric_times_inverse_is_identity() {
    let rn = ReissnerNordstrom::new(1.0, 0.8).unwrap();
    for x in sample_points() {
        let prod = rn.g(&x) * rn.g_inv(&x);
        assert_mat_close(
            &prod,
            &Mat4::identity(),
            1e-12,
            &format!("reissner-nordstrom g*g_inv at {x:?}"),
        );
    }
}

#[test]
fn reissner_nordstrom_fd_matches_analytic_dg_inv() {
    let rn = ReissnerNordstrom::new(1.0, 0.8).unwrap();
    for x in sample_points() {
        let analytic = rn.dg_inv(&x);
        let fd = fd_partials(|y| rn.g_inv(y), &x);
        let scale = analytic.iter().map(|m| m.amax()).fold(1.0_f64, f64::max);
        for (mu, fd_mu) in fd.iter().enumerate() {
            assert_mat_close(
                fd_mu,
                &analytic[mu],
                1e-7 * scale,
                &format!("rn dg_inv[{mu}] at {x:?}"),
            );
        }
    }
}

#[test]
fn reissner_nordstrom_with_zero_charge_matches_schwarzschild() {
    let rn = ReissnerNordstrom::new(1.0, 0.0).unwrap();
    let schw = Schwarzschild::new(1.0).unwrap();
    for x in sample_points() {
        assert_mat_close(&rn.g(&x), &schw.g(&x), 1e-12, &format!("rn g at {x:?}"));
        assert_mat_close(
            &rn.g_inv(&x),
            &schw.g_inv(&x),
            1e-12,
            &format!("rn g_inv at {x:?}"),
        );
    }
}

#[test]
fn reissner_nordstrom_rejects_overcharged_black_hole() {
    assert!(ReissnerNordstrom::new(1.0, 1.2).is_err());
    assert!(ReissnerNordstrom::new(1.0, -1.2).is_err());
    assert!(ReissnerNordstrom::new(1.0, 1.0).is_ok());
}

#[test]
fn schwarzschild_fd_matches_analytic_dg_inv() {
    let schw = Schwarzschild::new(1.0).unwrap();
    for x in sample_points() {
        let analytic = schw.dg_inv(&x);
        let fd = fd_partials(|y| schw.g_inv(y), &x);
        let scale = analytic.iter().map(|m| m.amax()).fold(1.0_f64, f64::max);
        for (mu, fd_mu) in fd.iter().enumerate() {
            assert_mat_close(
                fd_mu,
                &analytic[mu],
                1e-7 * scale,
                &format!("dg_inv[{mu}] at {x:?}"),
            );
        }
    }
}

#[test]
fn kerr_fd_matches_analytic_dg_inv() {
    for spin in [0.9, -0.6] {
        let kerr = Kerr::new(1.0, spin).unwrap();
        for x in sample_points() {
            let analytic = kerr.dg_inv(&x);
            let fd = fd_partials(|y| kerr.g_inv(y), &x);
            let scale = analytic.iter().map(|m| m.amax()).fold(1.0_f64, f64::max);
            for (mu, fd_mu) in fd.iter().enumerate() {
                assert_mat_close(
                    fd_mu,
                    &analytic[mu],
                    1e-7 * scale,
                    &format!("kerr(a={spin}) dg_inv[{mu}] at {x:?}"),
                );
            }
        }
    }
}

#[test]
fn kerr_dg_inv_with_zero_spin_matches_schwarzschild() {
    let kerr = Kerr::new(1.0, 0.0).unwrap();
    let schw = Schwarzschild::new(1.0).unwrap();
    for x in sample_points() {
        let dk = kerr.dg_inv(&x);
        let ds = schw.dg_inv(&x);
        for mu in 0..4 {
            assert_mat_close(
                &dk[mu],
                &ds[mu],
                1e-12,
                &format!("kerr a=0 dg_inv[{mu}] at {x:?}"),
            );
        }
    }
}

#[test]
fn minkowski_fd_dg_inv_is_zero() {
    let mink = Minkowski;
    for x in sample_points() {
        let fd = fd_partials(|y| mink.g_inv(y), &x);
        for fd_mu in &fd {
            assert_mat_close(fd_mu, &Mat4::zeros(), 1e-14, "minkowski fd");
        }
    }
}

#[test]
fn schwarzschild_is_ingoing_kerr_schild() {
    let m = 1.0;
    let schw = Schwarzschild::new(m).unwrap();
    let x = Vec4::new(0.0, 4.0, 0.0, 0.0);
    let g = schw.g(&x);
    let h = m / 4.0;

    assert!((g[(0, 0)] - (-(1.0 - 2.0 * h))).abs() < 1e-14, "g_tt");
    assert!((g[(0, 1)] - 2.0 * h).abs() < 1e-14, "g_tx sign (ingoing)");
    assert!((g[(1, 1)] - (1.0 + 2.0 * h)).abs() < 1e-14, "g_xx");
    assert!(
        g[(0, 2)].abs() < 1e-14 && g[(0, 3)].abs() < 1e-14,
        "g_ty, g_tz"
    );
}

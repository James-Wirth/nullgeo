use nullgeo::integrator::{hamiltonian, rk45_step, rk4_step, Tolerances};
use nullgeo::metric::{State4, Vec4};
use nullgeo::metrics::schwarzschild::Schwarzschild;

fn deflected_ray_start() -> State4 {
    State4 {
        x: Vec4::new(0.0, -15.0, 0.0, 0.0),
        p: Vec4::new(
            -9.309_493_362_512_625e-1,
            8.506_570_025_217_01e-1,
            3.793_532_811_354_438e-1,
            0.0,
        ),
    }
}

fn integrate_rk4(m: &Schwarzschild, mut s: State4, dl: f64, total: f64) -> State4 {
    let n = (total / dl).round() as usize;
    for _ in 0..n {
        s = rk4_step(m, &s, dl);
    }
    s
}

fn integrate_rk45_fixed(m: &Schwarzschild, mut s: State4, dl: f64, total: f64) -> State4 {
    let tol = Tolerances {
        rtol: 1e30,
        atol: 1e30,
    };
    let n = (total / dl).round() as usize;
    for _ in 0..n {
        let res = rk45_step(m, &s, dl, &tol);
        assert!(res.accepted);
        s = res.state;
    }
    s
}

fn integrate_adaptive(
    m: &Schwarzschild,
    mut s: State4,
    tol: &Tolerances,
    total: f64,
) -> (State4, usize, usize) {
    let mut dl = 0.1_f64;
    let mut lambda = 0.0;
    let mut accepted = 0;
    let mut rejected = 0;
    while lambda < total {
        let res = rk45_step(m, &s, dl.min(total - lambda), tol);
        if res.accepted {
            s = res.state;
            lambda += res.dl_used;
            accepted += 1;
        } else {
            rejected += 1;
        }
        dl = res.dl_next.max(1e-9);
        assert!(accepted + rejected < 1_000_000);
    }
    (s, accepted, rejected)
}

fn state_distance(a: &State4, b: &State4) -> f64 {
    let dx = a.x - b.x;
    let dp = a.p - b.p;
    (dx.dot(&dx) + dp.dot(&dp)).sqrt()
}

#[test]
fn rk45_is_fifth_order() {
    let m = Schwarzschild::new(1.0).unwrap();
    let s0 = deflected_ray_start();
    let total = 16.0;

    let reference = integrate_rk4(&m, s0, 1e-3, total);
    let err_coarse = state_distance(&integrate_rk45_fixed(&m, s0, 0.4, total), &reference);
    let err_fine = state_distance(&integrate_rk45_fixed(&m, s0, 0.2, total), &reference);

    let order = (err_coarse / err_fine).log2();
    assert!(
        (4.0..6.5).contains(&order),
        "observed order {order}, coarse {err_coarse:.3e}, fine {err_fine:.3e}"
    );
}

#[test]
fn adaptive_matches_fine_reference() {
    let m = Schwarzschild::new(1.0).unwrap();
    let s0 = deflected_ray_start();
    let total = 30.0;

    let reference = integrate_rk4(&m, s0, 1e-3, total);
    let (adaptive, accepted, _) = integrate_adaptive(&m, s0, &Tolerances::default(), total);

    let dist = state_distance(&adaptive, &reference);
    assert!(dist < 1e-6, "distance {dist:.3e} after {accepted} steps");
    assert!(accepted < 5000, "took {accepted} accepted steps");
}

#[test]
fn adaptive_conserves_null_constraint_and_energy() {
    let m = Schwarzschild::new(1.0).unwrap();
    let s0 = State4 {
        x: Vec4::new(0.0, -15.0, 0.0, 0.0),
        p: Vec4::new(
            -9.309493362512625e-1,
            8.723291675492424e-1,
            3.258_322_676_879_42e-1,
            0.0,
        ),
    };
    let h0 = hamiltonian(&m, &s0);
    assert!(h0.abs() < 1e-9, "initial H = {h0:.3e}");

    let tol = Tolerances {
        rtol: 1e-11,
        atol: 1e-11,
    };
    let mut s = s0;
    let mut dl = 0.1;
    let mut lambda = 0.0;
    let mut max_h = 0.0_f64;
    let mut max_pt_drift = 0.0_f64;
    while lambda < 60.0 {
        let res = rk45_step(&m, &s, dl, &tol);
        if res.accepted {
            s = res.state;
            lambda += res.dl_used;
            max_h = max_h.max(hamiltonian(&m, &s).abs());
            max_pt_drift = max_pt_drift.max((s.p[0] - s0.p[0]).abs());
        }
        dl = res.dl_next.max(1e-9);
    }

    assert!(max_h < 1e-9, "max |H| = {max_h:.3e}");
    assert!(max_pt_drift < 1e-9, "max p_t drift = {max_pt_drift:.3e}");
}

#[test]
fn oversized_step_is_rejected_and_shrunk() {
    let m = Schwarzschild::new(1.0).unwrap();
    let s0 = deflected_ray_start();
    let tol = Tolerances {
        rtol: 1e-12,
        atol: 1e-12,
    };

    let res = rk45_step(&m, &s0, 50.0, &tol);
    assert!(!res.accepted);
    assert!(res.err > 1.0);
    assert!(res.dl_next < 50.0);
    assert_eq!(res.state.x, s0.x);
    assert_eq!(res.state.p, s0.p);
}

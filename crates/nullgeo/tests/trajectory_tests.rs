use nullgeo::integrator::rk4_step;
use nullgeo::metric::{State4, Vec4};
use nullgeo::metrics::schwarzschild::Schwarzschild;

#[test]
#[allow(clippy::excessive_precision)]
fn b_6p11_ray_escapes() {
    let m = Schwarzschild::new(1.0).unwrap();
    let mut s = State4 {
        x: Vec4::new(0.0, -15.0, 0.0, 0.0),
        p: Vec4::new(
            -9.30949336251262527e-1,
            8.50657002521700956e-1,
            3.79353281135443809e-1,
            0.0,
        ),
    };
    let dl = 0.01;
    let mut rmin = f64::INFINITY;
    let mut outcome = "max_steps";
    for _ in 0..5000 {
        s = rk4_step(&m, &s, dl);
        let r = (s.x[1] * s.x[1] + s.x[2] * s.x[2] + s.x[3] * s.x[3]).sqrt();
        rmin = rmin.min(r);
        if r < 2.1 {
            outcome = "captured";
            break;
        }
        if s.x[1] > 10.0 || r > 1.0e6 {
            outcome = "escaped";
            break;
        }
    }
    assert_eq!(outcome, "escaped", "rmin = {rmin}");
    assert!((rmin - 4.59).abs() < 0.05, "rmin = {rmin}, expected ~4.59");
}

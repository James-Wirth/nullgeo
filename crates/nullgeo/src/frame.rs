use crate::metric::{Mat4, Metric, Vec4};
use crate::{Error, Result};

pub fn metric_dot(g: &Mat4, a: &Vec4, b: &Vec4) -> f64 {
    a.dot(&(g * b))
}

fn project_out(g: &Mat4, v: &Vec4, onto: &Vec4) -> Vec4 {
    v - onto * (metric_dot(g, v, onto) / metric_dot(g, onto, onto))
}

fn unit_spacelike(g: &Mat4, v: &Vec4, x: &Vec4) -> Result<Vec4> {
    let len_sq = metric_dot(g, v, v);
    if len_sq <= 1e-24 {
        return Err(Error::DegenerateFrame(*x));
    }
    Ok(v / len_sq.sqrt())
}

pub fn build_coframe<M: Metric + ?Sized>(m: &M, x: &Vec4) -> Result<[Vec4; 4]> {
    build_coframe_seeded(
        m,
        x,
        [
            Vec4::new(0.0, 1.0, 0.0, 0.0),
            Vec4::new(0.0, 0.0, 1.0, 0.0),
            Vec4::new(0.0, 0.0, 0.0, 1.0),
        ],
    )
}

pub fn build_coframe_seeded<M: Metric + ?Sized>(
    m: &M,
    x: &Vec4,
    seeds: [Vec4; 3],
) -> Result<[Vec4; 4]> {
    build_coframe_for(m, x, &Vec4::new(1.0, 0.0, 0.0, 0.0), seeds)
}

pub fn build_coframe_for<M: Metric + ?Sized>(
    m: &M,
    x: &Vec4,
    observer: &Vec4,
    seeds: [Vec4; 3],
) -> Result<[Vec4; 4]> {
    let g = m.g(x);

    let u_len_sq = metric_dot(&g, observer, observer);
    if u_len_sq >= -1e-12 {
        return Err(Error::NonTimelikeObserver(*x));
    }
    let e0 = observer / (-u_len_sq).sqrt();

    let e1 = unit_spacelike(&g, &project_out(&g, &seeds[0], &e0), x)?;

    let mut s2 = project_out(&g, &seeds[1], &e0);
    s2 = project_out(&g, &s2, &e1);
    let e2 = unit_spacelike(&g, &s2, x)?;

    let mut s3 = project_out(&g, &seeds[2], &e0);
    s3 = project_out(&g, &s3, &e1);
    s3 = project_out(&g, &s3, &e2);
    let e3 = unit_spacelike(&g, &s3, x)?;

    Ok([g * e0, g * e1, g * e2, g * e3])
}

pub fn make_null_covector(coframe: &[Vec4; 4], n: [f64; 3], energy: f64) -> Vec4 {
    let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-300);
    coframe[0] * energy
        + coframe[1] * (energy * n[0] / norm)
        + coframe[2] * (energy * n[1] / norm)
        + coframe[3] * (energy * n[2] / norm)
}

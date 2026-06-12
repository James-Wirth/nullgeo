mod ellis;
mod kerr;
pub mod kerr_schild;
mod minkowski;
mod pullback;
mod reissner_nordstrom;
mod schwarzschild;

pub use ellis::Ellis;
pub use kerr::Kerr;
pub use minkowski::Minkowski;
pub use pullback::{Pullback, Rotation, Transition};
pub use reissner_nordstrom::ReissnerNordstrom;
pub use schwarzschild::Schwarzschild;

use crate::geometry::{inner, Chart, Metric, Vec4};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkySide {
    Primary,
    Secondary,
}

pub trait Spacetime: Metric + Chart {
    fn is_captured(&self, _x: &Vec4) -> bool {
        false
    }

    fn sky_side(&self, _x: &Vec4) -> SkySide {
        SkySide::Primary
    }

    fn circular_orbits(&self) -> Option<&dyn CircularOrbits> {
        None
    }
}

pub trait CircularOrbits {
    fn isco_radius(&self) -> f64;

    fn four_velocity(&self, x: &Vec4) -> Option<Vec4>;
}

pub fn cartesian_circular_four_velocity<M: Metric + ?Sized>(
    metric: &M,
    x: &Vec4,
    omega: f64,
) -> Option<Vec4> {
    let u = Vec4::new(1.0, -omega * x[2], omega * x[1], 0.0);
    let norm_sq = inner(&metric.g(x), &u, &u);
    (norm_sq < 0.0).then(|| u / (-norm_sq).sqrt())
}

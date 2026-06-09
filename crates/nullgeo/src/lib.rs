#![forbid(unsafe_code)]

//! nullgeo conventions!
//! - The general relativist's metric signature: (-, +, +, +)
//! - Geometric units (G = c = 1), because life is too short for anything else.
//! - Coordinates: for flat space we use (t, x, y, z). For general metrics, there are no rules.

pub mod camera;
pub mod integrator;
pub mod metric;
pub mod metrics;
pub mod ray;

pub use camera::{Camera, CameraSpec};
pub use metric::{Mat4, Metric, State4, Vec4};
pub use ray::RayBundle;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub type Result<T> = std::result::Result<T, Error>;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("invalid argument: {0}")]
    InvalidArg(String),
}

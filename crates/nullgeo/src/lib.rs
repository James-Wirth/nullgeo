#![forbid(unsafe_code)]

//! nullgeo conventions!
//! - The general relativist's metric signature: (-, +, +, +)
//! - Geometric units (G = c = 1), because life is too short for anything else.
//! - Coordinates: for flat space we use (t, x, y, z). For general metrics, there are no rules.

pub mod camera;
pub mod frame;
pub mod integrator;
pub mod metric;
pub mod metrics;
pub mod ray;
pub mod spacetime;
pub mod tracer;

pub use camera::{Camera, CameraPose, CameraSpec};
pub use metric::{Mat4, Metric, State4, Vec4};
pub use ray::RayBundle;
pub use spacetime::{RayAlignment, SkySide, Spacetime};
pub use tracer::{trace, Termination, TraceConfig};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub type Result<T> = std::result::Result<T, Error>;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("invalid argument: {0}")]
    InvalidArg(String),
    #[error("coordinate time axis is not timelike at x = {0:?}")]
    NonTimelikeObserver(Vec4),
    #[error("frame seeds are degenerate at x = {0:?}")]
    DegenerateFrame(Vec4),
}

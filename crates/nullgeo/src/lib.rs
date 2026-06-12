#![forbid(unsafe_code)]

//! nullgeo conventions!
//! - The general relativist's metric signature: (-, +, +, +)
//! - Geometric units (G = c = 1), because life is too short for anything else.
//! - Coordinates: for flat space we use (t, x, y, z). For general metrics, there are no rules.

pub mod geometry;
pub mod integrator;
pub mod render;
pub mod spacetimes;
pub mod tracer;

pub use geometry::{Chart, Mat4, Metric, PhasePoint, RayAlignment, Vec4};
pub use render::{
    class_color, colorize, render, shade_beauty, shade_map, tone_map, trace_geometry, Camera,
    CameraPose, CameraSpec, Colormap, Disk, EquirectImage, GeometryBuffer, ImageF32, MapField,
    MapQuantity, RayClass, RayInfo, RayOutcome, Scene, SkyMap,
};
pub use spacetimes::{CircularOrbits, SkySide, Spacetime};
pub use tracer::{
    trace, trace_with_stats, EquatorialAnnulus, Termination, TraceConfig, TraceStats,
};

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

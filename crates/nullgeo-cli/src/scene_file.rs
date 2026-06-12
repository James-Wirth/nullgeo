use std::path::{Path, PathBuf};

use nullgeo::geometry::Vec4;
use nullgeo::integrator::Tolerances;
use nullgeo::spacetimes::{Ellis, Kerr, Minkowski, ReissnerNordstrom, Schwarzschild};
use nullgeo::{
    Camera, CameraPose, CameraSpec, Colormap, Disk, EquirectImage, MapQuantity, SkyMap, Spacetime,
    TraceConfig,
};
use serde::Deserialize;

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum MetricKind {
    Minkowski,
    Schwarzschild,
    ReissnerNordstrom,
    Kerr,
    Ellis,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SceneFile {
    pub metric: MetricSection,
    pub camera: CameraSection,
    pub disk: Option<DiskSection>,
    pub sky: SkySection,
    pub sky_secondary: Option<SkySection>,
    #[serde(default)]
    pub integrator: IntegratorSection,
    #[serde(rename = "output")]
    pub outputs: Vec<OutputSection>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MetricSection {
    pub kind: MetricKind,
    #[serde(default = "one")]
    pub mass: f64,
    #[serde(default)]
    pub spin: f64,
    #[serde(default)]
    pub charge: f64,
    #[serde(default = "one")]
    pub b0: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CameraSection {
    pub position: [f64; 3],
    #[serde(default)]
    pub look_at: [f64; 3],
    #[serde(default = "default_up")]
    pub up: [f64; 3],
    #[serde(default)]
    pub velocity: [f64; 3],
    #[serde(default = "default_fov")]
    pub fov_deg: f64,
    pub width: usize,
    pub height: usize,
    #[serde(default = "one_usize")]
    pub supersample: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiskSection {
    #[serde(default)]
    pub r_in: f64,
    pub r_out: f64,
    #[serde(default = "default_emissivity")]
    pub emissivity_index: f64,
    #[serde(default = "default_g_power")]
    pub g_power: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SkySection {
    pub checker_deg: Option<f64>,
    pub image: Option<PathBuf>,
    pub uniform: Option<[f32; 3]>,
    pub graticule_deg: Option<f64>,
    pub graticule_width_deg: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IntegratorSection {
    #[serde(default = "default_tol")]
    pub tol: f64,
    #[serde(default = "default_max_steps")]
    pub max_steps: usize,
    pub escape_radius: Option<f64>,
}

impl Default for IntegratorSection {
    fn default() -> Self {
        Self {
            tol: default_tol(),
            max_steps: default_max_steps(),
            escape_radius: None,
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputSection {
    pub path: PathBuf,
    #[serde(default)]
    pub kind: OutputKind,
    pub format: Option<OutputFormat>,
    pub colormap: Option<ColormapChoice>,
    #[serde(default = "default_exposure")]
    pub exposure: f32,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum OutputKind {
    #[default]
    Beauty,
    Classification,
    Redshift,
    ImageOrder,
    MinRadius,
    CoordTime,
    AffineLength,
    Steps,
    EscapeTheta,
    EscapePhi,
}

impl OutputKind {
    pub fn map_quantity(self) -> Option<MapQuantity> {
        match self {
            OutputKind::Beauty => None,
            OutputKind::Classification => Some(MapQuantity::Classification),
            OutputKind::Redshift => Some(MapQuantity::Redshift),
            OutputKind::ImageOrder => Some(MapQuantity::ImageOrder),
            OutputKind::MinRadius => Some(MapQuantity::MinRadius),
            OutputKind::CoordTime => Some(MapQuantity::CoordTime),
            OutputKind::AffineLength => Some(MapQuantity::AffineLength),
            OutputKind::Steps => Some(MapQuantity::Steps),
            OutputKind::EscapeTheta => Some(MapQuantity::EscapeTheta),
            OutputKind::EscapePhi => Some(MapQuantity::EscapePhi),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ColormapChoice {
    Viridis,
    Diverging,
}

impl ColormapChoice {
    pub fn to_colormap(self) -> Colormap {
        match self {
            ColormapChoice::Viridis => Colormap::Viridis,
            ColormapChoice::Diverging => Colormap::Diverging,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum OutputFormat {
    Png,
    Ppm,
    Pfm,
    Csv,
}

impl OutputSection {
    pub fn resolved_format(&self) -> OutputFormat {
        self.format
            .unwrap_or(match self.path.extension().and_then(|e| e.to_str()) {
                Some("ppm") => OutputFormat::Ppm,
                Some("pfm") => OutputFormat::Pfm,
                Some("csv") => OutputFormat::Csv,
                _ => OutputFormat::Png,
            })
    }
}

fn one() -> f64 {
    1.0
}

fn one_usize() -> usize {
    1
}

fn default_up() -> [f64; 3] {
    [0.0, 0.0, 1.0]
}

fn default_fov() -> f64 {
    60.0
}

fn default_emissivity() -> f64 {
    2.0
}

fn default_g_power() -> f64 {
    3.0
}

fn default_tol() -> f64 {
    1e-9
}

fn default_max_steps() -> usize {
    100_000
}

fn default_exposure() -> f32 {
    1.0
}

fn default_graticule_width() -> f64 {
    0.25
}

pub fn build_spacetime(m: &MetricSection) -> Result<Box<dyn Spacetime + Sync>, String> {
    let invalid = |e: nullgeo::Error| format!("invalid metric: {e}");
    Ok(match m.kind {
        MetricKind::Minkowski => Box::new(Minkowski),
        MetricKind::Schwarzschild => Box::new(Schwarzschild::new(m.mass).map_err(invalid)?),
        MetricKind::ReissnerNordstrom => {
            Box::new(ReissnerNordstrom::new(m.mass, m.charge).map_err(invalid)?)
        }
        MetricKind::Kerr => Box::new(Kerr::new(m.mass, m.spin).map_err(invalid)?),
        MetricKind::Ellis => Box::new(Ellis::new(m.b0).map_err(invalid)?),
    })
}

pub fn build_camera(c: &CameraSection) -> Result<Camera, String> {
    Camera::new(
        CameraSpec {
            fov_deg: c.fov_deg,
            res: (c.width, c.height),
            energy: 1.0,
            supersample: c.supersample,
        },
        CameraPose {
            position: Vec4::new(0.0, c.position[0], c.position[1], c.position[2]),
            look_at: c.look_at,
            up: c.up,
            velocity: c.velocity,
        },
    )
    .map_err(|e| format!("invalid camera: {e}"))
}

pub fn build_sky(s: &SkySection, base: &Path) -> Result<SkyMap, String> {
    let background = match (s.checker_deg, &s.image, s.uniform) {
        (Some(angular_size_deg), None, None) => {
            Some(SkyMap::checker(angular_size_deg).map_err(|e| e.to_string())?)
        }
        (None, Some(path), None) => {
            let full = if path.is_absolute() {
                path.clone()
            } else {
                base.join(path)
            };
            let img = image::open(&full)
                .map_err(|e| format!("cannot load sky image {}: {e}", full.display()))?
                .to_rgb32f();
            let (w, h) = img.dimensions();
            let data = img
                .pixels()
                .map(|p| [p.0[0].powf(2.2), p.0[1].powf(2.2), p.0[2].powf(2.2)])
                .collect();
            Some(
                EquirectImage::new(w as usize, h as usize, data)
                    .map(SkyMap::Equirect)
                    .map_err(|e| e.to_string())?,
            )
        }
        (None, None, Some(color)) => Some(SkyMap::Uniform(color)),
        (None, None, None) => None,
        _ => return Err("sky must set at most one of checker_deg, image, uniform".into()),
    };

    match (s.graticule_deg, background) {
        (Some(spacing_deg), background) => SkyMap::graticule(
            spacing_deg,
            s.graticule_width_deg.unwrap_or(default_graticule_width()),
            background,
        )
        .map_err(|e| e.to_string()),
        (None, _) if s.graticule_width_deg.is_some() => {
            Err("graticule_width_deg requires graticule_deg".into())
        }
        (None, Some(sky)) => Ok(sky),
        (None, None) => {
            Err("sky must set one of checker_deg, image, uniform, graticule_deg".into())
        }
    }
}

pub fn build_disk(d: &DiskSection) -> Disk {
    Disk {
        r_in: d.r_in,
        r_out: d.r_out,
        emissivity_index: d.emissivity_index,
        g_power: d.g_power,
    }
}

pub fn build_trace_config(i: &IntegratorSection, camera_position: [f64; 3]) -> TraceConfig {
    let [x, y, z] = camera_position;
    let r = (x * x + y * y + z * z).sqrt();
    TraceConfig {
        tol: Tolerances {
            rtol: i.tol,
            atol: i.tol,
        },
        max_steps: i.max_steps,
        escape_radius: i.escape_radius.unwrap_or((4.0 * r).max(100.0)),
        ..TraceConfig::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &str) -> SceneFile {
        toml::from_str(text).unwrap()
    }

    #[test]
    fn example_scenes_parse_and_build() {
        let examples = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples");
        for name in [
            "kerr_disk.toml",
            "ellis_checker.toml",
            "kerr_diagnostics.toml",
        ] {
            let text = std::fs::read_to_string(examples.join(name)).unwrap();
            let file = parse(&text);
            build_spacetime(&file.metric).unwrap();
            build_camera(&file.camera).unwrap();
            build_sky(&file.sky, &examples).unwrap();
            if let Some(sky) = &file.sky_secondary {
                build_sky(sky, &examples).unwrap();
            }
            build_trace_config(&file.integrator, file.camera.position);
            assert!(!file.outputs.is_empty());
        }
    }

    #[test]
    fn full_scene_round_trips() {
        let file = parse(
            r#"
            [metric]
            kind = "kerr"
            mass = 2.0
            spin = 2.5

            [camera]
            position = [-50.0, 0.0, 6.0]
            look_at = [0.0, 0.0, 1.0]
            up = [0.0, 1.0, 0.0]
            fov_deg = 25.0
            width = 64
            height = 48
            supersample = 2

            [disk]
            r_in = 7.0
            r_out = 30.0
            emissivity_index = 2.5
            g_power = 4.0

            [sky]
            checker_deg = 15.0

            [sky_secondary]
            uniform = [0.1, 0.2, 0.3]

            [integrator]
            tol = 1e-10
            max_steps = 50000
            escape_radius = 300.0

            [[output]]
            path = "out.ppm"
            format = "ppm"
            exposure = 2.5

            [[output]]
            path = "redshift.png"
            kind = "redshift"
            colormap = "diverging"

            [[output]]
            path = "order.csv"
            kind = "image-order"
            "#,
        );

        assert_eq!(file.metric.kind, MetricKind::Kerr);
        assert_eq!(file.metric.mass, 2.0);
        assert_eq!(file.metric.spin, 2.5);
        assert_eq!(file.camera.fov_deg, 25.0);
        assert_eq!(file.camera.width, 64);
        assert_eq!(file.camera.supersample, 2);
        let disk = file.disk.unwrap();
        assert_eq!(disk.r_in, 7.0);
        assert_eq!(disk.g_power, 4.0);
        assert!(matches!(
            file.sky_secondary,
            Some(SkySection {
                uniform: Some(_),
                ..
            })
        ));
        assert_eq!(file.integrator.tol, 1e-10);
        assert_eq!(file.integrator.escape_radius, Some(300.0));
        assert_eq!(file.outputs.len(), 3);
        assert_eq!(file.outputs[0].kind, OutputKind::Beauty);
        assert_eq!(file.outputs[0].resolved_format(), OutputFormat::Ppm);
        assert_eq!(file.outputs[0].exposure, 2.5);
        assert_eq!(file.outputs[1].kind, OutputKind::Redshift);
        assert_eq!(file.outputs[1].colormap, Some(ColormapChoice::Diverging));
        assert_eq!(file.outputs[1].resolved_format(), OutputFormat::Png);
        assert_eq!(file.outputs[2].kind, OutputKind::ImageOrder);
        assert_eq!(file.outputs[2].resolved_format(), OutputFormat::Csv);

        let cfg = build_trace_config(&file.integrator, file.camera.position);
        assert_eq!(cfg.escape_radius, 300.0);
        assert_eq!(cfg.max_steps, 50000);

        assert!(build_spacetime(&file.metric).is_err());
    }

    #[test]
    fn defaults_fill_optional_sections() {
        let file = parse(
            r#"
            [metric]
            kind = "schwarzschild"

            [camera]
            position = [-15.0, 0.0, 0.0]
            width = 8
            height = 8

            [sky]
            checker_deg = 20.0

            [[output]]
            path = "out.png"
            "#,
        );
        assert_eq!(file.metric.mass, 1.0);
        assert_eq!(file.camera.fov_deg, 60.0);
        assert_eq!(file.camera.up, [0.0, 0.0, 1.0]);
        assert_eq!(file.camera.velocity, [0.0; 3]);
        assert_eq!(file.camera.supersample, 1);
        assert_eq!(file.integrator.tol, 1e-9);
        assert_eq!(file.integrator.max_steps, 100_000);
        assert_eq!(file.outputs[0].kind, OutputKind::Beauty);
        assert_eq!(file.outputs[0].colormap, None);
        assert_eq!(file.outputs[0].resolved_format(), OutputFormat::Png);
        assert_eq!(file.outputs[0].exposure, 1.0);
        let cfg = build_trace_config(&file.integrator, file.camera.position);
        assert_eq!(cfg.escape_radius, 100.0);
    }

    fn sky_section() -> SkySection {
        SkySection {
            checker_deg: None,
            image: None,
            uniform: None,
            graticule_deg: None,
            graticule_width_deg: None,
        }
    }

    #[test]
    fn seamed_checker_rejected() {
        let sky = |deg| SkySection {
            checker_deg: Some(deg),
            ..sky_section()
        };
        assert!(build_sky(&sky(7.0), Path::new(".")).is_err());
        assert!(build_sky(&sky(40.0), Path::new(".")).is_err());
        assert!(build_sky(&sky(-15.0), Path::new(".")).is_err());
        assert!(build_sky(&sky(12.0), Path::new(".")).is_ok());
    }

    #[test]
    fn sky_requires_exactly_one_source() {
        let sky = SkySection {
            checker_deg: Some(10.0),
            uniform: Some([1.0, 1.0, 1.0]),
            ..sky_section()
        };
        assert!(build_sky(&sky, Path::new(".")).is_err());
        assert!(build_sky(&sky_section(), Path::new(".")).is_err());
    }

    #[test]
    fn graticule_sky_builds_alone_or_over_background() {
        let alone = SkySection {
            graticule_deg: Some(10.0),
            ..sky_section()
        };
        assert!(matches!(
            build_sky(&alone, Path::new(".")).unwrap(),
            SkyMap::Graticule { .. }
        ));

        let over_uniform = SkySection {
            graticule_deg: Some(15.0),
            graticule_width_deg: Some(0.5),
            uniform: Some([0.1, 0.1, 0.2]),
            ..sky_section()
        };
        assert!(matches!(
            build_sky(&over_uniform, Path::new(".")).unwrap(),
            SkyMap::Graticule { .. }
        ));

        let seamed = SkySection {
            graticule_deg: Some(7.0),
            ..sky_section()
        };
        assert!(build_sky(&seamed, Path::new(".")).is_err());

        let width_without_lines = SkySection {
            graticule_width_deg: Some(0.5),
            uniform: Some([0.0; 3]),
            ..sky_section()
        };
        assert!(build_sky(&width_without_lines, Path::new(".")).is_err());
    }
}

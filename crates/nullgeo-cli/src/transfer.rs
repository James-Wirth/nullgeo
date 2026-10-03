//! CLI container for lossless, contributing-sample transfer data. No extra tracing.
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};

use crate::scene_file::{MetricKind, MetricSection, SceneFile};
use nullgeo::{
    Camera, DiskModel, GeometryBuffer, RayInfo, RayOutcome, Scene, SkyMap, SkySide, TraceConfig,
};
use serde::Serialize;

// Packed structured NPY records; write every scalar explicitly (no native struct padding).
const DESCR: &str = "[('pixel_index', '<u8'), ('sample_index', '<u8'), ('offset_x', '<f8'), ('offset_y', '<f8'), ('screen_u', '<f8'), ('screen_v', '<f8'), ('weight', '<f8'), ('render_weight', '<f4'), ('radius', '<f8'), ('g', '<f8'), ('radiance', '<f4'), ('has_intersection', '|u1'), ('radius_valid', '|u1'), ('g_available', '|u1'), ('g_valid', '|u1'), ('radiance_valid', '|u1'), ('outcome', '|u1'), ('finished', '|u1'), ('steps_accepted', '<u8'), ('steps_rejected', '<u8')]";

pub fn validate(
    scene: &Scene,
    camera: &Camera,
    cfg: &TraceConfig,
    metric: &MetricSection,
) -> Result<(), String> {
    if !matches!(metric.kind, MetricKind::Kerr | MetricKind::Schwarzschild) {
        return Err("transfer export requires Kerr or Schwarzschild circular disk orbits".into());
    }
    if !metric.mass.is_finite() || metric.mass <= 0.0 || !metric.spin.is_finite() {
        return Err("transfer export requires positive finite mass and finite spin".into());
    }
    let Some(disk) = scene.disk else {
        return Err("transfer export requires a stylized thin opaque disk".into());
    };
    let DiskModel::Stylized {
        emissivity_index,
        g_power,
    } = disk.model
    else {
        return Err("transfer export supports only the stylized thin opaque disk, not blackbody/volume emission".into());
    };
    let black = |sky: &SkyMap| matches!(sky, SkyMap::Uniform(c) if *c == [0.0; 3]);
    if !black(&scene.sky) || scene.sky_secondary.as_ref().is_some_and(|s| !black(s)) {
        return Err(
            "transfer export requires a uniform black sky (including secondary sky)".into(),
        );
    }
    if ![
        disk.r_in,
        disk.r_out,
        emissivity_index,
        g_power,
        camera.spec.energy,
        cfg.tol.rtol,
        cfg.tol.atol,
        cfg.escape_radius,
    ]
    .iter()
    .all(|x| x.is_finite())
        || !camera
            .pose
            .position
            .iter()
            .chain(camera.pose.look_at.iter())
            .chain(camera.pose.up.iter())
            .all(|x| x.is_finite())
        || cfg.tol.rtol <= 0.0
        || cfg.tol.atol <= 0.0
        || cfg.escape_radius <= 0.0
    {
        return Err("transfer export requires finite scene parameters and positive tolerances/escape radius".into());
    }
    // All current camera modes (regular, deterministic jitter, adaptive replacement) are supported.
    Ok(())
}

fn outcome(info: &RayInfo) -> (u8, bool) {
    match info.outcome {
        RayOutcome::Captured => (0, true),
        RayOutcome::Escaped {
            side: SkySide::Primary,
            ..
        } => (1, true),
        RayOutcome::Escaped {
            side: SkySide::Secondary,
            ..
        } => (2, true),
        RayOutcome::Saturated => (3, true),
        RayOutcome::MaxSteps => (4, false),
        RayOutcome::Stalled => (5, false),
    }
}

fn write_record(
    w: &mut impl Write,
    camera: &Camera,
    pixel: usize,
    sample: usize,
    offset: (f64, f64),
    count: usize,
    info: &RayInfo,
) -> std::io::Result<()> {
    for value in [pixel, sample] {
        w.write_all(&(value as u64).to_le_bytes())?;
    }
    let (u, v) = camera.screen_coordinates(pixel, offset);
    for value in [offset.0, offset.1, u, v, 1.0 / count as f64] {
        w.write_all(&value.to_le_bytes())?;
    }
    w.write_all(&(1.0 / count as f32).to_le_bytes())?;
    let radius = info.first_crossing.map(|c| c.radius);
    let g = info.first_crossing.and_then(|c| c.g);
    for value in [radius.unwrap_or(f64::NAN), g.unwrap_or(f64::NAN)] {
        w.write_all(&value.to_le_bytes())?;
    }
    w.write_all(&info.disk_radiance[0].to_le_bytes())?;
    let (code, finished) = outcome(info);
    w.write_all(&[
        info.first_crossing.is_some() as u8,
        radius.is_some_and(|r| r.is_finite() && r > 0.0) as u8,
        g.is_some() as u8,
        g.is_some_and(|g| g.is_finite() && g > 0.0) as u8,
        info.disk_radiance[0].is_finite() as u8,
        code,
        finished as u8,
    ])?;
    for value in [info.stats.steps_accepted, info.stats.steps_rejected] {
        w.write_all(&(value as u64).to_le_bytes())?;
    }
    Ok(())
}

#[derive(Serialize)]
struct Metadata {
    schema: &'static str,
    schema_version: u32,
    width: usize,
    height: usize,
    sample_count: usize,
    record_bytes: usize,
    array: &'static str,
    source_scene: &'static str,
    sampling: &'static str,
    effective_r_in: f64,
    effective_r_out: f64,
    camera_energy: f64,
    camera_time: f64,
    rtol: f64,
    atol: f64,
    dl_init: f64,
    dl_min: f64,
    dl_max: f64,
    escape_radius: f64,
    max_steps: usize,
    package_version: &'static str,
    build_id: &'static str,
    git_revision: &'static str,
    git_state_at_build: &'static str,
    source_digest: &'static str,
    rustc: &'static str,
    target: &'static str,
    profile: &'static str,
    parallel: bool,
    resolved_scene: toml::Value,
}

#[allow(clippy::too_many_arguments)]
pub fn write(
    path: &Path,
    buffer: &GeometryBuffer,
    camera: &Camera,
    scene: &Scene,
    cfg: &TraceConfig,
    file: &SceneFile,
    original: &str,
) -> Result<(), String> {
    validate(scene, camera, cfg, &file.metric)?;
    let annulus = buffer
        .annulus
        .ok_or("transfer export missing effective annulus")?;
    // Boundary/degenerate intersections cannot be represented as an opaque first hit.
    for pixel in 0..buffer.pixel_count() {
        let mut unsupported = false;
        buffer.for_each_sample(pixel, |info| {
            if info.first_crossing.is_some()
                && (info.disk_transmission != 0.0 || !matches!(info.outcome, RayOutcome::Saturated))
            {
                unsupported = true;
            }
        });
        if unsupported {
            return Err(
                "transfer export encountered a non-opaque/degenerate first intersection".into(),
            );
        }
    }
    let count = (0..buffer.pixel_count())
        .map(|p| buffer.sample_count(p))
        .sum();
    let mut resolved = toml::Value::try_from(file).map_err(|e| e.to_string())?;
    resolved["camera"].as_table_mut().unwrap().insert(
        "supersample_max".into(),
        (camera.spec.supersample_max as i64).into(),
    );
    if let Some(disk) = scene.disk {
        if let DiskModel::Stylized {
            emissivity_index,
            g_power,
        } = disk.model
        {
            let table = resolved["disk"].as_table_mut().unwrap();
            table.insert("emissivity_index".into(), emissivity_index.into());
            table.insert("g_power".into(), g_power.into());
        }
    }
    let metadata = Metadata {
        schema: "nullgeo.thin-disk-transfer",
        schema_version: 1,
        width: buffer.width,
        height: buffer.height,
        sample_count: count,
        record_bytes: 103,
        array: "samples.npy",
        source_scene: "scene.toml",
        sampling: "camera-grid-v1-adaptive-replacement",
        effective_r_in: annulus.r_in,
        effective_r_out: annulus.r_out,
        camera_energy: camera.spec.energy,
        camera_time: camera.pose.position[0],
        rtol: cfg.tol.rtol,
        atol: cfg.tol.atol,
        dl_init: cfg.dl_init,
        dl_min: cfg.dl_min,
        dl_max: cfg.dl_max,
        escape_radius: cfg.escape_radius,
        max_steps: cfg.max_steps,
        package_version: env!("CARGO_PKG_VERSION"),
        build_id: env!("NULLGEO_BUILD_ID"),
        git_revision: env!("NULLGEO_REVISION"),
        git_state_at_build: env!("NULLGEO_SOURCE_STATE"),
        source_digest: "unknown",
        rustc: env!("NULLGEO_RUSTC"),
        target: env!("NULLGEO_TARGET"),
        profile: env!("NULLGEO_PROFILE"),
        parallel: cfg!(feature = "parallel"),
        resolved_scene: resolved,
    };
    let metadata = toml::to_string_pretty(&metadata).map_err(|e| e.to_string())?;
    let io = || -> std::io::Result<()> {
        fs::create_dir(path)?; // Refuse overwrites, including concurrent creators.
        let mut writer = BufWriter::new(File::create(path.join("samples.npy"))?);
        let mut header =
            format!("{{'descr': {DESCR}, 'fortran_order': False, 'shape': ({count},), }}");
        let padding = (64 - (10 + header.len() + 1) % 64) % 64;
        header.push_str(&" ".repeat(padding));
        header.push('\n');
        writer.write_all(b"\x93NUMPY\x01\x00")?;
        writer.write_all(&(header.len() as u16).to_le_bytes())?;
        writer.write_all(header.as_bytes())?;
        let base = camera.subpixel_offsets();
        let refined = if camera.spec.supersample_max > camera.spec.supersample {
            camera.subpixel_offsets_for(camera.spec.supersample_max)
        } else {
            Vec::new()
        };
        for pixel in 0..buffer.pixel_count() {
            let offsets = if buffer.refined[pixel].is_empty() {
                &base
            } else {
                &refined
            };
            let mut sample = 0;
            let mut result = Ok(());
            buffer.for_each_sample(pixel, |info| {
                if result.is_ok() {
                    result = write_record(
                        &mut writer,
                        camera,
                        pixel,
                        sample,
                        offsets[sample],
                        offsets.len(),
                        info,
                    );
                }
                sample += 1;
            });
            result?;
        }
        writer.flush()?;
        fs::write(path.join("scene.toml"), original)?;
        // Written last: a directory without metadata.toml is an incomplete export.
        fs::write(path.join("metadata.toml"), metadata)?;
        Ok(())
    };
    io().map_err(|e| format!("cannot write transfer export {}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nullgeo::{CameraPose, CameraSpec, DiskCrossing, TraceStats, Vec4};

    #[test]
    fn packed_records_preserve_masks_and_all_outcomes() {
        let camera = Camera::new(
            CameraSpec {
                fov_deg: 40.0,
                res: (2, 2),
                energy: 1.0,
                supersample: 1,
                supersample_max: 1,
                jitter: false,
            },
            CameraPose {
                position: Vec4::new(0.0, -20.0, 0.0, 5.0),
                look_at: [0.0; 3],
                up: [0.0, 0.0, 1.0],
                velocity: [0.0; 3],
            },
        )
        .unwrap();
        let outcomes = [
            RayOutcome::Captured,
            RayOutcome::Escaped {
                side: SkySide::Primary,
                dir: [1.0, 0.0, 0.0],
                g: None,
            },
            RayOutcome::Escaped {
                side: SkySide::Secondary,
                dir: [1.0, 0.0, 0.0],
                g: None,
            },
            RayOutcome::Saturated,
            RayOutcome::MaxSteps,
            RayOutcome::Stalled,
        ];
        for (code, outcome) in outcomes.into_iter().enumerate() {
            for crossing in [
                None,
                Some(DiskCrossing {
                    radius: 8.0,
                    g: None,
                }),
                Some(DiskCrossing {
                    radius: 8.0,
                    g: Some(0.75),
                }),
                Some(DiskCrossing {
                    radius: f64::NAN,
                    g: Some(f64::INFINITY),
                }),
            ] {
                let info = RayInfo {
                    outcome,
                    first_crossing: crossing,
                    disk_radiance: [0.5; 3],
                    disk_transmission: 0.0,
                    stats: TraceStats {
                        affine_length: 0.0,
                        coord_time: 0.0,
                        equatorial_crossings: 0,
                        min_radius: 8.0,
                        steps_accepted: 123,
                        steps_rejected: 7,
                    },
                };
                let mut bytes = Vec::new();
                write_record(&mut bytes, &camera, 0, 2, (0.25, 0.75), 9, &info).unwrap();
                assert_eq!(bytes.len(), 103);
                assert_eq!(
                    &bytes[80..87],
                    &[
                        crossing.is_some() as u8,
                        crossing.is_some_and(|c| c.radius.is_finite()) as u8,
                        crossing.is_some_and(|c| c.g.is_some()) as u8,
                        crossing.is_some_and(|c| c.g.is_some_and(f64::is_finite)) as u8,
                        1,
                        code as u8,
                        (code < 4) as u8,
                    ]
                );
                assert_eq!(u64::from_le_bytes(bytes[87..95].try_into().unwrap()), 123);
                assert_eq!(u64::from_le_bytes(bytes[95..103].try_into().unwrap()), 7);
                assert!(f64::from_le_bytes(bytes[32..40].try_into().unwrap()) < 0.0);
                assert!(f64::from_le_bytes(bytes[40..48].try_into().unwrap()) > 0.0);
            }
        }
    }
}

use nullgeo::geometry::{build_coframe, make_null_covector, Metric, PhasePoint, Vec4};
use nullgeo::integrator::Tolerances;
use nullgeo::spacetimes::Schwarzschild;
use nullgeo::{
    colorize, shade_beauty, shade_map, trace_geometry, trace_with_stats, Camera, CameraPose,
    CameraSpec, Colormap, Disk, MapField, MapQuantity, RayClass, RayOutcome, Scene, SkyMap,
    Termination, TraceConfig,
};

const B_CRIT: f64 = 5.196152422706632;

fn backward_ray<M: Metric>(m: &M, x: Vec4, dir: [f64; 3], energy: f64) -> PhasePoint {
    let coframe = build_coframe(m, &x).unwrap();
    let arriving = make_null_covector(&coframe, [-dir[0], -dir[1], -dir[2]], energy);
    PhasePoint { x, p: -arriving }
}

fn face_on_ray(m: &Schwarzschild, z0: f64, b_aim: f64) -> PhasePoint {
    let alpha = (b_aim / z0).atan();
    backward_ray(
        m,
        Vec4::new(0.0, 0.0, 0.0, z0),
        [alpha.sin(), 0.0, -alpha.cos()],
        1.0,
    )
}

fn precise_config() -> TraceConfig {
    TraceConfig {
        tol: Tolerances {
            rtol: 1e-10,
            atol: 1e-10,
        },
        escape_radius: 4000.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    }
}

#[test]
fn min_radius_near_critical_b_approaches_photon_sphere() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cfg = precise_config();

    let (term, stats) = trace_with_stats(&m, face_on_ray(&m, 1000.0, B_CRIT * 1.001), &cfg);
    assert!(matches!(term, Termination::Escaped { .. }));
    assert!(
        stats.min_radius > 3.0 && stats.min_radius < 3.25,
        "closest approach {} should hug the photon sphere r = 3M",
        stats.min_radius
    );

    let (term, stats) = trace_with_stats(&m, face_on_ray(&m, 1000.0, B_CRIT * 0.999), &cfg);
    assert!(matches!(term, Termination::Captured { .. }));
    assert!(
        stats.min_radius < 2.1,
        "captured ray should reach the horizon, min r = {}",
        stats.min_radius
    );
}

#[test]
fn image_order_increments_across_the_photon_ring() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cfg = precise_config();

    let order = |b: f64| {
        let (term, stats) = trace_with_stats(&m, face_on_ray(&m, 1000.0, b), &cfg);
        assert!(
            matches!(term, Termination::Escaped { .. }),
            "ray at b = {b} should escape"
        );
        stats.equatorial_crossings
    };

    assert_eq!(order(12.0), 1);
    assert_eq!(order(5.6), 2);
    assert_eq!(order(5.215), 3);
}

#[test]
fn trace_stats_record_shapiro_consistent_times() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cfg = precise_config();

    let (_, near) = trace_with_stats(&m, face_on_ray(&m, 1000.0, 8.0), &cfg);
    let (_, far) = trace_with_stats(&m, face_on_ray(&m, 1000.0, 60.0), &cfg);
    assert!(near.coord_time > 0.0 && far.coord_time > 0.0);
    assert!(
        near.coord_time - near.affine_length > far.coord_time - far.affine_length,
        "the closer passage should accumulate more gravitational delay"
    );
    assert!(near.steps_accepted > 0);
}

fn shadow_setup() -> (Schwarzschild, Camera, Scene, TraceConfig) {
    let m = Schwarzschild::new(1.0).unwrap();
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 60.0,
            res: (24, 24),
            energy: 1.0,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, -15.0, 0.0, 0.0),
            look_at: [0.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap();
    let scene = Scene {
        sky: SkyMap::Uniform([1.0; 3]),
        sky_secondary: None,
        disk: None,
    };
    let cfg = TraceConfig {
        escape_radius: 100.0,
        ..TraceConfig::default()
    };
    (m, camera, scene, cfg)
}

#[test]
fn classification_map_degenerates_to_the_shadow() {
    let (m, camera, scene, cfg) = shadow_setup();
    let buffer = trace_geometry(&m, &camera, &scene, &cfg).unwrap();
    let image = shade_beauty(&buffer, &scene);
    let field = shade_map(&buffer, MapQuantity::Classification);

    let mut captured = 0;
    for (pixel, (color, &class_value)) in image.data.iter().zip(&field.values).enumerate() {
        let class = buffer.primary(pixel).class();
        match class {
            RayClass::Captured => {
                assert_eq!(color[0], 0.0, "shadow pixel {pixel} should be black");
                assert_eq!(class_value, 0.0);
                captured += 1;
            }
            RayClass::EscapedPrimary => {
                let RayOutcome::Escaped { g: Some(g), .. } = buffer.primary(pixel).outcome else {
                    panic!("escaped pixel {pixel} should carry a sky g-factor");
                };
                assert!(
                    g > 1.0,
                    "a static camera at r = 15M sees the sky gravitationally blueshifted"
                );
                let boosted_white = (g * g * g * g) as f32;
                assert_eq!(
                    color[0], boosted_white,
                    "sky pixel {pixel} should be white times g^4"
                );
                assert_eq!(class_value, 1.0);
            }
            other => panic!("unexpected class {other:?} at pixel {pixel}"),
        }
    }
    assert!(captured > 0, "the shadow should appear in frame");
    assert!(
        captured < image.data.len(),
        "the sky should appear in frame"
    );
}

#[test]
fn redshift_map_matches_face_on_formula() {
    let m = Schwarzschild::new(1.0).unwrap();
    let z_cam = 1000.0;
    let camera = Camera::new(
        CameraSpec {
            fov_deg: 2.0,
            res: (3, 3),
            energy: 1.0,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, 0.0, 0.0, z_cam),
            look_at: [12.0, 0.0, 0.0],
            up: [0.0, 1.0, 0.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap();
    let scene = Scene {
        sky: SkyMap::Uniform([0.0; 3]),
        sky_secondary: None,
        disk: Some(Disk::stylized(40.0)),
    };
    let cfg = precise_config();

    let buffer = trace_geometry(&m, &camera, &scene, &cfg).unwrap();
    let center = 4;
    let info = buffer.primary(center);
    let Some(crossing) = info.first_crossing else {
        panic!("center pixel should hit the disk, got {:?}", info.outcome);
    };
    let (radius, g) = (crossing.radius, crossing.g.unwrap());
    assert!((6.0..=40.0).contains(&radius));

    let predicted = (1.0 - 3.0 / radius).sqrt() / (1.0 - 2.0 / z_cam).sqrt();
    assert!(
        (g - predicted).abs() < 1e-5 * predicted,
        "g = {g:.9} vs predicted {predicted:.9} at r = {radius}"
    );

    let field = shade_map(&buffer, MapQuantity::Redshift);
    assert_eq!(field.values[center], g);
    assert!(shade_map(&buffer, MapQuantity::EscapeTheta).values[center].is_nan());
}

#[test]
fn diverging_colormap_centers_on_unit_redshift() {
    let field = MapField {
        width: 4,
        height: 1,
        quantity: MapQuantity::Redshift,
        values: vec![0.8, 1.0, 1.3, f64::NAN],
    };
    assert_eq!(field.display_range(), Some((0.7, 1.3)));

    let colors = colorize(&field, Colormap::Diverging);
    assert_eq!(colors[1], [220, 220, 220], "g = 1 should map to the center");
    assert_eq!(colors[2], [180, 4, 38], "the range top should saturate red");
    assert_eq!(colors[3], [0, 0, 0], "undefined pixels should be black");
}

#[test]
fn continuous_colormap_spans_viridis_endpoints() {
    let field = MapField {
        width: 3,
        height: 1,
        quantity: MapQuantity::MinRadius,
        values: vec![2.0, 5.0, 8.0],
    };
    let colors = colorize(&field, Colormap::Viridis);
    assert_eq!(colors[0], [68, 1, 84]);
    assert_eq!(colors[2], [253, 231, 37]);
}

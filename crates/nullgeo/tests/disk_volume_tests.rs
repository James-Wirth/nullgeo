use nullgeo::geometry::Vec4;
use nullgeo::integrator::Tolerances;
use nullgeo::spacetimes::{Kerr, Schwarzschild};
use nullgeo::{
    render, shade_beauty, shakura_sunyaev_peak_radius, tone_map, trace_disk, trace_geometry,
    Camera, CameraPose, CameraSpec, Disk, DiskModel, DiskVolume, EquatorialAnnulus, RayClass,
    Scene, SkyMap, Termination, TraceConfig,
};

const GOLDEN: [u8; 588] = [
    103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110,
    120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120,
    103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110,
    120, 103, 110, 120, 103, 110, 120, 14, 7, 1, 13, 7, 1, 10, 5, 0, 103, 110, 120, 103, 110, 120,
    103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103,
    110, 120, 22, 13, 4, 45, 29, 13, 57, 39, 20, 48, 32, 15, 31, 19, 7, 17, 9, 2, 7, 4, 0, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 28, 17, 6,
    76, 54, 31, 131, 108, 78, 110, 86, 57, 27, 16, 6, 20, 12, 3, 24, 14, 5, 15, 8, 2, 6, 3, 0, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 30, 18, 7, 77, 55, 32, 164, 144, 117,
    206, 196, 185, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 16, 9, 2, 10, 5, 0, 5, 2, 0, 103, 110, 120,
    53, 35, 17, 81, 58, 34, 123, 98, 69, 178, 161, 138, 226, 222, 219, 103, 110, 120, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 11, 6, 0, 12, 6, 0, 8, 4, 0, 5, 2, 0, 49, 32, 15, 58, 39, 20, 65, 46,
    25, 70, 49, 27, 69, 48, 27, 62, 43, 23, 51, 34, 16, 38, 24, 10, 27, 16, 6, 19, 11, 3, 13, 7, 1,
    9, 5, 0, 6, 3, 0, 5, 2, 0, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 83, 60,
    36, 201, 190, 176, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 17, 9, 2, 103, 110, 120, 103, 110, 120,
    103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 26, 16, 5, 104, 80,
    52, 110, 86, 58, 25, 15, 5, 20, 12, 3, 25, 15, 5, 9, 5, 0, 103, 110, 120, 103, 110, 120, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 20, 11, 3,
    40, 25, 11, 39, 25, 11, 24, 14, 5, 10, 5, 0, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110,
    120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120,
    103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110,
    120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120,
    103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103,
    110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120, 103, 110, 120,
];

fn camera(position: [f64; 3], look_at: [f64; 3], fov_deg: f64, res: (usize, usize)) -> Camera {
    Camera::new(
        CameraSpec {
            fov_deg,
            res,
            energy: 1.0,
            supersample: 1,
            supersample_max: 1,
            jitter: false,
        },
        CameraPose {
            position: Vec4::new(0.0, position[0], position[1], position[2]),
            look_at,
            up: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
        },
    )
    .unwrap()
}

fn blackbody(
    r_out: f64,
    t_in: f64,
    optical_depth: f64,
    aspect_ratio: f64,
    edge_taper: f64,
) -> Disk {
    Disk {
        r_in: 0.0,
        r_out,
        model: DiskModel::Blackbody {
            t_in,
            doppler_beaming: true,
            redshift_color: true,
            optical_depth,
            aspect_ratio,
            density_index: 3.0,
            edge_taper,
        },
    }
}

#[test]
fn thin_opaque_disk_render_is_bit_identical_to_the_legacy_disk() {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let cam = camera([-90.0, 0.0, 9.0], [0.0, 0.0, 0.0], 18.0, (14, 14));
    let scene = Scene {
        sky: SkyMap::Uniform([0.1, 0.12, 0.15]),
        sky_secondary: None,
        disk: Some(Disk::blackbody(16.0, 10_000.0)),
    };
    let cfg = TraceConfig {
        escape_radius: 400.0,
        max_steps: 200_000,
        ..TraceConfig::default()
    };
    let img = render(&kerr, &cam, &scene, &cfg).unwrap();
    let bytes: Vec<u8> = tone_map(&img, 1.5)
        .iter()
        .flat_map(|p| p.iter().copied())
        .collect();
    assert_eq!(bytes, GOLDEN);
}

#[test]
fn vertical_column_integrates_to_tau_perp() {
    let volume = DiskVolume {
        r_in: 6.0,
        r_out: 20.0,
        tau0: 7.0,
        aspect_ratio: 0.1,
        density_index: 3.0,
        w_out: 0.0,
    };
    for r in [9.0, 13.0, 17.0] {
        let h = volume.scale_height(r);
        let dz = h / 400.0;
        let mut integral = 0.0;
        let mut z = -12.0 * h;
        while z <= 12.0 * h {
            integral += volume.density_alpha(r, z) * dz;
            z += dz;
        }
        let tau_perp = volume.tau_perp(r);
        assert!(
            (integral - tau_perp).abs() < 1e-3 * tau_perp,
            "column at r = {r}: {integral} vs tau_perp {tau_perp}"
        );
    }
}

#[test]
fn tau_eff_increases_toward_grazing_incidence() {
    let volume = DiskVolume {
        r_in: 6.0,
        r_out: 20.0,
        tau0: 5.0,
        aspect_ratio: 0.0,
        density_index: 3.0,
        w_out: 0.0,
    };
    let mut last = 0.0;
    for mu in [1.0, 0.7, 0.4, 0.2, 0.05] {
        let tau = volume.tau_eff(12.0, mu);
        assert!(tau > last, "tau_eff at mu = {mu} did not increase: {tau}");
        last = tau;
    }
}

#[test]
fn surface_density_peaks_near_the_inner_edge_then_declines_to_zero() {
    let r_in = 6.0;
    let density_index = 3.0;
    let volume = DiskVolume {
        r_in,
        r_out: 20.0,
        tau0: 9.0,
        aspect_ratio: 0.1,
        density_index,
        w_out: 3.0,
    };
    assert_eq!(volume.tau_perp(6.0), 0.0);
    assert_eq!(volume.tau_perp(20.0), 0.0);

    let r_peak = r_in * ((2.0 * density_index + 1.0) / (2.0 * density_index)).powi(2);
    assert!((volume.tau_perp(r_peak) - 9.0).abs() < 1e-9);

    let mut last = f64::INFINITY;
    let mut r = r_peak + 0.5;
    while r < 20.0 {
        let here = volume.tau_perp(r);
        assert!(here <= last, "tau_perp not monotone declining at r = {r}");
        assert!(here < 9.0, "tau_perp exceeds the peak at r = {r}");
        last = here;
        r += 0.5;
    }
}

#[test]
fn occultation_fades_outward_in_step_with_emission() {
    let r_in = 6.0;
    let volume = DiskVolume {
        r_in,
        r_out: 20.0,
        tau0: 2.0,
        aspect_ratio: 0.0,
        density_index: 3.0,
        w_out: 4.0,
    };
    let mu = 0.3;
    let transmission = |r: f64| (-volume.tau_eff(r, mu)).exp();

    let r_peak = shakura_sunyaev_peak_radius(r_in);
    let mut last = transmission(r_peak);
    let mut r = r_peak + 0.5;
    while r < 20.0 {
        let here = transmission(r);
        assert!(
            here >= last - 1e-12,
            "transmission does not increase outward at r = {r}: {here} < {last}"
        );
        last = here;
        r += 0.5;
    }

    assert!(
        transmission(r_peak) < 0.05,
        "the hot inner disk should still occult: {}",
        transmission(r_peak)
    );
    assert!(
        transmission(19.9) > 0.9,
        "the cool outer rim should be nearly transparent: {}",
        transmission(19.9)
    );
}

#[test]
fn thin_translucent_disk_composites_radiance_over_background() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cam = camera([0.0, 0.0, 1000.0], [14.0, 0.0, 0.0], 4.0, (9, 9));
    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-10,
            atol: 1e-10,
        },
        escape_radius: 4000.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    };
    let sky = SkyMap::Uniform([0.3, 0.4, 0.5]);
    let translucent = Scene {
        sky: sky.clone(),
        sky_secondary: None,
        disk: Some(blackbody(40.0, 1.0e4, 1.5, 0.0, 0.0)),
    };
    let bare = Scene {
        sky: sky.clone(),
        sky_secondary: None,
        disk: None,
    };

    let buffer = trace_geometry(&m, &cam, &translucent, &cfg).unwrap();
    let img = shade_beauty(&buffer, &translucent);
    let img_bg = shade_beauty(&trace_geometry(&m, &cam, &bare, &cfg).unwrap(), &bare);

    let pixel = (0..buffer.pixel_count())
        .find(|&p| {
            buffer.primary(p).first_crossing.is_some()
                && (0.05..0.95).contains(&buffer.primary(p).disk_transmission)
        })
        .expect("a translucent disk pixel");
    let info = buffer.primary(pixel);
    for c in 0..3 {
        let expected = info.disk_radiance[c] + info.disk_transmission * img_bg.data[pixel][c];
        assert!(
            (img.data[pixel][c] - expected).abs() < 1e-5 * expected.max(1e-3),
            "channel {c}: {} vs disk_radiance + T*I_bg = {expected}",
            img.data[pixel][c]
        );
    }
    assert!(
        img.data[pixel]
            .iter()
            .zip(&info.disk_radiance)
            .any(|(p, r)| *p > *r + 1e-6),
        "the lensed background should brighten a translucent pixel above its bare emission"
    );
}

#[test]
fn translucent_edge_lets_the_background_bleed_through_where_opaque_does_not() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cam = camera([0.0, 0.0, 1000.0], [16.0, 0.0, 0.0], 4.0, (9, 9));
    let cfg = TraceConfig {
        tol: Tolerances {
            rtol: 1e-10,
            atol: 1e-10,
        },
        escape_radius: 4000.0,
        max_steps: 1_000_000,
        ..TraceConfig::default()
    };
    let sky = SkyMap::Uniform([1.0; 3]);
    let make = |optical_depth| Scene {
        sky: sky.clone(),
        sky_secondary: None,
        disk: Some(blackbody(40.0, 1.0e4, optical_depth, 0.0, 0.0)),
    };

    let translucent = make(1.0);
    let opaque = make(1.0e9);
    let buf_t = trace_geometry(&m, &cam, &translucent, &cfg).unwrap();
    let buf_o = trace_geometry(&m, &cam, &opaque, &cfg).unwrap();
    let img_t = shade_beauty(&buf_t, &translucent);
    let img_o = shade_beauty(&buf_o, &opaque);

    let pixel = (0..buf_t.pixel_count())
        .find(|&p| {
            buf_t.primary(p).first_crossing.is_some() && buf_t.primary(p).disk_transmission > 0.1
        })
        .expect("a translucent disk pixel");

    assert!(buf_t.primary(pixel).disk_transmission > 0.1);
    assert!(buf_o.primary(pixel).disk_transmission < 1e-4);

    let bled = img_t.data[pixel]
        .iter()
        .zip(&buf_t.primary(pixel).disk_radiance)
        .any(|(p, r)| *p > *r + 1e-3);
    assert!(bled, "stars should bleed through the translucent disk");

    for c in 0..3 {
        let radiance = buf_o.primary(pixel).disk_radiance[c];
        assert!(
            (img_o.data[pixel][c] - radiance).abs() < 1e-4,
            "the opaque disk pixel should be pure emission, no sky"
        );
    }
}

fn disk_class_count(aspect_ratio: f64) -> usize {
    let kerr = Kerr::new(1.0, 0.9).unwrap();
    let cam = camera([-60.0, 0.0, 1.5], [0.0, 0.0, 0.0], 16.0, (18, 18));
    let scene = Scene {
        sky: SkyMap::Uniform([0.0; 3]),
        sky_secondary: None,
        disk: Some(blackbody(12.0, 1.0e4, 30.0, aspect_ratio, 0.1)),
    };
    let cfg = TraceConfig {
        escape_radius: 300.0,
        max_steps: 100_000,
        ..TraceConfig::default()
    };
    let buffer = trace_geometry(&kerr, &cam, &scene, &cfg).unwrap();
    (0..buffer.pixel_count())
        .filter(|&p| buffer.primary(p).class() == RayClass::Disk)
        .count()
}

#[test]
fn thickness_widens_the_edge_on_silhouette() {
    let thin = disk_class_count(0.0);
    let thick = disk_class_count(0.18);
    assert!(
        thick > thin,
        "a thick disk should occult more pixels edge-on: thick {thick} vs thin {thin}"
    );
}

#[test]
fn opaque_core_ray_saturates_early_with_bounded_steps() {
    let m = Schwarzschild::new(1.0).unwrap();
    let cam = camera([-100.0, 0.0, 0.0], [0.0, 0.0, 0.0], 4.0, (1, 1));
    let ray = cam.pixel_rays(&m).unwrap()[0];
    let volume = DiskVolume {
        r_in: 6.0,
        r_out: 20.0,
        tau0: 40.0,
        aspect_ratio: 0.1,
        density_index: 3.0,
        w_out: 0.0,
    };
    let cfg = TraceConfig {
        disk: Some(EquatorialAnnulus {
            r_in: 6.0,
            r_out: 20.0,
        }),
        escape_radius: 500.0,
        max_steps: 200_000,
        ..TraceConfig::default()
    };
    let emission = |_r: f64, _g: Option<f64>| [1.0f32; 3];
    let result = trace_disk(&m, ray, &cfg, &volume, 1.0, &emission);
    assert!(
        matches!(result.termination, Termination::Saturated { .. }),
        "an in-plane ray through the dense core should saturate"
    );
    assert!(result.disk_transmission < 1e-4);
    assert!(
        result.stats.steps_accepted < 50_000,
        "saturation should bound the step count: {}",
        result.stats.steps_accepted
    );
}

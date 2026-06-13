use nullgeo::{
    planck_xyz, quantize16, quantize8, shakura_sunyaev_peak_radius, shakura_sunyaev_temperature,
    tone_map, tone_map_curve, xyz_to_linear_srgb, ImageF32, ToneCurve,
};

fn planck_srgb_normalized(t: f64) -> [f64; 3] {
    let xyz = planck_xyz(t);
    xyz_to_linear_srgb(xyz.map(|c| c / xyz[1]))
}

#[test]
fn blackbody_at_6500k_is_close_to_the_srgb_white_point() {
    let [r, g, b] = planck_srgb_normalized(6500.0);
    assert!((0.85..1.15).contains(&(r / g)), "r/g = {}", r / g);
    assert!((0.85..1.15).contains(&(b / g)), "b/g = {}", b / g);
}

#[test]
fn blackbody_color_orders_with_temperature() {
    let [r3, g3, b3] = planck_srgb_normalized(3000.0);
    assert!(
        r3 > g3 && g3 > b3,
        "3000 K should glow warm: {r3} {g3} {b3}"
    );

    let [r15, _, b15] = planck_srgb_normalized(15_000.0);
    assert!(b15 > r15, "15000 K should look blue: {r15} vs {b15}");

    let ratio = |t: f64| {
        let [r, _, b] = planck_srgb_normalized(t);
        b / r
    };
    assert!(ratio(4000.0) < ratio(6000.0));
    assert!(ratio(6000.0) < ratio(10_000.0));
}

#[test]
fn zero_temperature_emits_nothing() {
    assert_eq!(planck_xyz(0.0), [0.0; 3]);
    assert_eq!(planck_xyz(-100.0), [0.0; 3]);
}

#[test]
fn shakura_sunyaev_profile_vanishes_at_the_inner_edge_and_peaks_at_49_36() {
    let (t_in, r_in) = (1.0e4, 6.0);
    let temperature = |r: f64| shakura_sunyaev_temperature(t_in, r_in, r);

    assert_eq!(temperature(r_in), 0.0);
    assert_eq!(temperature(0.5 * r_in), 0.0);

    let r_peak = shakura_sunyaev_peak_radius(r_in);
    assert!((r_peak - 6.0 * 49.0 / 36.0).abs() < 1e-12);
    assert!(temperature(r_peak) > temperature(r_peak * 1.01));
    assert!(temperature(r_peak) > temperature(r_peak * 0.99));

    let r = 2.0 * r_in;
    let expected = t_in * 2.0f64.powf(-0.75) * (1.0 - 0.5f64.sqrt()).powf(0.25);
    assert!((temperature(r) - expected).abs() < 1e-9 * expected);
}

#[test]
fn aces_curve_is_monotone_fixes_black_and_caps_at_one() {
    let img = ImageF32 {
        width: 6,
        height: 1,
        data: vec![
            [0.0; 3], [0.05; 3], [0.3; 3], [1.0; 3], [5.0; 3], [100.0; 3],
        ],
    };
    let display = tone_map_curve(&img, 1.0, ToneCurve::Aces);
    assert_eq!(display[0], [0.0; 3]);
    for pair in display.windows(2) {
        assert!(pair[0][0] < pair[1][0], "ACES not monotone: {display:?}");
    }
    assert!(display[5][0] <= 1.0);
    assert!(
        display[5][0] > 0.99,
        "deep highlights should approach white"
    );
}

#[test]
fn tone_map_is_reinhard_curve_plus_8_bit_quantization() {
    let img = ImageF32 {
        width: 4,
        height: 1,
        data: vec![[0.0; 3], [0.123; 3], [1.7; 3], [42.0; 3]],
    };
    let via_stages = quantize8(&tone_map_curve(&img, 2.5, ToneCurve::Reinhard));
    assert_eq!(tone_map(&img, 2.5), via_stages);
}

#[test]
fn quantizers_hit_their_endpoints() {
    let display = [[0.0, 0.5, 1.0]];
    assert_eq!(quantize8(&display), vec![[0, 128, 255]]);
    assert_eq!(quantize16(&display), vec![[0, 32768, 65535]]);
}

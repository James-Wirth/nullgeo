use assert_cmd::Command;

fn render_shadow(metric: &str, path: &std::path::Path) -> Vec<u8> {
    Command::cargo_bin("nullgeo")
        .unwrap()
        .args([
            "shadow",
            "--metric",
            metric,
            "--width",
            "16",
            "--height",
            "16",
            "--fov-deg",
            "60",
            "--out",
            path.to_str().unwrap(),
        ])
        .assert()
        .success();

    let data = std::fs::read(path).unwrap();
    let mut offset = 0;
    for _ in 0..3 {
        offset = data.iter().skip(offset).position(|&b| b == b'\n').unwrap() + offset + 1;
    }
    let header = std::str::from_utf8(&data[..offset]).unwrap();
    assert!(header.starts_with("P5\n16 16\n255\n"), "header: {header:?}");
    data[offset..].to_vec()
}

#[test]
fn schwarzschild_shadow_has_dark_center_and_bright_corners() {
    let path = std::env::temp_dir().join("nullgeo_cli_test_schwarzschild.ppm");
    let pixels = render_shadow("schwarzschild", &path);
    assert_eq!(pixels.len(), 256);
    assert_eq!(pixels[8 * 16 + 8], 0, "center pixel should be captured");
    assert_eq!(pixels[0], 255, "corner pixel should escape");
    assert_eq!(pixels[255], 255, "corner pixel should escape");
    std::fs::remove_file(path).ok();
}

#[test]
fn minkowski_shadow_is_empty() {
    let path = std::env::temp_dir().join("nullgeo_cli_test_minkowski.ppm");
    let pixels = render_shadow("minkowski", &path);
    assert!(pixels.iter().all(|&p| p == 255));
    std::fs::remove_file(path).ok();
}

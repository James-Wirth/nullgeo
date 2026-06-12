use assert_cmd::Command;

#[test]
fn render_command_writes_a_valid_png() {
    let dir = std::env::temp_dir().join(format!("nullgeo_render_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("out.png");
    let scene = format!(
        r#"
[metric]
kind = "schwarzschild"

[camera]
position = [-15.0, 0.0, 0.0]
fov_deg = 60.0
width = 12
height = 8

[sky]
checker_deg = 20.0

[integrator]
tol = 1e-8
max_steps = 20000

[[output]]
path = "{}"
"#,
        out.display()
    );
    let scene_path = dir.join("scene.toml");
    std::fs::write(&scene_path, scene).unwrap();

    Command::cargo_bin("nullgeo")
        .unwrap()
        .args(["render", scene_path.to_str().unwrap()])
        .assert()
        .success();

    let img = image::open(&out).unwrap();
    assert_eq!((img.width(), img.height()), (12, 8));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn render_command_writes_maps_from_one_trace() {
    let dir = std::env::temp_dir().join(format!("nullgeo_maps_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let scene = format!(
        r#"
[metric]
kind = "schwarzschild"

[camera]
position = [-15.0, 0.0, 0.0]
fov_deg = 60.0
width = 10
height = 6

[sky]
graticule_deg = 15.0

[integrator]
tol = 1e-8
max_steps = 20000

[[output]]
path = "{dir}/beauty.png"

[[output]]
path = "{dir}/class.png"
kind = "classification"

[[output]]
path = "{dir}/min_r.csv"
kind = "min-radius"

[[output]]
path = "{dir}/min_r.pfm"
kind = "min-radius"

[[output]]
path = "{dir}/order.ppm"
kind = "image-order"
"#,
        dir = dir.display()
    );
    let scene_path = dir.join("scene.toml");
    std::fs::write(&scene_path, scene).unwrap();

    Command::cargo_bin("nullgeo")
        .unwrap()
        .args(["render", scene_path.to_str().unwrap()])
        .assert()
        .success();

    let img = image::open(dir.join("class.png")).unwrap();
    assert_eq!((img.width(), img.height()), (10, 6));

    let csv = std::fs::read_to_string(dir.join("min_r.csv")).unwrap();
    let rows: Vec<&str> = csv.lines().collect();
    assert_eq!(rows.len(), 6);
    assert!(rows.iter().all(|r| r.split(',').count() == 10));

    let pfm = std::fs::read(dir.join("min_r.pfm")).unwrap();
    assert!(pfm.starts_with(b"Pf\n10 6\n-1.0\n"));
    assert_eq!(pfm.len(), b"Pf\n10 6\n-1.0\n".len() + 10 * 6 * 4);

    let ppm = std::fs::read(dir.join("order.ppm")).unwrap();
    assert!(ppm.starts_with(b"P6\n10 6\n255\n"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn render_command_rejects_invalid_scene() {
    let dir = std::env::temp_dir().join(format!("nullgeo_badscene_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let scene_path = dir.join("scene.toml");
    std::fs::write(&scene_path, "[metric]\nkind = \"warp-drive\"\n").unwrap();

    Command::cargo_bin("nullgeo")
        .unwrap()
        .args(["render", scene_path.to_str().unwrap()])
        .assert()
        .failure();
    std::fs::remove_dir_all(&dir).ok();
}

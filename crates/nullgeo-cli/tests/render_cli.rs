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

[output]
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

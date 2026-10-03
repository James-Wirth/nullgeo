use assert_cmd::Command;

fn scene(output: &str, model: &str, sky: &str, sampling: &str) -> String {
    format!(
        r#"
[metric]
kind = "kerr"
spin = 0.7
[camera]
position = [-30.0, 0.0, 20.0]
width = 5
height = 4
fov_deg = 55.0
{sampling}
[disk]
model = "{model}"
r_in = 0.0
r_out = 18.0
[sky]
{sky}
[integrator]
tol = 1e-7
max_steps = 4000
[[output]]
path = "{output}"
"#
    )
}

#[test]
fn export_preserves_render_and_rejects_unsupported_scenes() {
    let dir = std::env::temp_dir().join(format!("nullgeo_transfer_cli_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let config = dir.join("scene.toml");
    let pfm = dir.join("image.pfm");
    for (i, sampling) in [
        "",
        "supersample = 2",
        "supersample = 2\njitter = true",
        "supersample = 1\nsupersample_max = 3",
        "supersample = 2\nsupersample_max = 3\njitter = true",
    ]
    .iter()
    .enumerate()
    {
        std::fs::write(
            &config,
            scene(
                pfm.to_str().unwrap(),
                "stylized",
                "uniform = [0.0, 0.0, 0.0]",
                sampling,
            ),
        )
        .unwrap();
        Command::cargo_bin("nullgeo")
            .unwrap()
            .arg("render")
            .arg(&config)
            .assert()
            .success();
        let before = std::fs::read(&pfm).unwrap();
        let export = dir.join(format!("export-{i}"));
        Command::cargo_bin("nullgeo")
            .unwrap()
            .arg("render")
            .arg(&config)
            .arg("--transfer-export")
            .arg(&export)
            .assert()
            .success();
        assert_eq!(before, std::fs::read(&pfm).unwrap());
        let meta: toml::Value =
            toml::from_str(&std::fs::read_to_string(export.join("metadata.toml")).unwrap())
                .unwrap();
        assert_eq!(meta["schema_version"].as_integer(), Some(1));
        assert!(meta["effective_r_in"].as_float().unwrap() > 0.0); // ISCO clamp retained.
        let npy = std::fs::read(export.join("samples.npy")).unwrap();
        assert_eq!(&npy[..8], b"\x93NUMPY\x01\x00");
        let start = 10 + u16::from_le_bytes([npy[8], npy[9]]) as usize;
        assert_eq!(
            npy.len() - start,
            meta["sample_count"].as_integer().unwrap() as usize * 103
        );
        Command::cargo_bin("nullgeo")
            .unwrap()
            .arg("render")
            .arg(&config)
            .arg("--transfer-export")
            .arg(&export)
            .assert()
            .failure();
        assert_eq!(before, std::fs::read(&pfm).unwrap());
    }
    for (model, sky, sampling) in [
        ("blackbody", "uniform = [0.0, 0.0, 0.0]", ""),
        ("stylized", "checker_deg = 15.0", ""),
        ("stylized", "uniform = [1.0, 0.0, 0.0]", ""),
        ("stylized", "uniform = [0.0, 0.0, 0.0]", "supersample = 0"),
        (
            "stylized",
            "uniform = [0.0, 0.0, 0.0]",
            "supersample = 3\nsupersample_max = 2",
        ),
    ] {
        std::fs::write(&config, scene(pfm.to_str().unwrap(), model, sky, sampling)).unwrap();
        let export = dir.join("rejected");
        Command::cargo_bin("nullgeo")
            .unwrap()
            .arg("render")
            .arg(&config)
            .arg("--transfer-export")
            .arg(&export)
            .assert()
            .failure();
        assert!(!export.exists());
    }
    std::fs::remove_dir_all(dir).unwrap();
}

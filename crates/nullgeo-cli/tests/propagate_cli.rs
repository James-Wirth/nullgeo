use assert_cmd::Command;

#[test]
fn propagate_emits_monotone_null_ray_csv() {
    let assert = Command::cargo_bin("nullgeo")
        .unwrap()
        .args([
            "propagate",
            "--metric",
            "schwarzschild",
            "--pos",
            "-20,0,0",
            "--dir",
            "-1,0,0",
            "--max-steps",
            "10000",
        ])
        .assert()
        .success();

    let stdout = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    let mut lines = stdout.lines();
    assert_eq!(lines.next(), Some("lambda,t,x,y,z,H"));
    let rows: Vec<Vec<f64>> = lines
        .map(|line| line.split(',').map(|v| v.parse().unwrap()).collect())
        .collect();

    assert!(rows.len() > 10, "expected many steps, got {}", rows.len());
    for row in &rows {
        assert_eq!(row.len(), 6);
        assert!(row[5].abs() < 1e-9, "H drift {}", row[5]);
    }
    assert_eq!(rows[0][0], 0.0);
    assert!(
        rows.windows(2).all(|w| w[1][0] > w[0][0]),
        "lambda monotone"
    );
    assert!(
        rows.windows(2).all(|w| w[1][1] < w[0][1]),
        "past-directed t"
    );

    let last = rows.last().unwrap();
    assert!(last[2] < -100.0, "outgoing ray should reach escape radius");

    let stderr = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(stderr.contains("escaped"), "stderr: {stderr}");
}

#[test]
fn propagate_toward_horizon_is_captured() {
    let assert = Command::cargo_bin("nullgeo")
        .unwrap()
        .args([
            "propagate",
            "--metric",
            "schwarzschild",
            "--pos",
            "-20,0,0",
            "--dir",
            "1,0,0",
            "--max-steps",
            "10000",
        ])
        .assert()
        .success();

    let stderr = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(stderr.contains("captured"), "stderr: {stderr}");
}

use std::{
    env,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

fn output(program: &str, args: &[&str]) -> Option<String> {
    let result = Command::new(program).args(args).output().ok()?;
    result
        .status
        .success()
        .then(|| String::from_utf8_lossy(&result.stdout).trim().to_owned())
}

fn main() {
    // Source archives may have no Git metadata. Never infer a clean revision in that case.
    for path in [
        "src",
        "build.rs",
        "Cargo.toml",
        "../nullgeo/src",
        "../nullgeo/Cargo.toml",
        "../../Cargo.lock",
        "../../Cargo.toml",
        "../../.git/HEAD",
        "../../.git/index",
        "../../.git/refs",
        "../../.git/packed-refs",
    ] {
        println!("cargo:rerun-if-changed={path}");
    }
    let revision = output("git", &["rev-parse", "HEAD"]).unwrap_or("unknown".into());
    let state = output(
        "git",
        &["status", "--porcelain", "--untracked-files=normal"],
    )
    .map(|s| if s.is_empty() { "clean" } else { "dirty" })
    .unwrap_or("unknown");
    let rustc = output(&env::var("RUSTC").unwrap(), &["--version"]).unwrap_or("unknown".into());
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    for (key, value) in [
        ("NULLGEO_REVISION", revision),
        ("NULLGEO_SOURCE_STATE", state.into()),
        ("NULLGEO_RUSTC", rustc),
        ("NULLGEO_TARGET", env::var("TARGET").unwrap()),
        ("NULLGEO_PROFILE", env::var("PROFILE").unwrap()),
        ("NULLGEO_BUILD_ID", format!("transfer-v1-{stamp}")),
    ] {
        println!("cargo:rustc-env={key}={value}");
    }
}

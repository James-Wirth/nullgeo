fmt:; cargo fmt --all
lint:; cargo clippy --all-targets --all-features -- -D warnings
test:; cargo test --workspace
bench:; cargo bench
doc:; RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --workspace

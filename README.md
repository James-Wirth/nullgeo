<p align="center">
  <img src="https://github.com/James-Wirth/nullgeo/releases/download/assets/kerr_starfield_hdr.png" alt="Kerr black hole" width="800">
</p>


## nullgeo

This is a general relativistic ray tracing library written in Rust. 

Currently, we have implemented the Minkowski, Schwarzschild, Reissner–Nordström, Kerr, and the Ellis wormhole spacetimes.

All quantities are in geometrized units ($G = c = 1$), and we have adopted the signature $(-,+,+,+)$.

## Installation

### CLI

```
cargo install nullgeo-cli
```

### Library

```toml
[dependencies]
nullgeo = "0.2"
```

## Example Usage (with CLI)

A scene can be defined with a TOML file (see below). To render the scene, run:

```
nullgeo render scene.toml
```

### e.g. Kerr black hole, $a = 0.9M$, with accretion disk

```toml
[metric]
kind = "kerr"
mass = 1.0
spin = 0.9

[camera]
position = [-85.0, 0.0, 9.0]
fov_deg  = 24.0
width    = 640
height   = 360
supersample = 3
supersample_max = 4              

[disk]                           
r_out = 18.0
t_in  = 10000.0          
optical_depth = 2.5
aspect_ratio  = 0.05

[sky]
checker_deg = 15.0

[[output]]
path = "kerr_output.png"
kind = "beauty"
exposure = 4.0
tone = "aces"
```

### Trace a single ray (for debugging)

```
nullgeo propagate --metric kerr --spin 0.9 --pos=-20,0,5 --dir=1,0,0 --out ray.csv
```

creates a file `ray.csv` containing $(\lambda, t, x, y, z, H)$.

## License

MIT or Apache-2.0
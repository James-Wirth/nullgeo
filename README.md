<p align="center">
  <img src="https://raw.githubusercontent.com/James-Wirth/nullgeo/main/assets/logo.svg" alt="nullgeo" width="420">
</p>

## About

This is a fast general relativistic ray-tracing engine built with Rust. The long-term goal is to support relativistic visualisation (e.g. black hole shadows, gravitational lensing, and eventually radiative transfer effects like attenuation) in arbitrary spacetime geometries.  

All quantities are in geometrized units ($G=c=1$).

## Installation

### CLI

```
cargo install nullgeo-cli
```

This installs the `nullgeo` binary. Build with `--features parallel` to enable Rayon-based parallel ray tracing:

```
cargo install nullgeo-cli --features parallel
```

### Library

```
cargo add nullgeo
```

Or in `Cargo.toml`:

```toml
[dependencies]
nullgeo = "0.1"
```

Enable the optional `parallel` feature for Rayon support: `nullgeo = { version = "0.1", features = ["parallel"] }`. 

## Example Usage (with CLI)

### e.g. Schwarzschild shadow, $512 \times 512$, camera at $x=-30$

```
nullgeo shadow \
    --metric=schwarzschild \
    --mass=1.0 \
    --width=512 --height=512 \
    --fov-deg=30.0 \
    --cam-x=-30.0 \
    --dl=0.005 \
    --max-steps=20000 \
    --out=shadow.pgm
```

## Theory

### Hamiltonian Formulation

The Hamiltonian for photons in GR is

$H(x, p) = \frac{1}{2} g^{\mu \nu}(x) p_{\mu} p_{\nu}$

together with the null constraint $H=0$. Hamilton's equations give the trajectory:

$$
\dot{x}^{\mu} = \frac{\partial H}{\partial p_{\mu}} = g^{\mu \nu} p_{\nu}, \quad 
\dot{p}_{\mu} = -\frac{\partial H}{\partial x^{\mu}} = -\frac{1}{2} \partial_{\mu} g^{\rho \sigma} p_{\rho} p_{\sigma}
$$
where the dot signifies differentiation with respect to an affine parameter $\lambda$. 

### Pinhole Camera and the Local Orthonormal Frame

A local tetrad $\{ e^{\mu}_{a} \}$ is constructed at the camera's position. Note that the Latin $a$ in this notation labels the vector in the basis, and is not an index. 

The spatial unit vectors $e^{\mu}_i$ are constructed by projecting the coordinate basis $\partial_i$ onto the 3-plane orthogonal to the timelike unit vector $e^{\mu}_0$ ($\sim u$) and carrying out Gram-Schmidt orthogonalization with the induced metric $h:= g + u \otimes u$ (see e.g. Wald). 

A photon with spatial direction $\mathbf{n}$ and energy $E$ in the local tetrad has four momentum $p^a = E(1, n^i)$. The covariant components in the coordinate basis can then be obtained via $p_{\mu} = g_{\mu \nu} {e_a}^{\nu} p^a$.

## Numerics

We use Rayon's `par_iter_mut()` over the image buffer. For each pixel, a worker:

1. Initializes the ray state by computing the photon's null covector at the camera position.
2. Integrates the equations of motion using a 4th-order Runge-Kutta scheme to obtain the trajectory $(x^{\mu}(\lambda), p_{\mu}(\lambda))$.
3. Applies termination checks, e.g. stopping if the ray falls inside the horizon.
4. Writes the pixel value to the grayscale image buffer.
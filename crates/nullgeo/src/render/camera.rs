use crate::geometry::{build_coframe_for, inner, make_null_covector, Metric, PhasePoint, Vec4};
use crate::spacetimes::Spacetime;
use crate::{Error, Result};

#[derive(Clone, Copy, Debug)]
pub struct CameraSpec {
    pub fov_deg: f64,
    pub res: (usize, usize),
    pub energy: f64,
    pub supersample: usize,
    pub supersample_max: usize,
    pub jitter: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct CameraPose {
    pub position: Vec4,
    pub look_at: [f64; 3],
    pub up: [f64; 3],
    pub velocity: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct Camera {
    pub spec: CameraSpec,
    pub pose: CameraPose,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalized(a: [f64; 3]) -> Option<[f64; 3]> {
    let len = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    if len < 1e-12 {
        return None;
    }
    Some([a[0] / len, a[1] / len, a[2] / len])
}

impl Camera {
    pub fn new(spec: CameraSpec, pose: CameraPose) -> Result<Self> {
        let (w, h) = spec.res;
        if w == 0 || h == 0 {
            return Err(Error::InvalidArg("resolution must be positive".into()));
        }
        if !(0.0 < spec.fov_deg && spec.fov_deg < 180.0) {
            return Err(Error::InvalidArg("fov_deg must be in (0, 180)".into()));
        }
        if spec.energy <= 0.0 {
            return Err(Error::InvalidArg("energy must be positive".into()));
        }
        if spec.supersample == 0 {
            return Err(Error::InvalidArg("supersample must be at least 1".into()));
        }
        if spec.supersample_max < spec.supersample {
            return Err(Error::InvalidArg(
                "supersample_max must be at least supersample".into(),
            ));
        }
        if pose.velocity.iter().any(|c| !c.is_finite()) {
            return Err(Error::InvalidArg("velocity must be finite".into()));
        }
        let camera = Self { spec, pose };
        camera.view_basis([pose.position[1], pose.position[2], pose.position[3]])?;
        Ok(camera)
    }

    fn view_basis(&self, position: [f64; 3]) -> Result<[[f64; 3]; 3]> {
        let forward = normalized(sub(self.pose.look_at, position))
            .ok_or_else(|| Error::InvalidArg("look_at coincides with position".into()))?;
        let right = normalized(cross(forward, self.pose.up))
            .ok_or_else(|| Error::InvalidArg("up is parallel to the view direction".into()))?;
        let up = cross(right, forward);
        Ok([forward, right, up])
    }

    fn scales(&self) -> (f64, f64) {
        let (w, h) = self.spec.res;
        let aspect = w as f64 / h as f64;
        let scale_u = (0.5 * self.spec.fov_deg.to_radians()).tan();
        (scale_u, scale_u / aspect)
    }

    pub fn pixel_directions(&self) -> Vec<[f64; 3]> {
        self.pixel_directions_at((0.5, 0.5))
    }

    pub fn pixel_directions_at(&self, subpixel: (f64, f64)) -> Vec<[f64; 3]> {
        let (w, h) = self.spec.res;
        let (scale_u, scale_v) = self.scales();
        (0..w * h)
            .map(|p| pixel_direction(scale_u, scale_v, w, h, p % w, p / w, subpixel))
            .collect()
    }

    pub fn subpixel_offsets(&self) -> Vec<(f64, f64)> {
        self.subpixel_offsets_for(self.spec.supersample)
    }

    pub fn subpixel_offsets_for(&self, n: usize) -> Vec<(f64, f64)> {
        (0..n * n)
            .map(|k| {
                let (cell_x, cell_y) = (k % n, k / n);
                let (jx, jy) = if self.spec.jitter {
                    (
                        radical_inverse(k + 1, 2) - 0.5,
                        radical_inverse(k + 1, 3) - 0.5,
                    )
                } else {
                    (0.0, 0.0)
                };
                (
                    (cell_x as f64 + 0.5 + jx) / n as f64,
                    (cell_y as f64 + 0.5 + jy) / n as f64,
                )
            })
            .collect()
    }

    pub fn observer_four_velocity<M: Metric + ?Sized>(&self, m: &M) -> Result<Vec4> {
        let [vx, vy, vz] = self.pose.velocity;
        let u = Vec4::new(1.0, vx, vy, vz);
        let len_sq = inner(&m.g(&self.pose.position), &u, &u);
        if len_sq >= -1e-12 {
            return Err(Error::NonTimelikeObserver(self.pose.position));
        }
        Ok(u / (-len_sq).sqrt())
    }

    pub fn pixel_rays<S: Spacetime + ?Sized>(&self, s: &S) -> Result<Vec<PhasePoint>> {
        self.pixel_rays_at(s, (0.5, 0.5))
    }

    pub fn pixel_rays_at<S: Spacetime + ?Sized>(
        &self,
        s: &S,
        subpixel: (f64, f64),
    ) -> Result<Vec<PhasePoint>> {
        let (w, h) = self.spec.res;
        let generator = self.ray_generator(s)?;
        Ok((0..w * h).map(|p| generator(p, subpixel)).collect())
    }

    pub fn ray_generator<S: Spacetime + ?Sized>(
        &self,
        s: &S,
    ) -> Result<impl Fn(usize, (f64, f64)) -> PhasePoint> {
        let [forward, right, up] = self.view_basis(s.embed(&self.pose.position))?;
        let observer = self.observer_four_velocity(s)?;
        let seed = |d: [f64; 3]| s.lift_direction(&self.pose.position, d);
        let coframe = build_coframe_for(
            s,
            &self.pose.position,
            &observer,
            [seed(forward), seed(right), seed(up)],
        )?;

        let position = self.pose.position;
        let energy = self.spec.energy;
        let (w, h) = self.spec.res;
        let (scale_u, scale_v) = self.scales();
        Ok(move |pixel: usize, subpixel: (f64, f64)| {
            let [f, r, u] = pixel_direction(scale_u, scale_v, w, h, pixel % w, pixel / w, subpixel);
            let arriving = make_null_covector(&coframe, [-f, -r, -u], energy);
            PhasePoint {
                x: position,
                p: -arriving,
            }
        })
    }
}

fn pixel_direction(
    scale_u: f64,
    scale_v: f64,
    w: usize,
    h: usize,
    i: usize,
    j: usize,
    subpixel: (f64, f64),
) -> [f64; 3] {
    let u = (2.0 * ((i as f64 + subpixel.0) / w as f64) - 1.0) * scale_u;
    let v = (1.0 - 2.0 * ((j as f64 + subpixel.1) / h as f64)) * scale_v;
    let inv_norm = 1.0 / (1.0 + u * u + v * v).sqrt();
    [inv_norm, u * inv_norm, v * inv_norm]
}

fn radical_inverse(mut index: usize, base: usize) -> f64 {
    let mut result = 0.0;
    let mut denom = 1.0;
    while index > 0 {
        denom *= base as f64;
        result += (index % base) as f64 / denom;
        index /= base;
    }
    result
}

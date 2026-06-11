use crate::frame::{build_coframe_for, make_null_covector, metric_dot};
use crate::metric::{Metric, State4, Vec4};
use crate::spacetime::Spacetime;
use crate::{Error, Result};

#[derive(Clone, Copy, Debug)]
pub struct CameraSpec {
    pub fov_deg: f64,
    pub res: (usize, usize),
    pub energy: f64,
    pub supersample: usize,
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

    pub fn pixel_directions(&self) -> Vec<[f64; 3]> {
        self.pixel_directions_at((0.5, 0.5))
    }

    pub fn pixel_directions_at(&self, subpixel: (f64, f64)) -> Vec<[f64; 3]> {
        let (w, h) = self.spec.res;
        let aspect = w as f64 / h as f64;
        let scale_u = (0.5 * self.spec.fov_deg.to_radians()).tan();
        let scale_v = scale_u / aspect;

        let mut dirs = Vec::with_capacity(w * h);
        for j in 0..h {
            let v = (1.0 - 2.0 * ((j as f64 + subpixel.1) / h as f64)) * scale_v;
            for i in 0..w {
                let u = (2.0 * ((i as f64 + subpixel.0) / w as f64) - 1.0) * scale_u;
                let inv_norm = 1.0 / (1.0 + u * u + v * v).sqrt();
                dirs.push([inv_norm, u * inv_norm, v * inv_norm]);
            }
        }
        dirs
    }

    pub fn subpixel_offsets(&self) -> Vec<(f64, f64)> {
        let n = self.spec.supersample;
        let centered = |k: usize| (k as f64 + 0.5) / n as f64;
        (0..n * n)
            .map(|k| (centered(k % n), centered(k / n)))
            .collect()
    }

    pub fn observer_four_velocity<M: Metric + ?Sized>(&self, m: &M) -> Result<Vec4> {
        let [vx, vy, vz] = self.pose.velocity;
        let u = Vec4::new(1.0, vx, vy, vz);
        let len_sq = metric_dot(&m.g(&self.pose.position), &u, &u);
        if len_sq >= -1e-12 {
            return Err(Error::NonTimelikeObserver(self.pose.position));
        }
        Ok(u / (-len_sq).sqrt())
    }

    pub fn pixel_rays<S: Spacetime + ?Sized>(&self, s: &S) -> Result<Vec<State4>> {
        self.pixel_rays_at(s, (0.5, 0.5))
    }

    pub fn pixel_rays_at<S: Spacetime + ?Sized>(
        &self,
        s: &S,
        subpixel: (f64, f64),
    ) -> Result<Vec<State4>> {
        let [forward, right, up] = self.view_basis(s.cartesian_position(&self.pose.position))?;
        let observer = self.observer_four_velocity(s)?;
        let seed = |d: [f64; 3]| s.chart_direction(&self.pose.position, d);
        let coframe = build_coframe_for(
            s,
            &self.pose.position,
            &observer,
            [seed(forward), seed(right), seed(up)],
        )?;

        let energy = self.spec.energy;
        let rays = self
            .pixel_directions_at(subpixel)
            .into_iter()
            .map(|[f, r, u]| {
                let arriving = make_null_covector(&coframe, [-f, -r, -u], energy);
                State4 {
                    x: self.pose.position,
                    p: -arriving,
                }
            })
            .collect();
        Ok(rays)
    }
}

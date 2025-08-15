use super::metric::Vec4;

#[derive(Clone, Copy, Debug)]
pub struct CameraSpec {
    pub fov_deg: f64,
    pub res: (usize, usize),
    pub energy: f64,
}

#[derive(Clone, Debug)]
pub struct Camera {
    pub spec: CameraSpec,
}

impl Camera {
    pub fn from_spec(spec: CameraSpec) -> Self { Self { spec } }

    pub fn generate_rays(&self) -> Vec<Vec4> {
        let (w, h) = self.spec.res;
        let mut ps = Vec::with_capacity(w * h);
        let fov_rad = self.spec.fov_deg.to_radians();
        let half = (fov_rad * 0.5).tan();
        let e = self.spec.energy.max(1e-12);

        for j in 0..h {
            for i in 0..w {
                let u = (2.0 * (i as f64 + 0.5) / w as f64 - 1.0) * half;   
                let v = (1.0 - 2.0 * (j as f64 + 0.5) / h as f64) * half;   

                let sx = 1.0;
                let sy = u;
                let sz = v;
                let norm = (sx*sx + sy*sy + sz*sz).sqrt();
                let nx = sx / norm;
                let ny = sy / norm;
                let nz = sz / norm;

                ps.push(Vec4::new(-e, e*nx, e*ny, e*nz));
            }
        }
        ps
    }
}

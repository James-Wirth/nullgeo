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

    pub fn generate_directions(&self) -> Vec<[f64; 3]> {
        let (w, h) = self.spec.res;
        assert!(w > 0 && h > 0, "resolution must be > 0");

        let fov = self.spec.fov_deg;
        assert!((0.0..180.0).contains(&fov), "fov_deg must be in (0, 180)");
        let half = (0.5 * fov.to_radians()).tan();

        let aspect = w as f64 / h as f64;
        let inv_w = 1.0 / w as f64;
        let inv_h = 1.0 / h as f64;

        let scale_u = half;
        let scale_v = half / aspect;

        let mut dirs = Vec::with_capacity(w * h);
        for j in 0..h {
            let ndc_v = 1.0 - 2.0 * ((j as f64 + 0.5) * inv_h);
            for i in 0..w {
                let ndc_u = 2.0 * ((i as f64 + 0.5) * inv_w) - 1.0;

                let u = ndc_u * scale_u;
                let v = ndc_v * scale_v;

                let sx = 1.0;
                let sy = u;
                let sz = v;

                let inv_norm = 1.0 / (sx*sx + sy*sy + sz*sz).sqrt();
                dirs.push([sx * inv_norm, sy * inv_norm, sz * inv_norm]);
            }
        }
        dirs
    }

    pub fn generate_rays(&self) -> Vec<Vec4> {
        let e = self.spec.energy.max(1e-12);
        self.generate_directions()
            .into_iter()
            .map(|[nx, ny, nz]| Vec4::new(-e, e*nx, e*ny, e*nz))
            .collect()
    }
}

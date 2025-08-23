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
    
        let fov = self.spec.fov_deg;
        assert!((0.0..180.0).contains(&fov), "fov_deg must be in (0, 180)");
        let fov_rad = fov.to_radians();
        let half = (0.5 * fov_rad).tan();
    
        let aspect = w as f64 / h as f64;
        let inv_w = 1.0 / w as f64;
        let inv_h = 1.0 / h as f64;
    
        let e = self.spec.energy.max(1e-12);
    
        for j in 0..h {
            let ndc_v = 1.0 - 2.0 * ((j as f64 + 0.5) * inv_h);
            for i in 0..w {
                let ndc_u = 2.0 * ((i as f64 + 0.5) * inv_w) - 1.0;
    
                let u = ndc_u * half;
                let v = ndc_v * (half / aspect);
    
                let sx = 1.0;
                let sy = u;
                let sz = v;
    
                let inv_norm = 1.0 / (sx*sx + sy*sy + sz*sz).sqrt();
                let nx = sx * inv_norm;
                let ny = sy * inv_norm;
                let nz = sz * inv_norm;
    
                ps.push(Vec4::new(-e, e*nx, e*ny, e*nz));
            }
        }
        ps
    }    
}

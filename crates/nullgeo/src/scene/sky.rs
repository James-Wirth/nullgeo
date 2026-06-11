use std::f64::consts::PI;

use crate::{Error, Result};

const CHECKER_LIGHT: [f32; 3] = [0.85, 0.85, 0.85];
const CHECKER_DARK: [f32; 3] = [0.05, 0.05, 0.05];

#[derive(Debug, Clone)]
pub enum SkyMap {
    Uniform([f32; 3]),
    Checker { angular_size_deg: f64 },
    Equirect(EquirectImage),
}

impl SkyMap {
    pub fn checker(angular_size_deg: f64) -> Result<Self> {
        let half_cells = 180.0 / angular_size_deg;
        let seamless = angular_size_deg > 0.0
            && half_cells >= 1.0
            && (half_cells - half_cells.round()).abs() < 1e-9;
        if !seamless {
            return Err(Error::InvalidArg(format!(
                "checker cells must tile the sphere seamlessly: \
                 360/{angular_size_deg} is not an even integer"
            )));
        }
        Ok(SkyMap::Checker { angular_size_deg })
    }

    pub fn sample(&self, dir: [f64; 3]) -> [f32; 3] {
        let theta = dir[2].clamp(-1.0, 1.0).acos();
        let phi = dir[1].atan2(dir[0]);
        match self {
            SkyMap::Uniform(color) => *color,
            SkyMap::Checker { angular_size_deg } => {
                let cell = angular_size_deg.to_radians();
                let index = (theta / cell).floor() as i64 + (phi / cell).floor() as i64;
                if index.rem_euclid(2) == 0 {
                    CHECKER_LIGHT
                } else {
                    CHECKER_DARK
                }
            }
            SkyMap::Equirect(image) => image.sample((phi + PI) / (2.0 * PI), theta / PI),
        }
    }
}

#[derive(Debug, Clone)]
pub struct EquirectImage {
    width: usize,
    height: usize,
    data: Vec<[f32; 3]>,
}

impl EquirectImage {
    pub fn new(width: usize, height: usize, data: Vec<[f32; 3]>) -> Result<Self> {
        if width == 0 || height == 0 || data.len() != width * height {
            return Err(Error::InvalidArg(format!(
                "equirect image of {}x{} needs {} texels, got {}",
                width,
                height,
                width * height,
                data.len()
            )));
        }
        Ok(Self {
            width,
            height,
            data,
        })
    }

    pub fn sample(&self, u: f64, v: f64) -> [f32; 3] {
        let x = u * self.width as f64 - 0.5;
        let y = (v * self.height as f64 - 0.5).clamp(0.0, (self.height - 1) as f64);

        let x0 = x.floor();
        let fx = (x - x0) as f32;
        let wrap = |k: i64| k.rem_euclid(self.width as i64) as usize;
        let (xa, xb) = (wrap(x0 as i64), wrap(x0 as i64 + 1));

        let y0 = y.floor() as usize;
        let y1 = (y0 + 1).min(self.height - 1);
        let fy = (y - y0 as f64) as f32;

        let texel = |xi: usize, yi: usize| self.data[yi * self.width + xi];
        let lerp = |a: [f32; 3], b: [f32; 3], t: f32| {
            [
                a[0] + (b[0] - a[0]) * t,
                a[1] + (b[1] - a[1]) * t,
                a[2] + (b[2] - a[2]) * t,
            ]
        };
        let top = lerp(texel(xa, y0), texel(xb, y0), fx);
        let bottom = lerp(texel(xa, y1), texel(xb, y1), fx);
        lerp(top, bottom, fy)
    }
}

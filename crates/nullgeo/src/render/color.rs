const LAMBDA_MIN_NM: f64 = 380.0;
const LAMBDA_STEP_NM: f64 = 5.0;
const LAMBDA_SAMPLES: usize = 81;
const PLANCK_C2_NM_K: f64 = 1.438777e7;

fn cmf_lobe(lambda: f64, mu: f64, sigma_lo: f64, sigma_hi: f64) -> f64 {
    let sigma = if lambda < mu { sigma_lo } else { sigma_hi };
    let t = (lambda - mu) / sigma;
    (-0.5 * t * t).exp()
}

fn cie_xyz_bar(lambda: f64) -> [f64; 3] {
    let x = 1.056 * cmf_lobe(lambda, 599.8, 37.9, 31.0)
        + 0.362 * cmf_lobe(lambda, 442.0, 16.0, 26.7)
        - 0.065 * cmf_lobe(lambda, 501.1, 20.4, 26.2);
    let y =
        0.821 * cmf_lobe(lambda, 568.8, 46.9, 40.5) + 0.286 * cmf_lobe(lambda, 530.9, 16.3, 31.1);
    let z =
        1.217 * cmf_lobe(lambda, 437.0, 11.8, 36.0) + 0.681 * cmf_lobe(lambda, 459.0, 26.0, 13.8);
    [x, y, z]
}

fn planck_spectral_radiance(lambda_nm: f64, temperature: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    let x = (PLANCK_C2_NM_K / (lambda_nm * temperature)).exp() - 1.0;
    1.0 / (lambda_nm.powi(5) * x)
}

pub fn planck_xyz(temperature: f64) -> [f64; 3] {
    let mut xyz = [0.0; 3];
    for i in 0..LAMBDA_SAMPLES {
        let lambda = LAMBDA_MIN_NM + LAMBDA_STEP_NM * i as f64;
        let radiance = planck_spectral_radiance(lambda, temperature);
        let bar = cie_xyz_bar(lambda);
        for (acc, b) in xyz.iter_mut().zip(bar) {
            *acc += radiance * b;
        }
    }
    let scale = LAMBDA_STEP_NM * 1e15;
    xyz.map(|c| c * scale)
}

pub fn xyz_to_linear_srgb(xyz: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = xyz;
    let r = 3.2404542 * x - 1.5371385 * y - 0.4985314 * z;
    let g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z;
    let b = 0.0556434 * x - 0.2040259 * y + 1.0572252 * z;
    [r.max(0.0), g.max(0.0), b.max(0.0)]
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ToneCurve {
    #[default]
    Reinhard,
    Aces,
}

impl ToneCurve {
    fn apply(self, v: f32) -> f32 {
        match self {
            ToneCurve::Reinhard => v / (1.0 + v),
            ToneCurve::Aces => {
                ((v * (2.51 * v + 0.03)) / (v * (2.43 * v + 0.59) + 0.14)).clamp(0.0, 1.0)
            }
        }
    }
}

pub fn tone_map_curve(image: &super::ImageF32, exposure: f32, curve: ToneCurve) -> Vec<[f32; 3]> {
    image
        .data
        .iter()
        .map(|c| {
            let mut out = [0.0f32; 3];
            for (display, &channel) in out.iter_mut().zip(c) {
                let v = (channel * exposure).max(0.0);
                *display = curve.apply(v).powf(1.0 / 2.2);
            }
            out
        })
        .collect()
}

pub fn quantize8(display: &[[f32; 3]]) -> Vec<[u8; 3]> {
    display
        .iter()
        .map(|c| c.map(|v| (v * 255.0 + 0.5).clamp(0.0, 255.0) as u8))
        .collect()
}

pub fn quantize16(display: &[[f32; 3]]) -> Vec<[u16; 3]> {
    display
        .iter()
        .map(|c| c.map(|v| (v * 65535.0 + 0.5).clamp(0.0, 65535.0) as u16))
        .collect()
}

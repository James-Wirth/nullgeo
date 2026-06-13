#[derive(Debug, Clone, Copy)]
pub struct Disk {
    pub r_in: f64,
    pub r_out: f64,
    pub model: DiskModel,
}

#[derive(Debug, Clone, Copy)]
pub enum DiskModel {
    Stylized {
        emissivity_index: f64,
        g_power: f64,
    },
    Blackbody {
        t_in: f64,
        doppler_beaming: bool,
        redshift_color: bool,
        optical_depth: f64,
        aspect_ratio: f64,
        density_index: f64,
        edge_taper: f64,
    },
}

impl Disk {
    pub fn stylized(r_out: f64) -> Self {
        Self {
            r_in: 0.0,
            r_out,
            model: DiskModel::Stylized {
                emissivity_index: 2.0,
                g_power: 3.0,
            },
        }
    }

    pub fn blackbody(r_out: f64, t_in: f64) -> Self {
        Self {
            r_in: 0.0,
            r_out,
            model: DiskModel::Blackbody {
                t_in,
                doppler_beaming: true,
                redshift_color: true,
                optical_depth: f64::INFINITY,
                aspect_ratio: 0.0,
                density_index: 3.0,
                edge_taper: 0.0,
            },
        }
    }

    pub fn volume(&self, r_in: f64) -> DiskVolume {
        let (tau0, aspect_ratio, density_index, edge_taper) = match self.model {
            DiskModel::Stylized { .. } => (f64::INFINITY, 0.0, 3.0, 0.0),
            DiskModel::Blackbody {
                optical_depth,
                aspect_ratio,
                density_index,
                edge_taper,
                ..
            } => (optical_depth, aspect_ratio, density_index, edge_taper),
        };
        DiskVolume {
            r_in,
            r_out: self.r_out,
            tau0,
            aspect_ratio,
            density_index,
            w_out: edge_taper * (self.r_out - r_in),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct DiskVolume {
    pub r_in: f64,
    pub r_out: f64,
    pub tau0: f64,
    pub aspect_ratio: f64,
    pub density_index: f64,
    pub w_out: f64,
}

impl DiskVolume {
    pub fn scale_height(&self, r: f64) -> f64 {
        self.aspect_ratio * r
    }

    pub fn surface_density(&self, r: f64) -> f64 {
        if r <= self.r_in || r >= self.r_out {
            return 0.0;
        }
        let peak = self.surface_density_peak();
        if peak <= 0.0 {
            return 0.0;
        }
        surface_density_shape(self.r_in, self.density_index, r) * ramp(self.r_out - r, self.w_out)
            / peak
    }

    fn surface_density_peak(&self) -> f64 {
        let p = self.density_index;
        if p <= 0.0 {
            return surface_density_shape(self.r_in, p, self.r_out);
        }
        let r_peak = self.r_in * ((2.0 * p + 1.0) / (2.0 * p)).powi(2);
        surface_density_shape(self.r_in, p, r_peak)
    }

    pub fn tau_perp(&self, r: f64) -> f64 {
        let sigma = self.surface_density(r);
        if sigma == 0.0 {
            0.0
        } else {
            self.tau0 * sigma
        }
    }

    pub fn density_alpha(&self, r: f64, z: f64) -> f64 {
        let h = self.scale_height(r);
        let tau_perp = self.tau_perp(r);
        if h <= 0.0 || tau_perp == 0.0 {
            return 0.0;
        }
        let norm = 1.0 / ((2.0 * std::f64::consts::PI).sqrt() * h);
        tau_perp * norm * (-0.5 * (z / h) * (z / h)).exp()
    }

    pub fn tau_eff(&self, r: f64, mu: f64) -> f64 {
        self.tau_perp(r) / mu.abs()
    }
}

fn surface_density_shape(r_in: f64, density_index: f64, r: f64) -> f64 {
    if r <= r_in {
        return 0.0;
    }
    (1.0 - (r_in / r).sqrt()) * (r / r_in).powf(-density_index)
}

fn ramp(edge_distance: f64, width: f64) -> f64 {
    if width <= 0.0 {
        return if edge_distance > 0.0 { 1.0 } else { 0.0 };
    }
    let t = (edge_distance / width).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

pub fn shakura_sunyaev_temperature(t_in: f64, r_in: f64, r: f64) -> f64 {
    if r <= r_in {
        return 0.0;
    }
    t_in * (r / r_in).powf(-0.75) * (1.0 - (r_in / r).sqrt()).powf(0.25)
}

pub fn shakura_sunyaev_peak_radius(r_in: f64) -> f64 {
    r_in * 49.0 / 36.0
}

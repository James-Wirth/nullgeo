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
                edge_taper: 0.0,
            },
        }
    }

    pub fn volume(&self, r_in: f64) -> DiskVolume {
        let (tau0, aspect_ratio, edge_taper) = match self.model {
            DiskModel::Stylized { .. } => (f64::INFINITY, 0.0, 0.0),
            DiskModel::Blackbody {
                optical_depth,
                aspect_ratio,
                edge_taper,
                ..
            } => (optical_depth, aspect_ratio, edge_taper),
        };
        let width = edge_taper * (self.r_out - r_in);
        DiskVolume {
            r_in,
            r_out: self.r_out,
            tau0,
            aspect_ratio,
            w_out: width,
            w_in: width,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct DiskVolume {
    pub r_in: f64,
    pub r_out: f64,
    pub tau0: f64,
    pub aspect_ratio: f64,
    pub w_out: f64,
    pub w_in: f64,
}

impl DiskVolume {
    pub fn scale_height(&self, r: f64) -> f64 {
        self.aspect_ratio * r
    }

    pub fn taper(&self, r: f64) -> f64 {
        ramp(self.r_out - r, self.w_out) * ramp(r - self.r_in, self.w_in)
    }

    pub fn tau_perp(&self, r: f64) -> f64 {
        let taper = self.taper(r);
        if taper == 0.0 {
            0.0
        } else {
            self.tau0 * taper
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

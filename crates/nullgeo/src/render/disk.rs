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
            },
        }
    }
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

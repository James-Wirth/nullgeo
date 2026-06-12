#[derive(Debug, Clone, Copy)]
pub struct Disk {
    pub r_in: f64,
    pub r_out: f64,
    pub emissivity_index: f64,
    pub g_power: f64,
}

impl Disk {
    pub fn new(r_out: f64) -> Self {
        Self {
            r_in: 0.0,
            r_out,
            emissivity_index: 2.0,
            g_power: 3.0,
        }
    }
}

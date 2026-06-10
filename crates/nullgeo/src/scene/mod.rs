pub mod disk;
pub mod sky;

use crate::spacetime::SkySide;
use disk::Disk;
use sky::SkyMap;

#[derive(Debug, Clone)]
pub struct Scene {
    pub sky: SkyMap,
    pub sky_secondary: Option<SkyMap>,
    pub disk: Option<Disk>,
}

impl Scene {
    pub fn sky_for(&self, side: SkySide) -> &SkyMap {
        match side {
            SkySide::Primary => &self.sky,
            SkySide::Secondary => self.sky_secondary.as_ref().unwrap_or(&self.sky),
        }
    }
}

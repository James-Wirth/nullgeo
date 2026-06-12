use super::geometry::{GeometryBuffer, RayClass, RayInfo, RayOutcome};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MapQuantity {
    Classification,
    Redshift,
    ImageOrder,
    MinRadius,
    CoordTime,
    AffineLength,
    Steps,
    EscapeTheta,
    EscapePhi,
}

impl MapQuantity {
    pub fn is_categorical(self) -> bool {
        matches!(self, MapQuantity::Classification | MapQuantity::ImageOrder)
    }

    fn extract(self, info: &RayInfo) -> f64 {
        match self {
            MapQuantity::Classification => class_index(info.class()) as f64,
            MapQuantity::Redshift => match info.outcome {
                RayOutcome::DiskHit { g: Some(g), .. } => g,
                _ => f64::NAN,
            },
            MapQuantity::ImageOrder => info.stats.equatorial_crossings as f64,
            MapQuantity::MinRadius => info.stats.min_radius,
            MapQuantity::CoordTime => info.stats.coord_time,
            MapQuantity::AffineLength => info.stats.affine_length,
            MapQuantity::Steps => (info.stats.steps_accepted + info.stats.steps_rejected) as f64,
            MapQuantity::EscapeTheta => match info.outcome {
                RayOutcome::Escaped { dir, .. } => dir[2].clamp(-1.0, 1.0).acos(),
                _ => f64::NAN,
            },
            MapQuantity::EscapePhi => match info.outcome {
                RayOutcome::Escaped { dir, .. } => dir[1].atan2(dir[0]),
                _ => f64::NAN,
            },
        }
    }
}

#[derive(Debug, Clone)]
pub struct MapField {
    pub width: usize,
    pub height: usize,
    pub quantity: MapQuantity,
    pub values: Vec<f64>,
}

impl MapField {
    pub fn finite_range(&self) -> Option<(f64, f64)> {
        let mut range: Option<(f64, f64)> = None;
        for &v in &self.values {
            if v.is_finite() {
                let (lo, hi) = range.unwrap_or((v, v));
                range = Some((lo.min(v), hi.max(v)));
            }
        }
        range
    }

    pub fn display_range(&self) -> Option<(f64, f64)> {
        let (lo, hi) = self.finite_range()?;
        match self.quantity {
            MapQuantity::Redshift => {
                let spread = (hi - 1.0).abs().max((lo - 1.0).abs()).max(f64::EPSILON);
                Some((1.0 - spread, 1.0 + spread))
            }
            _ => Some((lo, hi)),
        }
    }
}

pub fn shade_map(buffer: &GeometryBuffer, quantity: MapQuantity) -> MapField {
    let values = (0..buffer.width * buffer.height)
        .map(|pixel| quantity.extract(buffer.primary(pixel)))
        .collect();
    MapField {
        width: buffer.width,
        height: buffer.height,
        quantity,
        values,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Colormap {
    Viridis,
    Diverging,
}

impl Colormap {
    pub fn default_for(quantity: MapQuantity) -> Self {
        match quantity {
            MapQuantity::Redshift => Colormap::Diverging,
            _ => Colormap::Viridis,
        }
    }

    fn lookup(self, t: f64) -> [u8; 3] {
        let anchors: &[[f64; 3]] = match self {
            Colormap::Viridis => &VIRIDIS,
            Colormap::Diverging => &DIVERGING,
        };
        let scaled = t.clamp(0.0, 1.0) * (anchors.len() - 1) as f64;
        let lo = scaled.floor() as usize;
        let hi = (lo + 1).min(anchors.len() - 1);
        let frac = scaled - lo as f64;
        let mut out = [0u8; 3];
        for (byte, (a, b)) in out.iter_mut().zip(anchors[lo].iter().zip(&anchors[hi])) {
            *byte = (a + (b - a) * frac).round().clamp(0.0, 255.0) as u8;
        }
        out
    }
}

const VIRIDIS: [[f64; 3]; 10] = [
    [68.0, 1.0, 84.0],
    [72.0, 40.0, 120.0],
    [62.0, 73.0, 137.0],
    [49.0, 104.0, 142.0],
    [38.0, 130.0, 142.0],
    [31.0, 158.0, 137.0],
    [53.0, 183.0, 121.0],
    [109.0, 205.0, 89.0],
    [180.0, 222.0, 44.0],
    [253.0, 231.0, 37.0],
];

const DIVERGING: [[f64; 3]; 5] = [
    [59.0, 76.0, 192.0],
    [124.0, 159.0, 249.0],
    [220.0, 220.0, 220.0],
    [245.0, 156.0, 125.0],
    [180.0, 4.0, 38.0],
];

const UNDEFINED_COLOR: [u8; 3] = [0, 0, 0];

const CLASS_COLORS: [[u8; 3]; 6] = [
    [0, 0, 0],
    [70, 130, 180],
    [230, 140, 30],
    [255, 255, 255],
    [220, 20, 60],
    [255, 0, 255],
];

const ORDER_COLORS: [[u8; 3]; 8] = [
    [25, 25, 35],
    [230, 159, 0],
    [86, 180, 233],
    [0, 158, 115],
    [240, 228, 66],
    [0, 114, 178],
    [213, 94, 0],
    [204, 121, 167],
];

fn class_index(class: RayClass) -> usize {
    match class {
        RayClass::Captured => 0,
        RayClass::EscapedPrimary => 1,
        RayClass::EscapedSecondary => 2,
        RayClass::Disk => 3,
        RayClass::MaxSteps => 4,
        RayClass::Stalled => 5,
    }
}

pub fn class_color(class: RayClass) -> [u8; 3] {
    CLASS_COLORS[class_index(class)]
}

pub fn colorize(field: &MapField, colormap: Colormap) -> Vec<[u8; 3]> {
    let categorical_palette: Option<&[[u8; 3]]> = match field.quantity {
        MapQuantity::Classification => Some(&CLASS_COLORS),
        MapQuantity::ImageOrder => Some(&ORDER_COLORS),
        _ => None,
    };

    if let Some(palette) = categorical_palette {
        return field
            .values
            .iter()
            .map(|&v| {
                if v.is_finite() && v >= 0.0 {
                    palette[v as usize % palette.len()]
                } else {
                    UNDEFINED_COLOR
                }
            })
            .collect();
    }

    let range = field.display_range();
    field
        .values
        .iter()
        .map(|&v| match (v.is_finite(), range) {
            (true, Some((lo, hi))) => {
                let span = (hi - lo).max(f64::EPSILON);
                colormap.lookup((v - lo) / span)
            }
            _ => UNDEFINED_COLOR,
        })
        .collect()
}

use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_ppm_gray(
    path: &str,
    width: usize,
    height: usize,
    pixels: &[u8],
) -> std::io::Result<()> {
    if pixels.len() != width * height {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "pixel buffer size mismatch",
        ));
    }
    let mut f = BufWriter::new(File::create(path)?);
    writeln!(f, "P5")?;
    writeln!(f, "{} {}", width, height)?;
    writeln!(f, "255")?;
    f.write_all(pixels)?;
    Ok(())
}

pub fn write_ppm_rgb(
    path: &str,
    width: usize,
    height: usize,
    pixels: &[[u8; 3]],
) -> std::io::Result<()> {
    if pixels.len() != width * height {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "pixel buffer size mismatch",
        ));
    }
    let mut f = BufWriter::new(File::create(path)?);
    writeln!(f, "P6")?;
    writeln!(f, "{} {}", width, height)?;
    writeln!(f, "255")?;
    for px in pixels {
        f.write_all(px)?;
    }
    Ok(())
}

pub fn write_pfm_gray(
    path: &str,
    width: usize,
    height: usize,
    values: &[f64],
) -> std::io::Result<()> {
    if values.len() != width * height {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "value buffer size mismatch",
        ));
    }
    let mut f = BufWriter::new(File::create(path)?);
    write!(f, "Pf\n{} {}\n-1.0\n", width, height)?;
    for row in (0..height).rev() {
        for &v in &values[row * width..(row + 1) * width] {
            f.write_all(&(v as f32).to_le_bytes())?;
        }
    }
    Ok(())
}

pub fn write_pfm_rgb(
    path: &str,
    width: usize,
    height: usize,
    pixels: &[[f32; 3]],
) -> std::io::Result<()> {
    if pixels.len() != width * height {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "pixel buffer size mismatch",
        ));
    }
    let mut f = BufWriter::new(File::create(path)?);
    write!(f, "PF\n{} {}\n-1.0\n", width, height)?;
    for row in (0..height).rev() {
        for px in &pixels[row * width..(row + 1) * width] {
            for channel in px {
                f.write_all(&channel.to_le_bytes())?;
            }
        }
    }
    Ok(())
}

pub fn write_csv_matrix(
    path: &str,
    width: usize,
    height: usize,
    values: &[f64],
) -> std::io::Result<()> {
    if values.len() != width * height {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "value buffer size mismatch",
        ));
    }
    let mut f = BufWriter::new(File::create(path)?);
    for row in values.chunks(width) {
        let line: Vec<String> = row.iter().map(|v| v.to_string()).collect();
        writeln!(f, "{}", line.join(","))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pfm_gray_layout_is_bottom_up_little_endian() {
        let dir = std::env::temp_dir().join("nullgeo_io_pfm_gray");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("map.pfm");
        let path = path.to_str().unwrap();
        write_pfm_gray(path, 2, 2, &[1.0, 2.0, 3.0, 4.0]).unwrap();
        let bytes = std::fs::read(path).unwrap();
        let header = b"Pf\n2 2\n-1.0\n";
        assert_eq!(&bytes[..header.len()], header);
        let mut floats = bytes[header.len()..]
            .chunks(4)
            .map(|c| f32::from_le_bytes(c.try_into().unwrap()));
        assert_eq!(floats.next(), Some(3.0));
        assert_eq!(floats.next(), Some(4.0));
        assert_eq!(floats.next(), Some(1.0));
        assert_eq!(floats.next(), Some(2.0));
        assert_eq!(floats.next(), None);
    }

    #[test]
    fn csv_matrix_rows_match_image_rows() {
        let dir = std::env::temp_dir().join("nullgeo_io_csv");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("map.csv");
        let path = path.to_str().unwrap();
        write_csv_matrix(path, 3, 2, &[1.0, 2.5, f64::NAN, 4.0, 5.0, 6.0]).unwrap();
        let text = std::fs::read_to_string(path).unwrap();
        assert_eq!(text, "1,2.5,NaN\n4,5,6\n");
    }
}

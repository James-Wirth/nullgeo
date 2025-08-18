use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_ppm_gray(path: &str, width: usize, height: usize, pixels: &[u8]) -> std::io::Result<()> {
    if pixels.len() != width*height {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "pixel buffer size mismatch"));
    }
    let mut f = BufWriter::new(File::create(path)?);
    writeln!(f, "P5")?;
    writeln!(f, "{} {}", width, height)?;
    writeln!(f, "255")?;
    f.write_all(pixels)?;
    Ok(())
}

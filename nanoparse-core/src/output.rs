//! Output serialization module

use anyhow::{Context, Result};
use serde::Serialize;
use std::fs::File;
use std::io::{self, Write};

/// Write output to file or stdout
pub fn write_output<T: Serialize>(data: &T, path: &str) -> Result<()> {
    let json = serde_json::to_string_pretty(data).context("Failed to serialize JSON")?;

    if path == "-" {
        io::stdout()
            .write_all(json.as_bytes())
            .context("Failed to write to stdout")?;
        println!(); // Newline at end
    } else {
        let mut file = File::create(path).context("Failed to create output file")?;
        file.write_all(json.as_bytes())
            .context("Failed to write to file")?;
        file.write_all(b"\n")?;
    }

    Ok(())
}

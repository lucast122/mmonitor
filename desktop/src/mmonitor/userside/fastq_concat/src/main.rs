use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use clap::Parser;

#[derive(Parser)]
#[command(about = "Fast parallel FASTQ file concatenation")]
struct Args {
    /// Input files to concatenate
    #[arg(required = true)]
    input_files: Vec<PathBuf>,

    /// Output file
    #[arg(short, long)]
    output: PathBuf,

    /// Buffer size in MB
    #[arg(short, long, default_value = "32")]
    buffer_size: usize,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let buffer_size = args.buffer_size * 1024 * 1024; // Convert MB to bytes

    // Create output file with gzip compression
    let output_file = File::create(&args.output)?;
    let mut encoder = GzEncoder::new(output_file, Compression::default());

    // Process files in parallel and collect their contents
    args.input_files.par_iter().try_for_each(|input_path| {
        let mut buffer = vec![0; buffer_size];
        let input_file = File::open(input_path)?;
        let mut reader = if input_path.extension().map_or(false, |ext| ext == "gz") {
            Box::new(GzDecoder::new(BufReader::new(input_file))) as Box<dyn Read>
        } else {
            Box::new(BufReader::new(input_file)) as Box<dyn Read>
        };

        loop {
            let bytes_read = reader.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }
            encoder.write_all(&buffer[..bytes_read])?;
        }
        Ok(()) as io::Result<()>
    })?;

    encoder.finish()?;
    Ok(())
}

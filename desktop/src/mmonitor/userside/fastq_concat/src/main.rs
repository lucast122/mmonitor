use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use clap::Parser;
use crossbeam::channel;

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

fn process_file(path: &PathBuf, buffer_size: usize) -> io::Result<Vec<u8>> {
    let mut buffer = vec![0; buffer_size];
    let input_file = File::open(path)?;
    let mut reader = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(GzDecoder::new(BufReader::new(input_file))) as Box<dyn Read>
    } else {
        Box::new(BufReader::new(input_file)) as Box<dyn Read>
    };

    let mut contents = Vec::new();
    loop {
        let bytes_read = reader.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        contents.extend_from_slice(&buffer[..bytes_read]);
    }
    Ok(contents)
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let buffer_size = args.buffer_size * 1024 * 1024; // Convert MB to bytes

    // Create output file with gzip compression
    let output_file = File::create(&args.output)?;
    let mut encoder = GzEncoder::new(output_file, Compression::default());

    // Create a channel to send file contents
    let (tx, rx) = channel::unbounded();

    // Process files in parallel
    args.input_files.par_iter().try_for_each(|input_path| {
        let contents = process_file(input_path, buffer_size)?;
        tx.send(contents).unwrap();
        Ok(()) as io::Result<()>
    })?;

    // Drop the sender to signal we're done
    drop(tx);

    // Write all received contents to the output file
    for contents in rx {
        encoder.write_all(&contents)?;
    }

    encoder.finish()?;
    Ok(())
}

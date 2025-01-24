use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use clap::Parser;
use std::thread;

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

    /// Number of threads to use (defaults to number of CPU cores)
    #[arg(short, long)]
    threads: Option<usize>,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let buffer_size = args.buffer_size * 1024 * 1024; // Convert MB to bytes

    // Set number of threads if specified, otherwise use number of CPU cores
    let thread_count = args.threads.unwrap_or_else(|| thread::available_parallelism().map_or(1, |p| p.get()));
    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .unwrap();
    println!("Using {} threads for parallel processing", thread_count);

    // First collect all file contents in parallel
    let all_contents: Vec<_> = args.input_files
        .par_iter()
        .map(|input_path| -> io::Result<Vec<u8>> {
            println!("Processing file: {:?}", input_path);
            let mut buffer = vec![0; buffer_size];
            let input_file = File::open(input_path)?;
            let mut reader = if input_path.extension().map_or(false, |ext| ext == "gz") {
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
            println!("Finished processing file: {:?}", input_path);
            Ok(contents)
        })
        .collect::<io::Result<Vec<_>>>()?;

    println!("Writing concatenated output to {:?}", args.output);

    // Write all contents sequentially
    let output_file = File::create(&args.output)?;
    let mut encoder = GzEncoder::new(output_file, Compression::default());
    for contents in all_contents {
        encoder.write_all(&contents)?;
    }

    encoder.finish()?;
    println!("Successfully finished concatenation");
    Ok(())
}

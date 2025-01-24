use std::fs::File;
use std::io::{self, BufReader, Read, Write, BufWriter};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use clap::Parser;
use std::thread;
use std::sync::mpsc;

const CHUNK_SIZE: usize = 1024 * 1024; // 1MB chunks

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

fn process_file_chunks(path: &PathBuf, sender: mpsc::Sender<Vec<u8>>) -> io::Result<()> {
    let input_file = File::open(path)?;
    let mut reader = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(GzDecoder::new(BufReader::new(input_file))) as Box<dyn Read>
    } else {
        Box::new(BufReader::new(input_file)) as Box<dyn Read>
    };

    let mut buffer = vec![0; CHUNK_SIZE];
    loop {
        let bytes_read = reader.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        sender.send(buffer[..bytes_read].to_vec()).unwrap();
    }
    Ok(())
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

    // Create output file
    let output_file = File::create(&args.output)?;
    let mut writer = BufWriter::new(GzEncoder::new(output_file, Compression::default()));

    // Create channel for sending chunks
    let (tx, rx) = mpsc::channel();

    // Process files in parallel
    args.input_files.par_iter().try_for_each(|input_path| {
        println!("Processing file: {:?}", input_path);
        let tx = tx.clone();
        process_file_chunks(input_path, tx)
    })?;

    // Drop original sender so receiver knows when we're done
    drop(tx);

    // Write chunks as they come in
    let mut total_bytes = 0;
    for chunk in rx {
        writer.write_all(&chunk)?;
        total_bytes += chunk.len();
        if total_bytes % (100 * 1024 * 1024) == 0 {  // Print progress every 100MB
            println!("Processed {} MB", total_bytes / (1024 * 1024));
        }
    }

    writer.flush()?;
    println!("Successfully finished concatenation. Total size: {} MB", total_bytes / (1024 * 1024));
    Ok(())
}

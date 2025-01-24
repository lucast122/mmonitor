use std::fs::File;
use std::io::{self, BufReader, Read, Write, BufWriter, Seek, SeekFrom};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use clap::Parser;
use std::thread;
use std::sync::{Arc, Mutex};
use std::time::Instant;

const CHUNK_SIZE: usize = 64 * 1024 * 1024; // 64MB chunks for better parallelization

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

// Read a chunk of data from a file at a specific offset
fn read_chunk(path: &PathBuf, offset: u64, chunk_size: usize) -> io::Result<Vec<u8>> {
    let mut file = File::open(path)?;
    file.seek(SeekFrom::Start(offset))?;
    let mut reader = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(GzDecoder::new(BufReader::new(file))) as Box<dyn Read>
    } else {
        Box::new(BufReader::new(file)) as Box<dyn Read>
    };

    let mut chunk = Vec::with_capacity(chunk_size);
    reader.take(chunk_size as u64).read_to_end(&mut chunk)?;
    Ok(chunk)
}

fn main() -> io::Result<()> {
    let start_time = Instant::now();
    let args = Args::parse();

    // Set number of threads if specified, otherwise use number of CPU cores
    let thread_count = args.threads.unwrap_or_else(|| thread::available_parallelism().map_or(1, |p| p.get()));
    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .unwrap();
    println!("Using {} threads for parallel processing", thread_count);

    // Create output file and wrap in Arc<Mutex> for thread-safe access
    let output_file = File::create(&args.output)?;
    let writer = Arc::new(Mutex::new(BufWriter::new(GzEncoder::new(output_file, Compression::default()))));
    let total_bytes = Arc::new(Mutex::new(0u64));

    // Process each file
    args.input_files.par_iter().try_for_each(|input_path| -> io::Result<()> {
        let file_size = std::fs::metadata(input_path)?.len();
        let num_chunks = (file_size + CHUNK_SIZE as u64 - 1) / CHUNK_SIZE as u64;
        println!("Processing file {:?} ({} MB) in {} chunks", 
                input_path, file_size / (1024 * 1024), num_chunks);

        // Process chunks in parallel
        (0..num_chunks).into_par_iter().try_for_each(|chunk_idx| -> io::Result<()> {
            let offset = chunk_idx as u64 * CHUNK_SIZE as u64;
            let chunk = read_chunk(input_path, offset, CHUNK_SIZE)?;
            
            // Write chunk to output file
            let writer = Arc::clone(&writer);
            let mut writer = writer.lock().unwrap();
            writer.write_all(&chunk)?;
            
            // Update progress
            let total_bytes = Arc::clone(&total_bytes);
            let mut total = total_bytes.lock().unwrap();
            *total += chunk.len() as u64;
            if *total % (500 * 1024 * 1024) == 0 {
                println!("Processed {} GB in {:.1} seconds", 
                        *total / (1024 * 1024 * 1024),
                        start_time.elapsed().as_secs_f32());
            }
            Ok(())
        })?;
        Ok(())
    })?;

    // Finish writing and print final statistics
    let writer = Arc::try_unwrap(writer).unwrap().into_inner().unwrap();
    writer.into_inner()?.finish()?;
    
    let total_bytes = *total_bytes.lock().unwrap();
    let duration = start_time.elapsed();
    println!("Successfully finished concatenation:");
    println!("Total size: {:.2} GB", total_bytes as f64 / (1024.0 * 1024.0 * 1024.0));
    println!("Time taken: {:.1} seconds", duration.as_secs_f32());
    println!("Speed: {:.1} MB/s", 
             (total_bytes as f64 / 1024.0 / 1024.0) / duration.as_secs_f64());
    Ok(())
}

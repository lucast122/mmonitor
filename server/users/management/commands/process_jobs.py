"""
Management command to process analysis jobs in the queue
"""
import time
import os
import subprocess
import logging
from django.core.management.base import BaseCommand
from django.utils import timezone
from users.models import AnalysisJob, SystemStatus

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'Process analysis jobs in the queue'
    
    def add_arguments(self, parser):
        parser.add_argument(
            '--daemon',
            action='store_true',
            help='Run as daemon, continuously processing jobs',
        )
        parser.add_argument(
            '--interval',
            type=int,
            default=10,
            help='Polling interval in seconds (default: 10)',
        )
    
    def handle(self, *args, **options):
        self.daemon = options['daemon']
        self.interval = options['interval']
        
        self.stdout.write(
            self.style.SUCCESS(f'Starting job processor (daemon: {self.daemon})')
        )
        
        if self.daemon:
            self.run_daemon()
        else:
            self.process_pending_jobs()
    
    def run_daemon(self):
        """Run continuously as a daemon"""
        self.stdout.write('Running in daemon mode...')
        
        while True:
            try:
                self.process_pending_jobs()
                time.sleep(self.interval)
            except KeyboardInterrupt:
                self.stdout.write('\nShutting down job processor...')
                break
            except Exception as e:
                logger.error(f"Error in job processor: {str(e)}")
                time.sleep(self.interval)
    
    def process_pending_jobs(self):
        """Process all pending jobs"""
        system_status = SystemStatus.get_current()
        
        # Update job counts
        system_status.update_job_counts()
        
        if system_status.maintenance_mode:
            return
        
        # Get pending jobs
        pending_jobs = AnalysisJob.objects.filter(
            status__in=['pending', 'queued']
        ).order_by('created_at')
        
        # Process jobs up to system capacity
        available_slots = system_status.max_concurrent_jobs - system_status.active_jobs
        
        for job in pending_jobs[:available_slots]:
            try:
                self.process_job(job)
            except Exception as e:
                logger.error(f"Error processing job {job.job_id}: {str(e)}")
                job.mark_failed(f"Job processing error: {str(e)}")
    
    def process_job(self, job):
        """Process a single analysis job"""
        self.stdout.write(f'Processing job: {job.job_id} - {job.sample_name}')
        
        # Mark job as running
        job.mark_started()
        
        try:
            if job.analysis_type == 'centrifuge':
                self.run_centrifuge_analysis(job)
            elif job.analysis_type == 'full_pipeline':
                self.run_full_pipeline_analysis(job)
            else:
                raise ValueError(f"Unknown analysis type: {job.analysis_type}")
                
        except Exception as e:
            logger.error(f"Job {job.job_id} failed: {str(e)}")
            job.mark_failed(str(e))
    
    def run_centrifuge_analysis(self, job):
        """Run Centrifuge taxonomic classification"""
        job.update_progress(10, "Preparing input files...")
        
        # Create working directory
        work_dir = self.create_work_directory(job)
        
        try:
            # Prepare input files
            job.update_progress(20, "Copying input files...")
            input_files = self.prepare_input_files(job, work_dir)
            
            # Run Centrifuge
            job.update_progress(30, "Running Centrifuge classification...")
            results = self.run_centrifuge(input_files, work_dir, job)
            
            # Process results
            job.update_progress(80, "Processing results...")
            self.process_centrifuge_results(results, work_dir, job)
            
            # Package results
            job.update_progress(90, "Packaging results...")
            results_path = self.package_results(work_dir, job)
            
            # Mark as completed
            job.mark_completed(results_path)
            
        except Exception as e:
            raise e
        finally:
            # Cleanup temporary files if needed
            pass
    
    def run_full_pipeline_analysis(self, job):
        """Run full MMonitor pipeline"""
        job.update_progress(5, "Starting full pipeline analysis...")
        
        # Create working directory
        work_dir = self.create_work_directory(job)
        
        try:
            # Prepare input files
            job.update_progress(10, "Preparing input files...")
            input_files = self.prepare_input_files(job, work_dir)
            
            # Quality control
            job.update_progress(20, "Running quality control...")
            self.run_quality_control(input_files, work_dir, job)
            
            # Assembly
            job.update_progress(40, "Running genome assembly...")
            assembly_results = self.run_assembly(input_files, work_dir, job)
            
            # Annotation
            job.update_progress(60, "Running gene annotation...")
            self.run_annotation(assembly_results, work_dir, job)
            
            # Taxonomic classification
            job.update_progress(80, "Running taxonomic classification...")
            self.run_centrifuge(input_files, work_dir, job)
            
            # Package results
            job.update_progress(95, "Packaging results...")
            results_path = self.package_results(work_dir, job)
            
            # Mark as completed
            job.mark_completed(results_path)
            
        except Exception as e:
            raise e
    
    def create_work_directory(self, job):
        """Create working directory for job"""
        from django.conf import settings
        
        work_dir = os.path.join(settings.MEDIA_ROOT, 'analysis_jobs', job.job_id)
        os.makedirs(work_dir, exist_ok=True)
        return work_dir
    
    def prepare_input_files(self, job, work_dir):
        """Copy uploaded files to working directory"""
        input_files = []
        
        for uploaded_file in job.uploaded_files.all():
            if not uploaded_file.is_valid:
                continue
                
            # Copy file to work directory
            input_path = os.path.join(work_dir, uploaded_file.original_filename)
            
            with open(uploaded_file.file.path, 'rb') as src:
                with open(input_path, 'wb') as dst:
                    dst.write(src.read())
            
            input_files.append(input_path)
        
        if not input_files:
            raise ValueError("No valid input files found")
        
        return input_files
    
    def run_centrifuge(self, input_files, work_dir, job):
        """Run Centrifuge classification (placeholder)"""
        # This is a placeholder - you would implement actual Centrifuge integration here
        # For now, create dummy results
        
        output_file = os.path.join(work_dir, 'centrifuge_results.txt')
        
        # Simulate processing time
        time.sleep(2)
        
        # Create dummy results
        with open(output_file, 'w') as f:
            f.write("# Centrifuge classification results\n")
            f.write("readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n")
            f.write("read1\tNC_000913.3\t511145\t150\t0\t100\t150\t1\n")
            f.write("read2\tNC_002695.2\t272561\t140\t0\t95\t140\t1\n")
        
        return output_file
    
    def run_quality_control(self, input_files, work_dir, job):
        """Run quality control analysis"""
        # Placeholder for QC analysis
        time.sleep(1)
        
        qc_file = os.path.join(work_dir, 'qc_report.txt')
        with open(qc_file, 'w') as f:
            f.write("Quality Control Report\n")
            f.write("Total reads: 10000\n")
            f.write("Mean quality: 35\n")
    
    def run_assembly(self, input_files, work_dir, job):
        """Run genome assembly"""
        # Placeholder for assembly
        time.sleep(3)
        
        assembly_file = os.path.join(work_dir, 'assembly.fasta')
        with open(assembly_file, 'w') as f:
            f.write(">contig1\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        
        return assembly_file
    
    def run_annotation(self, assembly_file, work_dir, job):
        """Run gene annotation"""
        # Placeholder for annotation
        time.sleep(2)
        
        annotation_file = os.path.join(work_dir, 'annotation.gff')
        with open(annotation_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("contig1\tprokka\tgene\t1\t20\t.\t+\t.\tID=gene1\n")
    
    def process_centrifuge_results(self, results_file, work_dir, job):
        """Process Centrifuge results"""
        # Parse results and create summary
        summary_file = os.path.join(work_dir, 'taxonomy_summary.txt')
        
        with open(summary_file, 'w') as f:
            f.write("Taxonomic Classification Summary\n")
            f.write("================================\n")
            f.write("E. coli: 60%\n")
            f.write("Salmonella: 30%\n")
            f.write("Unclassified: 10%\n")
    
    def package_results(self, work_dir, job):
        """Package results into a ZIP file"""
        import zipfile
        
        results_zip = os.path.join(work_dir, f'{job.sample_name}_results.zip')
        
        with zipfile.ZipFile(results_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(work_dir):
                for file in files:
                    if file.endswith('.zip'):
                        continue  # Don't include the zip file itself
                    
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, work_dir)
                    zipf.write(file_path, arcname)
        
        # Calculate results size
        job.results_size = os.path.getsize(results_zip)
        job.save()
        
        return work_dir

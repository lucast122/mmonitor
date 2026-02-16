from datetime import date
from django.utils import timezone
import hashlib
import os

from django.contrib.auth.models import User
from django.db import models
from django.db.models import JSONField
from django.conf import settings

class Feedback(models.Model):
    user = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE)
    subject = models.CharField(max_length=100)
    message = models.TextField()
    submitted_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"Feedback from {self.user.username} - {self.subject}"




class UserProfile(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    some_field = models.CharField(max_length=200)

    class Meta:
        app_label = 'users'  # If this is the Django app for your MySQL db


class NanoporeRecord(models.Model):
    read_id = models.AutoField(primary_key=True)
    project_id = models.IntegerField()
    sample_id = models.IntegerField()
    taxonomy = models.TextField()
    abundance = models.IntegerField()
    user_id = models.IntegerField()
    
    class Meta:
        db_table = 'users_nanoporerecord'  # Map to your existing mmonitor table
        managed = False  # Don't let Django manage this table
    
    # Add properties for fields that the dashboard expects but don't exist in the table
    @property
    def count(self):
        return float(self.abundance)  # Use abundance as count
    
    @property
    def subproject(self):
        return 'default'
    
    @property
    def date(self):
        return date.today()
    
    @property
    def tax_id(self):
        return 0
    
    @property
    def tax_genus(self):
        return 'Unknown'
    
    @property
    def tax_family(self):
        return 'Unknown'
    
    @property
    def tax_order(self):
        return 'Unknown'
    
    @property
    def tax_class(self):
        return 'Unknown'
    
    @property
    def tax_phylum(self):
        return 'Unknown'
    
    @property
    def tax_superkingdom(self):
        return 'Unknown'
    
    @property
    def tax_clade(self):
        return 'Unknown'
    
    @property
    def tax_subspecies(self):
        return 'Unknown'
    
    @property
    def tax_species_subgroup(self):
        return 'Unknown'
    
    @property
    def tax_species_group(self):
        return 'Unknown'



# class NanoporeRecord(models.Model):
#     read_id = models.AutoField(primary_key=True)
#     taxonomy = models.TextField()
#     abundance = models.FloatField()
#     sample_id = models.CharField(max_length=255)
#     project_id = models.CharField(max_length=255)
#     subproject = models.CharField(max_length=255)
#     date = models.CharField(max_length=255)

# model that gets basic statistics from the sequencing files
class SequencingStatistics(models.Model):
    sample_name = models.CharField(max_length=255)
    project_id = models.CharField(max_length=255)
    subproject_id = models.CharField(max_length=255)
    date = models.CharField(max_length=255)
    mean_read_length = models.FloatField(null=True, blank=True)
    median_read_length = models.FloatField(null=True, blank=True)
    mean_quality_score = models.FloatField(null=True, blank=True)
    mean_gc_content = models.FloatField(null=True, blank=True)
    read_lengths = models.TextField(null=True, blank=True)  # Serialized list of read lengths
    avg_qualities = models.TextField(null=True, blank=True)  # Serialized list of average qualities per length
    number_of_reads = models.IntegerField(null=True, blank=True)
    total_bases_sequenced = models.IntegerField(null=True, blank=True)
    q20_score = models.FloatField(null=True, blank=True)
    q30_score = models.FloatField(null=True, blank=True)
    avg_quality_per_read = models.TextField(null=True, blank=True)  # Serialized list
    base_quality_avg = models.TextField(null=True, blank=True)  # Serialized dictionary
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    gc_contents_per_sequence = models.TextField(null=True, blank=True)

    class Meta:
        verbose_name_plural = "Sequencing Statistics"


class Metadata(models.Model):
    user_id = models.IntegerField(default=0)
    sample_id = models.CharField(max_length=255, default='empty')
    data = JSONField()  # This field will store the dynamic metadata

    class Meta:
        verbose_name_plural = "Metadata"
    @classmethod
    def create_metadata(cls,sample_id,data,user_id):
        metadata = cls(sample_id=sample_id,data=data,user_id=user_id)
        return metadata




class MAG(models.Model):
    name = models.CharField(max_length=255)
    taxonomy = models.JSONField(null=True, blank=True)  
    sample_name = models.CharField(max_length=255, null=True, blank=True)
    gff_file = models.TextField(null=True, blank=True)  
    fasta_file = models.TextField(null=True, blank=True)  
    fai_file = models.TextField(null=True, blank=True)  
    protein_file = models.TextField(null=True, blank=True)  
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    created_at = models.DateTimeField(auto_now_add=True)
    
    # Quality metrics
    completeness = models.FloatField(null=True, blank=True)
    contamination = models.FloatField(null=True, blank=True)
    strain_heterogeneity = models.FloatField(null=True, blank=True)
    gene_count = models.IntegerField(null=True, blank=True)
    trna_count = models.IntegerField(null=True, blank=True)
    rrna_count = models.IntegerField(null=True, blank=True)
    cds_count = models.IntegerField(null=True, blank=True)
    annotation_summary = models.TextField(null=True, blank=True)
    quality_data = models.JSONField(null=True, blank=True)  
    annotations_data = models.JSONField(null=True, blank=True)  # Renamed from annotations to annotations_data

    class Meta:
        app_label = 'users'
        unique_together = ('name', 'user')

    def __str__(self):
        completeness_str = f"{self.completeness:.1f}%" if self.completeness is not None else "N/A"
        return f"{self.name} (Completeness: {completeness_str})"

class MAGAnnotation(models.Model):
    sample = models.ForeignKey(MAG, on_delete=models.CASCADE, related_name='mag_annotations')
    gene_id = models.CharField(max_length=255)
    contig = models.CharField(max_length=255)
    start = models.IntegerField()
    stop = models.IntegerField()
    strand = models.CharField(max_length=1)  # '+' or '-'
    product = models.CharField(max_length=1000)
    kegg = models.CharField(max_length=255, null=True, blank=True)
    ec_numbers = models.JSONField(default=list)  # List of strings
    uniref = models.JSONField(default=list)  # List of strings
    sequence_nt = models.TextField()
    sequence_aa = models.TextField()
    uploaded_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        app_label = 'users'
        indexes = [
            models.Index(fields=['sample', 'gene_id']),
            models.Index(fields=['contig']),
        ]

    def __str__(self):
        return f"{self.gene_id} ({self.product[:50]}...)" if len(self.product) > 50 else f"{self.gene_id} ({self.product})"


# Upload System Models
class UserQuota(models.Model):
    """Track daily user quotas for uploads and job submissions"""
    user = models.OneToOneField(User, on_delete=models.CASCADE, related_name='quota')
    jobs_used_today = models.IntegerField(default=0)
    storage_used_today = models.BigIntegerField(default=0)  # bytes
    last_reset = models.DateField(auto_now_add=True)
    
    # Quota limits
    DAILY_JOB_LIMIT = 3
    DAILY_STORAGE_LIMIT = 2 * 1024 * 1024 * 1024  # 2GB in bytes
    
    class Meta:
        app_label = 'users'
    
    def reset_if_new_day(self):
        """Reset quotas if it's a new day"""
        today = date.today()
        if self.last_reset < today:
            self.jobs_used_today = 0
            self.storage_used_today = 0
            self.last_reset = today
            self.save()
    
    def can_submit_job(self):
        """Check if user can submit another job today"""
        self.reset_if_new_day()
        return self.jobs_used_today < self.DAILY_JOB_LIMIT
    
    def can_upload_size(self, size_bytes):
        """Check if user can upload file of given size"""
        self.reset_if_new_day()
        return (self.storage_used_today + size_bytes) <= self.DAILY_STORAGE_LIMIT
    
    def use_job_quota(self):
        """Increment job quota usage"""
        self.jobs_used_today += 1
        self.save()
    
    def use_storage_quota(self, size_bytes):
        """Increment storage quota usage"""
        self.storage_used_today += size_bytes
        self.save()
    
    @property
    def remaining_jobs(self):
        self.reset_if_new_day()
        return max(0, self.DAILY_JOB_LIMIT - self.jobs_used_today)
    
    @property
    def remaining_storage_mb(self):
        self.reset_if_new_day()
        remaining_bytes = max(0, self.DAILY_STORAGE_LIMIT - self.storage_used_today)
        return remaining_bytes / (1024 * 1024)  # Convert to MB
    
    def __str__(self):
        return f"{self.user.username} - Jobs: {self.jobs_used_today}/{self.DAILY_JOB_LIMIT}, Storage: {self.storage_used_today//1024//1024}MB/{self.DAILY_STORAGE_LIMIT//1024//1024}MB"


def upload_file_path(instance, filename):
    """Generate upload path for user files"""
    # Create hash of filename for uniqueness
    file_hash = hashlib.md5(f"{instance.user.id}_{filename}_{timezone.now().timestamp()}".encode()).hexdigest()[:8]
    return f"uploads/{instance.user.id}/{file_hash}_{filename}"


class UploadedFile(models.Model):
    """Store information about uploaded FASTQ files"""
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='uploaded_files')
    original_filename = models.CharField(max_length=255)
    file = models.FileField(upload_to=upload_file_path)
    file_size = models.BigIntegerField()  # bytes
    file_hash = models.CharField(max_length=64, unique=True)  # SHA256 hash for deduplication
    uploaded_at = models.DateTimeField(auto_now_add=True)
    is_valid = models.BooleanField(default=False)  # Set after validation
    validation_errors = models.TextField(blank=True, null=True)
    
    # File type validation
    ALLOWED_EXTENSIONS = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    MAX_FILE_SIZE = 500 * 1024 * 1024  # 500MB
    
    class Meta:
        app_label = 'users'
        ordering = ['-uploaded_at']
    
    def clean(self):
        """Validate file extension and size"""
        from django.core.exceptions import ValidationError
        
        # Check file extension
        filename_lower = self.original_filename.lower()
        if not any(filename_lower.endswith(ext) for ext in self.ALLOWED_EXTENSIONS):
            raise ValidationError(f"File must be a FASTQ file with extension: {', '.join(self.ALLOWED_EXTENSIONS)}")
        
        # Check file size
        if self.file_size > self.MAX_FILE_SIZE:
            raise ValidationError(f"File size must be less than {self.MAX_FILE_SIZE // 1024 // 1024}MB")
    
    def save(self, *args, **kwargs):
        # Generate file hash if not set
        if not self.file_hash and self.file:
            self.file_hash = self.calculate_file_hash()
        super().save(*args, **kwargs)
    
    def calculate_file_hash(self):
        """Calculate SHA256 hash of file content"""
        if not self.file:
            return ""
        
        hasher = hashlib.sha256()
        try:
            self.file.seek(0)
            for chunk in iter(lambda: self.file.read(4096), b""):
                hasher.update(chunk)
            self.file.seek(0)
            return hasher.hexdigest()
        except:
            return ""
    
    @property
    def file_size_mb(self):
        return self.file_size / (1024 * 1024)
    
    def __str__(self):
        return f"{self.original_filename} ({self.file_size_mb:.1f}MB) - {self.user.username}"


class AnalysisJob(models.Model):
    """Track analysis jobs submitted by users"""
    
    STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('queued', 'Queued'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
        ('cancelled', 'Cancelled'),
    ]
    
    ANALYSIS_TYPES = [
        ('centrifuge', 'Taxonomic Classification (Centrifuge)'),
        ('full_pipeline', 'Full MMonitor Pipeline'),
    ]
    
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='analysis_jobs')
    job_id = models.CharField(max_length=32, unique=True)  # UUID for job tracking
    
    # Job configuration
    sample_name = models.CharField(max_length=255)
    analysis_type = models.CharField(max_length=50, choices=ANALYSIS_TYPES, default='centrifuge')
    uploaded_files = models.ManyToManyField(UploadedFile, related_name='analysis_jobs')
    
    # Job status and timing
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    progress = models.IntegerField(default=0)  # 0-100
    progress_message = models.CharField(max_length=500, blank=True)
    
    created_at = models.DateTimeField(auto_now_add=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    
    # Error handling
    error_message = models.TextField(blank=True, null=True)
    retry_count = models.IntegerField(default=0)
    max_retries = models.IntegerField(default=3)
    
    # Results
    results_path = models.CharField(max_length=500, blank=True, null=True)
    results_size = models.BigIntegerField(default=0)  # bytes
    results_available = models.BooleanField(default=False)
    
    # System information
    worker_id = models.CharField(max_length=100, blank=True, null=True)
    estimated_runtime = models.IntegerField(default=300)  # seconds
    
    class Meta:
        app_label = 'users'
        ordering = ['-created_at']
    
    def save(self, *args, **kwargs):
        # Generate job ID if not set
        if not self.job_id:
            import uuid
            self.job_id = uuid.uuid4().hex
        super().save(*args, **kwargs)
    
    def can_retry(self):
        """Check if job can be retried"""
        return self.status == 'failed' and self.retry_count < self.max_retries
    
    def mark_started(self):
        """Mark job as started"""
        self.status = 'running'
        self.started_at = timezone.now()
        self.save()
    
    def mark_completed(self, results_path=None):
        """Mark job as completed"""
        self.status = 'completed'
        self.completed_at = timezone.now()
        self.progress = 100
        self.progress_message = "Analysis completed successfully"
        if results_path:
            self.results_path = results_path
            self.results_available = True
        self.save()
    
    def mark_failed(self, error_message=None):
        """Mark job as failed"""
        self.status = 'failed'
        self.completed_at = timezone.now()
        if error_message:
            self.error_message = error_message
        self.save()
    
    def update_progress(self, progress, message=""):
        """Update job progress"""
        self.progress = min(100, max(0, progress))
        self.progress_message = message
        self.save()
    
    @property
    def runtime_seconds(self):
        """Get job runtime in seconds"""
        if self.started_at:
            end_time = self.completed_at or timezone.now()
            return (end_time - self.started_at).total_seconds()
        return 0
    
    @property
    def runtime_formatted(self):
        """Get formatted runtime string"""
        seconds = self.runtime_seconds
        if seconds < 60:
            return f"{int(seconds)}s"
        elif seconds < 3600:
            return f"{int(seconds//60)}m {int(seconds%60)}s"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            return f"{hours}h {minutes}m"
    
    def __str__(self):
        return f"{self.sample_name} - {self.get_analysis_type_display()} ({self.status})"


class SystemStatus(models.Model):
    """Track system-wide status and limits"""
    
    # System limits
    max_concurrent_jobs = models.IntegerField(default=5)
    max_queue_size = models.IntegerField(default=50)
    maintenance_mode = models.BooleanField(default=False)
    maintenance_message = models.TextField(blank=True, null=True)
    
    # Current status
    active_jobs = models.IntegerField(default=0)
    queued_jobs = models.IntegerField(default=0)
    total_users = models.IntegerField(default=0)
    total_jobs_processed = models.IntegerField(default=0)
    
    # System health
    last_updated = models.DateTimeField(auto_now=True)
    system_load = models.FloatField(default=0.0)  # 0.0 to 1.0
    disk_usage_percent = models.FloatField(default=0.0)  # 0.0 to 100.0
    
    class Meta:
        app_label = 'users'
        verbose_name_plural = "System Status"
    
    @classmethod
    def get_current(cls):
        """Get or create the current system status"""
        status, created = cls.objects.get_or_create(pk=1)
        return status
    
    def can_accept_job(self):
        """Check if system can accept a new job"""
        if self.maintenance_mode:
            return False
        return (self.active_jobs < self.max_concurrent_jobs and 
                self.queued_jobs < self.max_queue_size)
    
    def update_job_counts(self):
        """Update job counts from database"""
        from django.db.models import Count
        
        job_counts = AnalysisJob.objects.aggregate(
            active=Count('id', filter=models.Q(status__in=['running', 'queued'])),
            queued=Count('id', filter=models.Q(status='queued')),
        )
        
        self.active_jobs = job_counts['active'] or 0
        self.queued_jobs = job_counts['queued'] or 0
        self.save()
    
    def __str__(self):
        return f"System Status - Active: {self.active_jobs}/{self.max_concurrent_jobs}, Queue: {self.queued_jobs}/{self.max_queue_size}"

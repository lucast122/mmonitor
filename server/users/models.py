from datetime import date
from django.utils import timezone

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
    taxonomy = models.TextField(default='empty')
    abundance = models.FloatField(default=0.0)
    count = models.FloatField(default=0.0)
    sample_id = models.CharField(max_length=255, default='empty')
    project_id = models.CharField(max_length=255, default='empty')
    subproject = models.CharField(max_length=255, default='empty')
    date = models.DateField(default=date.today)
    tax_id = models.IntegerField(default=0)
    tax_genus = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_family = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_order = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_class = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_phylum = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_superkingdom = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_clade = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_subspecies = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_species_subgroup = models.CharField(max_length=255, null=True, blank=True, default='empty')
    tax_species_group = models.CharField(max_length=255, null=True, blank=True, default='empty')
    user = models.ForeignKey(User, on_delete=models.CASCADE)



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
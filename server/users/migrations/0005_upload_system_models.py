# Generated manually for upload system models

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion
import users.models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('users', '0004_add_protein_file'),
    ]

    operations = [
        migrations.CreateModel(
            name='SystemStatus',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('max_concurrent_jobs', models.IntegerField(default=5)),
                ('max_queue_size', models.IntegerField(default=50)),
                ('maintenance_mode', models.BooleanField(default=False)),
                ('maintenance_message', models.TextField(blank=True, null=True)),
                ('active_jobs', models.IntegerField(default=0)),
                ('queued_jobs', models.IntegerField(default=0)),
                ('total_users', models.IntegerField(default=0)),
                ('total_jobs_processed', models.IntegerField(default=0)),
                ('last_updated', models.DateTimeField(auto_now=True)),
                ('system_load', models.FloatField(default=0.0)),
                ('disk_usage_percent', models.FloatField(default=0.0)),
            ],
            options={
                'verbose_name_plural': 'System Status',
                'app_label': 'users',
            },
        ),
        migrations.CreateModel(
            name='UserQuota',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('jobs_used_today', models.IntegerField(default=0)),
                ('storage_used_today', models.BigIntegerField(default=0)),
                ('last_reset', models.DateField(auto_now_add=True)),
                ('user', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='quota', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'app_label': 'users',
            },
        ),
        migrations.CreateModel(
            name='UploadedFile',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('original_filename', models.CharField(max_length=255)),
                ('file', models.FileField(upload_to=users.models.upload_file_path)),
                ('file_size', models.BigIntegerField()),
                ('file_hash', models.CharField(max_length=64, unique=True)),
                ('uploaded_at', models.DateTimeField(auto_now_add=True)),
                ('is_valid', models.BooleanField(default=False)),
                ('validation_errors', models.TextField(blank=True, null=True)),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='uploaded_files', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-uploaded_at'],
                'app_label': 'users',
            },
        ),
        migrations.CreateModel(
            name='AnalysisJob',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('job_id', models.CharField(max_length=32, unique=True)),
                ('sample_name', models.CharField(max_length=255)),
                ('analysis_type', models.CharField(choices=[('centrifuge', 'Taxonomic Classification (Centrifuge)'), ('full_pipeline', 'Full MMonitor Pipeline')], default='centrifuge', max_length=50)),
                ('status', models.CharField(choices=[('pending', 'Pending'), ('queued', 'Queued'), ('running', 'Running'), ('completed', 'Completed'), ('failed', 'Failed'), ('cancelled', 'Cancelled')], default='pending', max_length=20)),
                ('progress', models.IntegerField(default=0)),
                ('progress_message', models.CharField(blank=True, max_length=500)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('started_at', models.DateTimeField(blank=True, null=True)),
                ('completed_at', models.DateTimeField(blank=True, null=True)),
                ('error_message', models.TextField(blank=True, null=True)),
                ('retry_count', models.IntegerField(default=0)),
                ('max_retries', models.IntegerField(default=3)),
                ('results_path', models.CharField(blank=True, max_length=500, null=True)),
                ('results_size', models.BigIntegerField(default=0)),
                ('results_available', models.BooleanField(default=False)),
                ('worker_id', models.CharField(blank=True, max_length=100, null=True)),
                ('estimated_runtime', models.IntegerField(default=300)),
                ('uploaded_files', models.ManyToManyField(related_name='analysis_jobs', to='users.uploadedfile')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='analysis_jobs', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-created_at'],
                'app_label': 'users',
            },
        ),
    ]

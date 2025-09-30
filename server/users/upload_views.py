"""
Upload system views for MMonitor Web Upload System
"""
import os
import zipfile
import gzip
import hashlib
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.contrib import messages
from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.core.paginator import Paginator

from .forms import FileUploadForm, AnalysisJobForm, JobFilterForm
from .models import UploadedFile, AnalysisJob, UserQuota, SystemStatus, NanoporeRecord, SequencingStatistics
from django.db.models import Q, Sum
from django.db import models
import json


@login_required
def upload_dashboard(request):
    """Main upload dashboard showing user's quota and recent jobs"""
    # Get or create user quota
    quota, created = UserQuota.objects.get_or_create(user=request.user)
    quota.reset_if_new_day()
    
    # Get recent jobs
    recent_jobs = AnalysisJob.objects.filter(user=request.user)[:5]
    
    # Get system status
    system_status = SystemStatus.get_current()
    
    context = {
        'quota': quota,
        'recent_jobs': recent_jobs,
        'system_status': system_status,
        'can_submit': quota.can_submit_job() and system_status.can_accept_job(),
    }
    
    return render(request, 'users/upload_dashboard.html', context)


@login_required
def upload_files(request):
    """Handle file upload and job creation"""
    quota, created = UserQuota.objects.get_or_create(user=request.user)
    quota.reset_if_new_day()
    
    if request.method == 'POST':
        upload_form = FileUploadForm(request.POST, request.FILES)
        job_form = AnalysisJobForm(request.POST, user=request.user)
        
        if upload_form.is_valid() and job_form.is_valid():
            files = request.FILES.getlist('files')
            
            # Check quota limits
            total_size = sum(f.size for f in files)
            if not quota.can_upload_size(total_size):
                messages.error(request, f"Upload exceeds daily storage limit. You have {quota.remaining_storage_mb:.1f}MB remaining.")
                return render(request, 'users/upload_files.html', {
                    'upload_form': upload_form,
                    'job_form': job_form,
                    'quota': quota
                })
            
            if not quota.can_submit_job():
                messages.error(request, f"Daily job limit reached. You have {quota.remaining_jobs} jobs remaining today.")
                return render(request, 'users/upload_files.html', {
                    'upload_form': upload_form,
                    'job_form': job_form,
                    'quota': quota
                })
            
            # Check system capacity
            system_status = SystemStatus.get_current()
            if not system_status.can_accept_job():
                if system_status.maintenance_mode:
                    messages.error(request, f"System is in maintenance mode: {system_status.maintenance_message}")
                else:
                    messages.error(request, "System is at capacity. Please try again later.")
                return render(request, 'users/upload_files.html', {
                    'upload_form': upload_form,
                    'job_form': job_form,
                    'quota': quota
                })
            
            try:
                # Create uploaded files
                uploaded_files = []
                for uploaded_file in files:
                    # Check for duplicate files
                    file_hash = calculate_file_hash(uploaded_file)
                    existing_file = UploadedFile.objects.filter(file_hash=file_hash).first()
                    
                    if existing_file:
                        # Use existing file if it belongs to the same user
                        if existing_file.user == request.user:
                            uploaded_files.append(existing_file)
                            continue
                        else:
                            # Different user has same file - create new record but don't store duplicate
                            pass
                    
                    # Create new uploaded file
                    new_file = UploadedFile.objects.create(
                        user=request.user,
                        original_filename=uploaded_file.name,
                        file=uploaded_file,
                        file_size=uploaded_file.size,
                        file_hash=file_hash
                    )
                    
                    # Validate file content (basic FASTQ check)
                    try:
                        validate_fastq_file(new_file)
                        new_file.is_valid = True
                        new_file.save()
                    except Exception as e:
                        new_file.is_valid = False
                        new_file.validation_errors = str(e)
                        new_file.save()
                        messages.warning(request, f"File {uploaded_file.name} failed validation: {str(e)}")
                    
                    uploaded_files.append(new_file)
                
                # Create analysis job
                job = job_form.save(commit=False)
                job.user = request.user
                job.save()
                
                # Add uploaded files to job
                job.uploaded_files.set(uploaded_files)
                
                # Update quotas
                quota.use_job_quota()
                quota.use_storage_quota(total_size)
                
                # Queue job for processing
                job.status = 'queued'
                job.save()
                
                messages.success(request, f"Job '{job.sample_name}' created successfully! Job ID: {job.job_id}")
                return redirect('users:job_status', job_id=job.job_id)
                
            except Exception as e:
                messages.error(request, f"Error creating job: {str(e)}")
    
    else:
        upload_form = FileUploadForm()
        job_form = AnalysisJobForm(user=request.user)
    
    context = {
        'upload_form': upload_form,
        'job_form': job_form,
        'quota': quota,
    }
    
    return render(request, 'users/upload_files.html', context)


@login_required
def job_list(request):
    """Display user's analysis jobs with filtering"""
    filter_form = JobFilterForm(request.GET)
    jobs = AnalysisJob.objects.filter(user=request.user)
    
    # Apply filters
    if filter_form.is_valid():
        if filter_form.cleaned_data.get('status'):
            jobs = jobs.filter(status=filter_form.cleaned_data['status'])
        
        if filter_form.cleaned_data.get('analysis_type'):
            jobs = jobs.filter(analysis_type=filter_form.cleaned_data['analysis_type'])
        
        if filter_form.cleaned_data.get('search'):
            search_term = filter_form.cleaned_data['search']
            jobs = jobs.filter(sample_name__icontains=search_term)
    
    # Pagination
    paginator = Paginator(jobs, 10)
    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)
    
    context = {
        'filter_form': filter_form,
        'page_obj': page_obj,
        'jobs': page_obj.object_list,
    }
    
    return render(request, 'users/job_list.html', context)


@login_required
def job_status(request, job_id):
    """Display detailed status for a specific job"""
    try:
        job = AnalysisJob.objects.get(job_id=job_id, user=request.user)
    except AnalysisJob.DoesNotExist:
        messages.error(request, "Job not found or you don't have permission to view it.")
        return redirect('users:job_list')
    
    # Handle job cancellation
    if request.method == 'POST' and 'cancel_job' in request.POST:
        if job.status in ['pending', 'queued']:
            job.status = 'cancelled'
            job.save()
            messages.success(request, f"Job {job.job_id} has been cancelled.")
        else:
            messages.error(request, "Job cannot be cancelled in its current state.")
        return redirect('users:job_status', job_id=job_id)
    
    # Handle job retry
    if request.method == 'POST' and 'retry_job' in request.POST:
        if job.can_retry():
            job.status = 'pending'
            job.retry_count += 1
            job.error_message = None
            job.save()
            messages.success(request, f"Job {job.job_id} has been queued for retry.")
        else:
            messages.error(request, "Job cannot be retried.")
        return redirect('users:job_status', job_id=job_id)
    
    context = {
        'job': job,
        'uploaded_files': job.uploaded_files.all(),
    }
    
    return render(request, 'users/job_status.html', context)


@login_required
def download_results(request, job_id):
    """Download job results as ZIP file"""
    try:
        job = AnalysisJob.objects.get(job_id=job_id, user=request.user)
    except AnalysisJob.DoesNotExist:
        messages.error(request, "Job not found or you don't have permission to access it.")
        return redirect('users:job_list')
    
    if not job.results_available or not job.results_path:
        messages.error(request, "Results are not available for this job.")
        return redirect('users:job_status', job_id=job_id)
    
    try:
        # Create ZIP file with results
        response = HttpResponse(content_type='application/zip')
        response['Content-Disposition'] = f'attachment; filename="{job.sample_name}_results.zip"'
        
        with zipfile.ZipFile(response, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            results_dir = job.results_path
            for root, dirs, files in os.walk(results_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, results_dir)
                    zip_file.write(file_path, arcname)
        
        return response
        
    except Exception as e:
        messages.error(request, f"Error downloading results: {str(e)}")
        return redirect('users:job_status', job_id=job_id)


def demo_page(request):
    """Demo page showcasing MMonitor capabilities"""
    # Get system statistics
    system_status = SystemStatus.get_current()
    system_status.update_job_counts()
    
    # Get some sample statistics
    total_users = User.objects.filter(is_active=True).count()
    total_jobs = AnalysisJob.objects.count()
    completed_jobs = AnalysisJob.objects.filter(status='completed').count()
    
    # Calculate success rate
    success_rate = (completed_jobs / total_jobs * 100) if total_jobs > 0 else 0
    
    context = {
        'system_status': system_status,
        'total_users': total_users,
        'total_jobs': total_jobs,
        'completed_jobs': completed_jobs,
        'success_rate': success_rate,
    }
    
    return render(request, 'users/demo_page.html', context)


# Utility functions
def calculate_file_hash(uploaded_file):
    """Calculate SHA256 hash of uploaded file"""
    hasher = hashlib.sha256()
    uploaded_file.seek(0)
    for chunk in iter(lambda: uploaded_file.read(4096), b""):
        hasher.update(chunk)
    uploaded_file.seek(0)
    return hasher.hexdigest()


def validate_fastq_file(uploaded_file):
    """Basic validation of FASTQ file format"""
    try:
        # Determine if file is gzipped
        is_gzipped = uploaded_file.original_filename.endswith('.gz')
        
        # Read first few lines
        uploaded_file.file.seek(0)
        if is_gzipped:
            with gzip.open(uploaded_file.file, 'rt') as f:
                lines = [f.readline().strip() for _ in range(4)]
        else:
            lines = [uploaded_file.file.readline().decode('utf-8').strip() for _ in range(4)]
        
        # Check FASTQ format (4-line structure)
        if len(lines) < 4:
            raise ValueError("File appears to be empty or too short")
        
        # Check first line starts with @
        if not lines[0].startswith('@'):
            raise ValueError("Invalid FASTQ format: first line should start with '@'")
        
        # Check third line starts with +
        if not lines[2].startswith('+'):
            raise ValueError("Invalid FASTQ format: third line should start with '+'")
        
        # Check sequence and quality lines have same length
        if len(lines[1]) != len(lines[3]):
            raise ValueError("Invalid FASTQ format: sequence and quality lines have different lengths")
        
        uploaded_file.file.seek(0)
        
    except Exception as e:
        uploaded_file.file.seek(0)
        raise ValueError(f"FASTQ validation failed: {str(e)}")


@login_required
def edit_samples(request):
    """Edit samples page - allows users to view, edit, and delete their samples"""
    user_id = request.user.id
    
    # Get all unique samples for the user
    samples = NanoporeRecord.objects.filter(user_id=user_id).values(
        'sample_id', 'project_id'
    ).distinct().order_by('sample_id')
    
    # Get sequencing statistics for additional sample info
    sequencing_stats = SequencingStatistics.objects.filter(user=request.user)
    
    # Create a mapping of sample_id to sequencing stats
    stats_by_sample = {}
    for stat in sequencing_stats:
        stats_by_sample[stat.sample_name] = stat
    
    # Add sequencing stats info to samples
    sample_list = []
    for sample in samples:
        sample_id = sample['sample_id']
        sample_info = {
            'sample_id': sample_id,
            'project_id': sample['project_id'],
            'record_count': NanoporeRecord.objects.filter(
                user_id=user_id, 
                sample_id=sample_id
            ).count(),
            'total_abundance': NanoporeRecord.objects.filter(
                user_id=user_id, 
                sample_id=sample_id
            ).aggregate(total=models.Sum('abundance'))['total'] or 0,
        }
        
        # Try to find matching sequencing stats
        for stat in sequencing_stats:
            if str(stat.sample_name) == str(sample_id):
                sample_info['sequencing_stats'] = stat
                break
        
        sample_list.append(sample_info)
    
    # Pagination
    paginator = Paginator(sample_list, 20)  # Show 20 samples per page
    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)
    
    context = {
        'page_obj': page_obj,
        'total_samples': len(sample_list),
    }
    
    return render(request, 'users/edit_samples.html', context)


@login_required
def edit_sample_detail(request, sample_id):
    """Edit specific sample details"""
    user_id = request.user.id
    
    # Get all records for this sample
    records = NanoporeRecord.objects.filter(user_id=user_id, sample_id=sample_id)
    
    if not records.exists():
        messages.error(request, f"Sample {sample_id} not found or you don't have permission to edit it.")
        return redirect('users:edit_samples')
    
    # Get sequencing stats if available
    try:
        sequencing_stat = SequencingStatistics.objects.get(
            user=request.user, 
            sample_name=str(sample_id)
        )
    except SequencingStatistics.DoesNotExist:
        sequencing_stat = None
    
    if request.method == 'POST':
        # Handle sample updates
        if 'update_sample' in request.POST:
            # Update sequencing statistics if they exist
            if sequencing_stat:
                sequencing_stat.sample_name = request.POST.get('sample_name', sequencing_stat.sample_name)
                sequencing_stat.project_id = request.POST.get('project_id', sequencing_stat.project_id)
                sequencing_stat.subproject_id = request.POST.get('subproject_id', sequencing_stat.subproject_id)
                sequencing_stat.date = request.POST.get('date', sequencing_stat.date)
                sequencing_stat.save()
                messages.success(request, f"Sample {sample_id} updated successfully!")
            else:
                messages.warning(request, f"No sequencing statistics found for sample {sample_id}.")
        
        elif 'delete_sample' in request.POST:
            # Delete all records for this sample
            deleted_count = records.count()
            records.delete()
            messages.success(request, f"Sample {sample_id} and {deleted_count} records deleted successfully!")
            return redirect('users:edit_samples')
        
        return redirect('users:edit_sample_detail', sample_id=sample_id)
    
    # Get unique taxonomies for this sample
    taxonomies = records.values('taxonomy').distinct()
    
    context = {
        'sample_id': sample_id,
        'records': records,
        'sequencing_stat': sequencing_stat,
        'taxonomies': taxonomies,
        'record_count': records.count(),
        'total_abundance': records.aggregate(total=models.Sum('abundance'))['total'] or 0,
    }
    
    return render(request, 'users/edit_sample_detail.html', context)


@login_required
def delete_sample_records(request, sample_id):
    """Delete specific records from a sample"""
    user_id = request.user.id
    
    if request.method == 'POST':
        record_ids = request.POST.getlist('record_ids')
        
        if record_ids:
            # Delete selected records
            deleted_count = NanoporeRecord.objects.filter(
                user_id=user_id,
                sample_id=sample_id,
                read_id__in=record_ids
            ).delete()[0]
            
            messages.success(request, f"Deleted {deleted_count} records from sample {sample_id}")
        else:
            messages.warning(request, "No records selected for deletion.")
    
    return redirect('users:edit_sample_detail', sample_id=sample_id)

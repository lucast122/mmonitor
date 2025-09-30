# from django.db import connections
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.contrib.admin.views.decorators import staff_member_required
from django.db import models
from users.models import NanoporeRecord, MAG, SequencingStatistics


def index(request):
    """Homepage"""
    # Get data statistics for authenticated users
    data_stats = None
    if request.user.is_authenticated:
        try:
            # Get user's data statistics
            user_id = request.user.id
            
            # Count samples (unique sample_ids)
            unique_samples = NanoporeRecord.objects.filter(user_id=user_id).values('sample_id').distinct().count()
            
            # Count unique taxa
            unique_taxa = NanoporeRecord.objects.filter(user_id=user_id).values('taxonomy').distinct().count()
            
            # Count total reads (sum of abundance)
            total_reads = NanoporeRecord.objects.filter(user_id=user_id).aggregate(
                total_reads=models.Sum('abundance')
            )['total_reads'] or 0
            
            # Count MAGs
            total_mags = MAG.objects.filter(user=request.user).count()
            
            # Count sequencing statistics records
            total_sequencing_runs = SequencingStatistics.objects.filter(user=request.user).count()
            
            data_stats = {
                'samples': unique_samples,
                'taxa': unique_taxa,
                'reads': total_reads,
                'mags': total_mags,
                'sequencing_runs': total_sequencing_runs
            }
        except Exception as e:
            # If there's an error getting stats, set to None
            data_stats = None
    
    return render(request, "main/index.html", {'data_stats': data_stats})

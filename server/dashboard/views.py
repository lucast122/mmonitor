from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, Http404
from django.shortcuts import render, get_object_or_404

from users.models import MAG


def get_user_id(request):
    """Get user id from request."""
    return request.user.id


def spa_view(request, path=''):
    """Serve the React SPA index.html for all dashboard routes.

    No @login_required — the React frontend handles its own
    authentication via JWT tokens against the /api/v1/ endpoints.
    """
    return render(request, 'index.html')


@login_required
def serve_gff_file(request, id):
    """Serve GFF3 file content for a MAG."""
    mag = get_object_or_404(MAG, pk=id)
    if mag.gff_file:
        response = HttpResponse(mag.gff_file, content_type='text/plain')
        response['Content-Disposition'] = f'attachment; filename="{mag.name}.gff3"'
        return response
    raise Http404("GFF3 file not found")


def extract_fasta_from_gff3(gff_content):
    """Extract FASTA sequences from GFF3 file content."""
    if not gff_content:
        return None

    # Find the FASTA section which starts with ##FASTA
    parts = gff_content.split('##FASTA')
    if len(parts) < 2:
        return None

    return parts[1].strip()


@login_required
def serve_fasta_from_gff(request, id):
    """Extract and serve FASTA content from GFF3 file."""
    mag = get_object_or_404(MAG, pk=id)
    if mag.gff_file:
        fasta_content = extract_fasta_from_gff3(mag.gff_file)
        if fasta_content:
            response = HttpResponse(fasta_content, content_type='text/plain')
            response['Content-Disposition'] = f'attachment; filename="{mag.name}.fasta"'
            return response
    raise Http404("FASTA content not found in GFF3 file")

from os.path import join, basename
import os

from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse, HttpResponse, Http404
from django.shortcuts import render, get_object_or_404
from django.views.generic import TemplateView
from django.contrib.auth.mixins import LoginRequiredMixin

from users.models import MAG
from .dashapp.diversity import Diversity
from .dashapp.taxonomy import Taxonomy
from .dashapp.qc import QC
from .dashapp.index import Index


def get_user_id(request):
    """Get user id from request."""
    return request.user.id


@login_required
def index(request):
    """Index view."""
    # Create and register the Index app with the new name
    app_name = 'dashboard_index'
    index_instance = Index(user_id=request.user.id)
    
    # Pass simple context to the template
    context = {
        'name': app_name,  # Changed to 'name' to match template parameter
        'user_id': request.user.id,
    }
    
    return render(request, 'dashboard/dashapp.html', context)


@login_required
def load_app(request, name):
    """Load different Dash apps by name."""
    user_id = request.user.id
    
    # Map old app name to new app name if needed
    app_name_mapping = {
        'Index': 'dashboard_index',
        'horizon': 'horizon'
        # Add other app name mappings here if needed
    }
    
    # Use the mapped name if available, otherwise use the original name
    display_name = app_name_mapping.get(name, name)
    
    # Initialize the appropriate app based on the name
    if name == 'Index':
        index_instance = Index(user_id=user_id)
    elif name == 'horizon':
        from .dashapp.horizon import Horizon
        horizon_instance = Horizon(user_id=user_id)
    
    # Return with the correct name for the template
    return render(request, 'dashboard/dashapp.html', context={'name': display_name})


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


@login_required
def horizon(request):
    """Render page for horizon plot."""
    if request.method == 'POST':
        try:
            from .dashapp.horizon import Horizon
            horizon = Horizon(user_id=request.user.id)
            data = request.POST
            project_id = data.get('project_id')
            taxa_count = int(data.get('taxa_count', 20))
            plot_data = horizon._get_cached_horizon_plot(project_id, taxa_count)
            return JsonResponse({'plot_data': plot_data})
        except Exception as e:
            return JsonResponse({'error': str(e)}, status=500)
    else:
        from .dashapp.horizon import Horizon
        horizon = Horizon(user_id=request.user.id)
        context = {
            'projects': horizon.projects,
            'plot_html': horizon._get_horizon_plot_html()
        }
        return render(request, 'dashboard/horizon.html', context)


class DashView(LoginRequiredMixin, TemplateView):
    template_name = 'dash_app.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        app_name = kwargs.get('app_name', '')
        context.update({
            'url_base_pathname': f'/dash/{app_name}/',
            'requests_pathname_prefix': f'/dash/{app_name}/',
        })
        return context
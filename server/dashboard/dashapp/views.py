from django.contrib.auth.decorators import login_required
from django.core.serializers import serialize
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render, get_object_or_404
from dashboard.dashapp.index import Index
from users.models import NanoporeRecord, MAG  # Adjust the import according to your model's location
from os.path import basename

@login_required
def chart_data_api(request):
    if request.user.is_authenticated:
        user_id = request.user.id

    print(f"user id in view: {user_id}")
    chart_data_queryset = NanoporeRecord.objects.filter(user_id=user_id)
    chart_data_json = serialize('json', chart_data_queryset)
    return HttpResponse(chart_data_json, content_type="application/json")


@login_required
def get_user_id(request):
    # Get the session ID from the URL parameters
    if request.user.is_authenticated:
        user_id = request.user.id
        # If a user ID was found, return it
        if user_id:
            return JsonResponse({'user_id': user_id})

    # If no session ID was provided, or no user ID was found, return an error
    return JsonResponse({'error': 'No session ID provided or no user found'}, status=400)


@login_required
def index(request):
    return render(request, 'dashboard/index.html')


@login_required
def load_app(request, name):
    user_id = request.user.id  # Assuming the user is authenticated
    index_instance = Index(user_id=user_id)
    return render(request, 'dashboard/dashapp.html', context={'name': name})

@login_required
def your_view(request):
    if request.user.is_authenticated:
        user_id = request.user.id

    return render(request, 'base.html', {'user_id': user_id})
def serve_gff_file(request, mag_id):
    mag = get_object_or_404(MAG, id=mag_id)
    response = HttpResponse(mag.gff_file, content_type='text/plain')
    # Remove the filename from the Content-Disposition header
    response['Content-Disposition'] = 'attachment'
    return response

def serve_fasta_file(request, mag_id):
    mag = get_object_or_404(MAG, id=mag_id)
    response = HttpResponse(mag.fasta_file, content_type='text/plain')
    # Remove the filename from the Content-Disposition header
    response['Content-Disposition'] = 'attachment'
    return response

def serve_fai_file(request, mag_id):
    mag = get_object_or_404(MAG, id=mag_id)
    response = HttpResponse(mag.fai_file, content_type='text/plain')
    # Remove the filename from the Content-Disposition header
    response['Content-Disposition'] = 'attachment'
    return response